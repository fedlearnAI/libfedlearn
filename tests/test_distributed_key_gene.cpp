/* Copyright 2020 The FedLearn Authors. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <iostream>
#include <vector>
#include "../include/NTL/ZZ.h"
#include "../include/gtest/gtest.h"
#include "../include/distributed_paillier.h"

class params : public ::testing::Test {
public:
    long length;
    long n, t;

    params() {
        length = 0;
        n = 0;
        t = 0;
    };

    void SetUp(long len, int num_party, int assessor) {
        this->length = len;
        this->n = num_party;
        this->t = assessor;
    }

    ~params() override = default;
};

void __dec_list__(std::vector<plain_text_t> &decrypted_lst,
                  const std::vector<cypher_text_t> &cypher_lst,
                  std::vector<dist_paillier_privkey_t> pri_key_lst,
                  int t,
                  int cypher_len) {
    int cnt = 0;
    decrypted_lst.reserve(cypher_lst.size());
    for (const auto &c : cypher_lst) {
        plain_text_t decrypted_p(cypher_len);
        std::vector<cypher_text_t> partial_cypher_lst;
        partial_cypher_lst.reserve(2 * t + 1);
        for (long i = 1; i <= 2 * t + 1; i++) {
            cypher_text_t partial_cypher(cypher_len);
            distributed_paillier::dec_partial(i, t, partial_cypher, c, pri_key_lst[i - 1]);
            partial_cypher_lst.push_back(partial_cypher);
        }
        distributed_paillier::dec_final(decrypted_p, c, partial_cypher_lst, t, to_ZZ(1 << 16), pri_key_lst[0]);
        decrypted_lst.push_back(decrypted_p);
        cnt += 1;
    }
}

 ZZ get_sumP_share(const std::vector<ZZ>& reorganized_pb, const LIST_OF_LONG& pa_lst, long i) {
    ZZ tmp = ZZ(0);
    for(const auto& pib: reorganized_pb) {
        tmp += pib;
    }
    for(const auto & pia: pa_lst) {
        assert(i == pia);
    }
    return tmp;
}

TEST_F(params, TestDistKeyGene) {
    long bit_len = 10;
    SetUp(bit_len, 3, 1);
    std::vector<NTL::ZZ> lambda_lst, beta_lst;
    std::vector<NTL::ZZ> pi_lst, qi_lst;
    std::vector<ShamirShares> pi_share_lst, qi_share_lst;
    ZZ P = NTL::GenPrime_ZZ(length*2);

    cout << "P = " << P << endl;
    distributed_paillier key_generator = distributed_paillier(P, n, t, length);

    // generate pi qi and pi share qi share
    ZZ sum_p, sum_q;
    sum_p = ZZ(0);
    sum_q = ZZ(0);
    for(auto i = 0; i < n; i++) {
        pi_lst.emplace_back(key_generator.gene_local_piqi()[0]);
        qi_lst.emplace_back(key_generator.gene_local_piqi()[1]);
        ASSERT_TRUE((pi_lst[i] > 0) && (qi_lst[i] > 0));

        ShamirSecretSharingScheme scheme = ShamirSecretSharingScheme(P, n, t);
        ShamirShares share_p = scheme.share_secret(pi_lst[i]);
        ShamirShares share_q = scheme.share_secret(qi_lst[i]);
        pi_share_lst.emplace_back(share_p);
        qi_share_lst.emplace_back(share_q);

        sum_p += pi_lst[i];
        sum_q += qi_lst[i];
        cout << pi_lst[i] << " " << qi_lst[i] << " " << endl;
    }
    cout << "sum_p = " << sum_p << " sum_q = " << sum_q << " sum_q*sum_p = "  << sum_q*sum_p << endl;

    NTL::Vec<NTL::Pair<long, ZZ>>  N_share, sum_p_share, sum_q_share;
    for (long i = 1; i <= n; i++) {
        LIST_OF_LONG pa_lst, qa_lst;
        std::vector<ZZ> reorganized_pb, reorganized_qb;
        distributed_paillier::__recv_and_reorganize__(pi_share_lst, reorganized_pb, pa_lst, n, i);
        distributed_paillier::__recv_and_reorganize__(qi_share_lst, reorganized_qb, qa_lst, n, i);

        ZZ tmp = NTL::MulMod(get_sumP_share(reorganized_qb, qa_lst, i), get_sumP_share(reorganized_pb, pa_lst, i), P);
        N_share.append( NTL::Pair<long, ZZ>(i,  tmp));

        sum_q_share.append(NTL::Pair<long, ZZ>(i, get_sumP_share(reorganized_qb, qa_lst, i)));
        sum_p_share.append(NTL::Pair<long, ZZ>(i, get_sumP_share(reorganized_pb, pa_lst, i)));
    }
    cout << ShamirShares(sum_q_share, P, n, t).reconstruct_secret(P) << " "
    << ShamirShares(sum_p_share, P, n, t).reconstruct_secret(P) << " "
    << ShamirShares(sum_q_share, P, n, t).__mult__(ShamirShares(sum_p_share, P, n, t)).reconstruct_secret(P)
    << endl;

    ASSERT_EQ(sum_q*sum_p, ShamirShares(sum_q_share, P, n, t).__mult__(ShamirShares(sum_p_share, P, n, t)).reconstruct_secret(P));
}



int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}

