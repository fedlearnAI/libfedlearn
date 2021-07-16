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
    NTL::ZZ P, Q;
    NTL::ZZ p, q, N;
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
        NTL::GenPrime(p, len);
        NTL::GenPrime(q, len);
        NTL::mul(N, p, q);
        // MAKE SURE that P is larger than secret
        NTL::GenPrime(P, len * 2 + 1);
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

void recv_and_reorganize_(const std::vector<LIST_OF_CHAR_LIST> &share_lst_b,
                          const LIST_OF_BOOL_LIST &share_lst_b_neg,
                          LIST_OF_CHAR_LIST &reorganized_b,
                          LIST_OF_BOOL &reorganized_b_neg,
                          LIST_OF_LONG &reorganized_a,
                          long n, long i) {
    reorganized_a.reserve(n);
    reorganized_b.reserve(n);
    reorganized_b_neg.reserve(n);

    for (long j = 1; j <= n; j++) {
        reorganized_b.push_back(share_lst_b[j - 1][i - 1]);
        reorganized_a.push_back(j);
        reorganized_b_neg.push_back(share_lst_b_neg[j - 1][i - 1]);
    }
}

void __pailler_dec__(NTL::ZZ &plain,
                     const NTL::ZZ &cypher_text,
                     const NTL::ZZ lambda,
                     const dist_paillier_privkey_t &priv_key) {
    NTL::ZZ lambda_invmod;
    NTL::InvMod(lambda_invmod, lambda, priv_key.n);
    NTL::PowerMod(plain, cypher_text, lambda, priv_key.n_squared);
    assert((plain - 1) % priv_key.n == 0);
    NTL::div(plain, plain - 1, priv_key.n);
    NTL::MulMod(plain, plain, lambda_invmod, priv_key.n);
}

bool test_char_list_ZZ_transformation_once(long length) {
    NTL::ZZ zz, zz_res, tmp;
    NTL::RandomLen(zz, length);
    NTL::RandomLen(tmp, length / 2);
    zz -= tmp;
    CHAR_LIST char_lst_zz, char_lst_zzp;
    bool is_negative;
    distributed_paillier::ZZ_2_byte(char_lst_zz, zz, is_negative);
    distributed_paillier::char_list_2_ZZ(zz_res, char_lst_zz, is_negative);

    return zz == zz_res;
}

TEST_F(params, test_char_list_ZZ_transformation_manyTimes) {
    SetUp(1024, 3, 1);
    std::vector<long> length_lst = {256, 512, 1024, 2048, 4096};
    for (auto length : length_lst) {
        NTL::ZZ P;
        NTL::GenPrime(P, length);
        NTL::ZZ_p::init(P);
        for (int n = 1; n < 10000; n++) {
            ASSERT_TRUE(test_char_list_ZZ_transformation_once(length));
        }
    }
}

TEST_F(params, test_standalone_key_gene) {
    SetUp(1024, 3, 1);
    std::vector<dist_paillier_privkey_t> pri_key_lst = distributed_paillier::gene_privkey_standalone(length, n, t);
    ASSERT_EQ(pri_key_lst.size(), 2 * t + 1);
}

// suppose the co-prime number N is known
TEST(test_enc_dec_arithmatic, test_generate_priv_key) {
    const long len = 512;
    const long n = 3;
    const long t = 1;
    const long plain_text_num = 100;
    // 2^N must be larger than sqr of the largest number in plain text
    // assert( (long)2<<len > (plain_text_num*plain_text_num));

    NTL::ZZ nfct = factorial_small(n);
    NTL::ZZ p, q, P;
    NTL::GenPrime(p, len);
    NTL::GenPrime(q, len);
    NTL::GenPrime(P, len);
    NTL::ZZ N = p * q;
    std::cout << "\np, q, N, P = " << p << " " << q << " " << N << " " << P << endl;

    NTL::Vec<NTL::ZZ> lambda_lst, beta_lst;
    NTL::Vec<NTL::ZZ> pi_lst, qi_lst;

    CHAR_LIST modulu_char;
    bool modulu_char_isneg = false;
    distributed_paillier::ZZ_2_byte(modulu_char, N, modulu_char_isneg);

    std::vector<LIST_OF_CHAR_LIST> lambda_share, beta_share;
    lambda_share.reserve(n);
    beta_share.reserve(n);
    LIST_OF_BOOL_LIST lambda_share_sign, beta_shares_sign;
    lambda_share_sign.reserve(n);
    beta_shares_sign.reserve(n);

    std::vector<distributed_paillier> paillier_class_lst;
    paillier_class_lst.reserve(n);

    ZZ p_sum(0), q_sum(0);
    for (long i = 1; i <= n; i++) {
        paillier_class_lst.emplace_back(modulu_char, n, t);
    }
    for (long i = 1; i < n; i++) {
        pi_lst.append(NTL::RandomBnd(p / n));
        qi_lst.append(NTL::RandomBnd(q / n));
        p_sum += pi_lst[i - 1];
        q_sum += qi_lst[i - 1];
    }
    pi_lst.append(p - p_sum);
    qi_lst.append(q - q_sum);

    // check qi, pi correctness
    p_sum = 0;
    q_sum = 0;
    for (long i = 0; i < n; i++) {
        p_sum += pi_lst[i];
        q_sum += qi_lst[i];
        cout << "pi = " << pi_lst[i] << " qi = " << qi_lst[i] << endl;
        assert((pi_lst[i] > 0) && (qi_lst[i] > 0));
    }
    ASSERT_TRUE(q_sum == q && p_sum == p);

    std::cout << "\np, q, N, P = " << p << " " << q << " " << N << " " << P << endl;
    std::cout << N << endl;

    for (long i = 1; i <= n; i++) {
        NTL::ZZ tmp;
        do {
            tmp = NTL::RandomBnd(N);
            tmp = 1;
        } while (tmp == 0);
        beta_lst.append(tmp);

        if (i == 1) {
            lambda_lst.append(N - pi_lst[i - 1] - qi_lst[i - 1] + 1);
        } else {
            lambda_lst.append(-pi_lst[i - 1] - qi_lst[i - 1]);
        }

        JAVA_CHAR lam_char, beta_char;
        std::vector<long> byte_arr_size_lambdai(n), byte_arr_size_betai((n));
        LIST_OF_BOOL is_negative_lst_lambda(n);
        LIST_OF_BOOL is_negative_lst_beta(n);
        distributed_paillier::ZZ_2_byte(lam_char.byte_arr,
                                        lambda_lst[i - 1],
                                        lam_char.is_negative);
        distributed_paillier::ZZ_2_byte(beta_char.byte_arr,
                                        beta_lst[i - 1],
                                        beta_char.is_negative);
        lambda_share.push_back(
                paillier_class_lst[i - 1].create_lambda_share(lam_char.byte_arr,
                                                              lam_char.is_negative,
                                                              byte_arr_size_lambdai,
                                                              is_negative_lst_lambda));
        beta_share.push_back(
                paillier_class_lst[i - 1].create_beta_share(beta_char.byte_arr,
                                                            beta_char.is_negative,
                                                            byte_arr_size_betai,
                                                            is_negative_lst_beta));

        lambda_share_sign.push_back(is_negative_lst_lambda);
        beta_shares_sign.push_back(is_negative_lst_beta);
    }

    // 1) simulate network transfer -- transferring shares to each party
    // 2) compute Lambda*beta
    // NB: be careful with the size of reconstruction party.
    NTL::Vec<ZZ> theta_lst, hi;
    for (long i = 1; i <= 2 * t + 1; i++) {
        LIST_OF_CHAR_LIST tmp_res(2);
        CHAR_LIST lamb_times_beta;
        CHAR_LIST theta;
        bool lamb_times_beta_isneg, theta_is_neg;
        vector<long> ret_char_size(n);
        // lambda recv
        LIST_OF_CHAR_LIST reorganized_lambda_b;
        LIST_OF_BOOL reorganized_lambda_b_neg;
        LIST_OF_LONG reorganized_a;
        recv_and_reorganize_(lambda_share, lambda_share_sign,
                             reorganized_lambda_b, reorganized_lambda_b_neg,
                             reorganized_a, n, i);

        // beta recv
        LIST_OF_CHAR_LIST reorganized_beta_b;
        LIST_OF_BOOL reorganized_beta_b_neg;
        LIST_OF_LONG reorganized_abeta;
        recv_and_reorganize_(beta_share, beta_shares_sign,
                             reorganized_beta_b, reorganized_beta_b_neg,
                             reorganized_abeta, n, i);

        long scalar = NTL::to_long(nfct);

        tmp_res = paillier_class_lst[i - 1].compute_lambda_times_beta_share(
                reorganized_lambda_b, reorganized_beta_b,
                reorganized_lambda_b_neg, reorganized_beta_b_neg,
                modulu_char, i, ret_char_size, lamb_times_beta_isneg,
                theta_is_neg, scalar);

        lamb_times_beta = tmp_res[0];

        // NOTE: the scaling factor is computed during secret sharing.
        //       this param should be passed out for computing theta here
        theta = tmp_res[1];

        // convert to zz type
        NTL::ZZ tmp1, tmp2;
        distributed_paillier::char_list_2_ZZ(tmp1, theta, theta_is_neg);
        theta_lst.append(tmp1);
        distributed_paillier::char_list_2_ZZ(tmp2, lamb_times_beta, lamb_times_beta_isneg);
        hi.append(tmp2);
    }

    // 1) compute theta
    // 2) check shares of lamb_times_beta can reconstruct real lamb_times_beta
    NTL::ZZ recon_lamb_times_beta(0);
    NTL::ZZ recon_theta(0);
    for (long i = 1; i <= 2 * t + 1; i++) {
        NTL::ZZ li, tmp;
        NTL::ZZ enume = nfct;
        NTL::ZZ denom = NTL::ZZ(1);
        for (long j = 1; j <= 2 * t + 1; j++) {
            if (i != j) {
                enume *= j;
                denom *= j - i;
            }
        }
        NTL::div(li, enume, denom);
        NTL::mul(tmp, li, hi[i - 1]);
        NTL::add(recon_lamb_times_beta, recon_lamb_times_beta, tmp);
        NTL::AddMod(recon_theta, recon_theta, li * theta_lst[i - 1], N);
    }
    NTL::ZZ scaler_inv;
    NTL::InvMod(scaler_inv, nfct * sqr(nfct), N);
    NTL::MulMod(recon_theta, recon_theta, scaler_inv, N);
    NTL::MulMod(recon_theta, recon_theta, nfct * sqr(nfct), N);

    NTL::ZZ lambda_lst_sum(0), beta_lst_sum(0);
    for (int i = 0; i < lambda_lst.length(); i++) {
        lambda_lst_sum += lambda_lst[i];
        beta_lst_sum += beta_lst[i];
    }

    ZZ real_lamb_times_beta = lambda_lst_sum * beta_lst_sum;
    std::cout << "\nrecon_lamb_times_beta  = " << recon_lamb_times_beta
              << "\n real_lamb_times_beta  = " << real_lamb_times_beta * nfct * sqr(nfct)
              << "\nrecon_theta = " << recon_theta
              << "\nrecon_lamb_times_beta % N = " << recon_lamb_times_beta % N
              << "\nN = " << N
              << endl;

    // NB: These two assertion is VERY IMPORTANT!
    ASSERT_TRUE(lambda_lst_sum * beta_lst_sum * nfct * sqr(nfct) == recon_lamb_times_beta);
    ASSERT_TRUE(recon_theta == recon_lamb_times_beta % N);
    std::cout << "key generation finished. Correct\n"
              << endl;

    // generate priv and pub key for each party
    dist_paillier_pubkey_t pub_key;
    std::vector<dist_paillier_privkey_t> pri_key_lst;
    for (long i = 1; i <= 2 * t + 1; i++) {
        dist_paillier_privkey_t pri_key;
        pri_key.n = N;
        pub_key.n = N;
        pri_key.n_squared = NTL::sqr(N);
        pub_key.n_squared = NTL::sqr(N);
        pub_key.bits = len;
        pub_key.n_plusone = N + 1;
        pri_key.hi = hi[i - 1];
        pri_key.t = t;
        pri_key.n_fat = NTL::to_long(nfct);
        NTL::InvMod(pri_key.theta_invmod, recon_theta, N);
        pri_key_lst.push_back(pri_key);
    }

    // generate numbers and do computation
    std::vector<plain_text_t> plain_lst, plain_lst_2;
    std::vector<cypher_text_t> cypher_lst, cypher_lst_2;
    plain_lst.reserve(plain_text_num);
    plain_lst_2.reserve(plain_text_num);
    cypher_lst.reserve(plain_text_num);
    cypher_lst_2.reserve(plain_text_num);

    std::cout << "doing enc..." << endl;
    time_t begin_enc = get_time_stamp();
    for (int i = 1; i <= plain_text_num; i++) {
        // cout << i << endl;
        plain_text_t tmp(len);
        cypher_text_t tmp2(len);
        distributed_paillier::ZZ_2_byte(tmp.plain, NTL::to_ZZ(i), tmp.is_negative);
        plain_lst.push_back(tmp);
        distributed_paillier::enc(tmp2, tmp, pub_key);
        cypher_lst.push_back(tmp2);
    }
    time_t end_enc = get_time_stamp();

    // test enc correctness
    for (long i = 1; i <= plain_text_num; i++) {
        NTL::ZZ c;
        distributed_paillier::char_list_2_ZZ(c, cypher_lst[i - 1].c, cypher_lst[i - 1].is_neg);
        NTL::ZZ test;
        NTL::PowerMod(test, c, recon_lamb_times_beta, NTL::sqr(N));
        if ((test - 1) % N != 0) {
            ZZ beta_sum(0);
            for (const auto &b : beta_lst)
                beta_sum += b;
            std::cout << "Enc is not correct, check pub_key. sum(pi)*sum(qi)*beta_sum may not be too large" << endl;
            std::cout << "test = " << test << " c = " << c << endl;
            std::cout << "(p-1)*(q-1)*beta_sum = " << (p - 1) * (q - 1) * beta_sum
                      << "(p-1)*(q-1) = " << (p - 1) * (q - 1)
                      << " recon_lamb_times_beta = " << recon_lamb_times_beta
                      << " original = " << i
                      << endl;

            NTL::ZZ test1;
            NTL::PowerMod(test1, c, (p - 1) * (q - 1), NTL::sqr(N));
            assert((test1 - 1) % N == 0);
            assert((test - 1) % N == 0);
        }
    }
    std::cout << "enc finished. Correct\n";
    std::cout << "Enc time : ";
    get_duration(begin_enc, end_enc, true);

    // test dec correctness
    int cnt = 0;
    time_t begin = get_time_stamp();
    double inner_dur = 0;
    for (const auto &c : cypher_lst) {
        plain_text_t decrypted_p(len);
        std::vector<cypher_text_t> partial_cypher_lst;
        partial_cypher_lst.reserve(2 * t + 1);
        time_t begin_inner = get_time_stamp();
        for (long i = 1; i <= 2 * t + 1; i++) {
            cypher_text_t partial_cypher(len);
            distributed_paillier::dec_partial(i, t, partial_cypher, c, pri_key_lst[i - 1]);
            partial_cypher_lst.push_back(partial_cypher);
        }

        NTL::ZZ decrypted_zz, original_zz;
        decrypted_zz = distributed_paillier::dec_final(c, partial_cypher_lst, t,
                                                       to_ZZ(1 << 16), pri_key_lst[0]);
        distributed_paillier::char_list_2_ZZ(original_zz, plain_lst[cnt].plain, plain_lst[cnt].is_negative);
        time_t end_inner = get_time_stamp();
        inner_dur += get_duration(begin_inner, end_inner, false);
        std::cout << "decrypted_p = " << decrypted_zz << "\toriginal_zz = " << original_zz << endl;
        assert(decrypted_zz == original_zz);
        cnt += 1;
    }
    time_t end = get_time_stamp();
    std::cout << "Dec finished. Correct. ";
    get_duration(begin, end, true);
    std::cout << "get rid of transferring: " << inner_dur << endl;

    /**
     * test adding time and correctness
     */
    std::vector<cypher_text_t> c_res_lst;
    std::vector<plain_text_t> plain_res;
    c_res_lst.reserve(cypher_lst.size());
    time_t add_start = get_time_stamp();
    for (const auto &c : cypher_lst) {
        cypher_text_t c_res(len);
        distributed_paillier::add(c_res, c, c, pub_key);
        c_res_lst.push_back(c_res);
    }
    time_t add_end = get_time_stamp();
    __dec_list__(plain_res, c_res_lst, pri_key_lst, t, len);
    for (int i = 0; i < plain_text_num; i++) {
        NTL::ZZ decrypted_zz, original_zz;
        distributed_paillier::char_list_2_ZZ(decrypted_zz, plain_res[i].plain, plain_res[i].is_negative);
        distributed_paillier::char_list_2_ZZ(original_zz, plain_lst[i].plain, plain_lst[i].is_negative);
        assert(decrypted_zz == original_zz + original_zz);
    }
    c_res_lst.clear();
    plain_res.clear();
    std::cout << "Add finished. Correct.";
    get_duration(add_start, add_end, true);

    /**
     * test mult time and correctness
     */
    c_res_lst.reserve(cypher_lst.size());
    time_t mul_start = get_time_stamp();
    cnt = 0;
    for (const auto &c : cypher_lst) {
        cypher_text_t c_res(len);
        distributed_paillier::mul(c_res, c, plain_lst[cnt], pub_key);
        c_res_lst.push_back(c_res);
        cnt++;
    }
    time_t mul_end = get_time_stamp();
    __dec_list__(plain_res, c_res_lst, pri_key_lst, t, len);
    for (int i = 0; i < plain_text_num; i++) {
        NTL::ZZ decrypted_zz, original_zz;
        distributed_paillier::char_list_2_ZZ(decrypted_zz, plain_res[i].plain, plain_res[i].is_negative);
        distributed_paillier::char_list_2_ZZ(original_zz, plain_lst[i].plain, plain_lst[i].is_negative);
        assert(decrypted_zz == original_zz * original_zz);
    }
    c_res_lst.clear();
    plain_res.clear();
    std::cout << "Mul finished. Correct. ";
    get_duration(mul_start, mul_end, true);

    /**
     * test inner-product time
     */
    c_res_lst.reserve(cypher_lst.size());
    time_t innerp_start = get_time_stamp();
    for (int i = 0; i < plain_text_num; i++) {
        cypher_text_t c_res(len);
        distributed_paillier::mul(c_res, cypher_lst[i], plain_lst[i], pub_key);
        c_res_lst.push_back(c_res);
        ZZ tmp(0);
        distributed_paillier::char_list_2_ZZ(tmp, c_res_lst[i].c, c_res_lst[i].is_neg);
    }
    cypher_text_t inp_res(len);
    ZZ tmp(1);
    distributed_paillier::ZZ_2_byte(inp_res.c, tmp, inp_res.is_neg);
    for (int i = 0; i < plain_text_num; i++) {
        distributed_paillier::add(inp_res, inp_res, c_res_lst[i], pub_key);
        ZZ tmpp, tmppp;
        distributed_paillier::char_list_2_ZZ(tmpp, c_res_lst[i].c, c_res_lst[i].is_neg);
        distributed_paillier::char_list_2_ZZ(tmppp, inp_res.c, inp_res.is_neg);
        // std::cout << tmpp << "\t" << tmppp << endl;
    }
    distributed_paillier::char_list_2_ZZ(tmp, inp_res.c, inp_res.is_neg);
    std::cout << tmp << endl;
    time_t innerp_end = get_time_stamp();
    std::cout << "Inner-product time : ";
    get_duration(innerp_start, innerp_end, true);
}

int main() {
    // test_char_list_ZZ_transformation();
    // test_N_generation(1024, 3, 1);
    //test_generate_priv_key(32, 3, 1, 100);

    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}

