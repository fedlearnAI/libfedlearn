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

#ifndef DIST_PAI_H
#define DIST_PAI_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <vector>
#include <chrono>
#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "../include/distributed_paillier.h"
#include "shamir_secret_sharing.h"
#include "shamir_secret_sharing_integer.h"

using namespace std;
using namespace NTL;

#define time_t std::chrono::steady_clock::time_point

typedef std::vector<unsigned char> CHAR_LIST;
typedef std::vector<long> LIST_OF_LONG;
typedef std::vector<bool> LIST_OF_BOOL;
typedef std::vector<LIST_OF_BOOL> LIST_OF_BOOL_LIST;
typedef std::vector<CHAR_LIST> LIST_OF_CHAR_LIST;

inline time_t get_time_stamp() {
    return std::chrono::steady_clock::now();
}

inline double get_duration(time_t begin, time_t end, bool show) {
    double dur = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    if (show)
        std::cout << "Time elapsed (sec) = " << dur << std::endl;
    return dur;
}

struct JAVA_CHAR {
    CHAR_LIST byte_arr;
    bool is_negative;

    JAVA_CHAR(CHAR_LIST byte_arr, bool is_negative) {
        this->byte_arr = byte_arr;
        this->is_negative = is_negative;
    }

    JAVA_CHAR() {}
};

typedef struct {
    int bits;          /* e.g., 1024 */
    NTL::ZZ n;         /* public modulus n = p q */
    NTL::ZZ n_squared; /* cached to avoid recomputing */
    NTL::ZZ n_plusone; /* cached to avoid recomputing */
} dist_paillier_pubkey_t;

typedef struct {
    int t;      // dec threshold is 2 * t + 1
    NTL::ZZ hi; // <lambda * beta>_(2t+1) evaluated at point i;
    NTL::ZZ n;
    NTL::ZZ n_squared;
    NTL::ZZ theta_invmod;
    long n_fat;
} dist_paillier_privkey_t;

typedef struct plain_text {
    CHAR_LIST plain;
    bool is_negative;

    plain_text(int len) {
        plain.reserve(len);
    }
} plain_text_t;

typedef struct cypher_text {
    CHAR_LIST c;
    bool is_neg;

    cypher_text(int len) {
        c.reserve(len);
        is_neg = false;
    }
} cypher_text_t;


class distributed_paillier {

private:
    const long kappa;
    const ZZ P;
    const long n;
    const long t;
    const long bit_len;
    long integer_share_scaling_fac;

    static bool
    check_points_in_order(const LIST_OF_LONG &input);

    static CHAR_LIST
    reshare(const LIST_OF_CHAR_LIST& in_lst,
            const vector<bool> &is_negative_lst_in,
            bool &is_negative_out,
            long &a,
            vector<long> &ret_char_size,
            const CHAR_LIST& modulu_char);

    static NTL::Pair<long, NTL::ZZ>
    reshare(const LIST_OF_CHAR_LIST& in_lst,
            const vector<bool> &is_negative_lst_in,
            long &a,
            const NTL::ZZ& modulu);

    static NTL::Pair<long, NTL::ZZ>
    __reshare__(NTL::Vec<NTL::Pair<long, NTL::ZZ>> in_lst,
                const NTL::ZZ& modulu);

    LIST_OF_CHAR_LIST
    __create_share_integer__(CHAR_LIST &secret,
                             bool &is_neg_in,
                             vector<long> &ret_char_size,
                             vector<bool> &is_negative_lst);

    void static
    __recv_and_reorganize__(const std::vector<LIST_OF_CHAR_LIST> &share_lst_b,
                            const LIST_OF_BOOL_LIST &share_lst_b_neg,
                            LIST_OF_CHAR_LIST &reorganized_b,
                            LIST_OF_BOOL &reorganized_b_neg,
                            LIST_OF_LONG &reorganized_a,
                            long num_party, long i);

public:
    distributed_paillier(ZZ P, long n, long t, long bit_len) : kappa(40), P(std::move(P)), n(n), t(t), bit_len(bit_len) {};

    LIST_OF_CHAR_LIST
    create_prime_share(vector<long> &ret_char_size,
                       vector<bool> &is_negative_lst,
                       NTL::ZZ &s);

    static ZZ get_rand_4_biprimetest(const ZZ& N);

    static ZZ biprime_test_step1(int party_id, const ZZ& N, const ZZ& pi, const ZZ& qi,  const ZZ& g);

    // this can be done on only one party
    static bool biprime_test_step2(const vector<ZZ>& other_v, const ZZ& v, const ZZ& N, const long &num_party );


        LIST_OF_CHAR_LIST
    create_lambda_share(CHAR_LIST &secret,
                        bool &is_neg_in,
                        vector<long> &ret_char_size,
                        vector<bool> &is_negative_lst);

    LIST_OF_CHAR_LIST
    create_beta_share(CHAR_LIST &secret,
                      bool &is_neg_in,
                      vector<long> &ret_char_size,
                      vector<bool> &is_negative_lst);

    LIST_OF_CHAR_LIST
    compute_lambda_times_beta_share(const LIST_OF_CHAR_LIST& in_lst_lambda,
                                    const LIST_OF_CHAR_LIST& in_lst_beta,
                                    const LIST_OF_BOOL &is_negative_lst_in_lam,
                                    const LIST_OF_BOOL &is_negative_lst_in_beta,
                                    const CHAR_LIST &modulu_char,
                                    long &a,
                                    vector<long> &ret_char_size,
                                    bool &is_negative_out,
                                    bool &theta_is_negative_out,
                                    long const &scalar);

    inline static NTL::ZZ
    compute_theta_share(const NTL::ZZ &l_times_b,
                                              const NTL::ZZ &modulu) {
        return (l_times_b % modulu);
    }

    static std::vector<NTL::ZZ> gene_local_piqi_4_first_party(int bit_len, int n);
    static std::vector<NTL::ZZ> gene_local_piqi_4_other_party(int bit_len, int n);

    /**
     * \brief use all shares of  theta to reconstruct theta
     * @param theta_shares ALL shares of theta
     * @param N modulo
     * @return reconstructed theta
     */


    void static enc(cypher_text_t &cypher_text,
                    const plain_text_t &plain,
                    const dist_paillier_pubkey_t &pb_key);

    void static enc(cypher_text_t &cypher_text,
                    const ZZ &plain,
                    const dist_paillier_pubkey_t &pb_key);

    void static dec_partial(const long i,
                            const long t,
                            cypher_text_t &plain_text,
                            const cypher_text_t &cypher_text,
                            const dist_paillier_privkey_t &priv_key);

    void static dec_final(plain_text_t &plain,
                          const cypher_text_t &cypher_text,
                          const std::vector<cypher_text_t> &lihi_lst,
                          const long &t,
                          const NTL::ZZ &max_negative_abs,
                          const dist_paillier_privkey_t &priv_key);

    /** 
     * Final decryption
     *  client 1, 2, 3, ..., 2t+1 send partial decrypted data to the
     *  one who needs decryption
     * 
     * Parameters:
     *  @param lihi_lst
     *  @param t
     *  @param max_negative_abs The greatest value of abs(to_be_decrypted) when to_be_decrypted is negative.
     *   Some Notes about this parameter: 
     *    The primitive Paillier only supports non-negative value output because it takes mod.
     *    Its final decrypted result is always between 0 and n. When decrypting a negative value
     *    (its plaintext is negative), mod n is taken at the last step. A negative value, say x
     *    will become n+x in the final result.
     *    Since that n is often a very large value, by carefully choosing a "max_negative_abs", 
     *    we can make n+x and normal positive results very easy to distinguish. 
     */
    NTL::ZZ static dec_final(const cypher_text_t &cypher_text,
                             const std::vector<cypher_text_t> &lihi_lst,
                             const long &t,
                             const NTL::ZZ &max_negative_abs,
                             const dist_paillier_privkey_t &priv_key);

    void static add(cypher_text_t &res,
                    const cypher_text_t &a,
                    const cypher_text_t &b,
                    const dist_paillier_pubkey_t &pk);

    void static mul(cypher_text_t &res,
                    const cypher_text_t &a,
                    const plain_text_t &b,
                    const dist_paillier_pubkey_t &pk);

    void static mul(cypher_text_t &res,
                    const cypher_text_t &a,
                    const NTL::ZZ &b,
                    const dist_paillier_pubkey_t &pk);

    /**
     * div(a, b) is basically mul(a, InvMod(b, modulo)), nothing special.
     * b should never be 0.
     */
    void static div(cypher_text_t &res,
                    const cypher_text_t &a,
                    const NTL::ZZ &b,
                    const dist_paillier_pubkey_t &pk);


    static std::vector<dist_paillier_privkey_t>
    gene_privkey_standalone(const long &len,
                            const long &n,
                            const long &t);


    static long
    ZZ_2_byte(CHAR_LIST &out,
              const NTL::ZZ &in,
              bool &is_negative);

    static void
    char_list_2_ZZ(NTL::ZZ &out,
                   const CHAR_LIST &in,
                   const bool &is_negative);

    void static
    __recv_and_reorganize__(const std::vector<ShamirShares> &share_lst_b,
                            std::vector<ZZ> &reorganized_b,
                            LIST_OF_LONG &reorganized_a,
                            long num_party, long i);


    Vec<Pair<long, NTL::ZZ>> __create_share_integer__(ZZ &s);

    /**
     * \brief adding shares of lambda and beta, then compute one share of sum(lambda)*sum(beta)
     * @param in_lst_lambda shares of lambda_i for one party
     * @param in_lst_beta shares of beta_i for one party
     * @param modulu
     * @param scalar scalar generated when multiplying 2 shares
     * @return
     */
    static vector<ZZ> compute_lambda_times_beta_share(const Vec<NTL::Pair<long, NTL::ZZ>> &in_lst_lambda,
                                               const Vec<NTL::Pair<long, NTL::ZZ>> &in_lst_beta, const ZZ &modulu,
                                               const long &scalar);

    static ZZ reveal_theta(const vector<NTL::ZZ> &theta_shares, const ZZ &N, const int &threshold, const int &num_party);
};

#endif //DIST_PAI_H
