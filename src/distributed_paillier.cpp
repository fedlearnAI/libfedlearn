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

#include <distributed_paillier.h>

/***************************************
 * Distributed Paillier key generation
***************************************/
LIST_OF_CHAR_LIST
distributed_paillier::create_prime_share(vector<long> &ret_char_size,
                                         vector<bool> &is_negative_lst,
                                         NTL::ZZ &s) {
    NTL::GenPrime(s, bit_len);

    ShamirSecretSharingScheme scheme = ShamirSecretSharingScheme(P, n, t);
    ShamirShares share = scheme.share_secret(s);

    // return shared points (a,b) to JAVA interface
    // NB: a is always 1, 2, 3, ..., t+1
    LIST_OF_CHAR_LIST ret;
    ret.reserve(share.shares.length());
    int i = 0;
    for (const auto &each_point : share.shares) {
        CHAR_LIST one_point;
        bool is_negative;
        ret_char_size[i] = ZZ_2_byte(one_point, each_point.b, is_negative);
        is_negative_lst[i] = is_negative;
        i++;
        ret.push_back(one_point);
    }

    // check if share is constructed correctly
    assert(s == share.reconstruct_secret(P));
    std::cout << "s = " << s << " P = " << P << endl;

    return ret;
}

LIST_OF_CHAR_LIST
distributed_paillier::create_lambda_share(CHAR_LIST &secret,
                                          bool &is_neg_in,
                                          vector<long> &ret_char_size,
                                          vector<bool> &is_negative_lst) {
    return __create_share_integer__(secret, is_neg_in, ret_char_size, is_negative_lst);
}

LIST_OF_CHAR_LIST
distributed_paillier::create_beta_share(CHAR_LIST &secret,
                                        bool &is_neg_in,
                                        vector<long> &ret_char_size,
                                        vector<bool> &is_negative_lst) {
    return __create_share_integer__(secret, is_neg_in, ret_char_size, is_negative_lst);
}

LIST_OF_CHAR_LIST
distributed_paillier::__create_share_integer__(CHAR_LIST &secret,
                                               bool &is_neg_in,
                                               vector<long> &ret_char_size,
                                               vector<bool> &is_negative_lst) {
    NTL::ZZ s;
    distributed_paillier::char_list_2_ZZ(s, secret, is_neg_in);
    if (s > P) {
        std::cout << "secret is larger then P." << endl;
        assert(s <= P);
    }

    ShamirSecretSharingSchemeInteger scheme = ShamirSecretSharingSchemeInteger(P, n, t, kappa);
    ShamirShares_integer share = scheme.share_secret(NTL::to_ZZ(s));
    integer_share_scaling_fac = NTL::to_long(share.scaling);

    LIST_OF_CHAR_LIST ret;
    ret.reserve(share.shares.length());
    int i = 0;
    for (const auto &each_point : share.shares) {
        CHAR_LIST one_point;
        bool is_negative;
        ret_char_size[i] = ZZ_2_byte(one_point, each_point.b, is_negative);
        is_negative_lst[i] = is_negative;
        ret.push_back(one_point);
        i++;
    }
    return ret;
}

/**
 * Re-share additive Shamir shares.
 * 
 * Parameters:
 *  a: the first value of a point, i.e. a in point(a, b)
 *  in_lst: b_j from other shares, i.e.  b_1, b_2, ..., b_t+1 in
 *          point(a, b_1), (a, b_2), ..., (a, b_t+1)
 */
NTL::Pair<long, NTL::ZZ>
distributed_paillier::__reshare__(NTL::Vec<NTL::Pair<long, NTL::ZZ>> in_lst,
                                  const NTL::ZZ& modulu) {
    NTL::ZZ b(0);
    long a = in_lst[0].a;
    for (const auto &pair : in_lst) {
        NTL::add(b, b, pair.b);
        if (a != pair.a) {
            std::cout << "a in points(a, b) are not the same. a = " + to_string(a) + ", pair.a = " + to_string(pair.a)
                      << endl;
            assert(a == pair.a);
        }
    }
    return NTL::Pair<long, NTL::ZZ>(a, b);
}

CHAR_LIST
distributed_paillier::reshare(const LIST_OF_CHAR_LIST& in_lst,
                              const vector<bool> &is_negative_lst_in,
                              bool &is_negative_out,
                              long &a,
                              vector<long> &ret_char_size,
                              const CHAR_LIST& modulu_char) {
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> share_part;
    int i = 0;
    for (const auto &char_list : in_lst) {
        NTL::ZZ b;
        char_list_2_ZZ(b, char_list, is_negative_lst_in[i]);
        share_part.append(NTL::Pair<long, NTL::ZZ>(a, b));
        i++;
    }
    NTL::ZZ modulu;
    char_list_2_ZZ(modulu, modulu_char, false);
    NTL::Pair<long, NTL::ZZ> res = __reshare__(share_part, modulu);
    CHAR_LIST one_point;
    ZZ_2_byte(one_point, res.b, is_negative_out);
    return one_point;
}

NTL::Pair<long, NTL::ZZ>
distributed_paillier::reshare(const LIST_OF_CHAR_LIST& in_lst,
                              const vector<bool> &is_negative_lst_in,
                              long &a,
                              const NTL::ZZ& modulu) {
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> share_part;
    int i = 0;
    for (const auto &char_list : in_lst) {
        NTL::ZZ b;
        char_list_2_ZZ(b, char_list, is_negative_lst_in[i]);
        share_part.append(NTL::Pair<long, NTL::ZZ>(a, b));
        i++;
    }
    NTL::Pair<long, NTL::ZZ> res = __reshare__(share_part, modulu);
    return res;
}

/** Compute share part of lambda*beta at point a
 * Return:
 *  res[0] -- lambda_times_beta
 *  res[1] -- theta
 * NB: 1. this would make the threshold increase to 2*t, which means 
 *        that 2*t+1 points are needed for decryption. Original threshold is t.
 *     2. we note that lambda_times_beta and theta must have the same sign
 *        this requires the scale factor to be positive.
 */
LIST_OF_CHAR_LIST
distributed_paillier::compute_lambda_times_beta_share(LIST_OF_CHAR_LIST in_lst_lambda,
                                                      LIST_OF_CHAR_LIST in_lst_beta,
                                                      const LIST_OF_BOOL &is_negative_lst_in_lam,
                                                      const LIST_OF_BOOL &is_negative_lst_in_beta,
                                                      const CHAR_LIST &modulu_char,
                                                      long &a,
                                                      vector<long> &ret_char_size,
                                                      bool &is_negative_out,
                                                      bool &theta_is_negative_out,
                                                      long const &scalar) {
    NTL::ZZ modulu;
    char_list_2_ZZ(modulu, modulu_char, false);
    NTL::Pair<long, NTL::ZZ> lamd = reshare(in_lst_lambda, is_negative_lst_in_lam, a, modulu);
    NTL::Pair<long, NTL::ZZ> beta = reshare(in_lst_beta, is_negative_lst_in_beta, a, modulu);
    NTL::ZZ lambda_times_beta, theta_share;
    NTL::mul(lambda_times_beta, lamd.b, beta.b);

    theta_share = compute_theta_share(lambda_times_beta, modulu);

    LIST_OF_CHAR_LIST res(2);
    ZZ_2_byte(res[0], lambda_times_beta, is_negative_out);
    ZZ_2_byte(res[1], theta_share, theta_is_negative_out);
    // std::cout << lambda_times_beta << endl;
    return res;
}

NTL::ZZ
distributed_paillier::compute_theta_share(const NTL::ZZ &l_times_b,
                                          const NTL::ZZ &modulu) {
    return (l_times_b % modulu);
}

NTL::ZZ
distributed_paillier::reveal_theta(const vector<NTL::ZZ> &theta_shares,
                                   const NTL::ZZ &N) {
    assert(theta_shares.size() >= 2 * t + 1);
    NTL::ZZ nfct = factorial_small(n);
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
        NTL::AddMod(recon_theta, recon_theta, li * theta_shares[i - 1], N);
    }
    NTL::ZZ scalar_inv;
    NTL::InvMod(scalar_inv, nfct * sqr(nfct), N);
    NTL::MulMod(recon_theta, recon_theta, scalar_inv, N);
    NTL::MulMod(recon_theta, recon_theta, nfct * sqr(nfct), N);
    return recon_theta;
}

/****************************
 * Paillier arithmetic
****************************/

void distributed_paillier::add(cypher_text_t &res,
                               const cypher_text_t &a,
                               const cypher_text_t &b,
                               const dist_paillier_pubkey_t &pk) {
    NTL::ZZ azz, bzz, reszz;
    char_list_2_ZZ(azz, a.c, a.is_neg);
    char_list_2_ZZ(bzz, b.c, b.is_neg);
    NTL::MulMod(reszz, azz, bzz, pk.n_squared);
    ZZ_2_byte(res.c, reszz, res.is_neg);

    // cout << reszz <<" "
    //      << azz <<" "
    //      << bzz <<" "
    //      << endl;
}

void distributed_paillier::mul(cypher_text_t &res,
                               const cypher_text_t &a,
                               const plain_text_t &b,
                               const dist_paillier_pubkey_t &pk) {
    NTL::ZZ azz, bzz, reszz;
    char_list_2_ZZ(azz, a.c, a.is_neg);
    char_list_2_ZZ(bzz, b.plain, b.is_negative);

    // here azz should never be larger than pk.n_squared
    assert(azz < pk.n_squared);
    NTL::PowerMod(reszz, azz, bzz, pk.n_squared);
    ZZ_2_byte(res.c, reszz, res.is_neg);
}

void distributed_paillier::mul(cypher_text_t &res,
                               const cypher_text_t &a,
                               const NTL::ZZ &b,
                               const dist_paillier_pubkey_t &pk) {
    NTL::ZZ azz, reszz;
    char_list_2_ZZ(azz, a.c, a.is_neg);
    assert(azz < pk.n_squared);
    NTL::PowerMod(reszz, azz, b, pk.n_squared);
    ZZ_2_byte(res.c, reszz, res.is_neg);
}

void distributed_paillier::div(cypher_text_t &res,
                               const cypher_text_t &a,
                               const NTL::ZZ &b,
                               const dist_paillier_pubkey_t &pk) {
    NTL::ZZ azz, reszz, b_inv;
    char_list_2_ZZ(azz, a.c, a.is_neg);
    assert(b != 0);
    NTL::InvMod(b_inv, b, pk.n);
    assert(azz < pk.n_squared);
    NTL::PowerMod(reszz, azz, b_inv, pk.n_squared);
    ZZ_2_byte(res.c, reszz, res.is_neg);

    // cout << "reszz = " << reszz << endl;
}

/****************************
 * Paillier Enc & Two-step Dec
****************************/

void distributed_paillier::enc(cypher_text_t &cypher_text,
                               const ZZ &plain,
                               const dist_paillier_pubkey_t &pb_key) {
    NTL::ZZ c_zz, r;
    do {
        NTL::RandomBnd(r, pb_key.n);
        //  r = pb_key.n-1; // only for debug
    } while (r == 0 || r == pb_key.n || (NTL::GCD(r, pb_key.n) != 1));

    NTL::PowerMod(c_zz, pb_key.n_plusone, plain, pb_key.n_squared);
    NTL::PowerMod(r, r, pb_key.n, pb_key.n_squared);
    NTL::MulMod(c_zz, c_zz, r, pb_key.n_squared);

    ZZ_2_byte(cypher_text.c, c_zz, cypher_text.is_neg);
}

void distributed_paillier::enc(cypher_text_t &cypher_text,
                               const plain_text_t &plain,
                               const dist_paillier_pubkey_t &pb_key) {
    NTL::ZZ p_zz, c_zz, r;
    do {
        NTL::RandomBnd(r, pb_key.n);
        // r = 1; // only for debug
    } while (r == 0 || r == pb_key.n || (NTL::GCD(r, pb_key.n) != 1));
    char_list_2_ZZ(p_zz, plain.plain, plain.is_negative);

    NTL::PowerMod(c_zz, pb_key.n_plusone, p_zz, pb_key.n_squared);
    NTL::PowerMod(r, r, pb_key.n, pb_key.n_squared);
    NTL::MulMod(c_zz, c_zz, r, pb_key.n_squared);

    ZZ_2_byte(cypher_text.c, c_zz, cypher_text.is_neg);
}

void distributed_paillier::dec_partial(const long i,
                                       const long t,
                                       cypher_text_t &partial_res,
                                       const cypher_text_t &cypher_text,
                                       const dist_paillier_privkey_t &priv_key) {
    NTL::ZZ c, e;
    char_list_2_ZZ(c, cypher_text.c, cypher_text.is_neg);

    // compute li*hi
    NTL::ZZ li;
    long enume = priv_key.n_fat;
    long denom = 1;
    for (long j = 1; j <= 2 * t + 1; j++) {
        if (i != j) {
            enume *= j;
            denom *= j - i;
        }
    }
    li = enume / denom;
    NTL::mul(c, li, priv_key.hi);
    ZZ_2_byte(partial_res.c, c, partial_res.is_neg);

    // std::cout << "li = " << li
    //           << endl;
}

void distributed_paillier::dec_final(plain_text_t &plain,
                                     const cypher_text_t &cypher_text,
                                     const std::vector<cypher_text_t> &lihi_lst,
                                     const long &t,
                                     const NTL::ZZ &max_negative_abs,
                                     const dist_paillier_privkey_t &priv_key) {
    NTL::ZZ res = dec_final(cypher_text, lihi_lst, t, max_negative_abs, priv_key);
    ZZ_2_byte(plain.plain, res, plain.is_negative);

    // std::cout << "res = " << res
    //           << " lihi_sum = " << lihi_sum
    //           << endl;
}

NTL::ZZ distributed_paillier::dec_final(const cypher_text_t &cypher_text,
                                        const std::vector<cypher_text_t> &lihi_lst,
                                        const long &t,
                                        const NTL::ZZ &max_negative_abs,
                                        const dist_paillier_privkey_t &priv_key) {
    NTL::ZZ c, lihi_sum(0), res;
    assert(lihi_lst.size() >= 2 * t + 1);
    for (long i = 1; i <= 2 * t + 1; i++) {
        NTL::ZZ tmp;
        char_list_2_ZZ(tmp, lihi_lst[i - 1].c, lihi_lst[i - 1].is_neg);
        NTL::add(lihi_sum, lihi_sum, tmp);
    }
    char_list_2_ZZ(c, cypher_text.c, cypher_text.is_neg);
    NTL::PowerMod(res, c, lihi_sum, priv_key.n_squared);
    NTL::add(res, res, -1);

    if ((res % priv_key.n) != 0) {
        std::cout << "res = " << res
                  << " lihi_sum = " << lihi_sum
                  << " c = " << c
                  << " priv_key.n_squared = " << priv_key.n_squared
                  << " priv_key.n = " << priv_key.n
                  << endl;
        assert(res % priv_key.n == 0);
    }

    NTL::div(res, res, priv_key.n);
    NTL::MulMod(res, res, priv_key.theta_invmod, priv_key.n);

    if (res > max_negative_abs) {
        res = res - priv_key.n;
    }
    // std::cout << "res = " << res
    //           << " max_negative_abs = " << max_negative_abs
    //           << endl;
    return res;
}

/****************************
 * class utils 
****************************/

long distributed_paillier::ZZ_2_byte(CHAR_LIST &out,
                                     const NTL::ZZ &in,
                                     bool &is_negative) {
    // return shared points (a,b) to JAVA interface
    // NB : a is always 1, 2, 3, ..., t+1
    long char_size = NTL::NumBytes(in);
    auto *pp = new unsigned char[char_size];
    NTL::BytesFromZZ(pp, in, char_size); // pp = byte-representation
    out.assign(pp, pp + char_size);
    delete[] pp;
    is_negative = (in < 0) != 0;
    return char_size;
}

void distributed_paillier::char_list_2_ZZ(NTL::ZZ &out,
                                          const CHAR_LIST &in,
                                          const bool &is_negative) {
    long char_size = in.size();
    const unsigned char *in_pnt = in.data();
    NTL::ZZFromBytes(out, in_pnt, char_size);
    out = is_negative ? -1 * out : out;
}

bool distributed_paillier::check_points_in_order(const LIST_OF_LONG &input) {
    long i = 1;
    for (auto x : input) {
        if (x != i) {
            std::cout << "Input points are not in order of 1, 2, ..., t+1" << endl;
            return false;
        }
        i += 1;
    }
    return true;
}

/**
 * A standalone version priv-key generation. Should switch to
 * distributed version
 * ==========
 * Parameters
 *  ...
 * Return
 *  std::vector containing 2t+1 dist_paillier_privkey_t's.
 *
 * FIXME: NOTE that only 2t+1 keys are generated. This code does
 *  not work when n!=2t+1. Fixing is easy -- generating n shares
 *  of lambda*beta instead of 2t+1.
 */
std::vector<dist_paillier_privkey_t>
distributed_paillier::gene_privkey_standalone(const long &len,
                                              const long &n,
                                              const long &t) {
    std::cout << "Generating priv-keys for all parties..." << endl;

    // 2^N must be larger than sqr of the largest number in plain text
    // assert( (long)2<<len > (plain_text_num*plain_text_num));

    // NOTE that N is generated by one party here, thus bi-prime test is not done.
    NTL::ZZ nfct = factorial_small(n);
    NTL::ZZ p, q, P;
    NTL::GenPrime(p, len);
    NTL::GenPrime(q, len);
    NTL::GenPrime(P, len);
    const NTL::ZZ N = p * q;
    std::cout << "\np = " << p
              << "; q=" << q
              << "; N=" << N
              << "; P=" << P
              << ";"
              << endl;

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
        paillier_class_lst.emplace_back(N, n, t, len);
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
        cout << "pi_lst[" << i << "]  = " << pi_lst[i]
             << "; qi_lst[" << i << "]  = " << qi_lst[i] << ";"
             << endl;
        assert((pi_lst[i] > 0) && (qi_lst[i] > 0));
    }
    assert(q_sum == q && p_sum == p);

    for (long i = 1; i <= n; i++) {
        NTL::ZZ tmp;
        do {
            tmp = NTL::RandomBnd(N);
//                tmp = 1; // for debug
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
        __recv_and_reorganize__(lambda_share, lambda_share_sign,
                                reorganized_lambda_b, reorganized_lambda_b_neg,
                                reorganized_a, n, i);

        // beta recv
        LIST_OF_CHAR_LIST reorganized_beta_b;
        LIST_OF_BOOL reorganized_beta_b_neg;
        LIST_OF_LONG reorganized_abeta;
        __recv_and_reorganize__(beta_share, beta_shares_sign,
                                reorganized_beta_b, reorganized_beta_b_neg,
                                reorganized_abeta, n, i);

        long scalar = NTL::to_long(nfct);

        tmp_res =
                paillier_class_lst[i - 1].compute_lambda_times_beta_share(
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
    // 2) check shares of lamb_times_beta can reconstruct rea l lamb_times_beta
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
    NTL::ZZ scalar_inv;
    NTL::InvMod(scalar_inv, nfct * sqr(nfct), N);
    NTL::MulMod(recon_theta, recon_theta, scalar_inv, N);
    NTL::MulMod(recon_theta, recon_theta, nfct * sqr(nfct), N);

    NTL::ZZ lambda_lst_sum(0), beta_lst_sum(0);
    for (int i = 0; i < lambda_lst.length(); i++) {
        lambda_lst_sum += lambda_lst[i];
        beta_lst_sum += beta_lst[i];
    }

    auto real_lamb_times_beta = lambda_lst_sum * beta_lst_sum;
    // std::cout << "\nrecon_lamb_times_beta  = " << recon_lamb_times_beta
    //           << "\n real_lamb_times_beta  = " << real_lamb_times_beta * nfct * sqr(nfct)
    //           << "\nrecon_theta = " << recon_theta
    //           << "\nN = " << N
    //           << endl;

    // NB: These two assertions are VERY IMPORTANT!
    assert(lambda_lst_sum * beta_lst_sum * nfct * sqr(nfct) == recon_lamb_times_beta);
    assert(recon_theta == recon_lamb_times_beta % N);
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
    return pri_key_lst;
}

void
distributed_paillier::__recv_and_reorganize__(const std::vector<LIST_OF_CHAR_LIST> &share_lst_b,
                                              const LIST_OF_BOOL_LIST &share_lst_b_neg,
                                              LIST_OF_CHAR_LIST &reorganized_b,
                                              LIST_OF_BOOL &reorganized_b_neg,
                                              LIST_OF_LONG &reorganized_a,
                                              long num_party, long i) {
    reorganized_a.reserve(num_party);
    reorganized_b.reserve(num_party);
    reorganized_b_neg.reserve(num_party);

    for (long j = 1; j <= num_party; j++) {
        reorganized_b.push_back(share_lst_b[j - 1][i - 1]);
        reorganized_a.push_back(j);
        reorganized_b_neg.push_back(share_lst_b_neg[j - 1][i - 1]);
    }
}

// TODO: reconstruct __recv_and_reorganize__ functions
void
distributed_paillier::__recv_and_reorganize__(const std::vector<ShamirShares> &share_lst_b,
                                              std::vector<ZZ> &reorganized_b,
                                              LIST_OF_LONG &reorganized_a,
                                              long num_party, long i) {
    reorganized_a.reserve(num_party);
    reorganized_b.reserve(num_party);

    for (long j = 1; j <= num_party; j++) {
        reorganized_b.push_back(share_lst_b[j - 1].shares.get(i - 1).b);
        assert(share_lst_b[j - 1].shares.get(i - 1).a == i);
        reorganized_a.push_back(i);
    }
}

std::vector<NTL::ZZ>
distributed_paillier::gene_local_piqi() {
    std::vector<ZZ> ret;
    ret.reserve(2);
    long p_q_bit_len = bit_len;
    long n_copy = n;
    while(n_copy>1) {
        p_q_bit_len --;
        n_copy = n_copy >> 1;
    }
    p_q_bit_len--;
    NTL::ZZ tmp;
    do {
        tmp = NTL::RandomBits_ZZ(p_q_bit_len);
//                tmp = 1; // for debug
    } while (tmp == 0);
    ret.emplace_back(tmp);
    do {
        tmp = NTL::RandomBits_ZZ(p_q_bit_len);
//                tmp = 1; // for debug
    } while (tmp == 0);
    ret.emplace_back(tmp);
    return ret;
}

ZZ distributed_paillier::get_rand_4_biprimetest(const ZZ& N) {
    ZZ tmp = NTL::RandomBnd(N);
    while(NTL::Jacobi(tmp, N) != 1) {
        tmp = NTL::RandomBnd(N);
    }
    return tmp;
}

ZZ distributed_paillier::biprime_test_step1(int party_id, const ZZ& N, const ZZ& pi, const ZZ& qi,  const ZZ& g) {
    if(party_id==1) {
        ZZ tmp =  (N-pi-qi+1);
        assert(tmp % 4 == 0);
        return NTL::PowerMod(g, tmp, N);
    } else {
        ZZ tmp =  pi+qi;
        assert(tmp % 4 == 0);
        return NTL::PowerMod(g, tmp, N);
    }
}

bool distributed_paillier::biprime_test_step2(const vector<ZZ>& other_v, const ZZ& v, const long &num_party ) {
    ZZ tmp = ZZ(0);
    assert(other_v.size() == num_party-1);
    for(const auto& value : other_v) {
        tmp *= value;
    }
    return ( (v==tmp) || (v==-1*tmp) );
}

