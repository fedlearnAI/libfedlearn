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
#include "distributed_paillier.h"

using namespace std;

namespace c_api {

    /**
     *  Enc one data
     */
    inline void enc(cypher_text_t &res, const plain_text_t &plain, const dist_paillier_pubkey_t &pubkey) {
        distributed_paillier::enc(res, plain, pubkey);
    }

    void enc(cypher_text_t &res, const NTL::ZZ &plain, const dist_paillier_pubkey_t &pubkey) {
        distributed_paillier::enc(res, plain, pubkey);
    }

    /**
     * partial dec for one data
     */
    inline void partial_dec(const long i,
                            const long t,
                            cypher_text_t &plain_text,
                            const cypher_text_t &cypher_text,
                            const dist_paillier_privkey_t &priv_key) {
        distributed_paillier::dec_partial(i, t, plain_text, cypher_text, priv_key);
    }

    /**
     * final dec for one data
     */
    inline void final_dec(plain_text_t &plain,
                          const cypher_text_t &cypher_text,
                          const std::vector<cypher_text_t> &lihi_lst,
                          const long &t,
                          const NTL::ZZ &max_negative_abs,
                          const dist_paillier_privkey_t &priv_key) {
        distributed_paillier::dec_final(plain, cypher_text, lihi_lst, t, max_negative_abs, priv_key);
    }

    inline NTL::ZZ final_dec(const cypher_text_t &cypher_text,
                             const std::vector<cypher_text_t> &lihi_lst,
                             const long &t,
                             const NTL::ZZ &max_negative_abs,
                             const dist_paillier_privkey_t &priv_key) {
        return distributed_paillier::dec_final(cypher_text, lihi_lst, t, max_negative_abs, priv_key);
    }

    /**
     * add one data
     */
    inline void
    add(cypher_text_t &res, const cypher_text_t &a, const cypher_text_t &b, const dist_paillier_pubkey_t &pk) {
        distributed_paillier::add(res, a, b, pk);
    }

    /**
     * mul one data
     */
    inline void mul(cypher_text_t &res, const cypher_text_t &a, const ZZ &b, const dist_paillier_pubkey_t &pk) {
        distributed_paillier::mul(res, a, b, pk);
    }

    /**
     * devide one data
     */
    inline void div(cypher_text_t &res, const cypher_text_t &a, const ZZ &b, const dist_paillier_pubkey_t &pk) {
        distributed_paillier::div(res, a, b, pk);
    }

} // namespace c_api
