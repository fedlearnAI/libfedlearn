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
#include <assert.h>
#include <stdlib.h>
#include <stdexcept>
#include "../include/NTL/ZZ.h"
#include "../include/NTL/vector.h"
#include "../include/NTL/pair.h"
#include "../include/NTL/tools.h"

using namespace std;
using namespace NTL;

inline NTL::ZZ factorial_small(long x) {
    assert(x > 0);
    if (x == (long) 1) {
        return NTL::to_ZZ(x);
    } else {
        return x * factorial_small(x - 1);
    }
}

class ShamirShares_integer {
public:
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> shares;
    NTL::ZZ n_fac;
    const NTL::ZZ P;   // the large prime to define Zp
    const long n;      // number of party
    const long degree; // degree of the poly
    const NTL::ZZ scaling;

    ShamirShares_integer(NTL::Vec<NTL::ZZ> &shares_vec,
                         const NTL::ZZ &P,
                         const long &n,
                         const long &degree,
                         const NTL::ZZ &scaling) : P(P), n(n), degree(degree), scaling(scaling) {
        n_fac = factorial_small(n);
        long cnt = 1;
        for (auto x: shares_vec) {
            shares.append(NTL::Pair<long, NTL::ZZ>(cnt, x));
            cnt += 1;
        }
    };

    ShamirShares_integer(NTL::Vec<NTL::Pair<long, NTL::ZZ>> &shares,
                         const NTL::ZZ &P,
                         const long &n,
                         const long &degree,
                         const NTL::ZZ &scaling) : shares(shares), P(P), n(n), degree(degree), scaling(scaling) {
        n_fac = factorial_small(n);
    };

    NTL::ZZ
    reconstruct_secret(const NTL::ZZ& N); // the Paillier modulu;

    NTL::Vec<NTL::Pair<long, NTL::ZZ>>
    get_partial_recon(const NTL::ZZ& N);

    ShamirShares_integer
    __add__(ShamirShares_integer other);

    ShamirShares_integer
    __mult__(ShamirShares_integer other);

    ShamirShares_integer
    __mult__(long other);

private:
    bool
    __check_if_compatible__(ShamirShares_integer other);

    void
    __check_2_op_doable__(ShamirShares_integer other);

    bool
    __check_points_in_order__(NTL::Vec<NTL::Pair<long, NTL::ZZ>> &input);
};

class ShamirSecretSharingSchemeInteger {
public:
    ShamirSecretSharingSchemeInteger(const NTL::ZZ &P,
                                     const long &n,
                                     const long &t,
                                     const long &kappa) : P(P), n(n), t(t), kappa(kappa) {
        randomness_interval = NTL::sqr(factorial_small(n)) * NTL::power2_ZZ(kappa) * P;
    }

    ShamirSecretSharingSchemeInteger(const NTL::ZZ &P,
                                     const long &n,
                                     const long &t) : P(P), n(n), t(t), kappa((long) 40) {
        randomness_interval = NTL::sqr(factorial_small(n)) * NTL::power2_ZZ(kappa) * P;
    }

    const NTL::ZZ P;  // the large prime to define Zp
    const long n;     // number of party
    const long t;     // decryption threshold, i.e. need t+1 party to decrypt ciphertext, degree of corresponding polynomial is t
    const long kappa; // statistical security coeff
    NTL::ZZ randomness_interval;

    ShamirShares_integer
    share_secret(const NTL::ZZ &s); // do secret share and return shares for each party.

    inline bool operator==(ShamirSecretSharingSchemeInteger other) {
        return (P != other.P || n != other.n || t != other.t);
    }
};