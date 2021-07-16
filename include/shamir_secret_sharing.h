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
#include <stdlib.h>
#include <stdexcept>
#include <assert.h>
#include "../include/NTL/ZZ.h"
#include "../include/NTL/vector.h"
#include "../include/NTL/pair.h"
#include "../include/NTL/tools.h"

using namespace std;
using namespace NTL;

class ShamirShares {
public:
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> shares;
    const NTL::ZZ P;    // the large prime to define Zp
    const long n; // number of party
    const long degree;

    ShamirShares(NTL::Vec<NTL::Pair<long, NTL::ZZ>> &shares,
                 const NTL::ZZ &P,
                 const long &n,
                 const long &degree) : shares(shares), P(P), n(n), degree(degree) {};

    NTL::ZZ
    reconstruct_secret(const NTL::ZZ &modulu);

    NTL::Vec<NTL::Pair<long, NTL::ZZ>>
    get_partial_recon(const NTL::ZZ &N);

    ShamirShares
    __add__(ShamirShares other);

    ShamirShares
    __mult__(ShamirShares other);

private:
    bool __check_if_compatible__(const ShamirShares &other);

    void __check_2_op_correct__(ShamirShares other);

    static bool __check_points_in_order__(NTL::Vec<NTL::Pair<long, NTL::ZZ>> &input);
};

class ShamirSecretSharingScheme {
public:
    ShamirSecretSharingScheme(const NTL::ZZ &P, const long &n, const long &t) : P(P), n(n), t(t) {}

    const NTL::ZZ P;    // the large prime to define Zp
    const long n; // number of party
    const long t; // decryption threshold, i.e. need t+1 party to decrypt cyphertext, degree of corresponding polynomial is t

    ShamirShares share_secret(const NTL::ZZ &s); // do secret share and return shares for each party.

    inline bool operator==(ShamirSecretSharingScheme other) {
        return (P != other.P || n != other.n || t != other.t);
    }
};