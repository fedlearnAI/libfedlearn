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
#include "../include/shamir_secret_sharing.h"

ShamirShares ShamirSecretSharingScheme::share_secret(const NTL::ZZ &s) {
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> shares;
    NTL::ZZ f0(s);
    NTL::Vec<NTL::ZZ> poly_param;

    for (int i = 1; i < t + 1; i++) {
        poly_param.append(NTL::RandomBnd(P));
    }

    // get n points from the polynomial
    // n points from 1 to n
    // polynomial has degree of t
    for (int i = 1; i <= n; i++) {
        NTL::ZZ secret(f0);
        for (int j = 1; j < t + 1; j++) {
            NTL::ZZ tmp;
            NTL::MulMod(tmp, poly_param[j - 1], (long) pow(i, j), P);
            NTL::AddMod(secret, secret, tmp, P);
            // secret += tmp;
        }
        shares.append(NTL::Pair<long, NTL::ZZ>((long) (i), secret % P));
    }

    ShamirShares Sharing = ShamirShares(shares, P, n, t);
    return Sharing;
}

NTL::ZZ ShamirShares::reconstruct_secret(const NTL::ZZ &modulu) {
    if (shares.length() < degree + 1) {
        std::cout << "no enough points for recon. Need = " + to_string(degree + 1) + ", get = " +
                     to_string(shares.length()) << endl;
        assert(shares.length() >= degree + 1);
    }

    NTL::ZZ partial_recon_sum = NTL::ZZ(0);
    NTL::ZZ res;
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> partial_recon = get_partial_recon(modulu);
    for (long i = 1; i <= partial_recon.length(); i++) {
        NTL::AddMod(partial_recon_sum, partial_recon_sum, partial_recon[i - 1].b, P); // sum the partial recon
    }
    return partial_recon_sum;
}

NTL::Vec<NTL::Pair<long, NTL::ZZ>> ShamirShares::get_partial_recon(const NTL::ZZ &N) {
    if (shares.length() < degree + 1) {
        std::cout << "no enough points for recon. Need = " + to_string(degree + 1) + ", get = " +
                     to_string(shares.length()) << endl;
        assert(shares.length() >= degree + 1);
    }
    assert(__check_points_in_order__(shares));

    NTL::Vec<NTL::Pair<long, NTL::ZZ>> partial_recon;
    for (long i = 1; i <= degree + 1; i++) {
        NTL::ZZ enume = NTL::ZZ(1);
        NTL::ZZ denom = NTL::ZZ(1);
        for (long j = 1; j <= degree + 1; j++) {
            if (i != j) {
                NTL::MulMod(enume, enume, j, P);
                NTL::MulMod(denom, denom, j - i, P);
            }
        }
        NTL::ZZ tmp;
        NTL::MulMod(tmp, shares[i - 1].b, enume, P); // s * upsilon * \pi(j)
        NTL::InvMod(denom, denom, P);
        NTL::MulMod(tmp, tmp, denom,
                    P);             // l_i = s * upsilon * (\pi(j)) / (\pi(j-i)) -- at here we get partial reconstruction

        partial_recon.append(NTL::Pair<long, NTL::ZZ>(i, tmp));
    }
    return partial_recon;
}

ShamirShares ShamirShares::__add__(ShamirShares other) {
    __check_2_op_correct__(other);

    long new_degree = max(degree, other.degree);
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> new_shares;

    for (int i = 1; i <= min(new_degree + 1, shares.length()); i++) {
        NTL::ZZ new_eval;
        NTL::AddMod(new_eval, shares[i - 1].b, other.shares[i - 1].b, P);
        new_shares.append(NTL::Pair<long, NTL::ZZ>((long) (i), new_eval));
    }
    ShamirShares new_sharing = ShamirShares(new_shares, P, n, new_degree);
    return new_sharing;
}

ShamirShares ShamirShares::__mult__(ShamirShares other) {
    __check_2_op_correct__(other);

    long new_degree = degree + other.degree;
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> new_shares;

    for (int i = 1; i <= min(new_degree + 1, shares.length()); i++) {
        NTL::ZZ new_eval;
        NTL::MulMod(new_eval, shares[i - 1].b, other.shares[i - 1].b, P);
        new_shares.append(NTL::Pair<long, NTL::ZZ>((long) (i), new_eval));
    }
    ShamirShares new_sharing = ShamirShares(new_shares, P, n, new_degree);
    return new_sharing;
}

bool ShamirShares::__check_points_in_order__(NTL::Vec<NTL::Pair<long, NTL::ZZ>> &input) {
    long i = 1;
    for (const auto &x : input) {
        if (x.a != i) {
            std::cout << "Input points are not in order of 1, 2, ..., t+1" << endl;
            return false;
        }
        i += 1;
    }
    return true;
}

bool ShamirShares::__check_if_compatible__(const ShamirShares &other) {
    if (!(P == other.P && n == other.n && degree == other.degree)) {
        std::cout << "two shares are not compatible" << endl;
        return false;
    }
    return true;
};

void ShamirShares::__check_2_op_correct__(ShamirShares other) {
    bool flag = true;
    flag = __check_if_compatible__(other);
    if (shares.length() != other.shares.length()) {
        std::cout << "Length of share do not match. a.length = " + to_string(shares.length())
                     + ", b.length = " + to_string(other.shares.length()) << endl;
    }
    flag = __check_points_in_order__(shares);
    flag = __check_points_in_order__(other.shares);
    assert(flag);
}
