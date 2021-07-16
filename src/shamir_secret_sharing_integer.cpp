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
#include <math.h>
#include "../include/shamir_secret_sharing_integer.h"

ShamirShares_integer
ShamirSecretSharingSchemeInteger::share_secret(const NTL::ZZ &s) {
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> shares;
    NTL::ZZ f0(factorial_small(n) * s);
    NTL::Vec<NTL::ZZ> poly_param;

    NTL::ZZ poly_param_bound = 2 * randomness_interval + 1;
    for (int i = 1; i < t + 1; i++) {
        poly_param.append(NTL::RandomBnd(poly_param_bound));
    }

    // get n points from the polynomial
    // n points from 1 to n
    // polynomial has degree of t
    for (int i = 1; i <= n; i++) {
        NTL::ZZ secret(f0);
        for (int j = 1; j < t + 1; j++) {
            NTL::ZZ tmp;
            NTL::mul(tmp, poly_param[j - 1], (long) pow(i, j));
            NTL::add(secret, secret, tmp);
        }
        shares.append(NTL::Pair<long, NTL::ZZ>((long) (i), secret));
    }

    NTL::ZZ scaling = factorial_small(n);
    ShamirShares_integer Sharing = ShamirShares_integer(shares, P, n, t, scaling);
    return Sharing;
}


NTL::ZZ
ShamirShares_integer::reconstruct_secret(NTL::ZZ modulu) {
    if (shares.length() < degree + 1) {
        std::cout << "no enough points for recon. Need = " +
                     to_string(degree + 1) + ", get = " + to_string(shares.length())
                  << endl;
        assert(shares.length() >= degree + 1);
    }

    NTL::ZZ partial_recon_sum = NTL::ZZ(0);
    NTL::ZZ res;
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> partial_recon = get_partial_recon(modulu);
    for (long i = 1; i <= partial_recon.length(); i++) {
        NTL::AddMod(partial_recon_sum, partial_recon_sum,
                    partial_recon[i - 1].b, modulu);
        // sum the partital recon
    }
    NTL::ZZ tmp;
    NTL::InvMod(tmp, scaling * n_fac, modulu);
    NTL::MulMod(res, partial_recon_sum, tmp, P);
    return res;
}

NTL::Vec<NTL::Pair<long, NTL::ZZ>>
ShamirShares_integer::get_partial_recon(NTL::ZZ modulu) {
    if (shares.length() < degree + 1) {
        std::cout << "no enough points for recon. Need = " + to_string(degree + 1) + ", get = " +
                     to_string(shares.length()) << endl;
        assert(shares.length() >= degree + 1);
    }
    if (scaling % modulu == 0) {
        std::cout << "scaling factor should not be devided by N." << endl;
        assert(scaling % modulu != 0);
    }
    assert(__check_points_in_order__(shares));

    NTL::Vec<NTL::Pair<long, NTL::ZZ>> partial_recon;
    for (long i = 1; i <= degree + 1; i++) {
        NTL::ZZ enume = NTL::ZZ(1);
        NTL::ZZ denom = NTL::ZZ(1);
        for (long j = 1; j <= degree + 1; j++) {
            if (i != j) {
                enume *= j;
                denom *= j - i;
            }
        }
        NTL::ZZ tmp;
        NTL::mul(enume, enume, n_fac);
        assert(enume % denom == 0);
        NTL::div(tmp, enume, denom);
        NTL::MulMod(tmp, tmp, shares[i - 1].b, modulu);

        partial_recon.append(NTL::Pair<long, NTL::ZZ>(i, tmp));
    }
    return partial_recon;
}

ShamirShares_integer
ShamirShares_integer::__add__(ShamirShares_integer other) {
    __check_2_op_doable__(other);

    long new_degree = NTL::max(degree, other.degree);
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> new_shares;
    for (int i = 1; i <= min(new_degree + 1, shares.length()); i++) {
        NTL::ZZ new_eval;
        NTL::add(new_eval, shares[i - 1].b, other.shares[i - 1].b);
        new_shares.append(NTL::Pair<long, NTL::ZZ>((long) (i), new_eval));
    }
    ShamirShares_integer new_sharing = ShamirShares_integer(new_shares, P, n, new_degree, scaling);
    return new_sharing;
}

ShamirShares_integer ShamirShares_integer::__mult__(ShamirShares_integer other) {
    // std::cout << "doing integer mult" << endl;
    __check_2_op_doable__(other);
    long new_degree = degree + other.degree;
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> new_shares;
    for (int i = 1; i <= min(new_degree + 1, shares.length()); i++) {
        NTL::ZZ new_eval;
        NTL::mul(new_eval, shares[i - 1].b, other.shares[i - 1].b);
        new_shares.append(NTL::Pair<long, NTL::ZZ>((long) (i), new_eval));
    }

    NTL::ZZ new_scaling = scaling * other.scaling;
    ShamirShares_integer new_sharing = ShamirShares_integer(new_shares, P, n, new_degree, new_scaling);
    return new_sharing;
}

/**
 * when a integer share times K
*/
ShamirShares_integer ShamirShares_integer::__mult__(long K) {
    long new_degree = degree;
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> new_shares;
    for (int i = 1; i <= min(new_degree + 1, shares.length()); i++) {
        NTL::ZZ new_eval;
        NTL::mul(new_eval, shares[i - 1].b, K);
        new_shares.append(NTL::Pair<long, NTL::ZZ>((long) (i), new_eval));
    }

    NTL::ZZ new_scaling = scaling;
    ShamirShares_integer new_sharing = ShamirShares_integer(new_shares, P, n, new_degree, new_scaling);
    return new_sharing;
}

bool ShamirShares_integer::__check_points_in_order__(NTL::Vec<NTL::Pair<long, NTL::ZZ>> &input) {
    long i = 1;
    for (auto x : input) {
        if (x.a != i) {
            std::cout << "Input points are not in order of 1, 2, ..., t+1" << endl;
            return false;
        }
        i += 1;
    }
    return true;
}

bool ShamirShares_integer::__check_if_compatible__(ShamirShares_integer other) {
    if (!(P == other.P && n == other.n && degree == other.degree)) {
        std::cout << "two shares are not compatible" << endl;
        return false;
    }
    return true;
};

void ShamirShares_integer::__check_2_op_doable__(ShamirShares_integer other) {
    bool flag = true;
    flag = __check_if_compatible__(other);
    if (shares.length() != other.shares.length()) {
        std::cout << "Length of share do not match. a.length = " + to_string(shares.length()) + ", b.length = " +
                     to_string(other.shares.length()) << endl;
    }
    flag = __check_points_in_order__(shares);
    flag = __check_points_in_order__(other.shares);
    assert(flag);
}
