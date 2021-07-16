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

#include <vector>
#include "../include/gtest/gtest.h"
#include "../include/shamir_secret_sharing.h"
#include "../include/shamir_secret_sharing_integer.h"

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

template<typename SSSST, typename ST, typename DT>
DT test_sharing(NTL::ZZ &P, long &n, long &t, DT s) {
    SSSST Scheme = SSSST(P, n, t);
    ST share = Scheme.share_secret(s);
    DT res = share.reconstruct_secret(P);
    return res;
}

template<typename T, typename S>
void test_mul(T a, T b, NTL::Vec<S> &res, NTL::ZZ modulo) {
    T new_share_1 = a.__mult__(b);
    T new_share_2 = b.__mult__(a);

    S res_1 = new_share_1.reconstruct_secret(modulo);
    S res_2 = new_share_2.reconstruct_secret(modulo);
    res.append(res_1);
    res.append(res_2);
}

template<typename T, typename S>
void test_add(T a, T b, NTL::Vec<S> &res, NTL::ZZ modulo) {
    T new_share_1 = a.__add__(b);
    T new_share_2 = b.__add__(a);

    S res_1 = new_share_1.reconstruct_secret(modulo);
    S res_2 = new_share_2.reconstruct_secret(modulo);
    res.append(res_1);
    res.append(res_2);
}


TEST_F(params, TestSharingInteger) {
    SetUp(1024, 3, 1);
    NTL::ZZ res = test_sharing<ShamirSecretSharingSchemeInteger, ShamirShares_integer, NTL::ZZ>(P, n, t,
                                                                                                to_ZZ(1000000));
    ASSERT_EQ(1000000, to_int(res));
}

TEST_F(params, TestSharing) {
    SetUp(1024, 3, 1);
    NTL::ZZ res = test_sharing<ShamirSecretSharingScheme, ShamirShares, NTL::ZZ>(P, n, t, to_ZZ(1000));
    ASSERT_EQ(1000, to_int(res));
}

TEST_F(params, TestAddMul) {
    SetUp(1024, 3, 1);
    NTL::ZZ s1, s2;
    NTL::RandomLen(s1, length);
    NTL::RandomLen(s2, length);
    ShamirSecretSharingScheme ssss = ShamirSecretSharingScheme(P, n, t);
    ShamirShares a = ssss.share_secret(s1);
    ShamirShares b = ssss.share_secret(s2);
    NTL::Vec<ZZ> res_2;
    test_add<ShamirShares, NTL::ZZ>(a, b, res_2, P);
    ASSERT_TRUE((res_2[0] == res_2[1]) && (s1 + s2 == res_2[1])); // <a> + <b> == <b> + <a> && <a+b> == <a> + <b>
    res_2.kill();
    test_mul<ShamirShares, NTL::ZZ>(a, b, res_2, P);
    ASSERT_TRUE((res_2[0] == res_2[1]) && (s1 * s2 == res_2[1]));
    res_2.kill();
}

TEST_F(params, TestAddMulInteger) {
    SetUp(1024, 3, 1);
    NTL::ZZ s1, s2;
    NTL::RandomLen(s1, length);
    NTL::RandomLen(s2, length);
    ShamirSecretSharingSchemeInteger ssss = ShamirSecretSharingSchemeInteger(P, n, t);
    ShamirShares_integer a = ssss.share_secret(s1);
    ShamirShares_integer b = ssss.share_secret(s2);
    NTL::Vec<ZZ> res_2;
    test_add<ShamirShares_integer, NTL::ZZ>(a, b, res_2, P);
    ASSERT_TRUE((res_2[0] == res_2[1]) && (s1 + s2 == res_2[1])); // <a> + <b> == <b> + <a> && <a+b> == <a> + <b>
    res_2.kill();
    test_mul<ShamirShares_integer, NTL::ZZ>(a, b, res_2, P);
    ASSERT_TRUE((res_2[0] == res_2[1]) && (s1 * s2 == res_2[1]));
    res_2.kill();
}

int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}