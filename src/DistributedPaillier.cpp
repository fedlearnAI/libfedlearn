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
#include <jni.h>
// #include <omp.h>
#include "../include/com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative.h"
#include "../include/cppWrapper.hpp"
#include "../include/DistributedPaillierHepler.h"

using namespace std;
using namespace c_api;

/**
 * ============================
 * Utils
 * ============================
 */

long witness(const ZZ &n, const ZZ &x) {
    ZZ m, y, z;
    long j, k;
    if (x == 0) return 0;
    // compute m, k such that n-1 = 2^k * m, m odd:
    k = 1;
    m = n / 2;
    while (m % 2 == 0) {
        k++;
        m /= 2;
    }
    z = PowerMod(x, m, n); // z = x^m % n
    if (z == 1) return 0;
    j = 0;
    do {
        y = z;
        z = (y * y) % n;
        j++;
    } while (j < k && z != 1);
    return z != 1 || y != n - 1;
}


long PrimeTest(const ZZ &n, long t) {
    if (n <= 1) return 0;
    // first, perform trial division by primes up to 2000
    PrimeSeq s;  // a class for quickly generating primes in sequence
    long p;
    p = s.next();  // first prime is always 2
    while (p && p < 2000) {
        if ((n % p) == 0) return (n == p);
        p = s.next();
    }
    // second, perform t Miller-Rabin tests
    ZZ x;
    long i;
    for (i = 0; i < t; i++) {
        x = RandomBnd(n); // random number between 0 and n-1
        if (witness(n, x))
            return 0;
    }
    return 1;
}

template<typename T>
void __print_vec_as_int__(const std::vector<T> &in) {
    for (auto x : in) {
        std::cout << (int) x << " ";
    }
    std::cout << std::endl;
}

void __print_zz_vec__(const std::vector<ZZ> &in) {
    for (const auto &x : in) {
        std::cout << x << " ";
    }
    std::cout << std::endl;
}

void __print_cypherVec_as_int__(const std::vector<cypher_text_t> &in) {
    std::cout << "cypher text len = " << in.size() << std::endl;
    for (const auto &x : in) {
        std::cout << "isNeg = " << x.is_neg << " ";
        __print_vec_as_int__(x.c);
    }
    std::cout << std::endl;
}

ZZ geneNoneZeroRandBnd(const ZZ &bnd) {
    ZZ out(0);
    while (out == 0) {
        out = NTL::RandomBnd(bnd);
    }
    return out;
}

jbyteArray as_jbyte_array(JNIEnv *env, const CHAR_LIST &in) {
    // copy the original data for safety reason
    CHAR_LIST cache(in);

    int len = (int) in.size();
    jbyteArray array = env->NewByteArray(len);
    if (nullptr == array)
        return nullptr;

    // convert c type into java type
    auto tmp = cache.data();
    auto *j_version = reinterpret_cast<jbyte *>(tmp);
    // for (int i = 0; i < len; i++)
    // std::cout << "char list after reinterpret_cast = " << (int)(j_version[i]) << endl;
    env->SetByteArrayRegion(array, 0, len, j_version);
    // std::cout << "char list after reinterpret_cast OK " << endl;
    return array;
}

void as_uncharlist(CHAR_LIST &out, JNIEnv *env, jobject array) {
    auto *arr = reinterpret_cast<jbyteArray *>(&array);
    const int len = env->GetArrayLength(*arr);
    if (len == 0) {
        std::cout << "Error, allocating and free space of size 0 " << endl;
        assert(len != 0);
    }
    auto *buf = new unsigned char[len];
    env->GetByteArrayRegion(*arr, 0, len, reinterpret_cast<jbyte *>(buf));
    out.reserve(len);
    for (int i = 0; i < len; i++) {
        out.push_back(buf[i]);
    }
    delete[] buf;
    buf = nullptr;
}

jobjectArray as_2djbyte_array(JNIEnv *env,
                              LIST_OF_CHAR_LIST &in) {
    unsigned int length1D = in.size();

    jclass byteArrayClass = env->FindClass("[B"); // FIXME
    jobjectArray out = env->NewObjectArray((jsize) length1D, byteArrayClass, nullptr);

    // Go through the firs dimension and add the second dimension arrays
    for (unsigned int i = 0; i < length1D; i++) {
        jbyteArray inner_arr = as_jbyte_array(env, in[i]);
        env->SetObjectArrayElement(out, (jsize) i, inner_arr);
        env->DeleteLocalRef(inner_arr);
    }
    return out;
}

dist_paillier_privkey_t
__privkey_fromJprivkey_helper__(JNIEnv *env,
                                jobjectArray jPrivkey,
                                jclass cls) {
    jint jlen = env->GetArrayLength(jPrivkey);

    assert(jlen == 3);
    NTL::ZZ hi, n, theta_invmod;
    CHAR_LIST out_hi, out_n, out_theta_invmod;
    bool hi_isNeg, n_isNeg, theta_invmod_isNeg;
    __get_arrOfSignedbyteArray_helper__(env, jPrivkey, out_hi, hi_isNeg, 0);
    __get_arrOfSignedbyteArray_helper__(env, jPrivkey, out_n, n_isNeg, 1);
    __get_arrOfSignedbyteArray_helper__(env, jPrivkey, out_theta_invmod, theta_invmod_isNeg, 2);

    distributed_paillier::char_list_2_ZZ(hi, out_hi, hi_isNeg);
    distributed_paillier::char_list_2_ZZ(n, out_n, n_isNeg);
    distributed_paillier::char_list_2_ZZ(theta_invmod, out_theta_invmod, theta_invmod_isNeg);

    dist_paillier_privkey_t priv_key;
    priv_key.hi = hi;
    priv_key.n = n;
    priv_key.n_squared = NTL::sqr(n);
    priv_key.theta_invmod = theta_invmod;

    // std::cout << "priv_key.hi = " << priv_key.hi << " "
    //           << "priv_key.n = " << priv_key.n << " "
    //           << "priv_key.theta_invmod = " << priv_key.theta_invmod
    //           << endl;
    return priv_key;
}

dist_paillier_pubkey_t
__pubkey_fromJpubkey_helper__(JNIEnv *env,
                              jobject jpubkey,
                              jclass cls) {
    NTL::ZZ n;
    CHAR_LIST out_n;
    bool n_isNeg;

    __get_SignedbyteArray_helper__(env, jpubkey, out_n, n_isNeg);

    distributed_paillier::char_list_2_ZZ(n, out_n, n_isNeg);
    dist_paillier_pubkey_t pub_key;

    pub_key.n = n;
    pub_key.n_plusone = n + 1;
    pub_key.n_squared = NTL::sqr(n);
    return pub_key;
}

void __get_SignedbyteArray_helper__(JNIEnv *env,
                                    const jobject &SbyteArrayObj,
                                    CHAR_LIST &target_lst,
                                    bool &target_isNeg) {
    jclass jcls_jobj_sbyteArray = env->GetObjectClass(SbyteArrayObj);
    if (jcls_jobj_sbyteArray == nullptr) {
        std::cout << "Error in Native: jcls_jobj_sbyteArray = env->GetObjectClass(SbyteArrayObj) is NULL" << endl;
        return;
    }

    jfieldID fid_byteArr = env->GetFieldID(jcls_jobj_sbyteArray, "byteArr", "[B");
    if (fid_byteArr == nullptr) {
        std::cout << "Error in Native: fid_byteArr = env->GetFieldID(jcls_jobj_sbyteArray... is NULL" << endl;
        return;
    }

    jfieldID fid_isNeg = env->GetFieldID(jcls_jobj_sbyteArray, "isNeg", "Z");
    if (fid_isNeg == nullptr) {
        std::cout << "Error in Native: fid_isNeg = env->GetFieldID(jcls_jobj_sbyteArray, is NULL" << endl;
        return;
    }

    jobject obj_byteArr = env->GetObjectField(SbyteArrayObj, fid_byteArr);
    if (obj_byteArr == nullptr) {
        std::cout << "Error in Native: obj_byteArr = env->GetObjectField(SbyteArrayObj, fid_byteArr) is NULL" << endl;
        return;
    }

    jboolean isNeg = env->GetBooleanField(SbyteArrayObj, fid_isNeg);
    as_uncharlist(target_lst, env, obj_byteArr);
    target_isNeg = isNeg;

    // std::cout << "target_isNeg = " << target_isNeg << " isNeg = " << (int)isNeg << endl;

    env->DeleteLocalRef(jcls_jobj_sbyteArray);
    env->DeleteLocalRef(obj_byteArr);
}

void __get_arrOfSignedbyteArray_helper__(JNIEnv *env,
                                         const jobjectArray &arrOfByteArrayObj,
                                         CHAR_LIST &target_lst,
                                         bool &target_isNeg,
                                         const int idx) {
    jobject jobj_signedByteArray = env->GetObjectArrayElement(arrOfByteArrayObj, idx);
    __get_SignedbyteArray_helper__(env, jobj_signedByteArray, target_lst, target_isNeg);
    env->DeleteLocalRef(jobj_signedByteArray);
}

void __fill_SByteArray_helper__(JNIEnv *env,
                                jobject &sByteArrayObj,
                                const CHAR_LIST &src_lst,
                                const bool &src_isNeg) {
    jbyteArray byte_arr = as_jbyte_array(env, src_lst);
    jclass cls_sbyte = env->GetObjectClass(sByteArrayObj);
    if (cls_sbyte == nullptr) {
        std::cout << "Error in Native: env->GetObjectClass(sByteArrayObj) is NULL" << endl;
        return;
    }

    jfieldID fid_byteArray = env->GetFieldID(cls_sbyte, "byteArr", "[B");
    env->SetObjectField(sByteArrayObj, fid_byteArray, byte_arr);
    jfieldID fid_isNeg = env->GetFieldID(cls_sbyte, "isNeg", "Z");
    env->SetBooleanField(sByteArrayObj, fid_isNeg, src_isNeg);

    env->DeleteLocalRef(byte_arr);
    env->DeleteLocalRef(cls_sbyte);
}

void __fill_arrOfByteArray_helper__(JNIEnv *env,
                                    jobjectArray arrOfByteArrayObj,
                                    const CHAR_LIST &src_lst,
                                    const bool &src_isNeg,
                                    const int &idx) {
    jbyteArray byte_arr = as_jbyte_array(env, src_lst);

    jobject innerobj = env->GetObjectArrayElement(arrOfByteArrayObj, idx);
    if (innerobj == nullptr) {
        std::cout << "Error in Native: env->GetObjectArrayElement(arrObj, cnt) is NULL" << endl;
        return;
    }

    jclass cls_innerobj = env->GetObjectClass(innerobj);
    if (cls_innerobj == nullptr) {
        std::cout << "Error in Native: env->GetObjectClass(innerobj) is NULL" << endl;
        return;
    }

    jfieldID fid_byteArray = env->GetFieldID(cls_innerobj, "byteArr", "[B");
    env->SetObjectField(innerobj, fid_byteArray, byte_arr);
    jfieldID fid_isNeg = env->GetFieldID(cls_innerobj, "isNeg", "Z");
    env->SetBooleanField(innerobj, fid_isNeg, src_isNeg);

    env->DeleteLocalRef(byte_arr);
    env->DeleteLocalRef(innerobj);
    env->DeleteLocalRef(cls_innerobj);
}

ZZ signedbyteArray_2_ZZ(JNIEnv *env,
                        const jobject &SbyteArrayObj) {
    NTL::ZZ value;
    CHAR_LIST out_value;
    bool value_isNeg;
    __get_SignedbyteArray_helper__(env, SbyteArrayObj, out_value, value_isNeg);
    distributed_paillier::char_list_2_ZZ(value, out_value, value_isNeg);
    return value;
}

std::vector<ZZ> arrOfSignedbyteArray_2_zzList(JNIEnv *env, const jobjectArray &arrOfByteArrayObj) {
    int array_len = env->GetArrayLength(arrOfByteArrayObj);
    std::vector<ZZ> ret;
    ret.reserve(array_len);
    for (int i = 0; i < array_len; i++) {
        jobject jobj_signedByteArray = env->GetObjectArrayElement(arrOfByteArrayObj, i);
        NTL::ZZ value = signedbyteArray_2_ZZ(env, jobj_signedByteArray);
        env->DeleteLocalRef(jobj_signedByteArray);
        ret.push_back(value);
    }
    return ret;
}

void fill_SByteArray_with_zz(JNIEnv *env,
                             jobject &sByteArrayObj,
                             const ZZ &zz_number) {
    CHAR_LIST zz_number_in_char;
    bool zz_number_isNeg;
    distributed_paillier::ZZ_2_byte(zz_number_in_char, zz_number, zz_number_isNeg);
    __fill_SByteArray_helper__(env, sByteArrayObj, zz_number_in_char, zz_number_isNeg);
}

void fill_arrOfSByteArray_with_zz(JNIEnv *env,
                                  const jobjectArray &arrOfByteArrayObj,
                                  const ZZ &zz_number,
                                  const int &idx) {
    CHAR_LIST zz_number_in_char;
    bool zz_number_isNeg;
    distributed_paillier::ZZ_2_byte(zz_number_in_char, zz_number, zz_number_isNeg);
    __fill_arrOfByteArray_helper__(env, arrOfByteArrayObj, zz_number_in_char, zz_number_isNeg, idx);
}

/**
 * ============================
 * JNI layer functions
 * ============================
 */


/**
 * Generate Priv and PubKey for all parties
 * ==========
 * Parameters
 *  arrObj: JAVA obj to be filled by native code. 
 *          MUST be initialized before passing to native !
 * Return
 *  void. JAVA object fields are filled in native layer.
 */
JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1generate_1privpub_1key_1_1(
        JNIEnv *env,
        jclass cls,
        jobjectArray arrObj,
        jint len,
        jint t,
        jint n,
        jboolean is_dealer) {
    std::vector<dist_paillier_privkey_t> ret = distributed_paillier::gene_privkey_standalone(len, n, t);
    jint list_size = env->GetArrayLength(arrObj);
    assert(list_size == (3 * n)); // 每个party有3个bytearray要返回
    assert(ret.size() == n);
    // std::cout << list_size << "\t" << n << endl;
    // std::cout << "ret[0].n = " << ret[0].n << endl;

    int cnt = 0;
    for (int i = 0; i < n; i++) {
        // hi
        {
            CHAR_LIST hi;
            bool hi_isNeg;
            distributed_paillier::ZZ_2_byte(hi, ret[i].hi, hi_isNeg);
            // std::cout << "ret[i].hi = " << ret[i].hi << " in vec : ";
            // __print_vec_as_int__(hi);
            __fill_arrOfByteArray_helper__(env, arrObj, hi, hi_isNeg, cnt++);
        }

        // N
        {
            CHAR_LIST N;
            bool N_isNeg;
            distributed_paillier::ZZ_2_byte(N, ret[i].n, N_isNeg);
            // std::cout << "ret[i].n = ";
            // __print_vec_as_int__(N);
            __fill_arrOfByteArray_helper__(env, arrObj, N, N_isNeg, cnt++);
        }

        // theta_invmod
        {
            CHAR_LIST theta_invmod;
            bool theta_invmod_isNeg;
            distributed_paillier::ZZ_2_byte(theta_invmod, ret[i].theta_invmod, theta_invmod_isNeg);
            // std::cout << "ret[i].theta_invmod = " << ret[i].theta_invmod << " in vec : ";
            // __print_vec_as_int__(theta_invmod);
            __fill_arrOfByteArray_helper__(env, arrObj, theta_invmod, theta_invmod_isNeg, cnt++);
        }
    }
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1enc_1_1(JNIEnv *env,
                                                                                                jclass cls,
                                                                                                jobject out,
                                                                                                jlong plain,
                                                                                                jint bits,
                                                                                                jobject jpubkey) {
    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;

    // std::cout << "pub_key.n = " << pub_key.n << endl;

    NTL::ZZ x(plain), xx;
    cypher_text_t ret(bits);
    c_api::enc(ret, x, pub_key);
    // distributed_paillier::char_list_2_ZZ(xx, ret.c, ret.is_neg);
    __fill_SByteArray_helper__(env, out, ret.c, ret.is_neg);
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1enc_1list_1_1(JNIEnv *env,
                                                                                                      jclass cls,
                                                                                                      jobjectArray out_lst,
                                                                                                      jlongArray plain_lst,
                                                                                                      jint bits,
                                                                                                      jobject jpubkey) {
    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    jsize const len = env->GetArrayLength(plain_lst);
    jsize const lenout = env->GetArrayLength(out_lst);
    assert(lenout == len);

    jlong *body = env->GetLongArrayElements(plain_lst, 0);
    std::vector<cypher_text_t> c_text_lst;
    c_text_lst.reserve(len);

    for (int i = 0; i < len; i++) {
        cypher_text_t ret(bits);
        c_text_lst.push_back(ret);
    }

    // #pragma omp parallel num_threads(4)
    {
        // int nthreads = omp_get_num_threads();
        // cout << "using parallel?" << " num thread = " << nthreads << endl;
//#pragma omp parallel for shared(c_text_lst)
        for (int i = 0; i < len; i++) {
            NTL::ZZ x(body[i]);
            // cypher_text_t ret(bits);
            // c_text_lst[i]  = ret;
            c_api::enc(c_text_lst[i], x, pub_key);
            // c_text_lst.push_back(ret);
            // c_text_lst[i] = ret;

            NTL::ZZ xx;
            distributed_paillier::char_list_2_ZZ(xx, c_text_lst[i].c, c_text_lst[i].is_neg);
            // cout << "body[i] = " << body[i]
            //      << " x = " << x
            //      << " xx = " << xx
            //      << endl;
        }
    }
    for (int i = 0; i < len; i++) {
        __fill_arrOfByteArray_helper__(env, out_lst, c_text_lst[i].c, c_text_lst[i].is_neg, i);
    }
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1partial_1dec_1_1(JNIEnv *env,
                                                                                                         jclass cls,
                                                                                                         jobject out,
                                                                                                         jobject in,
                                                                                                         jint bits,
                                                                                                         jint party_id,
                                                                                                         jint t,
                                                                                                         jlong n_fact,
                                                                                                         jobjectArray jprivkey) {
    dist_paillier_privkey_t priv_key = __privkey_fromJprivkey_helper__(env, jprivkey, cls);
    priv_key.t = t;
    priv_key.n_fat = n_fact;
    cypher_text_t partial_dec_text(bits), cypher_text(bits);

    // cout <<  priv_key.hi  << priv_key.n << priv_key.n_squared << priv_key.theta_invmod << priv_key.t << priv_key.n_fat << endl;

    __get_SignedbyteArray_helper__(env, in, cypher_text.c, cypher_text.is_neg);

    // cout << "in info -- cypher_text.c : "; __print_vec_as_int__(cypher_text.c);
    // cout << "cypher_text.is_neg : " << cypher_text.is_neg << endl;

    partial_dec(party_id, t, partial_dec_text, cypher_text, priv_key);

    __fill_SByteArray_helper__(env, out, partial_dec_text.c, partial_dec_text.is_neg);
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1partial_1dec_1lst_1_1(
        JNIEnv *env,
        jclass cls,
        jobjectArray out_lst,
        jobjectArray in_lst, jint bits,
        jint party_id,
        jint t,
        jlong n_fact,
        jobjectArray jprivkey) {
    const jsize len = env->GetArrayLength(in_lst);
    const jsize len2 = env->GetArrayLength(out_lst);
    assert(len == len2);

    dist_paillier_privkey_t priv_key = __privkey_fromJprivkey_helper__(env, jprivkey, cls);
    priv_key.t = t;
    priv_key.n_fat = n_fact;

    std::vector<cypher_text_t> c_text_lst, im_res_lst;
    c_text_lst.reserve(len);
    im_res_lst.reserve(len);

    // FIXME: 三个循环并在一起是否会效率更高?
    for (int i = 0; i < len; i++) {
        cypher_text_t tmp(bits);
        __get_arrOfSignedbyteArray_helper__(env, in_lst, tmp.c, tmp.is_neg, i);
        c_text_lst.push_back(tmp);
    }

    for (int i = 0; i < len; i++) {
        cypher_text_t tmp(bits);
        partial_dec(party_id, t, tmp, c_text_lst[i], priv_key);
        im_res_lst.push_back(tmp);
    }

    for (int i = 0; i < len; i++) {
        __fill_arrOfByteArray_helper__(env, out_lst, im_res_lst[i].c, im_res_lst[i].is_neg, i);
    }
}

JNIEXPORT jlong JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1final_1dec_1_1(JNIEnv *env,
                                                                                                       jclass cls,
                                                                                                       jobjectArray im_res,
                                                                                                       jobject jcypher_text,
                                                                                                       jint bits,
                                                                                                       jint t,
                                                                                                       jlong n_fact,
                                                                                                       jlong max_neg_abs,
                                                                                                       jobjectArray jprivkey) {
    const jint inlen = env->GetArrayLength(im_res);
    assert(inlen >= 2 * t + 1);

    dist_paillier_privkey_t priv_key = __privkey_fromJprivkey_helper__(env, jprivkey, cls);
    priv_key.t = t;
    priv_key.n_fat = n_fact;

    std::vector<cypher_text_t> partial_dec_lst;
    partial_dec_lst.reserve(inlen);
    for (int i = 0; i < inlen; i++) {
        cypher_text_t tmp(bits);
        __get_arrOfSignedbyteArray_helper__(env, im_res, tmp.c, tmp.is_neg, i);
        partial_dec_lst.push_back(tmp);
    }
    cypher_text_t cypher_text(bits);
    __get_SignedbyteArray_helper__(env, jcypher_text, cypher_text.c, cypher_text.is_neg);
    NTL::ZZ res = c_api::final_dec(cypher_text, partial_dec_lst, t, NTL::to_ZZ(max_neg_abs), priv_key);

    long long_res = NTL::to_long(res);
    return long_res;
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1add_1_1(JNIEnv *env,
                                                                                                jclass cls,
                                                                                                jobject out,
                                                                                                jobject a,
                                                                                                jobject b,
                                                                                                jint bits,
                                                                                                jobject jpubkey) {
    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    cypher_text_t a_c(bits), b_c(bits), res(bits);
    __get_SignedbyteArray_helper__(env, a, a_c.c, a_c.is_neg);
    __get_SignedbyteArray_helper__(env, b, b_c.c, b_c.is_neg);
    c_api::add(res, a_c, b_c, pub_key);
    __fill_SByteArray_helper__(env, out, res.c, res.is_neg);
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1mul_1_1(JNIEnv *env,
                                                                                                jclass cls,
                                                                                                jobject out,
                                                                                                jobject a,
                                                                                                jlong b,
                                                                                                jint bits,
                                                                                                jobject jpubkey) {
    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    cypher_text_t a_c(bits), res(bits);
    __get_SignedbyteArray_helper__(env, a, a_c.c, a_c.is_neg);
    c_api::mul(res, a_c, NTL::to_ZZ(b), pub_key);
    __fill_SByteArray_helper__(env, out, res.c, res.is_neg);
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1mul_1elementWise_1_1(
        JNIEnv *env,
        jclass cls,
        jobjectArray out_lst,
        jobjectArray a_lst,
        jlongArray b_lst,
        jint bits,
        jobject jpubkey) {
    jint len = env->GetArrayLength(a_lst);
    jint len2 = env->GetArrayLength(b_lst);
    jint len3 = env->GetArrayLength(out_lst);
    assert(len == len2 && len == len3 && len > 0);

    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    std::vector<cypher_text_t> a_text_lst, res_lst;
    a_text_lst.reserve(len);
    jlong *body = env->GetLongArrayElements(b_lst, 0);

    for (int i = 0; i < len; i++) {
        cypher_text_t tmp_a(bits);
        __get_arrOfSignedbyteArray_helper__(env, a_lst, tmp_a.c, tmp_a.is_neg, i);
        a_text_lst.push_back(tmp_a);

        ZZ tmpp;
        distributed_paillier::char_list_2_ZZ(tmpp, tmp_a.c, tmp_a.is_neg);
        // cout << "a_text_lst as ZZ: " << tmpp << endl;
        // cout << "in info -- b : " << NTL::to_ZZ(body[i]) << endl;
    }
    // cout << "in info -- a_text_lst: "<< endl; __print_cypherVec_as_int__(a_text_lst);

    for (int i = 0; i < len; i++) {
        cypher_text_t tmp(bits);
        c_api::mul(tmp, a_text_lst[i], NTL::to_ZZ(body[i]), pub_key);
        res_lst.push_back(tmp);
    }

    for (int i = 0; i < len; i++) {
        __fill_arrOfByteArray_helper__(env, out_lst, res_lst[i].c, res_lst[i].is_neg, i);
    }
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1div_1_1(JNIEnv *env,
                                                                                                jclass cls,
                                                                                                jobject out,
                                                                                                jobject a,
                                                                                                jlong b,
                                                                                                jint bits,
                                                                                                jobject jpubkey) {
    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    cypher_text_t a_c(bits), res(bits);
    __get_SignedbyteArray_helper__(env, a, a_c.c, a_c.is_neg);
    c_api::div(res, a_c, NTL::to_ZZ(b), pub_key);
    __fill_SByteArray_helper__(env, out, res.c, res.is_neg);
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1div_1elementWise_1_1(
        JNIEnv *env,
        jclass cls,
        jobjectArray out_lst,
        jobjectArray a_lst,
        jlongArray b_lst,
        jint bits,
        jobject jpubkey) {
    jint len = env->GetArrayLength(a_lst);
    jint len2 = env->GetArrayLength(b_lst);
    jint len3 = env->GetArrayLength(out_lst);
    assert(len == len2 && len == len3 && len > 0);

    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    std::vector<cypher_text_t> a_text_lst, res_lst;
    a_text_lst.reserve(len);
    jlong *body = env->GetLongArrayElements(b_lst, 0);

    for (int i = 0; i < len; i++) {
        cypher_text_t tmp_a(bits);
        __get_arrOfSignedbyteArray_helper__(env, a_lst, tmp_a.c, tmp_a.is_neg, i);
        a_text_lst.push_back(tmp_a);
    }
    for (int i = 0; i < len; i++) {
        cypher_text_t tmp(bits);
        c_api::div(tmp, a_text_lst[i], NTL::to_ZZ(body[i]), pub_key);
        res_lst.push_back(tmp);
    }

    for (int i = 0; i < len; i++) {
        __fill_arrOfByteArray_helper__(env, out_lst, res_lst[i].c, res_lst[i].is_neg, i);
    }
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1add_1vec_1_1(JNIEnv *env,
                                                                                                     jclass cls,
                                                                                                     jobjectArray out_lst,
                                                                                                     jobjectArray a_lst,
                                                                                                     jobjectArray b_lst,
                                                                                                     jint bits,
                                                                                                     jobject jpubkey) {
    jint len = env->GetArrayLength(a_lst);
    jint len2 = env->GetArrayLength(b_lst);
    jint len3 = env->GetArrayLength(out_lst);
    assert(len == len2 && len == len3 && len > 0);

    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    std::vector<cypher_text_t> a_text_lst, b_text_lst, res_lst;
    a_text_lst.reserve(len);
    b_text_lst.reserve(len);

    // FIXME: 三个循环并在一起是否会效率更高?
    for (int i = 0; i < len; i++) {
        cypher_text_t tmp_a(bits), tmp_b(bits);
        __get_arrOfSignedbyteArray_helper__(env, a_lst, tmp_a.c, tmp_a.is_neg, i);
        __get_arrOfSignedbyteArray_helper__(env, b_lst, tmp_b.c, tmp_b.is_neg, i);
        a_text_lst.push_back(tmp_a);
        b_text_lst.push_back(tmp_b);
    }
    for (int i = 0; i < len; i++) {
        cypher_text_t tmp(bits);
        c_api::add(tmp, a_text_lst[i], b_text_lst[i], pub_key);
        res_lst.push_back(tmp);
    }

    for (int i = 0; i < len; i++) {
        // cout << "in add_1vec, ii = " << i << endl;
        __fill_arrOfByteArray_helper__(env, out_lst, res_lst[i].c, res_lst[i].is_neg, i);
    }
}

JNIEXPORT void JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1sum_1vec_1_1(JNIEnv *env,
                                                                                                     jclass cls,
                                                                                                     jobject out,
                                                                                                     jobjectArray in_lst,
                                                                                                     jint bits,
                                                                                                     jobject jpubkey) {
    // cout << "in sum_vec" << endl;
    jint len = env->GetArrayLength(in_lst);
    assert(len >= 1);

    dist_paillier_pubkey_t pub_key = __pubkey_fromJpubkey_helper__(env, jpubkey, cls);
    pub_key.bits = bits;
    // std::cout << "pub_key.n = " << pub_key.n << endl;

    std::vector<cypher_text_t> in_text_lst;
    in_text_lst.reserve(len);

    cypher_text_t res(bits);
    c_api::enc(res, 0, pub_key);

    for (int i = 0; i < len; i++) {
        cypher_text_t tmp(bits);
        __get_arrOfSignedbyteArray_helper__(env, in_lst, tmp.c, tmp.is_neg, i);
        in_text_lst.push_back(tmp);
    }
    for (int i = 0; i < len; i++) {
        c_api::add(res, in_text_lst[i], res, pub_key);

        ZZ tmp1, tmp2;
        distributed_paillier::char_list_2_ZZ(tmp1, res.c, res.is_neg);
        distributed_paillier::char_list_2_ZZ(tmp2, in_text_lst[i].c, in_text_lst[i].is_neg);
        // cout << "in sum_vec, i = "  << i
        //      << "res = " << tmp1
        //      << "in_text_lst[i] = " << tmp2
        //      << endl;
    }
    __fill_SByteArray_helper__(env, out, res.c, res.is_neg);
}

JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_genePiQiShares(JNIEnv *env,
                                                                                                   jclass cls,
                                                                                                   jobjectArray piOut,
                                                                                                   jobject pi, 
                                                                                                   jobjectArray qiOut,
                                                                                                   jobject qi, 
                                                                                                   jint bit_len,
                                                                                                   jint t,
                                                                                                   jint num_party,
                                                                                                   jint partyID,
                                                                                                   jobject P_char) {
    jint len = num_party;
    ZZ P = signedbyteArray_2_ZZ(env, P_char);
    // distributed_paillier key_generator = distributed_paillier(P, num_party, t, bit_len / 2);
    std::vector<NTL::ZZ> piqi;
    if (partyID == 1) {
        piqi = distributed_paillier::gene_local_piqi_4_first_party(bit_len, num_party);
    } else {
        piqi = distributed_paillier::gene_local_piqi_4_other_party(bit_len, num_party);
    }

    // debug
    // piqi[0] = 1;
    // piqi[1] = 2;
    // cout << "pi, qi = ";
    // __print_zz_vec__(piqi);

    ShamirSecretSharingScheme scheme = ShamirSecretSharingScheme(P, num_party, t);
    ShamirShares share_p = scheme.share_secret(piqi[0]);
    ShamirShares share_q = scheme.share_secret(piqi[1]);

    assert(share_p.shares.length() == num_party);
    assert(share_q.shares.length() == num_party);

    // fill results in java side vectors
    fill_SByteArray_with_zz(env, pi, piqi[0]);
    fill_SByteArray_with_zz(env, qi, piqi[1]);
    for (int i = 0; i < len; i++) {
        // debug
        // cout << "pi(" << i + 1 << ") =  " << share_p.shares.get(i).b << "\t";
        // cout << "qi(" << i + 1 << ") =  " << share_q.shares.get(i).b << "\t";


        fill_arrOfSByteArray_with_zz(env, piOut, share_p.shares.get(i).b, i);
        fill_arrOfSByteArray_with_zz(env, qiOut, share_q.shares.get(i).b, i);
    }
    // cout << "\n" << endl;
}

ZZ get_sumP_share(const std::vector<ZZ> &reorganized_pb) {
    ZZ tmp = ZZ(0);
    for (const auto &pib: reorganized_pb) {
        tmp += pib;
    }

    // TODO: this check is done on Java Side
//    for(const auto & pia: pa_lst) {
//        assert(i == pia);
//    }

    return tmp;
}

JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneNShares(JNIEnv *env,
                                                                                                jclass cls,
                                                                                                jobjectArray allPiShares,
                                                                                                jobjectArray allQiShares,
                                                                                                jobject N_share_out,
                                                                                                jint bit_len,
                                                                                                jint t,
                                                                                                jint num_party,
                                                                                                jobject P_char) {
    jint pi_share_list_size = env->GetArrayLength(allPiShares);
    jint qi_share_list_size = env->GetArrayLength(allQiShares);

    assert((pi_share_list_size == num_party) && (qi_share_list_size == num_party));

    std::vector<ZZ> reorganized_pb, reorganized_qb;
    reorganized_pb.reserve(pi_share_list_size);
    reorganized_qb.reserve(qi_share_list_size);

    for (int i = 0; i < num_party; i++) {
        NTL::ZZ pi_share_zz, qi_share_zz;
        CHAR_LIST in_pi_shares, in_qi_shares;
        bool in_pi_isNeg, in_qi_isNeg;
        __get_arrOfSignedbyteArray_helper__(env, allPiShares, in_pi_shares, in_pi_isNeg, i);
        __get_arrOfSignedbyteArray_helper__(env, allQiShares, in_qi_shares, in_qi_isNeg, i);
        distributed_paillier::char_list_2_ZZ(pi_share_zz, in_pi_shares, in_pi_isNeg);
        distributed_paillier::char_list_2_ZZ(qi_share_zz, in_qi_shares, in_qi_isNeg);
        reorganized_pb.push_back(pi_share_zz);
        reorganized_qb.push_back(qi_share_zz);

        // debug
        // cout << "pi(" << i + 1 << ") =  " << pi_share_zz << "\t";
        // cout << "qi(" << i + 1 << ") =  " << qi_share_zz << "\t";
    }
    // cout << endl;

    ZZ P = signedbyteArray_2_ZZ(env, P_char);
    ZZ N_share_part = NTL::MulMod(get_sumP_share(reorganized_qb), get_sumP_share(reorganized_pb), P);
    fill_SByteArray_with_zz(env, N_share_out, N_share_part);

    // cout << "N_share = " << N_share_part
    //      << " get_sumP_share(reorganized_qb) = " << get_sumP_share(reorganized_qb)
    //      << " get_sumP_share(reorganized_pb) = " << get_sumP_share(reorganized_pb)
    //      << endl;

//    CHAR_LIST N_share_part_char;
//    bool N_share_part_isNeg;
//    distributed_paillier::ZZ_2_byte(N_share_part_char, N_share_part, N_share_part_isNeg);
//    __fill_SByteArray_helper__(env, N_share_out, N_share_part_char, N_share_part_isNeg);
}

JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_revealN(JNIEnv *env,
                                                                                            jclass cls,
                                                                                            jobjectArray allNShares,
                                                                                            jobject N_out,
                                                                                            jint bit_len,
                                                                                            jint t,
                                                                                            jint num_party,
                                                                                            jobject P_char) {
    jint allNShares_size = env->GetArrayLength(allNShares);

    assert((allNShares_size == num_party));

    std::vector<ZZ> N_shares_b;
    NTL::Vec<NTL::Pair<long, ZZ>> N_share;
    N_shares_b.reserve(allNShares_size);
    for (int i = 0; i < allNShares_size; i++) {
        NTL::ZZ N_share_zz;
        CHAR_LIST N_share_char;
        bool N_share_isNeg, in_qi_isNeg;
        __get_arrOfSignedbyteArray_helper__(env, allNShares, N_share_char, N_share_isNeg, i);
        distributed_paillier::char_list_2_ZZ(N_share_zz, N_share_char, N_share_isNeg);

        // note that i+1 should be consistent with the real A value of the reconstruction points.
        N_share.append(NTL::Pair<long, ZZ>(i + 1, N_share_zz));

        // cout << "N_share = " << N_share_zz << "\t";
    }
    // cout << endl;

    ZZ P = signedbyteArray_2_ZZ(env, P_char);

    // NOTE: degree of N's shares should be t*2 !
    ZZ N_zz = ShamirShares(N_share, P, num_party, t * 2).reconstruct_secret(P);

    CHAR_LIST N_char;
    bool N_isNeg;
    distributed_paillier::ZZ_2_byte(N_char, N_zz, N_isNeg);
    __fill_SByteArray_helper__(env, N_out, N_char, N_isNeg);

    // cout << "N = " << N_zz << endl;
}

JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getRand4Biprimetest(JNIEnv *env,
                                                                                                        jclass cls,
                                                                                                        jobject out,
                                                                                                        jobject N_char) {
    ZZ N = signedbyteArray_2_ZZ(env, N_char);
    // cout << "g _init = " << N << endl;
    fill_SByteArray_with_zz(env, out, distributed_paillier::get_rand_4_biprimetest(N));
}

/*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    biPrimeTestStage1
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;I)V
 */
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_biPrimeTestStage1(JNIEnv *env,
                                                                                                      jclass cls,
                                                                                                      jobject N_char,
                                                                                                      jobject pi_char,
                                                                                                      jobject qi_char,
                                                                                                      jobject g_char,
                                                                                                      jobject v_char,
                                                                                                      jint partyId) {
    ZZ N = signedbyteArray_2_ZZ(env, N_char);
    ZZ pi = signedbyteArray_2_ZZ(env, pi_char);
    ZZ qi = signedbyteArray_2_ZZ(env, qi_char);
    ZZ g = signedbyteArray_2_ZZ(env, g_char);
    ZZ v = distributed_paillier::biprime_test_step1(partyId, N, pi, qi, g);
    // cout << "v = " << v << endl;
    fill_SByteArray_with_zz(env, v_char, v);
}

/*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    biPrimeTestStage2
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IJ)V
 */
JNIEXPORT jlong  JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_biPrimeTestStage2(JNIEnv *env,
                                                                                                      jclass cls,
                                                                                                      jobjectArray other_v_char,
                                                                                                      jobject v_char,
                                                                                                      jobject N_char,
                                                                                                      jlong num_party) {
    ZZ N = signedbyteArray_2_ZZ(env, N_char);
    std::vector<ZZ> other_v = arrOfSignedbyteArray_2_zzList(env, other_v_char);
    ZZ v = signedbyteArray_2_ZZ(env, v_char);
    if(distributed_paillier::biprime_test_step2(other_v, v, N, num_party)){
        return 1;
    } else {
        return 0;
    }
}


/*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    geneLambdaBetaShares
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;III)V
 */
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneLambdaBetaShares(JNIEnv *env,
                                                                                                         jclass cls,
                                                                                                         jobjectArray lambdaiOut,
                                                                                                         jobjectArray betaiOut,
                                                                                                         jobject N_char,
                                                                                                         jobject pi_char,
                                                                                                         jobject qi_char,
                                                                                                         jint party_id,
                                                                                                         jint bit_len,
                                                                                                         jint t,
                                                                                                         jint num_party,
                                                                                                         jobject P_char) {
    ZZ lambda, beta;
    ZZ pi = signedbyteArray_2_ZZ(env, pi_char);
    ZZ qi = signedbyteArray_2_ZZ(env, qi_char);
    ZZ N = signedbyteArray_2_ZZ(env, N_char);
    if (party_id == 1) {
        lambda = N - pi - qi + 1;
    } else {
        lambda = -pi - qi;
    }
    do {
//        beta = NTL::RandomBnd(N);
        beta = 1; // for debug
    } while (beta == 0);
    ZZ P = signedbyteArray_2_ZZ(env, P_char);

    // debug
    // cout << "geneLambdaBetaShares\t"
    //      << "P = " << P << "\t"
    //      << "lambda = " << lambda << "\t"
    //      << "pi = " << pi << "\t"
    //      << "qi = " << qi << "\t"
    //      << endl;

    ShamirSecretSharingSchemeInteger scheme = ShamirSecretSharingSchemeInteger(P, num_party, t);
    ShamirShares_integer share_lambda = scheme.share_secret(lambda);
    ShamirShares_integer share_beta = scheme.share_secret(beta);
    assert(share_lambda.shares.length() == num_party);
    assert(share_beta.shares.length() == num_party);

    // fill results in java side vectors
    for (int i = 0; i < num_party; i++) {
        assert(share_lambda.shares.get(i).a == i + 1);
        assert(share_beta.shares.get(i).a == i + 1);

        // cout << "share_lambda = " << share_lambda.shares.get(i).b << "\t"
        //      << "share_beta = " << share_beta.shares.get(i).b << "\t";

        fill_arrOfSByteArray_with_zz(env, lambdaiOut, share_lambda.shares.get(i).b, i);
        fill_arrOfSByteArray_with_zz(env, betaiOut, share_beta.shares.get(i).b, i);
    }
    // cout << "\n" << endl;
}

/*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    geneLambdaTimesBetaShares
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;III)V
 */
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneLambdaTimesBetaShares(
        JNIEnv *env,
        jclass cls,
        jobjectArray lambdaSumShares,
        jintArray lambda_a,
        jobjectArray betaSumShares,
        jintArray beta_a,
        jobject N_char,
        jobject lambdaTimesBetaShareOut,
        jobject thetaOut,
        jint num_party) {

    // get data from JAVA
    std::vector<ZZ> in_lst_lambda_b = arrOfSignedbyteArray_2_zzList(env, lambdaSumShares);
    std::vector<ZZ> in_lst_beta_b = arrOfSignedbyteArray_2_zzList(env, betaSumShares);

    jint *body_lambda_a = env->GetIntArrayElements(lambda_a, 0);
    jint *body_beta_a = env->GetIntArrayElements(beta_a, 0);
    assert((in_lst_lambda_b.size() == in_lst_beta_b.size())
           && (in_lst_beta_b.size() == env->GetArrayLength(lambda_a))
           && (in_lst_beta_b.size() == env->GetArrayLength(beta_a))
           && (in_lst_beta_b.size() == num_party));

    // prepare args
    NTL::Vec<NTL::Pair<long, NTL::ZZ>> in_lst_lambda, in_lst_beta;
    for (int i = 0; i < in_lst_lambda_b.size(); i++) {
        in_lst_lambda.append(NTL::Pair<long, ZZ>(body_lambda_a[i], in_lst_lambda_b[i]));
        in_lst_beta.append(NTL::Pair<long, ZZ>(body_beta_a[i], in_lst_beta_b[i]));

        // cout << "body_lambda_a[i] = " << body_lambda_a[i] << " in_lst_lambda_b[i] = " << in_lst_lambda_b[i] << "\n"
        //      << "body_beta_a[i] = " << body_beta_a[i] << " in_lst_beta_b[i] = " << in_lst_beta_b[i]
        //      << endl;
    }
    ZZ N = signedbyteArray_2_ZZ(env, N_char);
    long scalar = NTL::to_long(factorial_small(num_party));
    std::vector<ZZ> ret = distributed_paillier::compute_lambda_times_beta_share(in_lst_lambda,
                                                                                in_lst_beta,
                                                                                N,
                                                                                scalar);
    // cout << "lambdaTimesBetaShareOut, thetaOut =  " << "\t";
    // __print_zz_vec__(ret);

    // write back
    fill_SByteArray_with_zz(env, lambdaTimesBetaShareOut, ret[0]);
    fill_SByteArray_with_zz(env, thetaOut, ret[1]);
}

/*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    revealThetaGeneKeys
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;III)V
 */
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_revealThetaGeneKeys(JNIEnv *env,
                                                                                                        jclass cls,
                                                                                                        jobjectArray allThetaShares,
                                                                                                        jintArray theta_a,
                                                                                                        jobject thetaOut,
                                                                                                        jobject N_char,
                                                                                                        jint t,
                                                                                                        jint n) {
    // get data from JAVA
    ZZ theta_invmod;
    ZZ N = signedbyteArray_2_ZZ(env, N_char);
    jint *body_theta_a = env->GetIntArrayElements(theta_a, 0);
    jint array_size = env->GetArrayLength(theta_a);
    std::vector<NTL::ZZ> theta_shares_raw = arrOfSignedbyteArray_2_zzList(env, allThetaShares);
    std::vector<NTL::ZZ> theta_shares;
    theta_shares.resize(array_size);
    assert(theta_shares.size() == array_size);

    // prepare args
    for (int i = 0; i < array_size; i++) {
        assert((body_theta_a[i] > 0) && (body_theta_a[i] <= array_size));
        theta_shares[body_theta_a[i] - 1] = theta_shares_raw[i];
    }

    ZZ theta = distributed_paillier::reveal_theta(theta_shares, N, t, n);
    NTL::InvMod(theta_invmod, theta, N);

    // write back
    fill_SByteArray_with_zz(env, thetaOut, theta_invmod);
}

JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargePrime(JNIEnv *env,
                                                                                                  jclass cls,
                                                                                                  jobject out,
                                                                                                  jint in) {
    ZZ P;
    NTL::GenPrime(P, in);
    fill_SByteArray_with_zz(env, out, P);
    // cout << "P = " << P << endl;
}

JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargeComposite(JNIEnv *env,
                                                                                                      jclass cls,
                                                                                                      jobject out,
                                                                                                      jint in) {
    fill_SByteArray_with_zz(env, out, NTL::RandomLen_ZZ(in / 2) * NTL::RandomLen_ZZ(in / 2));
}

// test function
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_testIO
        (JNIEnv *env, jclass cls, jobject out, jobject in) {
    ZZ data = signedbyteArray_2_ZZ(env, in);
    fill_SByteArray_with_zz(env, out, data);
}

// test function
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_testIOVec
        (JNIEnv *env, jclass cls, jobjectArray out, jobjectArray in) {
    std::vector<ZZ> data = arrOfSignedbyteArray_2_zzList(env, in);
    int i = 0;
    for (const auto &elem : data) {
        fill_arrOfSByteArray_with_zz(env, out, elem, i);
        i++;
    }
}

// test function
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneNFromPQProduct
        (JNIEnv *env, jclass cls, jint num_party, jobject n_out, jobjectArray piOut, jobjectArray qiOut, jobject pIn,
         jobject qIn) {
    NTL::Vec<NTL::ZZ> pi_lst, qi_lst;
    ZZ p = signedbyteArray_2_ZZ(env, pIn);
    ZZ q = signedbyteArray_2_ZZ(env, qIn);
    ZZ p_sum(0), q_sum(0);
    for (long i = 1; i < num_party; i++) {
        pi_lst.append(geneNoneZeroRandBnd(p / num_party));
        qi_lst.append(geneNoneZeroRandBnd(q / num_party));
        p_sum += pi_lst[i - 1];
        q_sum += qi_lst[i - 1];
    }
    pi_lst.append(p - p_sum);
    qi_lst.append(q - q_sum);

    // check qi, pi correctness
    p_sum = 0;
    q_sum = 0;

    for (int i = 0; i < num_party; i++) {
        fill_arrOfSByteArray_with_zz(env, piOut, pi_lst[i], i);
        fill_arrOfSByteArray_with_zz(env, qiOut, qi_lst[i], i);
    }
    fill_SByteArray_with_zz(env, n_out, p * q);

    for (long i = 0; i < num_party; i++) {
        p_sum += pi_lst[i];
        q_sum += qi_lst[i];
        // cout << "pi = " << pi_lst[i] << " qi = " << qi_lst[i] << endl;
        assert((pi_lst[i] > 0) && (qi_lst[i] > 0));
    }
    // cout << "p*q = " << p * q << endl;
    assert(q_sum == q && p_sum == p);
}

// test function
// beta is fixed to be 1
JNIEXPORT jboolean JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_checkCorrectness
        (JNIEnv *env, jclass cls, jobject pIn, jobject qIn, jobject N_char,
         jobjectArray lambdaTimesBetaShareChar, jobject thetaChar, jint n, jint t) {
    ZZ p = signedbyteArray_2_ZZ(env, pIn);
    ZZ q = signedbyteArray_2_ZZ(env, qIn);
    ZZ N = signedbyteArray_2_ZZ(env, N_char);
    ZZ nfct = factorial_small(n);
    ZZ theta = signedbyteArray_2_ZZ(env, thetaChar);
    vector<ZZ> lambdaTimesBetaShareVector = arrOfSignedbyteArray_2_zzList(env, lambdaTimesBetaShareChar);
    NTL::Vec<ZZ> lambdaTimesBetaShare;
    for (const auto &elem: lambdaTimesBetaShareVector) {
        lambdaTimesBetaShare.append(elem);
    }

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
        NTL::mul(tmp, li, lambdaTimesBetaShare[i - 1]);
        NTL::add(recon_lamb_times_beta, recon_lamb_times_beta, tmp); // this reconstruction does not mod N
    }

    NTL::InvMod(theta, theta, N);

    // cout << "\nreconLambdaTimesBeta = " << recon_lamb_times_beta << "\t"
    //      << "LambdaTimesBeta = " << (N - p - q + 1) * 3 * nfct * sqr(nfct) << "\n"
    //      << "theta = " << theta << "\t"
    //      << "reconLambdaTimesBeta % N = " << recon_lamb_times_beta % N << "\n"
    //      << "N = " << N
    //      << endl;
    return ((N - p - q + 1) * 3 * nfct * sqr(nfct) == recon_lamb_times_beta) && (theta == recon_lamb_times_beta % N);
}

// test function
//generate p,a without pum
JNIEXPORT void  JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargePrime4Test
        (JNIEnv *env, jclass cls, jobject pOut, jobjectArray piOut, jint bitLen, jint numParty) {
    std::vector<ZZ> pi_list;
    pi_list.resize(numParty);
    ZZ p_sum(0);
    while((p_sum==0) || (PrimeTest(p_sum, 40)==0)) {
        p_sum = ZZ(0);
        for (int i = 0; i < numParty; i++) {
            NTL::ZZ tmp;
            if (i == 0) {
                do {
                    tmp = NTL::RandomBits_ZZ(bitLen);
                } while (tmp == 0 || ( (tmp%4)!=3) );
            } else {
                do {
                    tmp = NTL::RandomBits_ZZ(bitLen);
                } while (tmp == 0 || ( (tmp%4)!=0) );
            }
            pi_list[i] = tmp;
            p_sum += tmp;
        }
    }
    fill_SByteArray_with_zz(env, pOut, p_sum);
    for(int i = 0; i < numParty; i++) {
        fill_arrOfSByteArray_with_zz(env, piOut, pi_list[i], i) ;
    }
    // cout << "p_sum = " << p_sum << endl;
    // cout << "pi List = ";
    // __print_zz_vec__(pi_list);
}


// test function
JNIEXPORT void  JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargeComposite4Test
        (JNIEnv *env, jclass cls, jobject pOut, jobjectArray piOut, jint bitLen, jint numParty) {
    std::vector<ZZ> pi_list;
    pi_list.resize(numParty);
    ZZ p_sum(0);
    while((p_sum==0) || (PrimeTest(p_sum, 40)==1)) {
        p_sum = ZZ(0);
        for (int i = 0; i < numParty; i++) {
            NTL::ZZ tmp;
            if (i == 0) {
                do {
                    tmp = NTL::RandomBits_ZZ(bitLen);
                } while ((tmp == 0 )|| ( (tmp%4)!=3) );
            } else {
                do {
                    tmp = NTL::RandomBits_ZZ(bitLen);
                } while ((tmp == 0) || ( (tmp%4)!=0) );
            }
            pi_list[i] = tmp;
            p_sum += tmp;
        }
    }
    fill_SByteArray_with_zz(env, pOut, p_sum);
    for(int i = 0; i < numParty; i++) {
        fill_arrOfSByteArray_with_zz(env, piOut, pi_list[i], i) ;
    }
    // cout << "pi = " << p_sum << "   ========" << endl;
}

JNIEXPORT void  JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_zzaTimeszzb
        (JNIEnv *env, jclass cls, jobject a, jobject b, jobject out) {
    ZZ N = signedbyteArray_2_ZZ(env, a)* signedbyteArray_2_ZZ(env, b);
    // cout << "N = " << N << endl;
    fill_SByteArray_with_zz(env,out,N);
}

JNIEXPORT jlong  JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_checkPiSumPrime
(JNIEnv *env, jclass cls, jobjectArray listIn) {
    std::vector<ZZ> listIn_zz = arrOfSignedbyteArray_2_zzList(env, listIn);
    ZZ pi_sum(0);
    for(auto value : listIn_zz) {
        pi_sum += value;
    }

    // cout << "sum_value = " << pi_sum;
    return PrimeTest(pi_sum, 80);
}



