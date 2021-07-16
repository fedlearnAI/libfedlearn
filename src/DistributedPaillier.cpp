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
#include <omp.h>
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
template<typename T>
void __print_vec_as_int__(const std::vector<T> &in) {
    for (auto x : in) {
        std::cout << (int) x << " ";
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

/**
 * ============================
 * JNI layer functions
 * ============================
 */

JNIEXPORT jobjectArray JNICALL
JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1create_1share_1_1(JNIEnv *env,
                                                                                                          jclass cls,
                                                                                                          jbyteArray secret,
                                                                                                          jint, jint) {
    std::cout << "\nRuning from Cpp, not inplemented." << endl;
    assert(false); // FIXME, not inplemented.
}

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
    std::cout << "ret[0].n = " << ret[0].n << endl;

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
#pragma omp parallel for shared(c_text_lst)
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
                                                                                                   jobjectArray qiOut,
                                                                                                   jint bit_len,
                                                                                                   jint t,
                                                                                                   jint num_party) {
    jint len = num_party;
    distributed_paillier key_generator = distributed_paillier(P, num_party, t, bit_len);
    std::vector<NTL::ZZ> piqi = key_generator.gene_local_piqi();

    ShamirSecretSharingScheme scheme = ShamirSecretSharingScheme(P, num_party, t);
    ShamirShares share_p = scheme.share_secret(piqi[0]);
    ShamirShares share_q = scheme.share_secret(piqi[1]);

    assert(share_p.shares.length() == num_party);
    assert(share_q.shares.length() == num_party);

    // fill results in java side vectors
    for (int i = 0; i < len; i++) {
        CHAR_LIST pi_char, qi_char;
        bool pi_isNeg, qi_isNeg;
        assert(share_p.shares.get(i).a == i+1);
        assert(share_q.shares.get(i).a == i+1);
        distributed_paillier::ZZ_2_byte(pi_char, share_p.shares.get(i).b, pi_isNeg);
        distributed_paillier::ZZ_2_byte(qi_char, share_q.shares.get(i).b, qi_isNeg);
        __fill_arrOfByteArray_helper__(env, piOut, pi_char, pi_isNeg, i);
        __fill_arrOfByteArray_helper__(env, qiOut, qi_char, qi_isNeg, i);
    }
}

ZZ get_sumP_share(const std::vector<ZZ>& reorganized_pb) {
    ZZ tmp = ZZ(0);
    for(const auto& pib: reorganized_pb) {
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
                                                                                            jint num_party) {
    jint len = num_party;
    jint pi_share_list_size = env->GetArrayLength(allPiShares);
    jint qi_share_list_size = env->GetArrayLength(allQiShares);

    assert((pi_share_list_size == num_party ) && (qi_share_list_size == num_party ));

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
    }
    ZZ N_share_part = NTL::MulMod(get_sumP_share(reorganized_qb), get_sumP_share(reorganized_pb), P);

    CHAR_LIST N_share_part_char;
    bool N_share_part_isNeg;
    distributed_paillier::ZZ_2_byte(N_share_part_char, N_share_part, N_share_part_isNeg);
    __fill_SByteArray_helper__(env, N_share_out, N_share_part_char, N_share_part_isNeg);
}

JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_revealN(JNIEnv *env,
                                                                                            jclass cls,
                                                                                            jobjectArray allNShares,
                                                                                            jobject N_out,
                                                                                            jint bit_len,
                                                                                            jint t,
                                                                                            jint num_party) {
    jint len = num_party;
    jint allNShares_size = env->GetArrayLength(allNShares);

    assert((allNShares_size == num_party ));

    std::vector<ZZ> N_shares_b;
    NTL::Vec<NTL::Pair<long, ZZ>>  N_share;
    N_shares_b.reserve(allNShares_size);
    for (int i = 0; i < allNShares_size; i++) {
        NTL::ZZ N_share_zz;
        CHAR_LIST N_share_char;
        bool N_share_isNeg, in_qi_isNeg;
        __get_arrOfSignedbyteArray_helper__(env, allNShares, N_share_char, N_share_isNeg, i);
        distributed_paillier::char_list_2_ZZ(N_share_zz, N_share_char, N_share_isNeg);

        // note that i+1 should be consistent with the real A value of the reconstruction points.
        N_share.append( NTL::Pair<long, ZZ>(i+1,  N_share_zz));
    }
    ZZ N_zz = ShamirShares(N_share, P, num_party, t).reconstruct_secret(P);

    CHAR_LIST N_char;
    bool N_isNeg;
    distributed_paillier::ZZ_2_byte(N_char, N_zz, N_isNeg);
    __fill_SByteArray_helper__(env, N_out, N_char, N_isNeg);
}

/*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    biPrimeTestStage1
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;I)V
 */
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_biPrimeTestStage1(JNIEnv *env,
                                                                                                      jclass cls,
                                                                                                      jobject N,
                                                                                                      jobject pi,
                                                                                                      jobject qi,
                                                                                                      jobject g,
                                                                                                      jobject v,
                                                                                                      jint partyId) {
}

/*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    biPrimeTestStage2
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IJ)V
 */
JNIEXPORT void JNICALL
Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_biPrimeTestStage2 (JNIEnv *env,
                                                                                                       jclass cls,
                                                                                                       jobjectArray other_v,
                                                                                                       jobject v,
                                                                                                       jint out,
                                                                                                       jlong n){
}