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

#ifndef DIST_HELPER_H
#define DIST_HELPER_H

#include <vector>
#include <jni.h>
#include "distributed_paillier.h"

template<typename T>
void __print_vec_as_int__(const std::vector<T> &in);

jbyteArray as_jbyte_array(JNIEnv *env, const CHAR_LIST &in);

void as_uncharlist(CHAR_LIST &out, JNIEnv *env, jobject array);

jobjectArray as_2djbyte_array(JNIEnv *env,
                              LIST_OF_CHAR_LIST &in);

void as_listofcharlist(LIST_OF_CHAR_LIST &out,
                       JNIEnv *env,
                       jobjectArray array);

dist_paillier_privkey_t
__privkey_fromJprivkey_helper__(JNIEnv *env,
                                jobjectArray jPrivkey,
                                jclass cls);

dist_paillier_pubkey_t
__pubkey_fromJpubkey_helper__(JNIEnv *env,
                              jobject jpubkey,
                              jclass cls);

void __get_SignedbyteArray_helper__(JNIEnv *env,
                                    const jobject &SbyteArrayObj,
                                    CHAR_LIST &target_lst,
                                    bool &target_isNeg);

void __get_arrOfSignedbyteArray_helper__(JNIEnv *env,
                                         const jobjectArray &arrOfByteArrayObj,
                                         CHAR_LIST &target_lst,
                                         bool &target_isNeg,
                                         const int idx);

void __fill_arrOfByteArray_helper__(JNIEnv *env,
                                    jobjectArray arrOfByteArrayObj,
                                    const CHAR_LIST &src_lst,
                                    const bool &src_isNeg,
                                    const int &idx);

void __fill_SByteArray_helper__(JNIEnv *env,
                                jobject &sByteArrayObj,
                                const CHAR_LIST &src_lst,
                                const bool &src_isNeg);

#endif