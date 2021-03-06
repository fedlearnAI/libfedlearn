/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative */

#ifndef _Included_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
#define _Included_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
#ifdef __cplusplus
extern "C"
{
#endif
        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __create_share__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[BII)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1create_1share_1_1(JNIEnv *, jclass, jobjectArray, jbyteArray, jint, jint);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __generate_privpub_key__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIIZ)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1generate_1privpub_1key_1_1(JNIEnv *, jclass, jobjectArray, jint, jint, jint, jboolean);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __enc__
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;JILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1enc_1_1(JNIEnv *, jclass, jobject, jlong, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __enc_list__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[JILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1enc_1list_1_1(JNIEnv *, jclass, jobjectArray, jlongArray, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __partial_dec__
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIIJ[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1partial_1dec_1_1(JNIEnv *, jclass, jobject, jobject, jint, jint, jint, jlong, jobjectArray);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __partial_dec_lst__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIIJ[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1partial_1dec_1lst_1_1(JNIEnv *, jclass, jobjectArray, jobjectArray, jint, jint, jint, jlong, jobjectArray);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __final_dec__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIJJ[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)J
 */
        JNIEXPORT jlong JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1final_1dec_1_1(JNIEnv *, jclass, jobjectArray, jobject, jint, jint, jlong, jlong, jobjectArray);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __add__
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;ILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1add_1_1(JNIEnv *, jclass, jobject, jobject, jobject, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __add_vec__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;ILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1add_1vec_1_1(JNIEnv *, jclass, jobjectArray, jobjectArray, jobjectArray, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __sum_vec__
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;ILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1sum_1vec_1_1(JNIEnv *, jclass, jobject, jobjectArray, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __mul__
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;JILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1mul_1_1(JNIEnv *, jclass, jobject, jobject, jlong, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __mul_elementWise__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[JILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1mul_1elementWise_1_1(JNIEnv *, jclass, jobjectArray, jobjectArray, jlongArray, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __div__
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;JILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1div_1_1(JNIEnv *, jclass, jobject, jobject, jlong, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    __div_elementWise__
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[JILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative__1_1div_1elementWise_1_1(JNIEnv *, jclass, jobjectArray, jobjectArray, jlongArray, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    genePiQiShares
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIIILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_genePiQiShares(JNIEnv *env,
                                                                                                                                  jclass cls,
                                                                                                                                  jobjectArray piOut,
                                                                                                                                  jobject pi,
                                                                                                                                  jobjectArray qiOut,
                                                                                                                                  jobject qi,
                                                                                                                                  jint bit_len,
                                                                                                                                  jint t,
                                                                                                                                  jint num_party,
                                                                                                                                  jint partyID,
                                                                                                                                  jobject P_char);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    geneNShares
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneNShares(JNIEnv *, jclass, jobjectArray, jobjectArray, jobject, jint, jint, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    revealN
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT jlong JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_revealN(JNIEnv *, jclass, jobjectArray, jobject, jint, jint, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    biPrimeTestStage1
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;I)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_biPrimeTestStage1(JNIEnv *, jclass, jobject, jobject, jobject, jobject, jobject, jint);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    biPrimeTestStage2
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;J)J
 */
        JNIEXPORT jlong JNICALL
        Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_biPrimeTestStage2(JNIEnv *env,
                                                                                                              jclass cls,
                                                                                                              jobjectArray other_v_char,
                                                                                                              jobject v_char,
                                                                                                              jobject N_char,
                                                                                                              jlong num_party);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    geneLambdaBetaShares
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;IIIILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneLambdaBetaShares(JNIEnv *, jclass, jobjectArray, jobjectArray, jobject, jobject, jobject, jint, jint, jint, jint, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    geneLambdaTimesBetaShares
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[I[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[ILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;I)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneLambdaTimesBetaShares(JNIEnv *, jclass, jobjectArray, jintArray, jobjectArray, jintArray, jobject, jobject, jobject, jint);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    revealThetaGeneKeys
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[ILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;II)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_revealThetaGeneKeys(JNIEnv *, jclass, jobjectArray, jintArray, jobject, jobject, jint, jint);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    getLargePrime
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;I)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargePrime(JNIEnv *, jclass, jobject, jint);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    getLargeComposite
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;I)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargeComposite(JNIEnv *, jclass, jobject, jint);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    getRand4Biprimetest
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getRand4Biprimetest(JNIEnv *, jclass, jobject, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    testIO
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_testIO(JNIEnv *, jclass, jobject, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    testIOVec
 * Signature: ([Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_testIOVec(JNIEnv *, jclass, jobjectArray, jobjectArray);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    geneNFromPQProduct
 * Signature: (ILcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;)V
 */
        JNIEXPORT void JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_geneNFromPQProduct(JNIEnv *, jclass, jint, jobject, jobjectArray, jobjectArray, jobject, jobject);

        /*
 * Class:     com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative
 * Method:    checkCorrectness
 * Signature: (Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;[Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;Lcom/jdt/fedlearn/core/encryption/distributedPaillier/DistributedPaillierNative/signedByteArray;II)Z
 */
        JNIEXPORT jboolean JNICALL Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_checkCorrectness(JNIEnv *, jclass, jobject, jobject, jobject, jobjectArray, jobject, jint, jint);

        JNIEXPORT void JNICALL
        Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargePrime4Test(JNIEnv *env, jclass cls, jobject pOut, jobjectArray piOut, jint bitLen, jint numParty);

        JNIEXPORT void JNICALL
        Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_getLargeComposite4Test(JNIEnv *env, jclass cls, jobject pOut, jobjectArray piOut, jint bitLen, jint numParty);

        JNIEXPORT void JNICALL
        Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_zzaTimeszzb(JNIEnv *env, jclass cls, jobject a, jobject b, jobject out);

        JNIEXPORT jlong  JNICALL
        Java_com_jdt_fedlearn_core_encryption_distributedPaillier_DistributedPaillierNative_checkPiSumPrime(JNIEnv *env, jclass cls, jobjectArray listIn);

#ifdef __cplusplus
}
#endif
#endif
