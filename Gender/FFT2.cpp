#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <tfhe/tlwe.h>
#include <time.h>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

#define K  4
#define M  200
#define M_SLOT 200
#define ABSCALE 50

// Generate IntPolynomial of coefficients of result[N-1-i] = 2^K * sin(2*pi*k*i/N) *inverse packing
IntPolynomial *generate_sinpoly(const TLweParams* tlwe_param, int fft_scale, int k){
    const int32_t N = tlwe_param->N;
    int32_t SCALE = 1<<fft_scale;
    IntPolynomial* result = new_IntPolynomial(N);
    for(int i=0; i<N; i++){
        result->coefs[N-1-i] = (int) round(sin(2*M_PI*k*i/N)*SCALE);
    }
    return result;
}

// Generate IntPolynomial of coefficients of result[N-1-i] = 2^K * cos(2*pi*k*i/N) *inverse packing
IntPolynomial *generate_cospoly(const TLweParams* tlwe_param, int fft_scale, int k){
    const int32_t N = tlwe_param->N;
    int32_t SCALE = 1<<fft_scale;
    IntPolynomial* result = new_IntPolynomial(N);
    for(int i=0; i<N; i++){
        result->coefs[N-1-i] = (int) round(cos(2*M_PI*k*i/N)*SCALE);
    }
    return result;
}

// calculate IntPolynomial array of size 2, whose coefs of sinpoly, cospoly [N-1-i] = 2^K * sin/cos(2*pi*k*i/N) *inverse packing 
IntPolynomial *generate_fftpoly(const TLweParams* tlwe_param, int fft_scale, int k){
    const int32_t N = tlwe_param->N;
    IntPolynomial* fftpoly = new_IntPolynomial_array(2, N);
    int32_t SCALE = 1<<fft_scale;
    for(int i=0; i<N; i++){
        fftpoly->coefs[N-1-i] = (int) round(cos(2*M_PI*k*i/N)*SCALE);
        (fftpoly+1)->coefs[N-1-i] = (int) round(sin(2*M_PI*k*i/N)*SCALE);
        // printf("%d ",fftpoly->coefs[N-1-i]);
        // printf("%d\n",(fftpoly+1)->coefs[N-1-i]);
            // *(fftpoly+1)->coefs[0];
    }
    return fftpoly;
}

// return 2 N-LWE ciphertexts, encrypting k-th bin of DCT, DST result respectively 
LweSample* calculate_bin(TLweSample* ciphertext, const TLweParams* rparams, int k){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* result = new_LweSample_array(2, extract_params);
    const TLweSample* temp = new_TLweSample(rparams);
    IntPolynomial* fftpoly = generate_fftpoly(rparams, K, k);
    for(int i=0; i<2; i++){
        torusPolynomialMultFFT(temp->a, (fftpoly+i), ciphertext->a);
        torusPolynomialMultFFT(temp->b, (fftpoly+i), ciphertext->b);
        tLweExtractLweSampleIndex((result+i), temp, N-1, extract_params, rparams);
    }
    return result;
}

LweSample* calculate_bin_KS(TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample_array(2, extract_params);
    LweSample* result = new_LweSample_array(2, lweparam);
    const TLweSample* temp = new_TLweSample(rparams);
    IntPolynomial* fftpoly = generate_fftpoly(rparams, K, k);
    for(int i=0; i<2; i++){
        torusPolynomialMultFFT(temp->a, (fftpoly+i), ciphertext->a);
        torusPolynomialMultFFT(temp->b, (fftpoly+i), ciphertext->b);
        tLweExtractLweSampleIndex((Nresult+i), temp, N-1, extract_params, rparams);
        lweKeySwitch((result+i),bk->bk->ks, (Nresult+i));
    }
    return result;
}

TorusPolynomial* generate_abs(int32_t N, int32_t alph){
    int32_t halfN = N/2;
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    Torus32 mu = modSwitchToTorus32(1, M_SLOT);
    for(int i=0; i<halfN; i++){
        testvect->coefsT[i] = mu*round(i/alph + 1);
        testvect->coefsT[N-1-i] = -mu*round(i/alph + 1);  
    }
    return testvect;
}

TorusPolynomial* generate_logabs(int32_t N){
    int32_t halfN = N/2;
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    Torus32 mu = modSwitchToTorus32(1, M_SLOT);
    for(int i=0; i<halfN; i++){
        int32_t slot = floor(log2(i+1)/2+1);
        Torus32 mu_slot = modSwitchToTorus32(slot, M_SLOT);
        testvect->coefsT[i] = mu_slot;
        testvect->coefsT[N-1-i] = -mu_slot;  
    }
    return testvect;
}

void bootstrap_woKS_FFT(LweSample *result, const LweBootstrappingKeyFFT *bk, const LweSample *x) {

    const TGswParams *bk_params = bk->bk_params;
    const TLweParams *accum_params = bk->accum_params;
    const LweParams *in_params = bk->in_out_params;
    const int32_t N = accum_params->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = in_params->n;

    TorusPolynomial *testvect = generate_abs(N, ABSCALE);
    // TorusPolynomial *testvect = generate_logabs(N);
    int32_t *bara = new int32_t[N];

    // Modulus switching
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    for (int32_t i = 0; i < n; i++) {
        bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
    }
    // Bootstrapping rotation and extraction
    tfhe_blindRotateAndExtract_FFT(result, testvect, bk->bkFFT, barb, bara, n, bk_params);
}

void bootstrap_FFT(LweSample *result, const LweBootstrappingKeyFFT *bk, const LweSample *x) {

    LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

    bootstrap_woKS_FFT(u, bk, x);
    // Key switching
    lweKeySwitch(result, bk->ks, u);
}

using namespace std;

int main() {
    FILE* cloud_key = fopen("cloud.key", "rb");
    TFheGateBootstrappingCloudKeySet* bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    fclose(cloud_key);

    const TFheGateBootstrappingParameterSet* params = bk->params;
    const LweParams* lweparams = bk->params->in_out_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;
    const int N = tlwe_params->N;
    const LweParams* exparams = &tlwe_params->extracted_lweparams;
    TLweSample* ciphertext = new_TLweSample(tlwe_params);

    FILE* cloud_data = fopen("cloud.data", "rb");
    import_tlweSample_fromFile(cloud_data, ciphertext, tlwe_params);
    fclose(cloud_data);

    float time = -clock();
    LweSample *DCT = new_LweSample_array(M, lweparams);
    LweSample *DST = new_LweSample_array(M, lweparams);
    for(int i=0; i<M; i++){
        LweSample * temp = calculate_bin_KS(ciphertext, lweparams ,tlwe_params, i+100, bk);
        lweCopy(&DCT[i], &temp[0], lweparams);
        lweCopy(&DST[i], &temp[1], lweparams);
    }
    delete_TLweSample(ciphertext); 
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    cout << "-----------------  We calculated FFT for 100~300Hz bin in "<< time << " seconds!  ----------------\n" << endl;

    LweSample *Bootstrapped_DCT = new_LweSample_array(M, lweparams);
    LweSample *Bootstrapped_DST = new_LweSample_array(M, lweparams);
    time = -clock();

    for(int i=0; i<M; i++){
        bootstrap_FFT(&Bootstrapped_DCT[i], bk->bkFFT, &DCT[i]);
        bootstrap_FFT(&Bootstrapped_DST[i], bk->bkFFT, &DST[i]);
        lweAddTo(&Bootstrapped_DCT[i], &Bootstrapped_DST[i], lweparams);
        //bootstrap_FFT(&Bootstrapped_DCT[i], bk->bkFFT, &Bootstrapped_DCT[i]);
    }
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    cout << "---------  We bootstrapped "<< M*2 << " LWE samples into Z/"<< M_SLOT <<"Z, and it took " << time << " seconds!  --------\n" << endl;
    cout << " Hello!\n" <<endl;

    FILE* secret_key = fopen("secret.key", "rb");
    TFheGateBootstrappingSecretKeySet* key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
    fclose(secret_key);

    const TFheGateBootstrappingParameterSet* Param = key->params;
    const LweKey* lwekey = new_LweKey(lweparams);
    lwekey = key->lwe_key;
    LweSample* res = new_LweSample_array(2, lweparams);
    // int32_t halfM = M/2;
    // for(int i=0; i<halfM; i++){
    //     lweAddTo(&res[0], &Bootstrapped_DCT[i], lweparams);
    //     lweAddTo(&res[1], &Bootstrapped_DCT[i+100], lweparams);
    // }

    int32_t decr[M];

    for(int i=0; i<M; i++){
        Torus32 dec = lwePhase(&Bootstrapped_DCT[i], lwekey);
        // printf("%d\n", dec);
        decr[i] = modSwitchFromTorus32(dec, M_SLOT);
        //decr[i] = dec;
    }
    ofstream myfile ("example.txt");
    for(int count=0; count<M; count++){
        myfile << decr[count] << " ";
    }
    myfile.close();


    delete_LweSample_array(M, DCT);
    delete_LweSample_array(M, DST);


    printf("Hello world!");

}