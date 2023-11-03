#include <tfhe/lwe-functions.h>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

#define K  9
#define M  100
#define ABSCALE pow(2., 5)
#define M_SLOT 700
#define C_NUM 138


// Generate IntPolynomial of coefficients of result[N-1-i] = 2^K * sin(2*pi*k*i/N) *inverse packing
IntPolynomial *generate_sinpoly(const TLweParams* tlwe_param, int fft_scale, int k){
    const int32_t N = tlwe_param->N;
    int32_t SCALE = 1<<fft_scale;
    IntPolynomial* result = new_IntPolynomial(N);
    for(int i=0; i<N; i++){
        result->coefs[N-1-i] = (int) round(-sin(2*M_PI*k*i/N)*SCALE);
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
        (fftpoly+1)->coefs[N-1-i] = (int) round(-sin(2*M_PI*k*i/N)*SCALE);
    }
    return fftpoly;
}

// return 2 N-LWE ciphertexts, encrypting k-th bin of DCT, DST result respectively 
LweSample* calculate_bin(TLweSample* ciphertext, const TLweParams* rparams, int k){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* result = new_LweSample_array(2, extract_params);
    TLweSample* temp = new_TLweSample(rparams);
    IntPolynomial* fftpoly = generate_fftpoly(rparams, K, k);
    for(int i=0; i<2; i++){
        torusPolynomialMultFFT(temp->a, (fftpoly+i), ciphertext->a);
        torusPolynomialMultFFT(temp->b, (fftpoly+i), ciphertext->b);
        tLweExtractLweSampleIndex((result+i), temp, N-1, extract_params, rparams);
    }
    delete_IntPolynomial(fftpoly);
    delete_TLweSample(temp);
    return result;
}

LweSample* calculate_bin_KS(TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample_array(2, extract_params);
    LweSample* result = new_LweSample_array(2, lweparam);
    TLweSample* temp = new_TLweSample(rparams);
    IntPolynomial* fftpoly = generate_fftpoly(rparams, K, k);
    for(int i=0; i<2; i++){
        torusPolynomialMultFFT(temp->a, (fftpoly+i), ciphertext->a);
        torusPolynomialMultFFT(temp->b, (fftpoly+i), ciphertext->b);
        tLweExtractLweSampleIndex((Nresult+i), temp, N-1, extract_params, rparams);
        lweKeySwitch((result+i),bk->bk->ks, (Nresult+i));
        tLweClear(temp, rparams);
    }
    delete_LweSample_array(2, Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial_array(2, fftpoly);
    return result;
}

void calculate_bin_KS_v2(LweSample* result, TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample_array(2, extract_params);
    TLweSample* temp = new_TLweSample(rparams);
    IntPolynomial* fftpoly = generate_fftpoly(rparams, K, k);
    for(int i=0; i<2; i++){
        torusPolynomialMultFFT(temp->a, (fftpoly+i), ciphertext->a);
        torusPolynomialMultFFT(temp->b, (fftpoly+i), ciphertext->b);
        tLweExtractLweSampleIndex((Nresult+i), temp, N-1, extract_params, rparams);
        lweKeySwitch((result+i),bk->bk->ks, (Nresult+i));
        tLweClear(temp, rparams);
    }
    delete_LweSample_array(2, Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial_array(2, fftpoly);
}

LweSample* calculate_DCT_bin_KS(TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample(extract_params);
    LweSample* result = new_LweSample(lweparam);
    TLweSample* temp = new_TLweSample(rparams);

    IntPolynomial* cospoly = generate_cospoly(rparams, K, k);
    torusPolynomialMultFFT(temp->a, cospoly, ciphertext->a);
    torusPolynomialMultFFT(temp->b, cospoly, ciphertext->b);
    tLweExtractLweSampleIndex((Nresult), temp, N-1, extract_params, rparams);
    lweKeySwitch(result, bk->bk->ks, Nresult);

    delete_LweSample(Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial(cospoly);
    return result;
}

void calculate_DCT_bin_KS_v2(LweSample* result ,TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample(extract_params);
    // LweSample* result = new_LweSample(lweparam);
    TLweSample* temp = new_TLweSample(rparams);

    IntPolynomial* cospoly = generate_cospoly(rparams, K, k);
    torusPolynomialMultFFT(temp->a, cospoly, ciphertext->a);
    torusPolynomialMultFFT(temp->b, cospoly, ciphertext->b);
    tLweExtractLweSampleIndex((Nresult), temp, N-1, extract_params, rparams);
    lweKeySwitch(result, bk->bk->ks, Nresult);

    delete_LweSample(Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial(cospoly);
}

void DCT_bin_KS(LweSample* result ,TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample(extract_params);
    TLweSample* temp = new_TLweSample(rparams);

    IntPolynomial* cospoly = generate_cospoly(rparams, K, k);
    torusPolynomialMultFFT(temp->a, cospoly, ciphertext->a);
    torusPolynomialMultFFT(temp->b, cospoly, ciphertext->b);
    tLweExtractLweSampleIndex((Nresult), temp, N-1, extract_params, rparams);
    lweKeySwitch(result, bk->bk->ks, Nresult);

    delete_LweSample(Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial(cospoly);
}

LweSample* calculate_DST_bin_KS(TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample(extract_params);
    LweSample* result = new_LweSample(lweparam);
    TLweSample* temp = new_TLweSample(rparams);

    IntPolynomial* sinpoly = generate_sinpoly(rparams, K, k);
    torusPolynomialMultFFT(temp->a, sinpoly, ciphertext->a);
    torusPolynomialMultFFT(temp->b, sinpoly, ciphertext->b);
    tLweExtractLweSampleIndex((Nresult), temp, N-1, extract_params, rparams);
    lweKeySwitch(result, bk->bk->ks, Nresult);

    delete_LweSample(Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial(sinpoly);
    return result;
}

void calculate_DST_bin_KS(LweSample* result , TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample(extract_params);
    // LweSample* result = new_LweSample(lweparam);
    TLweSample* temp = new_TLweSample(rparams);

    IntPolynomial* sinpoly = generate_sinpoly(rparams, K, k);
    torusPolynomialMultFFT(temp->a, sinpoly, ciphertext->a);
    torusPolynomialMultFFT(temp->b, sinpoly, ciphertext->b);
    tLweExtractLweSampleIndex((Nresult), temp, N-1, extract_params, rparams);
    lweKeySwitch(result, bk->bk->ks, Nresult);

    delete_LweSample(Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial(sinpoly);
}

LweSample* calculate_DCT_bin(TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample(extract_params);
    // LweSample* result = new_LweSample(lweparam);
    TLweSample* temp = new_TLweSample(rparams);

    IntPolynomial* cospoly = generate_cospoly(rparams, K, k);
    torusPolynomialMultFFT(temp->a, cospoly, ciphertext->a);
    torusPolynomialMultFFT(temp->b, cospoly, ciphertext->b);
    tLweExtractLweSampleIndex((Nresult), temp, N-1, extract_params, rparams);
    // lweKeySwitch(result, bk->bk->ks, Nresult);

    // delete_LweSample(Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial(cospoly);
    return Nresult;
}

LweSample* calculate_DST_bin(TLweSample* ciphertext, const LweParams* lweparam ,const TLweParams* rparams, int k, TFheGateBootstrappingCloudKeySet* bk){
    const int32_t N = rparams->N;
    const LweParams *extract_params = &rparams->extracted_lweparams;
    LweSample* Nresult = new_LweSample(extract_params);
    // LweSample* result = new_LweSample(lweparam);
    TLweSample* temp = new_TLweSample(rparams);

    IntPolynomial* sinpoly = generate_sinpoly(rparams, K, k);
    torusPolynomialMultFFT(temp->a, sinpoly, ciphertext->a);
    torusPolynomialMultFFT(temp->b, sinpoly, ciphertext->b);
    tLweExtractLweSampleIndex((Nresult), temp, N-1, extract_params, rparams);
    // lweKeySwitch(result, bk->bk->ks, Nresult);

    // delete_LweSample(Nresult);
    delete_TLweSample(temp);
    delete_IntPolynomial(sinpoly);
    return Nresult;
}

TorusPolynomial* generate_abs(int32_t N, double alph){
    int32_t halfN = N/2;
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    Torus32 mu = modSwitchToTorus32(1, M_SLOT);
    for(int i=0; i<halfN; i++){ 
        testvect->coefsT[i] = modSwitchToTorus32(round((i - halfN)/alph), M_SLOT);
        testvect->coefsT[N-1-i] = -modSwitchToTorus32(round((i - halfN)/alph), M_SLOT); 
    }
    return testvect;
}

TorusPolynomial* generate_sign(int32_t N, double alph){
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    Torus32 mu = modSwitchToTorus32(1, M_SLOT);
    for(int i=0; i<N; i++){ 
        testvect->coefsT[i] = modSwitchToTorus32(1, M_SLOT);
    }
    return testvect;
}

TorusPolynomial* generate_logabs(int32_t N){
    int32_t halfN = N/2;
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    Torus32 mu = modSwitchToTorus32(1, M_SLOT);
    for(int i=0; i<halfN; i++){
        int32_t slot = floor(log2(i+1)+1);
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
    int32_t halfN = N/2;
    Torus32 mu = modSwitchToTorus32(halfN/ABSCALE, M_SLOT);

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
    result->b += mu;

    delete[] bara;
    delete_TorusPolynomial(testvect);
}

void Signbootstrap_woKS_FFT(LweSample *result, const LweBootstrappingKeyFFT *bk, const LweSample *x) {

    const TGswParams *bk_params = bk->bk_params;
    const TLweParams *accum_params = bk->accum_params;
    const LweParams *in_params = bk->in_out_params;
    const int32_t N = accum_params->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = in_params->n;
    int32_t halfN = N/2;

    TorusPolynomial *testvect = generate_sign(N, ABSCALE);
    // TorusPolynomial *testvect = generate_logabs(N);
    int32_t *bara = new int32_t[N];

    // Modulus switching
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    for (int32_t i = 0; i < n; i++) {
        bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
    }
    // Bootstrapping rotation and extraction
    tfhe_blindRotateAndExtract_FFT(result, testvect, bk->bkFFT, barb, bara, n, bk_params);

    delete[] bara;
    delete_TorusPolynomial(testvect);
}

void bootstrap_FFT(LweSample *result, const LweBootstrappingKeyFFT *bk, const LweSample *x) {

    LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

    bootstrap_woKS_FFT(u, bk, x);
    // Key switching
    lweKeySwitch(result, bk->ks, u);
    delete_LweSample(u);

}

void multi_FFT(LweSample* result, TLweSample* ciphertext, TFheGateBootstrappingCloudKeySet* bk, int32_t m){
    const LweParams* lweparams = bk->params->in_out_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;
    int32_t halfm = m/2;


    LweSample *DCT = new_LweSample_array(m, lweparams);
    LweSample *DST = new_LweSample_array(m, lweparams);
    for(int i=0; i<m; i++){
        LweSample * temp = calculate_bin_KS(ciphertext, lweparams ,tlwe_params, i+100, bk);
        bootstrap_FFT(&DCT[i], bk->bkFFT, &temp[0]);
        bootstrap_FFT(&DST[i], bk->bkFFT, &temp[1]);
        lweAddTo(&result[i], &DST[i], lweparams);
        lweAddTo(&result[i], &DCT[i], lweparams);
    }
    delete_LweSample_array(m, DST);
    delete_LweSample_array(m, DCT);
}

void binQFT_woKS(LweSample* result, TLweSample* ciphertext, int k, TFheGateBootstrappingCloudKeySet* bk){
    const LweParams* extract_param = bk->bkFFT->extract_params;
    const LweParams* lweparam = bk->params->in_out_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;

    LweSample *DCT = new_LweSample(extract_param);
    LweSample *DST = new_LweSample(extract_param);

    LweSample * temp = calculate_bin_KS(ciphertext, lweparam ,tlwe_params, k, bk);
    bootstrap_woKS_FFT(DCT, bk->bkFFT, &temp[0]);
    bootstrap_woKS_FFT(DST, bk->bkFFT, &temp[1]);
    lweAddTo(result, DST, extract_param);
    lweAddTo(result, DCT, extract_param);
    delete_LweSample(DST);
    delete_LweSample(DCT);
    delete_LweSample_array(2, temp);
}

LweSample* multi_FFT_exp(TLweSample* ciphertext, TFheGateBootstrappingCloudKeySet* bk, int32_t m){
    const LweParams* lweparams = bk->params->in_out_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;
    const LweParams* extracted_param = bk->bkFFT->extract_params;
    int32_t halfm = m/2;

    LweSample *DCT = new_LweSample_array(m, extracted_param);
    LweSample *DST = new_LweSample_array(m, extracted_param);
    for(int i=0; i<m; i++){
        LweSample * temp = calculate_bin_KS(ciphertext, lweparams ,tlwe_params, i+100, bk);
        bootstrap_woKS_FFT(&DCT[i], bk->bkFFT, &temp[0]);
        bootstrap_woKS_FFT(&DST[i], bk->bkFFT, &temp[1]);
        lweAddTo(&DCT[i], &DST[i], extracted_param);
    }
    delete_LweSample_array(m, DST);
    return DCT;
}

LweSample* multi_DCT_exp(TLweSample* ciphertext, TFheGateBootstrappingCloudKeySet* bk, int32_t m){
    const LweParams* lweparams = bk->params->in_out_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;
    int32_t halfm = m/2;

    LweSample *DCT = new_LweSample_array(m, lweparams);
    for(int i=0; i<m; i++){
        LweSample * temp = calculate_DCT_bin_KS(ciphertext, lweparams ,tlwe_params, i+100, bk);
        bootstrap_FFT(&DCT[i], bk->bkFFT, &temp[0]);
        delete_LweSample(temp);
    }
    return DCT;
}

void classify_vowel_QCT_30(LweSample* result, TLweSample* ciphertext, TFheGateBootstrappingCloudKeySet* bk, const LweKey* tlwe_key, std::string fname, std::ofstream& timeeval, std::ofstream& dinneval){
    const LweParams* lweparams = bk->params->in_out_params;
    const LweParams* extracted_params = bk->bkFFT->extract_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;

    FILE* fp;
    int32_t feat = 513;
    int32_t firstweight[feat][30];
    int32_t firstbias[30];
    int32_t secondweight[30][6];
    int32_t secondbias[6];
    
    fp = fopen(("weight/"+fname+"/firstweight.txt").c_str(), "r");
    for(int i=0; i<feat; i++){
        for(int j=0; j<30; j++){
            fscanf(fp, "%d ", &firstweight[i][j]);
        }
    }
    fclose(fp);

    fp = fopen(("weight/"+fname+"/firstbias.txt").c_str(), "r");
    for(int i=0; i<30; i++){
        fscanf(fp, "%d ", &firstbias[i]);
        firstbias[i] = modSwitchToTorus32(firstbias[i], M_SLOT); 
    }
    fclose(fp);

    fp = fopen(("weight/"+fname+"/secondweight.txt").c_str(), "r");
    for(int i=0; i<30; i++){
        for(int j=0; j<6; j++){
            fscanf(fp, "%d ", &secondweight[i][j]);
        }
    }
    fclose(fp);

    fp = fopen(("weight/"+fname+"/secondbias.txt").c_str(), "r");
    for(int i=0; i<6; i++){
        fscanf(fp, "%d ", &secondbias[i]);
        secondbias[i] = modSwitchToTorus32(secondbias[i], M_SLOT); 
    }
    fclose(fp);
    float time = -clock();
    // LweSample* QCT = new_LweSample_array(513, extracted_params);
    LweSample* QCT = new_LweSample_array(feat, extracted_params);
    for(int i=0; i<feat; i++){
        LweSample* temp = calculate_DCT_bin_KS(ciphertext, lweparams, tlwe_params, i, bk);
        bootstrap_woKS_FFT(&QCT[i], bk->bkFFT, temp);
        delete_LweSample(temp);
    }
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    std::cout <<"QCT took " << time << " seconds! "<< std::endl;
    timeeval << time << std::endl;
    LweSample* temp = new_LweSample(extracted_params);
    LweSample* KStemp = new_LweSample(lweparams);

    time = -clock();
    LweSample* FLayer = new_LweSample_array(30, extracted_params);
    for(int i=0; i<30; i++){
        for(int j=0; j<feat; j++){
            lweAddMulTo(temp, firstweight[j][i], &QCT[j], extracted_params);
        }
        lweKeySwitch(KStemp, bk->bkFFT->ks, temp);
        KStemp->b += firstbias[i];
        Signbootstrap_woKS_FFT(&FLayer[i], bk->bkFFT, KStemp);
        lweClear(temp, extracted_params);
        lweClear(KStemp, lweparams);
    }
    delete_LweSample(temp);
    delete_LweSample(KStemp);
    std::cout<<" "  << std::endl;
    

    for(int i=0; i<6; i++){
        for(int j=0; j<30; j++){
            lweAddMulTo(&result[i], secondweight[j][i], &FLayer[j], extracted_params);
        }
        (&result[i])->b += secondbias[i];
        Torus32 dec = lwePhase(&result[i], tlwe_key);
        dec = modSwitchFromTorus32(dec, M_SLOT);
        if(dec > M_SLOT/2){
        dec = (Torus32) -(M_SLOT - dec);
        }
        std::cout << dec << " ";
    }
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    std::cout << " " << std::endl;
    std::cout<< "Evaluation of NN took "<< time << " seconds! "<< std::endl;
    dinneval << time << std::endl;
    delete_LweSample_array(30, FLayer);
    delete_LweSample_array(feat, QCT);
}

void classify_vowel_QFT_30(LweSample* result, TLweSample* ciphertext, TFheGateBootstrappingCloudKeySet* bk, const LweKey* tlwe_key, const std::string fname, std::ofstream& timeeval, std::ofstream& dinneval){
    const LweParams* lweparams = bk->params->in_out_params;
    const LweParams* extracted_params = bk->bkFFT->extract_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;

    FILE* fp;
    int32_t feat = 513;
    int32_t firstweight[feat][30];
    int32_t firstbias[30];
    int32_t secondweight[30][6];
    int32_t secondbias[6];

    fp = fopen(("weight/"+fname+"/firstweight.txt").c_str(), "r");
    for(int i=0; i<feat; i++){
        for(int j=0; j<30; j++){
            fscanf(fp, "%d ", &firstweight[i][j]);
        }
    }
    fclose(fp);

    fp = fopen(("weight/"+fname+"/firstbias.txt").c_str(), "r");
    for(int i=0; i<30; i++){
        fscanf(fp, "%d ", &firstbias[i]);
        firstbias[i] = modSwitchToTorus32(firstbias[i], M_SLOT); 
    }
    fclose(fp);

    fp = fopen(("weight/"+fname+"/secondweight.txt").c_str(), "r");
    for(int i=0; i<30; i++){
        for(int j=0; j<6; j++){
            fscanf(fp, "%d ", &secondweight[i][j]);
        }
    }
    fclose(fp);

    fp = fopen(("weight/"+fname+"/secondbias.txt").c_str(), "r");
    for(int i=0; i<6; i++){
        fscanf(fp, "%d ", &secondbias[i]);
        secondbias[i] = modSwitchToTorus32(secondbias[i], M_SLOT); 
    }
    fclose(fp);
    float time = -clock();
    // LweSample* QCT = new_LweSample_array(513, extracted_params);
    LweSample* QFT = new_LweSample_array(feat, extracted_params);
    LweSample* ct = new_LweSample(lweparams);
    LweSample* st = new_LweSample(lweparams);
    LweSample *ctemp = new_LweSample(extracted_params);
    LweSample *stemp = new_LweSample(extracted_params);
    // LweSample* re = new_LweSample_array(2, lweparams);
    for(int i=0; i<feat; i++){
        lweClear(&QFT[i], extracted_params);
        // LweSample *kstemp = new_LweSample(lweparams);
        // binQFT_woKS(temp, ciphertext, bw[i][0], bk);
        // calculate_bin_KS_v2(re, ciphertext, lweparams ,tlwe_params, i, bk);
        calculate_DCT_bin_KS_v2(ct, ciphertext, lweparams ,tlwe_params, i, bk);
        calculate_DST_bin_KS(st, ciphertext, lweparams ,tlwe_params, i, bk);
        bootstrap_woKS_FFT(ctemp, bk->bkFFT, ct);
        bootstrap_woKS_FFT(stemp, bk->bkFFT, st);
        lweAddTo(&QFT[i], ctemp, extracted_params);
        lweAddTo(&QFT[i], stemp, extracted_params);
        // lweKeySwitch(kstemp, bk->bkFFT->ks, &QFT[i]);
        // bootstrap_woKS_FFT(&QFT[i], bk->bkFFT, kstemp);
        // delete_LweSample(kstemp);
        // for(int i=0; i<2; i++){
        //     lweClear(&re[i], lweparams);
        // }
        lweClear(ct, lweparams);
        lweClear(st, lweparams);
        lweClear(ctemp, extracted_params);
        lweClear(stemp, extracted_params);
    }
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    std::cout <<"QFT took " << time << " seconds! "<< std::endl;
    timeeval << time << std::endl;

    LweSample* temp = new_LweSample(extracted_params);
    LweSample* KStemp = new_LweSample(lweparams);

    time = -clock();
    LweSample* FLayer = new_LweSample_array(30, extracted_params);
    for(int i=0; i<30; i++){
        for(int j=0; j<feat; j++){
            lweAddMulTo(temp, firstweight[j][i], &QFT[j], extracted_params);
        }
        lweKeySwitch(KStemp, bk->bkFFT->ks, temp);
        KStemp->b += firstbias[i];
        Signbootstrap_woKS_FFT(&FLayer[i], bk->bkFFT, KStemp);
        lweClear(temp, extracted_params);
        lweClear(KStemp, lweparams);
    }
    std::cout<<" "  << std::endl;
    

    for(int i=0; i<6; i++){
        for(int j=0; j<30; j++){
            lweAddMulTo(&result[i], secondweight[j][i], &FLayer[j], extracted_params);
        }
        (&result[i])->b += secondbias[i];
        Torus32 dec = lwePhase(&result[i], tlwe_key);
        dec = modSwitchFromTorus32(dec, M_SLOT);
        if(dec > M_SLOT/2){
        dec = (Torus32) -(M_SLOT - dec);
        }
        std::cout << dec << " ";
    }
    std::cout << " " << std::endl;
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    std::cout << " " << std::endl;
    std::cout<< "Evaluation of NN took "<< time << " seconds! "<< std::endl;
    dinneval << time << std::endl;

    // delete_LweSample_array(2, re);
    delete_LweSample(ct);
    delete_LweSample(st);
    delete_LweSample(ctemp);
    delete_LweSample(stemp);
    delete_LweSample(temp);
    delete_LweSample(KStemp);
    delete_LweSample_array(30, FLayer);
    delete_LweSample_array(feat, QFT);
}

void classify_vowel_twolayer_QCT_30(LweSample* result, TLweSample* ciphertext, TFheGateBootstrappingCloudKeySet* bk, const LweKey* tlwe_key, std::string fname){
    const LweParams* lweparams = bk->params->in_out_params;
    const LweParams* extracted_params = bk->bkFFT->extract_params;
    const TLweParams* tlwe_params = bk->params->tgsw_params->tlwe_params;

    FILE* fp;
    int32_t feat = 513;
    int32_t firstweight[feat][30];
    int32_t firstbias[30];
    int32_t secondweight[30][30];
    int32_t secondbias[30];
    int32_t thirdweight[30][6];
    int32_t thirdbias[6];
    
    fp = fopen("twolayer/firstweight.txt", "r");
    for(int i=0; i<feat; i++){
        for(int j=0; j<30; j++){
            fscanf(fp, "%d ", &firstweight[i][j]);
        }
    }
    fclose(fp);

    fp = fopen("twolayer/firstbias.txt", "r");
    for(int i=0; i<30; i++){
        fscanf(fp, "%d ", &firstbias[i]);
        firstbias[i] = modSwitchToTorus32(firstbias[i], M_SLOT); 
    }
    fclose(fp);

    fp = fopen("twolayer/secondweight.txt", "r");
    for(int i=0; i<30; i++){
        for(int j=0; j<30; j++){
            fscanf(fp, "%d ", &secondweight[i][j]);
        }
    }
    fclose(fp);

    fp = fopen("twolayer/secondbias.txt", "r");
    for(int i=0; i<30; i++){
        fscanf(fp, "%d ", &secondbias[i]);
        secondbias[i] = modSwitchToTorus32(secondbias[i], M_SLOT); 
    }
    fclose(fp);

    fp = fopen("twolayer/thirdweight.txt", "r");
    for(int i=0; i<30; i++){
        for(int j=0; j<6; j++){
            fscanf(fp, "%d ", &thirdweight[i][j]);
        }
    }
    fclose(fp);

    fp = fopen("twolayer/thirdbias.txt", "r");
    for(int i=0; i<6; i++){
        fscanf(fp, "%d ", &thirdbias[i]);
        thirdbias[i] = modSwitchToTorus32(thirdbias[i], M_SLOT); 
    }
    fclose(fp);

    float time = -clock();
    // LweSample* QCT = new_LweSample_array(513, extracted_params);
    LweSample* QCT = new_LweSample_array(feat, extracted_params);
    for(int i=0; i<feat; i++){
        LweSample* temp = calculate_DCT_bin_KS(ciphertext, lweparams, tlwe_params, i, bk);
        bootstrap_woKS_FFT(&QCT[i], bk->bkFFT, temp);
        delete_LweSample(temp);
    }
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    std::cout <<"QCT took " << time << " seconds! "<< std::endl;

    LweSample* temp = new_LweSample(extracted_params);
    LweSample* KStemp = new_LweSample(lweparams);

    time = -clock();
    LweSample* FLayer = new_LweSample_array(30, extracted_params);
    for(int i=0; i<30; i++){
        for(int j=0; j<feat; j++){
            lweAddMulTo(temp, firstweight[j][i], &QCT[j], extracted_params);
        }
        lweKeySwitch(KStemp, bk->bkFFT->ks, temp);
        KStemp->b += firstbias[i];
        Signbootstrap_woKS_FFT(&FLayer[i], bk->bkFFT, KStemp);
        lweClear(temp, extracted_params);
        lweClear(KStemp, lweparams);
    }
    delete_LweSample(temp);
    delete_LweSample(KStemp);
    

    LweSample* SLayer = new_LweSample_array(30, extracted_params);
    LweSample* Stemp = new_LweSample(extracted_params);
    LweSample* SKStemp = new_LweSample(lweparams);
    for(int i=0; i<30; i++){
        for(int j=0; j<30; j++){
            lweAddMulTo(Stemp, secondweight[j][i], &FLayer[j], extracted_params);
        }
        lweKeySwitch(SKStemp, bk->bkFFT->ks, Stemp);
        SKStemp->b += secondbias[i];
        Signbootstrap_woKS_FFT(&SLayer[i], bk->bkFFT, SKStemp);
        lweClear(Stemp, extracted_params);
        lweClear(SKStemp, lweparams);
    }
    delete_LweSample(Stemp);
    delete_LweSample(SKStemp);

    for(int i=0; i<6; i++){
        for(int j=0; j<30; j++){
            lweAddMulTo(&result[i], thirdweight[j][i], &SLayer[j], extracted_params);
        }
        (&result[i])->b += thirdbias[i];
        Torus32 dec = lwePhase(&result[i], tlwe_key);
        dec = modSwitchFromTorus32(dec, M_SLOT);
        if(dec > M_SLOT/2){
        dec = (Torus32) -(M_SLOT - dec);
        }
        std::cout << dec << " ";
    }
    
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    std::cout << " " << std::endl;
    std::cout<< "Evaluation of NN took "<< time << " seconds! "<< std::endl;
    delete_LweSample_array(30, FLayer);
    delete_LweSample_array(30, SLayer);
    delete_LweSample_array(feat, QCT);
}

using namespace std;

int main() {

    const string fname = "12"; // filename : 01 ~ 12
    // const string fname = "12";
    const string feat = "QCT";
    // const string feat = "QFT"; // QCT or QFT
    // const string param = "old";
    const string param = "new"; // 'old' params, or 'new' params

    FILE* cloud_key = fopen(("cloud110b_"+param+".key").c_str(), "rb");
    TFheGateBootstrappingCloudKeySet* bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    fclose(cloud_key);

    const TFheGateBootstrappingParameterSet* params = bk->params;
    const LweParams* lweparams = params->in_out_params;
    const TLweParams* tlwe_params = params->tgsw_params->tlwe_params;
    const LweParams* extracted_params = bk->bkFFT->extract_params;
    const int N = tlwe_params->N;
    TLweSample* ciphertext = new_TLweSample_array(C_NUM, tlwe_params);

    FILE* cloud_data = fopen(("cloudMulti110b_"+param+".data").c_str(), "rb");
    for(int i=0; i<C_NUM; i++){
        import_tlweSample_fromFile(cloud_data, &ciphertext[i], tlwe_params);
    }
    fclose(cloud_data);

    FILE* secret_key = fopen(("secret110b_"+param+".key").c_str(), "rb");
    TFheGateBootstrappingSecretKeySet* key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
    fclose(secret_key);

    const TFheGateBootstrappingParameterSet* Param = key->params;
    const LweKey* lwekey = new_LweKey(lweparams);
    const TLweKey* tlwekey = &key->tgsw_key->tlwe_key;
    lwekey = key->lwe_key;
    LweKey* tlwe_key = new_LweKey(extracted_params);
    tLweExtractKey(tlwe_key, tlwekey);
    // LweSample* res = new_LweSample_array(2, lweparams);
    int32_t halfM = M/2;

    ofstream vowel ("./result/"+feat+"_vowel_110b_A"+to_string((int32_t) ABSCALE)+"_"+fname+"_"+param+".txt");
    ofstream timeeval ("./result/"+feat+"_vowel_timeeval_"+param+".txt", ios::app);
    ofstream dinneval ("./result/"+feat+"_vowel_dinneval_"+param+".txt", ios::app);
    LweSample* res = new_LweSample_array(6, extracted_params);
    float time = -clock();
    for(int i=0; i<C_NUM; i++){
        // classify_vowel_twolayer_QCT_30(res, &ciphertext[i], bk, tlwe_key);
        if(feat == "QCT"){
            classify_vowel_QCT_30(res, &ciphertext[i], bk, tlwe_key, fname, timeeval, dinneval);
        }
        if(feat == "QFT"){
            classify_vowel_QFT_30(res, &ciphertext[i], bk, tlwe_key, fname, timeeval, dinneval);
        }

        for(int j=0; j<6; j++){
            Torus32 dec = lwePhase(&res[j], tlwe_key);
            dec = modSwitchFromTorus32(dec, M_SLOT);
            if(dec > M_SLOT/2){
                dec = -(M_SLOT - dec);
            }
            if (j < 5){
                vowel << dec << " ";
            }
            else{
                vowel << dec << endl;
            }
            lweClear(&res[j], extracted_params);
        }
        cout << "-- " << i+1 << "-th vowel classification finished! --" << endl;
        cout << " " << endl;
    }
    time += clock();
    time = time/(CLOCKS_PER_SEC);
    cout << "---  We classified " << C_NUM << "voices in "<< time << " seconds!  ---" << endl;
    cout << "--- Average time per classification is " << time/C_NUM << " seconds ---" << endl;
    // int32_t decr[C_NUM][M];
    // ofstream cosine ("QCTResult110b_L15_K9.txt");
    // ofstream sine ("QSTResult0b_L15_K9.txt");
    // ofstream qft ("QFTResult110b_L15_K9.txt");
    // float time = -clock();
    // LweSample* res = new_LweSample(extracted_params);
    // uint32_t a = 0;
    // for(int i=0; i<C_NUM; i++){
    //     for(int j=0; j<513; j++){
    //         // binQFT_woKS(res, &ciphertext[i], j, bk);
    //         LweSample* rel = calculate_DCT_bin_KS(&ciphertext[i], lweparams, tlwe_params, j, bk);
    //         // LweSample* img = calculate_DST_bin(&ciphertext[i], lweparams, tlwe_params, j, bk);
    //         // LweSample* QFT = calculate_bin(&ciphertext[i], tlwe_params, j);
    //         Signbootstrap_woKS_FFT(res, bk->bkFFT, rel);
    //         Torus32 cosdec = lwePhase(res, tlwe_key);
    //         // Torus32 sindec = lwePhase(img, tlwe_key);
    //         // Torus32 qftdec = lwePhase(res, tlwe_key);
    //         cosdec = modSwitchFromTorus32(cosdec, M_SLOT);
    //         // qftdec = modSwitchFromTorus32(qftdec, M_SLOT);
    //         if (j < 512){
    //             cosine << cosdec << " ";
    //             // qft << qftdec << " ";
    //             // sine << sindec << " ";
    //             // myfile << sq << " ";
    //         }
    //         else{ 
    //             cosine << cosdec << endl;
    //             // qft << qftdec << endl;
    //             // sine << sindec << endl;
    //             // myfile << sq << endl;
    //         }
    //         // delete_LweSample(rel);
    //         // delete_LweSample(img);
    //         // delete_LweSample(img);
    //         // delete_LweSample_array(2, QFT);
    //         lweClear(res, extracted_params);
    //         // delete_LweSample(res);
    //     }
    //     // delete_LweSample_array(M, FFTR);
    //     cout << "-- " << i+1 << "-th DCT finished! --" << endl;
    // }
    // cosine.close();
    // // sine.close();
    // // qft.close();
    // time += clock();
    // time = time/(CLOCKS_PER_SEC);
    // cout << "---  We calculated " << C_NUM << "FFTs in "<< time << " seconds!  ---" << endl;
    // cout << "--- Average time per FFT is " << time/C_NUM << " seconds ---" << endl;

    // ofstream myfile ("NEW_QFT_110b_woKS_A32_svm.txt");
    // ofstream scmyfile ("NEW_QFT_110b_score_woKS_A32_svm.txt");
    // // // ofstream scmyfile ("Only_one.txt");
    // float time = -clock();
    // LweSample* res = new_LweSample(extracted_params);

    // for(int i=0; i<C_NUM; i++){
    //     classify_WM_woKS_ww_QFT(res ,&ciphertext[i], bk, 1);
    //     // int32_t ai = lwePhase(res, tlwe_key);
    //     // int32_t dec = ai > 0? 1 : 0;
    //     int32_t dec = lwePhase(res, tlwe_key);
    //     int32_t classif = dec > 0? 1 : 0;
    //     int32_t score = modSwitchFromTorus32(dec, M_SLOT);
    //     // cout << "decresult : " << dec << "\n" <<endl;
    //     // decr[i] = dec;
    //     myfile << classif << endl;
    //     scmyfile<< score << endl;
    //     lweClear(res, extracted_params);
    //     cout<< i+1 <<"-th ciphertext classification is done! " << (C_NUM - i -1)*100/C_NUM << "% left..."<< endl;
    // }
    // time += clock();
    // time = time/(CLOCKS_PER_SEC);
    // cout << "---  We classified " << C_NUM << " samples in "<< time << " seconds!  ---" << endl;
    // cout << "--- Average time per classification is " << time/C_NUM << " seconds ---" << endl;
    // myfile.close();
    // scmyfile.close();
    // printf("Hello World!");
    
    delete_TLweParams((TLweParams*)tlwe_params);
    delete_LweParams((LweParams*)lweparams);
    delete_gate_bootstrapping_cloud_keyset(bk);
    delete_TLweSample_array(C_NUM, ciphertext);
    delete_LweKey((LweKey*)lwekey);
    // delete_LweSample_array(2, res);
    
}