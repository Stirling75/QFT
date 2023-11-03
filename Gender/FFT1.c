#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define AS_l 21
#define AUDIO_LENGTH  1024

// Extract tLwe key(Actually a tGSW key) & tLwe parameters from secret key & tLwe parameters
TLweKey *ExtracttLweKey(TFheGateBootstrappingSecretKeySet *secret_key, const TLweParams *params) {
    TLweKey *tlwe_key = new_TLweKey(params);
    tlwe_key->key = secret_key->tgsw_key->key;
    return tlwe_key;
}

// Define struct of the audio file, (fs : frequency, length : length of our audio data,audio data)
typedef struct audiofile{
    int fs;
    int length; 
    float audio[AUDIO_LENGTH];
} audiofile;


// Read [-1, 1]-normalized audio file, which is in a csv format
audiofile *read_audio(char *filename){
    audiofile *Audio = malloc(sizeof(audiofile));
    FILE *fp;
    Audio->fs = 1000;
    Audio->length = AUDIO_LENGTH;
    float temp;
    fp = fopen(filename, "r");
    for(int i=0; i<AUDIO_LENGTH; i++){
        fscanf(fp, "%f ", &temp);
        Audio->audio[i] = temp;
    }
    fclose(fp);
    return Audio;
}


// Encode [-1, 1]-normalized audio file into TorusPolynomial & We don't use dtot32 function in TFHE
// It is actually a morphism from [-1, 1] to [-2^l, 2^l], where l is the rounded log2 value of our scaling 
TorusPolynomial *Encode_audio(audiofile *file, int scale, const TLweParams *tlwe_params){ 
    const int32_t N = tlwe_params->N;
    TorusPolynomial *result = new_TorusPolynomial(N);
    const int32_t SCALE = 1<<scale;
    for(int i=0; i<N; i++){
        result->coefsT[i] = (Torus32) round(file->audio[i] * SCALE);
    }
    return result;
}

int main() {
    TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(110);
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);

    FILE* secret_key = fopen("secret.key", "wb");
    export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
    fclose(secret_key);

    FILE* cloud_key = fopen("cloud.key", "wb");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
    fclose(cloud_key);

    // const LweBootstrappingKeyFFT *bs_key = key->cloud.bkFFT;
    const TLweParams* tlwe_param = params->tgsw_params->tlwe_params;
    TLweKey *TLWE_KEY = ExtracttLweKey(key, tlwe_param);
    audiofile *audio = read_audio("mat_4874_1024.csv");
    TorusPolynomial *encoded_message = Encode_audio(audio, AS_l, tlwe_param);
    TorusPolynomial *result = new_TorusPolynomial(tlwe_param->N);
    free(audio);
    TLweSample *ciphertext = new_TLweSample(tlwe_param);
    
    tLweSymEncrypt(ciphertext, encoded_message,tlwe_param->alpha_min, TLWE_KEY);
    FILE* cloud_data = fopen("cloud.data", "wb");
    export_tlweSample_toFile(cloud_data, ciphertext, tlwe_param);
    fclose(cloud_data);

    // tLwePhase(result, ciphertext, TLWE_KEY);
    printf("Hello world!\n");
    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);
    delete_TorusPolynomial(encoded_message);
    delete_TorusPolynomial(result);
    delete_TLweParams((TLweParams*)tlwe_param);
    delete_TLweSample(ciphertext);

}