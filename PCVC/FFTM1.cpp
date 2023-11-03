#include <tfhe/polynomials.h>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <tfhe/tfhe_garbage_collector.h>

#define AS_l 15
#define FILELENGTH 138
#define AUDIO_LENGTH  1024
#define KS_STDEV  pow(2., -24) //80bit/ 110bit -> -18
#define BK_STDEV  pow(2., -31)


// Extract tLwe key(Actually a tGSW key) & tLwe parameters from secret key & tLwe parameters
TLweKey *ExtracttLweKey(TFheGateBootstrappingSecretKeySet *secret_key, const TLweParams *params) {
    TLweKey *tlwe_key = new_TLweKey(params);
    tlwe_key->key = secret_key->tgsw_key->key;
    return tlwe_key;
}

// Define struct of the audio file, (fs : frequency, length : length of our audio data,audio data)
struct audiofile{
    int fs;
    int length; 
    float audio[FILELENGTH][AUDIO_LENGTH];
};


// Read [-1, 1]-normalized audio file, which is in a csv format
audiofile *read_audio(const char *filename){
    audiofile *Audio = (audiofile*)malloc(sizeof(audiofile));
    FILE *fp;
    Audio->fs = 1000;
    Audio->length = AUDIO_LENGTH;
    float temp;
    fp = fopen(filename, "r");
    for(int i=0; i<FILELENGTH; i++){
        for(int j=0; j<AUDIO_LENGTH; j++){
            fscanf(fp, "%f ", &temp);
        Audio->audio[i][j] = temp;
        }
    }
    fclose(fp);
    return Audio;
}


// Encode [-1, 1]-normalized audio file into TorusPolynomial & We don't use dtot32 function in TFHE
// It is actually a morphism from [-1, 1] to [-2^l, 2^l], where l is the rounded log2 value of our scaling 
TorusPolynomial *Encode_audio(audiofile *file, int scale, const TLweParams *tlwe_params){ 
    const int32_t N = tlwe_params->N;
    TorusPolynomial *result = new_TorusPolynomial_array(FILELENGTH, N);
    const int32_t SCALE = 1<<scale;
    for(int i=0; i<FILELENGTH; i++){
        for (int j=0; j<N; j++) {
            TorusPolynomial *temppolynomial = result + i;
            temppolynomial->coefsT[j] = (Torus32) round(file->audio[i][j] * SCALE);
        }
    }
    return result;
}

// Encode [-1, 1]-normalized audio file into TorusPolynomial & We don't use dtot32 function in TFHE
// It is actually a morphism from [-1, 1] to [-2^l, 2^l], where l is the rounded log2 value of our scaling 
TLweSample *Encode_audio_encrypt(audiofile *file, int scale, const TLweParams *tlwe_params, TLweKey *tlwe_key){ 
    const int32_t N = tlwe_params->N;
    TorusPolynomial *result = new_TorusPolynomial(N);
    TLweSample* encAudio = new_TLweSample_array(FILELENGTH, tlwe_params);
    const int32_t SCALE = 1<<scale;
    for(int i=0; i<FILELENGTH; i++){
        for (int j=0; j<N; j++) {
            result->coefsT[j] = (Torus32) round(file->audio[i][j] * SCALE);
            tLweSymEncrypt(&encAudio[i], result, tlwe_params->alpha_min, tlwe_key);
        }
    }
    return encAudio;
}

TFheGateBootstrappingParameterSet *default_80bit_gate_bootstrapping_parameters()
{
    static const int n = 630;
    static const int N = 1024;
    static const int k = 1;
    static const double max_stdev = pow(2., -11);

    static const int bk_Bgbit    = 6;  
    static const int bk_l        = 4;
    static const double bk_stdev = BK_STDEV; 

    static const int ks_basebit  = 2; 
    static const int ks_length   = 8;
    static const double ks_stdev = KS_STDEV;


    LweParams  *params_in    = new_LweParams (n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk    = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}

TFheGateBootstrappingParameterSet *default_110bit_gate_bootstrapping_parameters()
{
    static const int n = 630;
    static const int N = 1024;
    static const int k = 1;
    static const double max_stdev = pow(2., -11);

    static const int bk_Bgbit    = 6;  
    static const int bk_l        = 4;
    static const double bk_stdev = pow(2., -29); 

    static const int ks_basebit  = 2; 
    static const int ks_length   = 8;
    static const double ks_stdev = pow(2., -17);


    LweParams  *params_in    = new_LweParams (n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk    = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}

TFheGateBootstrappingParameterSet *default_0bit_gate_bootstrapping_parameters()
{
    static const int n = 630;
    static const int N = 1024;
    static const int k = 1;
    static const double max_stdev = 0;

    static const int bk_Bgbit    = 6;  
    static const int bk_l        = 4;
    static const double bk_stdev = 0; 

    static const int ks_basebit  = 2; 
    static const int ks_length   = 8;
    static const double ks_stdev = 0;


    LweParams  *params_in    = new_LweParams (n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk    = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}

TFheGateBootstrappingParameterSet *default_new110bit_gate_bootstrapping_parameters()
{
    static const int n = 750;
    static const int N = 1024;
    static const int k = 1;
    static const double max_stdev = pow(2., -11);

    static const int bk_Bgbit    = 3;  
    static const int bk_l        = 9;
    static const double bk_stdev = pow(2., -29); 

    static const int ks_basebit  = 8; 
    static const int ks_length   = 3;
    static const double ks_stdev = pow(2., -21);


    LweParams  *params_in    = new_LweParams (n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk    = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}

TFheGateBootstrappingParameterSet *default_new0bit_gate_bootstrapping_parameters()
{
    static const int n = 750;
    static const int N = 1024;
    static const int k = 1;
    static const double max_stdev = pow(2., -11);

    static const int bk_Bgbit    = 3;  
    static const int bk_l        = 9;
    static const double bk_stdev = 0; 

    static const int ks_basebit  = 8; 
    static const int ks_length   = 3;
    static const double ks_stdev = 0;


    LweParams  *params_in    = new_LweParams (n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk    = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}

using namespace std;

int main() {
    const string filename[13] = {"01", "01_2", "02", "03","04", "05", "06", "07", "08", "09","10","11","12"};
    const string param = "new";
    TFheGateBootstrappingParameterSet* params;
    
    if(param=="old"){
        params = default_110bit_gate_bootstrapping_parameters();
        // TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);
    }
    if(param=="new"){
        params = default_new110bit_gate_bootstrapping_parameters();
        // TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);
    }
    
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);

    FILE* secret_key = fopen(("secret110b_"+param+".key").c_str(), "wb");
    export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
    fclose(secret_key);

    FILE* cloud_key = fopen(("cloud110b_"+param+".key").c_str(), "wb");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
    fclose(cloud_key);

    // const LweBootstrappingKeyFFT *bs_key = key->cloud.bkFFT;
    const TLweParams* tlwe_param = params->tgsw_params->tlwe_params;
    TLweKey *TLWE_KEY = ExtracttLweKey(key, tlwe_param);
    // audiofile *audio = read_audio("voweldat.csv");
    // audiofile *audio = read_audio("12.txt");
    // TorusPolynomial *encoded_message = Encode_audio(audio, AS_l, tlwe_param);
    // free(audio);
    // TLweSample *ciphertext = new_TLweSample(tlwe_param);
    // tLweSymEncrypt(ciphertext, encoded_message,tlwe_param->alpha_min, TLWE_KEY);
    float time = -clock();
    printf("---------------  We Now Encode & Encrypt your Data!  ---------------\n");
    printf("\n");
    printf("---------------------------  Encoding...  --------------------------\n");
    printf("\n");
    for(int i=0; i < 13; i++){
        audiofile *audio = read_audio((filename[i]+".txt").c_str());
        TLweSample *EncAudio = Encode_audio_encrypt(audio, AS_l, tlwe_param, TLWE_KEY);
        FILE* cloud_data = fopen(("cloudMulti110b_"+param+"_"+filename[i]+".data").c_str(), "wb");
        for(int j=0; j<FILELENGTH; j++){
            export_tlweSample_toFile(cloud_data, &EncAudio[j], tlwe_param);
        }
        fclose(cloud_data);
    }

    time += clock();
    time = time/(CLOCKS_PER_SEC);
    printf("----------------  Finished Encoding & Encryption!  -----------------\n");
    printf("-----------------  We done it in %f seconds!  ----------------\n", time);
    printf("----------------  Per sample took %f seconds!  ---------------\n", time/FILELENGTH);
    printf("\n");
    printf("\n");

    printf("---------  We now send your Encrypted data to the cloud!  ----------\n");



    printf("\n");
    printf("------------------------  Sending Finished!  -----------------------\n");

    // tLwePhase(result, ciphertext, TLWE_KEY);
    printf("Hello world!\n");
    // free(audio);
    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);
    // delete_TorusPolynomial(encoded_message);
    delete_TLweParams((TLweParams*)tlwe_param);
    // delete_TLweSample_array(FILELENGTH, EncAudio);
    // delete_TLweSample(ciphertext);

}