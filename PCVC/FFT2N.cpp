#include "fft_functions.h"
#include "fft_functions.cpp"
// Generate IntPolynomial of coefficients of result[N-1-i] = 2^K * sin(2*pi*k*i/N) *inverse packing


int main() {

    const string filename[13] = {"01", "01_2", "02", "03","04", "05", "06", "07", "08", "09","10","11","12"};
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

    // FILE* cloud_data = fopen(("cloudMulti110b_"+param+".data").c_str(), "rb");
    // for(int i=0; i<C_NUM; i++){
    //     import_tlweSample_fromFile(cloud_data, &ciphertext[i], tlwe_params);
    // }
    // fclose(cloud_data);

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

    // ofstream vowel ("./result/"+feat+"_vowel_110b_A"+to_string((int32_t) ABSCALE)+"_"+fname+"_"+param+".txt");
    ofstream timeeval ("./result/"+feat+"_vowel_timeeval_"+param+".txt", ios::app);
    ofstream dinneval ("./result/"+feat+"_vowel_dinneval_"+param+".txt", ios::app);
    ofstream clastime ("./result/"+feat+"_vowel_classtime_"+param+".txt", ios::app);
    LweSample* res = new_LweSample_array(6, extracted_params);
    // float time = -clock();
    for(int k=0; k<13; k++){
        string fname = filename[k];
        FILE* cloud_data = fopen(("cloudMulti110b_"+param+"_"+fname+".data").c_str(), "rb");
        for(int i=0; i<C_NUM; i++){
            import_tlweSample_fromFile(cloud_data, &ciphertext[i], tlwe_params);
        }
        fclose(cloud_data);
        float time = -clock();
        ofstream vowel ("./result/"+feat+"_vowel_110b_A"+to_string((int32_t) ABSCALE)+"_"+fname+"_"+param+".txt");
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
            cout << "-- " << k << "-" << i+1 << "-th vowel classification finished! --" << endl;
            cout << " " << endl;
        }
        time += clock();
        time = time/(CLOCKS_PER_SEC);
        cout << "---  We classified " << C_NUM << "voices in "<< time << " seconds!  ---" << endl;
        cout << "--- Average time per classification is " << time/C_NUM << " seconds ---" << endl;
        clastime << time/C_NUM << endl;
        vowel.close();
    }
    // time += clock();
    // time = time/(CLOCKS_PER_SEC);
    // cout << "---  We classified " << C_NUM << "voices in "<< time << " seconds!  ---" << endl;
    // cout << "--- Average time per classification is " << time/C_NUM << " seconds ---" << endl;
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