#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/tfhe_io.h>

TLweKey *ExtracttLweKey(TFheGateBootstrappingSecretKeySet *secret_key, TLweParams *params) {
    TLweKey *tlwe_key = new_TLweKey(params);
    tlwe_key->key = secret_key->tgsw_key->key;
    return tlwe_key;
}