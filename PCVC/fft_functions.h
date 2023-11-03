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
#include <vector>
#include <thread>


#define K  9
#define M  100
#define ABSCALE pow(2., 5)
#define M_SLOT 700
#define C_NUM 138

using namespace std;