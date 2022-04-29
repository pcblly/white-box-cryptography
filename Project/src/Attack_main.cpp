#include <stdio.h>
#include <NTL/GF2X.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <stdlib.h>
#include <iostream>

#include "wbsm4_function.h"
#include "Collision_Attack.h"
#include "DualCiphers.h"
#include "sm4_function.h"
#include "stdlib.h"
#include "time.h"
#include "GF2E_function.h"

using namespace std;
using namespace NTL;



int main(void){

    
    huifuSBOX();

    huifuSboxNi();

    huifu_invEij();

   
    key1_function();
    key2_function();
    key3_function();
    key4_function();
    key0_function();
    Ki_found();
    key_end();
    


    /*
    int i = 1,j;
    uint8_t k;
    mat_GF2 x =  matGF2FromUint8(i);
    j = uint8FromMatGF2(x);
    k = uint8FromMatGF2(x);
    cout << transpose(x) <<endl;
    cout << j <<"aa" <<endl;
    printf("%d",k);

    uint8_t x = 0x01;
    uint8_t y;
    mat_GF2 output;
    output = Ni_sm4Sbox(x);
    //output = transpose(matGF2FromUint8(y)) ;
    cout << output;
*/

    return 0;
}