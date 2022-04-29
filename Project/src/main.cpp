#include <NTL/GF2X.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <stdlib.h>
#include <iostream>

#include "DualCiphers.h"
#include "sm4_function.h"
#include "stdlib.h"
#include "time.h"
#include "GF2E_function.h"
#include "wbsm4_function.h"
#include "Collision_Attack.h"

using namespace std;
using namespace NTL;

int main()
{
    mat_GF2 x,x1,x2,x3;
    vec_GF2 xx;
    uint32_t y;


    unsigned char key[16] =
    {0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
    unsigned char plain[16]=
    {0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
    cout << "输入明文:" << endl;
    printstate_sm4(plain);
    cout << "输入密钥:" << endl;
    printstate_sm4(key);
    cout << "原SM4分组密码加密结果:" <<endl;
    SM4_SelfCheck();
    mat_GF2 * Oi = oi_matrix(); 
    cout << "SM4白盒加密结果:" <<endl;
    wbsm4_gen(Oi,key);
    wbsm4_Encrypt(plain,Oi);
    cout << endl;

    return 0;
}