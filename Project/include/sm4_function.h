#pragma once
#include "stdint.h"

using namespace NTL;
//rotate n bits to the left in a 32bit buffer
#define SM4_Rotl32(buf, n) (((buf)<<n)|((buf)>>(32-n)))



void SM4_KeySchedule(unsigned char MK[], unsigned int rk[]);

void SM4_Encrypt(unsigned char MK[],unsigned char PlainText[],unsigned char CipherText[]);

void SM4_Decrypt(unsigned char MK[],unsigned char CipherText[], unsigned char PlainText[]);

int SM4_SelfCheck();

unsigned char sm4Sbox(unsigned char inch);

NTL::mat_GF2 Ni_sm4Sbox(mat_GF2 & input);

//NTL::mat_GF2 Construct_B(const uint8_t x[]);

//NTL::mat_GF2 *Lj_matrix() ; //创建白盒线性矩阵

NTL::mat_GF2 *Lj_matrix11() ;

NTL::mat_GF2 *oi_matrix(); //创建在白盒中随机生成的32x32的非奇异矩阵

NTL::mat_GF2 * Fi_matrix(const mat_GF2 * temp);

NTL::mat_GF2 * Gi_matrix(const mat_GF2 * temp);

NTL::mat_GF2 setup_matB1();

NTL::mat_GF2 setup_matB2();

NTL::mat_GF2 setup_matB3();

void printstate_sm4(unsigned char * in);

//mat_GF2 * oi_ceshi();

//uint32_t SM411_Encrypt(unsigned char MK[],unsigned char PlainText[],unsigned char CipherText[]);
