#include <stdio.h>
#include <string.h>
#include <iostream>
#include <NTL/GF2X.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>


#include "sm4_function.h"
#include "sm4_key.h"
#include "GF2E_function.h"

using namespace NTL;
using namespace std;

void SM4_KeySchedule(unsigned char MK[],unsigned int rk[])
{
    unsigned int tmp,buf,K[36];
    int i;
    for(i=0;i<4;i++)
{
    K[i]=SM4_FK[i]^( (MK[4*i]<<24) | (MK[4*i+1]<<16)
    |(MK[4*i+2]<<8) | (MK[4*i+3]) );
}
    for(i=0;i<32;i++)
    {
    tmp =K[i+1]^K[i+2]^K[i+3]^ SM4_CK[i];
//nonlinear operation
    buf= (SM4_Sbox[(tmp >> 24) & 0xFF]) << 24
    |(SM4_Sbox[(tmp >> 16) & 0xFF]) << 16
    |(SM4_Sbox[(tmp >> 8) & 0xFF]) << 8
    |(SM4_Sbox[tmp & 0xFF]);
//linear operation
    K[i+4]=K[i]^((buf)^(SM4_Rotl32((buf),13))^(SM4_Rotl32((buf),23)));
    rk[i]=K[i+4];
   // cout << transpose(matGF2FromUint32(rk[i]))  << endl;
}
}

void SM4_Encrypt(unsigned char MK[],unsigned char PlainText[],unsigned char CipherText[])
{
unsigned int rk[32],X[36],tmp,buf;
int i,j;
SM4_KeySchedule(MK,rk);
for(j=0;j<4;j++)
{
X[j]=(PlainText[j*4]<<24) |(PlainText[j*4+1]<<16)
|(PlainText[j*4+2]<<8)|(PlainText[j*4+3]);
}
for(i=0;i<32;i++)
{
tmp = X[i+1]^X[i+2]^X[i+3]^rk[i];
//nonlinear operation
buf= ( SM4_Sbox[(tmp >> 24) & 0xFF]) << 24
|(SM4_Sbox[(tmp >> 16) & 0xFF]) << 16
|(SM4_Sbox[(tmp >> 8) & 0xFF]) << 8
|(SM4_Sbox[tmp & 0xFF]);
//linear operation
X[i+4]=X[i]^(buf^SM4_Rotl32((buf),2)^ SM4_Rotl32((buf),10)
^ SM4_Rotl32((buf),18)^ SM4_Rotl32((buf),24));
}
for(j=0;j<4;j++)
{
CipherText[4*j]=(X[35-j]>> 24)& 0xFF;
CipherText[4*j+1]=(X[35-j]>> 16)& 0xFF;
CipherText[4*j+2]=(X[35-j]>> 8)& 0xFF;
CipherText[4*j+3]=(X[35-j])& 0xFF;
}
}

void SM4_Decrypt(unsigned char MK[],unsigned char CipherText[],unsigned char PlainText[])
{
unsigned int rk[32],X[36],tmp,buf;
int i,j;
SM4_KeySchedule(MK,rk);
for(j=0;j<4;j++)
{
X[j]=(CipherText[j*4]<<24) |(CipherText[j*4+1]<<16)|
(CipherText[j*4+2]<<8)|(CipherText[j*4+3]);
}
for(i=0;i<32;i++)
{
tmp = X[i+1]^X[i+2]^X[i+3]^rk[31-i];
//nonlinear operation
buf= (SM4_Sbox[(tmp >> 24) & 0xFF]) << 24
|(SM4_Sbox[(tmp >> 16) & 0xFF]) << 16
|(SM4_Sbox[(tmp >> 8) & 0xFF]) << 8
|(SM4_Sbox[tmp & 0xFF]);
//linear operation
X[i+4]=X[i]^(buf^SM4_Rotl32((buf),2)^ SM4_Rotl32((buf),10)
^ SM4_Rotl32((buf),18)^ SM4_Rotl32((buf),24));
}
for(j=0;j<4;j++)
{
PlainText[4*j]=(X[35-j]>> 24)& 0xFF;
PlainText[4*j+1]=(X[35-j]>>16)& 0xFF;
PlainText[4*j+2]=(X[35-j]>> 8)& 0xFF;
 PlainText[4*j+3]=(X[35-j])& 0xFF;
}
}

int SM4_SelfCheck()
{
int i;
//Standard data
unsigned char key[16] =
{0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
unsigned char plain[16]=
{0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
unsigned char
cipher[16]={0x68,0x1e,0xdf,0x34,0xd2,0x06,0x96,0x5e,0x86,0xb3,0xe9,0x4f,0x53,0x6e,0x42,0x46}
;
unsigned char En_output[16];
unsigned char De_output[16];
SM4_Encrypt(key,plain,En_output);
SM4_Decrypt(key,cipher,De_output);
for(i=0;i<16;i++)
{
if ( (En_output[i]!=cipher[i]) | (De_output[i]!=plain[i]) )
{
  printf("Self-check error");
 return 1;
}
}
    printstate_sm4(En_output);
  //printf("Self-check success");
return 0;
}

void printstate_sm4(unsigned char * in)
{
    int i;
    for(i = 0; i < 16; i++) 
    {
        printf("%.2X", in[i]);
    }
    printf("\n");
}

unsigned char sm4Sbox(unsigned char inch)
{
	unsigned char *pTable = (unsigned char *)SboxTable;
	unsigned char retVal = (unsigned char)(pTable[inch]);
	return retVal;
}

mat_GF2 Ni_sm4Sbox(mat_GF2 & input)//1x8
{
    mat_GF2 temp;
    temp.SetDims(1,8);
    vec_GF2 dop;
    uint8_t inch;
    dop = reverse(input[0]);
    VectorCopy(temp[0],dop,8);
    inch= uint8FromMatGF2(transpose(temp));

	unsigned char *pTable = (unsigned char *)SM4_NiSbox;
	unsigned char retVal = (unsigned char)(pTable[inch]);
    mat_GF2 x = transpose(matGF2FromUint8(retVal));
    vec_GF2 y = x[0];
    y = reverse(y);
    mat_GF2 output;
    output.SetDims(1,8);
    VectorCopy(output[0],y,8);
	return output;
}



/*
mat_GF2 *Lj_matrix() //构造线性变换的矩阵形式
{
   // mat_GF2 L_Mat;
    mat_GF2* Lj = new mat_GF2[4];
    vec_GF2 l[4],temp;
    int m=0,n=0;
    //L_Mat.SetDims(32,32);
    Lj[0].SetDims(8,32);
    Lj[1].SetDims(8,32);
    Lj[2].SetDims(8,32);
    Lj[3].SetDims(8,32);
    for(int i=0;i<32;i++){
        temp.SetLength(0);

            for(int k=0;k<4;k++)
            {
                l[k] = vecGF2FromUint8(LMat[i][k]);
                append(temp,l[k]);
                //cout << temp;
            }
            //l[j] = vecGF2FromUint8(LMat[i][j]);
            if(m == 8)
            {
                n++;
                m = 0;
            }
            //VectorCopy(L_Mat[i],temp,32);
            VectorCopy(Lj[n][m],temp,32);
            m ++ ;

        }
    //cout << Lj[1];


    return Lj;

}
*/

mat_GF2 *Lj_matrix11() //构造线性变换的矩阵形式
{
   // mat_GF2 L_Mat;
    mat_GF2* Lj = new mat_GF2[4];
    vec_GF2 l[4],temp;
    int m=0,n=0;
    //L_Mat.SetDims(32,32);
    Lj[0].SetDims(8,32);
    Lj[1].SetDims(8,32);
    Lj[2].SetDims(8,32);
    Lj[3].SetDims(8,32);
    for(int i=0;i<32;i++){
        temp.SetLength(0);

            for(int k=0;k<4;k++)
            {
                l[k] = vecGF2FromUint8(LMat2[i][k]);
                append(temp,l[k]);
                //cout << temp;
            }
            //l[j] = vecGF2FromUint8(LMat[i][j]);
            if(m == 8)
            {
                n++;
                m = 0;
            }
            //VectorCopy(L_Mat[i],temp,32);
            VectorCopy(Lj[n][m],temp,32);
            m ++ ;

        }
    //cout << Lj[1];


    return Lj;

}

mat_GF2 * oi_matrix() //32x32
{
    mat_GF2* oi = new mat_GF2[36];
    for(int i=0;i<36;i++)
    {
        oi[i] = random_invertible_mat_GF2(32);
    }
    return oi;

}

mat_GF2 * Fi_matrix(const mat_GF2 * temp) //32x32
{
    mat_GF2* Fi = new mat_GF2[4];
    //mat_GF2 * temp =  oi_matrix();
    for(int i=0;i<4;i++)
    {
        Fi[i].SetDims(32,32);
        Fi[i] = temp[i];
    }
    return Fi;

}

mat_GF2 * Gi_matrix(const mat_GF2 * temp) //32x32 ,这里已经逆过了！！
{
    mat_GF2* Gi = new mat_GF2[4];
    //mat_GF2 * temp = oi_matrix();
    for(int i=0;i<4;i++)
    {
        Gi[i].SetDims(32,32);
        Gi[i] = inv(temp [32+i]);
    }
    return Gi;
}

mat_GF2 setup_matB1()
{
    vec_GF2 temp;
    mat_GF2 Mat_B1;
    Mat_B1.SetDims(8,8);
    for(int i =0;i<8;i++)
    {
        temp = vecGF2FromUint8(B1_mat[i][0]);
        VectorCopy(Mat_B1[i],temp,8);   
    }
    //cout << Mat_B1;
    return Mat_B1;
    
}

mat_GF2 setup_matB2()
{
    vec_GF2 temp;
    mat_GF2 Mat_B2;
    Mat_B2.SetDims(8,8);
    for(int i =0;i<8;i++)
    {
        temp = vecGF2FromUint8(B2_mat[i][0]);
        VectorCopy(Mat_B2[i],temp,8);   
    }
    //cout << Mat_B1;
    return Mat_B2;
    
}

mat_GF2 setup_matB3()
{
    vec_GF2 temp;
    mat_GF2 Mat_B3;
    Mat_B3.SetDims(8,8);
    for(int i =0;i<8;i++)
    {
        temp = vecGF2FromUint8(B3_mat[i][0]);
        VectorCopy(Mat_B3[i],temp,8);   
    }
    //cout << Mat_B1;
    return Mat_B3;
    
}



