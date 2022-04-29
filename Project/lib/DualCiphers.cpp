#include <NTL/GF2X.h>
#include <NTL/mat_GF2.h>
#include <iostream>

#include "DualCiphers.h"


using namespace std;
using namespace NTL;



mat_GF2 ** DualMat()
{
    GF2X a,b,c,d,e,f,g,h,i,x1,x2,x3;
    SetCoeff(a,0,1);
    SetCoeff(b,1,1);
    SetCoeff(c,2,1);
    SetCoeff(d,3,1);
    SetCoeff(e,4,1);
    SetCoeff(f,5,1);
    SetCoeff(g,6,1);
    SetCoeff(h,7,1);
    SetCoeff(i,8,1);
    x1 = a + b + d+ e + i; //GF28上的不可约多项式; x8+x4+x3+x+1;
    mat_GF2** M = new mat_GF2* [3];
   // mat_GF2 R1,R2;
    M[0] = InitialDualToEight(x1);
    //构建不可约多项式的本原元;
    GF2X benyuan_x2 ,benyuan_x3;
    benyuan_x2 = a + b;//x+1
    x2 = a + c + d + e + i;//x8+x4+x3+x2+1
    //cout << R1;
    M[1] = DualToEight (x2 ,benyuan_x2);
    //cout << M[1][0];
    x3 = a + c + e + f + g + h + i;//x8+x7+x6+x5+x4+x2+1
    benyuan_x3 = b + c + d + e + f;//x+x2+x3+x4+x5
    M[2] = DualToEight(x3,benyuan_x3);
	return M;



}

mat_GF2 * InitialDualToEight(const GF2X& x) //第一个AES的不可约多项式8个对偶矩阵生成
{
    mat_GF2* M = new mat_GF2[8];
    long n = 8;
    GF2X q[8];
    for(int i=0;i<8;i++) //标准基生成
    {
        SetCoeff(q[i],i,1);
    }
    M[0].SetDims(n,n);
    for(int i =0;i<8;i++)//标准矩阵生成
    {
        VectorCopy(M[0][i],q[i],n);
    }
    for(int i=1;i<8;i++) //构建其余7个对偶矩阵
    {
        M[i].SetDims(n,n);
        for(int j=0;j<8;j++)
        {
            q[j] = SqrMod(q[j],x);
            VectorCopy(M[i][j],q[j],n);
        }
        M[i] = transpose(M[i]);
    }
    return M;
}

mat_GF2 * DualToEight(const GF2X& x,const GF2X benyuan) //一般的8个对偶矩阵生成
{
    mat_GF2* M = new mat_GF2[8];
    long n = 8;
    GF2X q[8],r[8];
    mat_GF2 R;
    R.SetDims(8,8);
    for(int i=0;i<8;i++)
    {
        r[i] = PowerMod(benyuan,i,x) ;
        VectorCopy(R[i],r[i],8);
    }
    M[0].SetDims(n,n);
    M[0] = transpose(R);
    for(int i=1;i<8;i++) //构建其余7个对偶矩阵
    {
        M[i].SetDims(n,n);
        for(int j=0;j<8;j++)
        {
            if(i==1){
                q[j] = SqrMod(r[j],x);
                VectorCopy(M[i][j],q[j],n);
            }
            else{
            q[j] = SqrMod(q[j],x);
            VectorCopy(M[i][j],q[j],n);}
        }
        M[i] = transpose(M[i]);
    }
    return M;
}

