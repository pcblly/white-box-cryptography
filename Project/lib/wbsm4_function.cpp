#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>
#include <iostream>

#include "DualCiphers.h"
#include "sm4_function.h"
#include "stdlib.h"
#include "time.h"
#include "GF2E_function.h"
#include "wbsm4_function.h"

using namespace std;
using namespace NTL;

mat_GF2 DiagMatrix(mat_GF2 * Eij) //8x8生成对角矩阵32x32
{
    mat_GF2 DiagMat;
    vec_GF2 n8,n16,n24;
    n8.SetLength(8);
    n16.SetLength(16);
    n24.SetLength(24);
    int flag = 0,temp = 0;
    DiagMat.SetDims(32,32);
    for(int k=0 ;k<32;k++)
    {
        if(flag == 0)
        {
            VectorCopy(DiagMat[k],Eij[flag][temp],32);
            temp++;
            if(temp == 8)
            {   
                temp = 0;
                flag = 1;
            }
        }
        else if(flag == 1)
        {
            append(n8,Eij[flag][temp]) ;
            VectorCopy(DiagMat[k],n8,32);
            temp++;
            clear(n8);
            n8.SetLength(8);
            if(temp == 8)
            {
                temp = 0;
                flag = 2;
            }
        }
        else if (flag == 2)
        {
            append(n16,Eij[flag][temp]);
            VectorCopy(DiagMat[k],n16,32);
            temp ++;
            clear(n16);
            n16.SetLength(16);
            if(temp == 8)
            {
                temp =0;
                flag = 3;
            }
        }
        else
        {
            append(n24,Eij[flag][temp]);
            VectorCopy(DiagMat[k],n24,32);
            temp ++ ;
            clear(n24);
            n24.SetLength(24);
        }

    }
    return DiagMat;

}

mat_GF2 Xi123 (const mat_GF2& x1,const mat_GF2& x2,const mat_GF2& x3,int i,const mat_GF2 & Mi,const mat_GF2 &InvDiag_Mat,const mat_GF2* Oi )
//传入的xi要是32x1的列向量，出1x32的xor向量
{
    mat_GF2 xi1,xi2,xi3,xi_xor;

    xi1 = JiLian(Mi,transpose(transpose(x1) * inv(Oi[i+1]))) * InvDiag_Mat; //出1x32
    //对xi1进行处理
    xi2 = JiLian(Mi,transpose(transpose(x2) * inv(Oi[i+2]))) * InvDiag_Mat;
    xi3 = JiLian(Mi,transpose(transpose(x3) * inv(Oi[i+3]))) * InvDiag_Mat;
    xi_xor = xi1 + xi2 +xi3;//三个值异或

    return xi_xor;//出1x32

}

mat_GF2 QiMat(const mat_GF2& xi0,int i,const mat_GF2 * Oi) //对xi0进行处理32x1进，
{

    mat_GF2 Qi_Mat;
    Qi_Mat = transpose(xi0) * inv(Oi[i]) * Oi[i+4];
    return Qi_Mat; //1x32

}




mat_GF2 Tbox_ij(const uint8_t x0,int i,int j, mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi,const mat_GF2& aij,const mat_GF2& M_xing,const mat_GF2* Oi)
{

    mat_GF2 N,Ki,temp,x,y1,y2,Lj;
    x = matGF2FromUint8(x0); //8x1
    N = inv(Mi); //对偶矩阵的逆
    Ki = transpose(rki) * Mi ; //轮密钥的8bit分密钥 1x8

    Lj = Lj_matrix11()[j]; 

    y1 = Sibox(Ki + transpose(x) * Eij,Mi) * N * Lj * Oi[i+4]; //1x32
    y2 = JiLian(N,transpose(aij * M_xing)) * Oi[i+4]; //1x32

    temp = y1 + y2;

    return temp;//1x32

}

mat_GF2 Sibox(const mat_GF2 & x, const mat_GF2 & Mi)//进1x8
//结果出来一个1 x 8列向量
{
    mat_GF2 Si_end,a,N = inv(Mi),pi;
    uint8_t tmp;
    vec_GF2 dop;
    pi = x * N;
    dop = reverse(pi[0]);
    VectorCopy(pi[0],dop,8);
    tmp = uint8FromMatGF2(transpose(pi));

    a = matGF2FromUint8(sm4Sbox(tmp));
    dop = reverse(transpose(a)[0]);
    VectorCopy(pi[0],dop,8);

    Si_end = pi * Mi;
    //cout << Si_end;
    return Si_end;

}

mat_GF2 LinearTrans(const mat_GF2 & Mi ,int j ,const mat_GF2& x) //在这篇文章里面的线性变换,返回一个8x32矩阵M*ij矩阵
//进1x8的x,出8x32的M*
{
    mat_GF2 Lj,temp,M_32,M_Xing;
    temp.SetDims(32,1);
    Lj.SetDims(8,32);
    M_Xing.SetDims(8,32);
    Lj = Lj_matrix11()[j]; 
    temp = x * inv(Mi) * Lj;
    M_32 = JiLian(Mi,transpose(temp));
    M_Xing = transpose(x) * M_32; //用8x1 和 1x32的矩阵生成8x32的矩阵
    //cout << M_Xing;
    return M_Xing; //出8x32
}

mat_GF2 JiLian(const mat_GF2& Mi,const mat_GF2&x) //可以直接进入级联操作，得1x32的矩阵
//进32x1，出1x32
{
    mat_GF2 M_32;
    M_32.SetDims(1,32);
    uint32_t in;
    uint8_t out[4];
    mat_GF2 po[4];
    vec_GF2 pi;
    in = uint32FromMatGF2(x);
    *(uint32_t *)&out = in;
    for(int i =0;i<4;i++) //进入4个级联的Mi矩阵
    {
        po[i] = transpose(matGF2FromUint8(out[i])) * Mi; //1x8
        append(pi,po[i][0]);

    }
    VectorCopy(M_32[0],pi,32);
    return M_32;


}

mat_GF2 aij_part(const mat_GF2& aij,const mat_GF2& M_xing,const mat_GF2& Mi,int i,const mat_GF2* Oi)//出1x32
{

    mat_GF2 aiMat;
    aiMat = JiLian(inv(Mi),transpose(aij * M_xing)) * Oi[i+4];
    return aiMat; 

}




mat_GF2 roundF(const mat_GF2& xi0,const mat_GF2& xi1,const mat_GF2& xi2,const mat_GF2& xi3,int i,const mat_GF2& rki,const mat_GF2 * Oi)
//轮函数
{
    mat_GF2 **M = DualMat();
    //mat_GF2 * Oi = oi_matrix();
    mat_GF2 Eij[4],InvEij[4],aij[4],temp[4],M_xing[4],rkey[4];
    mat_GF2 Mi,xi123_xor,InvDiag_Mat,Tbox_x0,Tbox_x1,Tbox_x2,Tbox_x3,Tbox_xor
    ,ai0,ai1,ai2,ai3,ai_xor,Qi,end_xor;
    uint32_t in,in_key;
    uint8_t out[4],out_key[4];  
    int flag = 0;

    in_key = uint32FromMatGF2(rki);
    *(uint32_t *)&out_key = in_key;


    //每轮Mi生成
    srand((unsigned)time(NULL));
    int m = rand()%3;
    int n = rand()%8;
    Mi = M[m][n];


    //生成Ei矩阵数组
    for(int k=0;k<4;k++) //保留的全部Eij的8x8矩阵，以便于后面约掉
    {
        Eij[k]=random_invertible_mat_GF2(8);
    }
    for(int k=0;k<4;k++)
    {
        InvEij[k] = inv(Eij[k]);
    }
    InvDiag_Mat = DiagMatrix(InvEij); //Ei(-1)的生成


    //生成随机的gf28元素aij,1x8
    for(int k=0;k<4;k++)
    {
        aij[k] = random_mat_GF2(1,8);
    }
    

    //xi123的xor完成，出来一个1x32的矩阵
    xi123_xor = Xi123(xi1,xi2,xi3,i,Mi,InvDiag_Mat,Oi); 
    //cout << xi123_xor;

    //将1x32bit的矩阵分为4份进入Tbox
    in = uint32FromMatGF2(transpose(xi123_xor));
    *(uint32_t *)&out = in;

    for(int k=0;k<4;k++)
    {
        temp[k] = matGF2FromUint8(out[k]);
        rkey[k] = matGF2FromUint8(out_key[k]);

    } //temp 8x1

    for(int k=0;k<4;k++)
    {
        M_xing[k] = LinearTrans(Mi,k,transpose(temp[k]));
        //cout << M_xing[k];
    }


    //进入Tbox的分box
    //这里的密钥需要分块进取，还没改，暂时是先写成这样
    Tbox_x0 = Tbox_ij(out[0],i,0,rkey[0],Eij[0],Mi,aij[0],M_xing[0],Oi);
    Tbox_x1 = Tbox_ij(out[1],i,1,rkey[1],Eij[1],Mi,aij[1],M_xing[1],Oi);
    Tbox_x2 = Tbox_ij(out[2],i,2,rkey[2],Eij[2],Mi,aij[2],M_xing[2],Oi);
    Tbox_x3 = Tbox_ij(out[3],i,3,rkey[3],Eij[3],Mi,aij[3],M_xing[3],Oi);
    //cout << Tbox_x0;

    Tbox_xor = Tbox_x0 + Tbox_x1 + Tbox_x2 + Tbox_x3;//出来1x32矩阵


    //ai生成
    ai0 = aij_part(aij[0],M_xing[0],Mi,i,Oi);
    ai1 = aij_part(aij[1],M_xing[1],Mi,i,Oi);
    ai2 = aij_part(aij[2],M_xing[2],Mi,i,Oi);
    ai3 = aij_part(aij[3],M_xing[3],Mi,i,Oi);

    ai_xor = ai0 + ai1 + ai2 + ai3; //1x32

    //Qi生成
    Qi = QiMat(xi0,i,Oi); //1x32

    //最后xor
    end_xor = Qi + ai_xor + Tbox_xor; //1x32

    return end_xor; //出来1x32

}

void printstate(unsigned char * in)
{
    int i;
    for(i = 3; i >= 0; i--) 
    {
        printf("%.2X", in[i]);
    }
    //printf("\n");
}



//重写有查找表的部分
mat_GF2 * Produce_x() //产生需要0-255的输入
{
    uint8_t a;
    mat_GF2 temp1;
    mat_GF2* array_a = new mat_GF2[256];
    for(int i=0;i<256;i++)
    {
        a = i;
        temp1 = matGF2FromUint8_reverse(a);
        array_a[i] = temp1;//1x8

    }
    return array_a;
}


mat_GF2 Tbox_ceshi(const mat_GF2& x,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi,const mat_GF2* Oi)
//x 1x8
{

    mat_GF2 N,Ki,temp,y1,Lj;
   // x = matGF2FromUint8(x0); //8x1
    N = inv(Mi); //对偶矩阵的逆
    Ki = transpose(rki) * Mi ; //轮密钥的8bit分密钥 1x8

    Lj = Lj_matrix11()[j]; 

    y1 = Sibox(Ki + x * Eij,Mi) * N * Lj * Oi[i+4]; //1x32
   // y2 = JiLian(N,transpose(aij * M_xing)) * Oi[i+4]; //1x32

 //   temp = y1 + y2;

    return y1;//1x32

}



static mat_GF2 Tbox1[32][4][256];
static mat_GF2 TBox_attack[32][4][256];
static mat_GF2 InvDiag_Mat[32];
static mat_GF2 Eij_inv[32][4];
static mat_GF2 Eij_xing[32][4];
static mat_GF2 invEij_xing32Mat[32];
static mat_GF2 Mi_box[32];
static mat_GF2 invEij_xing[32][4];
static mat_GF2 Sbox_xing[32][4][256];

void wbsm4_gen(const mat_GF2 * Oi,unsigned char MK[])
{
    mat_GF2 **M = DualMat();
   // cout << M[1][2]<<endl;
    mat_GF2 *x = Produce_x() ;
    uint32_t rk[32];
    vec_GF2 rktemp[32];
    mat_GF2 rki[32],rkii[32],Eij[32][4],key[32][4];
    //mat_GF2 InvDiag_Mat[32];
    //static mat_GF2 Tbox[32][4][256];
    

    uint32_t in_key;
    uint8_t out_key[4];  


    SM4_KeySchedule(MK,rk);
    //cout << transpose(matGF2FromUint32(rk[0])) <<endl;

    for(int k=0;k<32;k++)
    {
        rkii[k].SetDims(1,32);
        rktemp[k] = reverse(transpose(matGF2FromUint32(rk[k]))[0]);
        rkii[k][0] = rktemp[k];
        rki[k] = transpose(rkii[k]);
    }//rki 32x1

   // cout << transpose(rki[0]) << endl;

    for(int i=0;i<32;i++)
    {
        in_key = uint32FromMatGF2(rki[i]);
        *(uint32_t *)&out_key = in_key;
        for(int j=0;j<4;j++)
        {
            key[i][j] = matGF2FromUint8(out_key[j]) ;//8x1
        }
    }

    //cout << key[2][1] <<endl;
    

    for(int i =0;i<32;i++)
    {
        srand((unsigned)time(NULL));
        int m = rand()%3;
        int n = rand()%8;
        Mi_box[i] = M[m][n];

        for(int j=0;j<4;j++)
        {
            Eij[i][j] = random_invertible_mat_GF2(8);
            Eij_inv[i][j] = inv(Eij[i][j]);

            //aij_box[i][j] = random_mat_GF2(1,8);

        }
    }
   // cout << Mi_box[2] << endl;

   for(int i =0 ;i < 32 ; i++)
   {
       for(int j=0;j<4;j++)
       {
           Eij_xing[i][j] = Eij[i][j] * inv(Mi_box[i]);
           invEij_xing[i][j] = inv(Eij_xing[i][j]);
       }
       invEij_xing32Mat[i] = DiagMatrix(invEij_xing[i]);
       
   }


    for(int i=0;i<32;i++)
    {
        InvDiag_Mat[i] = DiagMatrix(Eij_inv[i]);
    }

    //cout << InvDiag_Mat[31] <<endl;

    for(int i=0;i<32;i++)
    {
        for(int j=0;j<4;j++)
        {
            for(int k=0;k<256;k++)
            {
                Tbox1[i][j][k] = Tbox_ceshi(x[k],i,j,key[i][j],Eij[i][j],Mi_box[i],Oi);
                TBox_attack[i][j][k] = Tbox_change(x[k],i,j,key[i][j],Eij[i][j],Mi_box[i]);
                Sbox_xing[i][j][k] = Sbox_change(x[k],i,j,key[i][j],Eij[i][j],Mi_box[i]);
                //cout <<TBox_attack[i][j][k] <<endl;

            }


        }
    }

    //cout << Mi_box[0] <<endl;

}


mat_GF2 roundF_ceshi(const mat_GF2& xi0,const mat_GF2& xi1,const mat_GF2& xi2,const mat_GF2& xi3,int i,const mat_GF2 * Oi)
{
    mat_GF2 xi123_xor,xi123_inv,Tbox_output0,Tbox_output1,Tbox_output2,Tbox_output3,Tbox_xor,Qi,end_xor,doi[4];
    vec_GF2 temp,dop[4];
    uint32_t in;
    uint8_t out[4],xiba[4];

   // cout << Mi_box[1][2] << endl;
    xi123_xor = Xi123(xi1,xi2,xi3,i,Mi_box[i],InvDiag_Mat[i],Oi); 

    //cout << Mi_box[0] <<endl;

   //cout << xi123_xor << endl;
    //temp = reverse(xi123_xor[0]);
    //VectorCopy(xi123_xor[0],temp,32);
    //cout << xi123_xor << endl;

    in = uint32FromMatGF2(transpose(xi123_xor));
    *(uint32_t *)&out = in;
    for(int k =0;k<4;k++)
    {
        dop[k] = transpose(matGF2FromUint8(out[k]))[0];
        dop[k] = reverse(dop[k]);
        doi[k].SetDims(1,8);
        VectorCopy(doi[k][0],dop[k],8); 
    }
    for(int k = 0;k<4;k++)
    {
        xiba[k] = uint8FromMatGF2(transpose(doi[k]));
    }



  //  cout << transpose(matGF2FromUint8(out[0]))<< endl;

    Tbox_output0 = Tbox1[i][0][xiba[0]];
    Tbox_output1 = Tbox1[i][1][xiba[1]];
    Tbox_output2 = Tbox1[i][2][xiba[2]];
    Tbox_output3 = Tbox1[i][3][xiba[3]];

    Tbox_xor = Tbox_output0 + Tbox_output1 + Tbox_output2 + Tbox_output3;

    Qi = QiMat(xi0,i,Oi); //1x32

    //最后xor
    end_xor = Qi + Tbox_xor; //1x32

    return end_xor;


}


void wbsm4_Encrypt(unsigned char CipherText[],const mat_GF2 * Oi)
{
    vec_GF2 orig[16],vec_pio,pio[4];
    mat_GF2 xx[4];
    mat_GF2 temp,x0,x1,x2,x3,input_Mat,y0,y1,y2,y3;
    int flag = 0;

    uint32_t out_ciph[4];
    uint8_t in_ciph[4][4];


    //把输入的明文串起来为1x128的矩阵
    //同时生成1x32的明文分矩阵
    for(int k =0;k<16;k++)
    {
        orig[k] = vecGF2FromUint8(CipherText[k]);
        
        if (k!=0 && k%4 == 0)
        {
            xx[flag].SetDims(1,32);
            VectorCopy(xx[flag][0],vec_pio,32);
            flag ++ ;
            clear(vec_pio);
            vec_pio.SetLength(0);
 
        }
            append(vec_pio,orig[k]);
            if(k == 15)
            {
                xx[3].SetDims(1,32);
                VectorCopy(xx[3][0],vec_pio,32);
            }
        
    }

    //cout << xx[0];

    x0 = transpose(xx[0] * Oi[0]); //32x1
    x1 = transpose(xx[1] * Oi[1]); 
    x2 = transpose(xx[2] * Oi[2]); 
    x3 = transpose(xx[3] * Oi[3]); 

    for(int k=0;k<32;k++)
    {
        temp = roundF_ceshi(x0,x1,x2,x3,k,Oi);
        x0 = x1;
        x1 = x2;
        x2 = x3;
        x3 = transpose(temp);
    }

    y0 = transpose(x3) * inv(Oi[35]);
    y1 = transpose(x2) * inv(Oi[34]);
    y2 = transpose(x1) * inv(Oi[33]);
    y3 = transpose(x0) * inv(Oi[32]);
    
//为了输出十六进制
    pio[0] = reverse(y0[0]);
    pio[1] = reverse(y1[0]);
    pio[2] = reverse(y2[0]);
    pio[3] = reverse(y3[0]);

    VectorCopy(y0[0],pio[0],32);
    VectorCopy(y1[0],pio[1],32);
    VectorCopy(y2[0],pio[2],32);
    VectorCopy(y3[0],pio[3],32);

    out_ciph[0] = uint32FromMatGF2(transpose(y0));
    out_ciph[1] = uint32FromMatGF2(transpose(y1));
    out_ciph[2] = uint32FromMatGF2(transpose(y2));
    out_ciph[3] = uint32FromMatGF2(transpose(y3));

    *(uint32_t *)&in_ciph[0] = out_ciph[0];
    *(uint32_t *)&in_ciph[1] = out_ciph[1];
    *(uint32_t *)&in_ciph[2] = out_ciph[2];
    *(uint32_t *)&in_ciph[3] = out_ciph[3];

    for(int k =0 ;k<4;k++)
    {
        printstate(in_ciph[k]);
    }

}






//碰撞攻击函数、测试等


mat_GF2 Tbox_change(const mat_GF2& x,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi)
//x 1x8
{

    mat_GF2 N,Ki,temp,y1,Lj;
   // x = matGF2FromUint8(x0); //8x1
    N = inv(Mi); //对偶矩阵的逆
    Ki = transpose(rki) * Mi ; //轮密钥的8bit分密钥 1x8

    Lj = Lj_matrix11()[j]; 

    y1 = Sibox(Ki + x * Eij,Mi) * N * Lj ; //1x32


    return y1;//1x32

}

mat_GF2 Sbox_change(const mat_GF2& x,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi)
{
    mat_GF2 N,Ki,temp,y1,Lj;
   // x = matGF2FromUint8(x0); //8x1
    N = inv(Mi); //对偶矩阵的逆
    Ki = transpose(rki) * Mi ; //轮密钥的8bit分密钥 1x8

    //Lj = Lj_matrix11()[j]; 

    y1 = Sibox(Ki + x * Eij,Mi) * N  ; //1x8


    return y1;//1x8

    
}



mat_GF2 * collision_function(int i,const mat_GF2& x0,const mat_GF2& x1,const mat_GF2& x2,const mat_GF2& x3)
//x0 1x8
{
    
    mat_GF2 Tbox_out0,Tbox_out1,Tbox_out2,Tbox_out3,dop[4],output0,output1,compare_jieguo[4],eij;
    vec_GF2 dot[4];
    uint8_t out11[4],out1[4];
    //int n=4;
    mat_GF2 *out_mat = new mat_GF2[4];
    uint32_t in;

    dot[0] = x0[0];
    dot[1] = x1[0];
    dot[2] = x2[0];
    dot[3] = x3[0];

    //cout << "shiji" << invEij_xing[1][0]<<endl;
/*
    mat_GF2 * x = Produce_x();
    eij = x[2] * invEij_xing[1][0];
    cout << "eij" << eij <<endl;
    */
    

    for(int k=0;k<4;k++)
    {
        dot[k] = reverse(dot[k]);
        dop[k].SetDims(1,8);
        VectorCopy(dop[k][0],dot[k],8);
       // cout << transpose(dop[k]) << endl;
        //cout << uint8FromMatGF2(transpose(dop[k])) <<endl;
        out11[k] = uint8FromMatGF2(transpose(dop[k])) ;
        //cout << out[k] << endl;
    }

    Tbox_out0 = TBox_attack[i][0][out11[0]];
    Tbox_out1 = TBox_attack[i][1][out11[1]];
    Tbox_out2 = TBox_attack[i][2][out11[2]];
    Tbox_out3 = TBox_attack[i][3][out11[3]];

    output0 = Tbox_out0 + Tbox_out1 + Tbox_out2 + Tbox_out3; //1x32

    //cout << "output0"<<output0 << endl;
    
    //cout << Sbox_xing[1][0][2] << endl;

    output1 = output0 * invEij_xing32Mat[i+1]; //1x32
   // cout << output1 << endl;

    in = uint32FromMatGF2(transpose(output1));
   // cout << in << endl;
    *(uint32_t *)&out1 = in;
    for(int k=0;k<4;k++)
    {
        out_mat[k] = transpose(matGF2FromUint8(out1[k])) ;
    }




    return out_mat;
    
 


    



}



mat_GF2 * ui0_xing(int i)
{
    mat_GF2* ui_ci = new mat_GF2[255];
    //mat_GF2 output[255];
    for(int t =1;t<256;t++)
    {
        ui_ci[t-1] = Sbox_xing[i][0][0] + Sbox_xing[i][0][t];

    }
 
    return ui_ci;
 
}

mat_GF2 * ui1_xing(int i)
{
    mat_GF2* ui_ci = new mat_GF2[255];
    //mat_GF2 output[255];
    for(int t =1;t<256;t++)
    {
        ui_ci[t-1] = Sbox_xing[i][1][0] + Sbox_xing[i][1][t];

    }
  
    return ui_ci;
  
}

mat_GF2 * ui2_xing(int i)
{
    mat_GF2* ui_ci = new mat_GF2[255];
    //mat_GF2 output[255];
    for(int t =1;t<256;t++)
    {
        ui_ci[t-1] = Sbox_xing[i][2][0] + Sbox_xing[i][2][t];

    }
 
    return ui_ci;
    
}

mat_GF2 * ui3_xing(int i)
{
    mat_GF2* ui_ci = new mat_GF2[255];
    for(int t =1;t<256;t++)
    {
        ui_ci[t-1] = Sbox_xing[i][3][0] + Sbox_xing[i][3][t];

    }
  
    return ui_ci;
 
}

mat_GF2 **  invEij_ceshi()
{
    mat_GF2** ga = new mat_GF2* [32];
   // mat_GF2 (*ga)[4] = new mat_GF2 [32][4];

   for(int k =0;k<32;k++)
   {
       ga[k] = new mat_GF2 [4];

   }

    for(int k=0;k<32;k++)
    {
        for(int t =0;t<4;t++)
        {
            ga[k][t] = invEij_xing[k][t];

        }
    }

   // cout << ga[0][0] <<endl;

    return ga;
}










