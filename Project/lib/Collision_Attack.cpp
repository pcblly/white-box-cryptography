#include <NTL/GF2X.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "DualCiphers.h"
#include "sm4_function.h"
#include "stdlib.h"
#include "time.h"
#include "GF2E_function.h"
#include "wbsm4_function.h"

using namespace std;
using namespace NTL;



mat_GF2 * Produce_temp() //产生需要1-255的输入
{
    uint8_t a;
    mat_GF2 temp1;
    mat_GF2* array_a = new mat_GF2[255];
    for(int i=1;i<256;i++)
    {
        a = i;
        temp1 = matGF2FromUint8_reverse(a);
        array_a[i-1] = temp1;//1x8

    }
    return array_a;
}

void coll()
{
    unsigned char key[16] =
    {0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
    unsigned char plain[16]=
    {0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
    cout << "输入密钥为:"<< endl ;
    printstate_sm4(key);
    mat_GF2 * Oi = oi_matrix(); 
    wbsm4_gen(Oi,key);
    mat_GF2* temp = Produce_temp(); 
   // cout << temp[0] << endl;
    uint8_t x = 0;
    uint8_t y = 1;
    mat_GF2 x1 = transpose(matGF2FromUint8(x));
    //cout << x1 << endl;
    mat_GF2 y1 = matGF2FromUint8_reverse(y);
    //cout << y1 << endl;
    /*
    mat_GF2 *out_array1,*out_array2;

   out_array1 = collision_function(0,y1,x1,x1,x1);
   //cout << out_array1[0] << endl;
   
   for(int k=0;k<255;k++)
   {
       out_array2 = collision_function(0,x1,temp[k],x1,x1);
       //cout << out_array2[0] << endl;
       if(out_array1[0] == out_array2[0])
       {
           cout << temp[k] << endl;
           break;
       }
       delete [] out_array2;
   }
   */
}

static mat_GF2 ui0_ci_test[32][255];
static mat_GF2 ui1_ci_test[32][255];
static mat_GF2 ui2_ci_test[32][255];
static mat_GF2 ui3_ci_test[32][255];



void coll_test()
{
    unsigned char key[16] =
    {0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
    unsigned char plain[16]=
    {0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10};
    cout << "输入密钥为:"<< endl ;
    printstate_sm4(key);
    mat_GF2 * Oi = oi_matrix(); 
    wbsm4_gen(Oi,key);
    mat_GF2* temp = Produce_temp(); 

    uint8_t x = 0;
    mat_GF2 x1 = transpose(matGF2FromUint8(x));//00000000
    //cout << x1 << endl;

    
    mat_GF2 *out_array1,*out_array2;
    mat_GF2 B1,B2,B3,dop;


    B1 = setup_matB1();
    B2 = setup_matB2();
    B3 = setup_matB3();



    for(int t=0;t<32;t++)
    {
        for(int m=0;m<255;m++)
        {
            out_array1 = collision_function_ceshi(t,temp[m],x1,x1,x1);
            out_array2 = collision_function_ceshi(t,x1,x1,x1,x1);
            dop = out_array1[0] + out_array2[0];//1x8
            ui0_ci_test[t][m] = dop * inv(B1);
        }
    }

    for(int t=0;t<32;t++)
    {
        for(int m=0;m<255;m++)
        {
            out_array1 = collision_function_ceshi(t,x1,temp[m],x1,x1);
            out_array2 = collision_function_ceshi(t,x1,x1,x1,x1);
            dop = out_array1[1] + out_array2[1];//1x8
            ui1_ci_test[t][m] = dop * inv(B1);
        }
    }

    for(int t=0;t<32;t++)
    {
        for(int m=0;m<255;m++)
        {
            out_array1 = collision_function_ceshi(t,x1,x1,temp[m],x1);
            out_array2 = collision_function_ceshi(t,x1,x1,x1,x1);
            dop = out_array1[2] + out_array2[2];//1x8
            ui2_ci_test[t][m] = dop * inv(B1);
        }
    }

    for(int t=0;t<32;t++)
    {
        for(int m=0;m<255;m++)
        {
            out_array1 = collision_function_ceshi(t,x1,x1,x1,temp[m]);
            out_array2 = collision_function_ceshi(t,x1,x1,x1,x1);
            dop = out_array1[3] + out_array2[3];//1x8
            ui3_ci_test[t][m] = dop * inv(B1);
        }
    }




   
}




static mat_GF2 Sbox_huifu[32][4][256];
static mat_GF2 NiSbox_huifu[32][4][256];


void huifuSBOX()
{

    mat_GF2 u1_ci,x1,dfp;
    static mat_GF2 um_ceshi0[32][256][255],um_ceshi1[32][256][255],um_ceshi2[32][256][255],um_ceshi3[32][256][255];
    vec_GF2 dop;

    mat_GF2* x =  Produce_x(); //产生0-255的输入，猜测u0可能等于0-255
    mat_GF2* temp = Produce_temp();//产生1-255的输入
    x1 = temp[0];
    dfp.SetDims(1,8);
    //cout << "ceshi" <<endl;

    //coll();
    coll_test();
   


    for(int m=0;m<32;m++)
    {
        for(int k=0;k<256;k++)
        {
            for(int t=0;t<255;t++)
            {
                um_ceshi0[m][k][t] = ui0_ci_test[m][t] + x[k];//1x8
                um_ceshi1[m][k][t] = ui1_ci_test[m][t] + x[k];
                um_ceshi2[m][k][t] = ui2_ci_test[m][t] + x[k];
                um_ceshi3[m][k][t] = ui3_ci_test[m][t] + x[k];

            }
        }
    }


    mat_GF2 n1,n2,n3,outi10,outi11,outi12,outi13,outi20,outi21,outi22,outi23,out3[255];
    uint8_t m1,m2;
    int flag,tt=-1,niuwa=-1;


    for(int m =0;m<32;m++)
    {
        for(int k=0;k<256;k++)
        {
            for(int t = 0 ;t<255;t++)
            {
                outi10 = Ni_sm4Sbox(um_ceshi0[m][k][t]);
                m1 = t + 1;
                n1 = matGF2FromUint8_reverse(m1);
                n1 = n1+x1;
                if(n1 == x[0])
                {
                    outi20 = Ni_sm4Sbox(x[k]);
                }
                else{
                    dop = reverse(n1[0]);
                    VectorCopy(dfp[0],dop,8);
                    m2 = uint8FromMatGF2(transpose(dfp));
                    flag = m2 - 1;
                    outi20 = Ni_sm4Sbox(um_ceshi0[m][k][flag]);

                }
                out3[t] = outi10 + outi20;
            
                if(t == 0 || t%2 == 1)
                {
                    continue;
                }
                else
                {
                    if(out3[t] == out3[t-2])
                    {
                    tt = t;
                    continue;
                    }
                    else{
                        break;
                    }
                }


            }
            if(tt>100)
            {
                niuwa = k;
                n3 = matGF2FromUint8_reverse(niuwa);
                Sbox_huifu[m][0][0] = n3;
                for(int k=0 ; k<255 ; k++)
                {
                    Sbox_huifu[m][0][k+1] = um_ceshi0[m][niuwa][k];
                }
               
                tt = -1;
                break;
            }
        }
    }

    for(int m =0;m<32;m++)
    {
        for(int k=0;k<256;k++)
        {
            for(int t = 0 ;t<255;t++)
            {
                outi11 = Ni_sm4Sbox(um_ceshi1[m][k][t]);
                m1 = t + 1;
                n1 = matGF2FromUint8_reverse(m1);
                n1 = n1+x1;
                if(n1 == x[0])
                {
                    outi21 = Ni_sm4Sbox(x[k]);
                }
                else{
                    dop = reverse(n1[0]);
                    VectorCopy(dfp[0],dop,8);
                    m2 = uint8FromMatGF2(transpose(dfp));
                    flag = m2 - 1;
                    outi21 = Ni_sm4Sbox(um_ceshi1[m][k][flag]);

                }
                out3[t] = outi11 + outi21;
            
                if(t == 0 || t%2 == 1)
                {
                    continue;
                }
                else
                {
                    if(out3[t] == out3[t-2])
                    {
                    tt = t;
                    continue;
                    }
                    else{
                        break;
                    }
                }


            }
            if(tt>100)
            {
                niuwa = k;
                n3 = matGF2FromUint8_reverse(niuwa);
                Sbox_huifu[m][1][0] = n3;
                for(int k=0 ; k<255 ; k++)
                {
                    Sbox_huifu[m][1][k+1] = um_ceshi1[m][niuwa][k];
                }
                tt = -1;
                break;
            }
        }
    }

     for(int m =0;m<32;m++)
    {
        for(int k=0;k<256;k++)
        {
            for(int t = 0 ;t<255;t++)
            {
                outi12 = Ni_sm4Sbox(um_ceshi2[m][k][t]);
                m1 = t + 1;
                n1 = matGF2FromUint8_reverse(m1);
                n1 = n1+x1;
                if(n1 == x[0])
                {
                    outi22 = Ni_sm4Sbox(x[k]);
    
                }
                else{
                    dop = reverse(n1[0]);
                    VectorCopy(dfp[0],dop,8);
                    m2 = uint8FromMatGF2(transpose(dfp));
                    flag = m2 - 1;
                    outi22 = Ni_sm4Sbox(um_ceshi2[m][k][flag]);

                }
                out3[t] = outi12 + outi22;

            
                if(t == 0 || t%2 == 1)
                {
                    continue;
                }
                else
                {
                    if(out3[t] == out3[t-2])
                    {
                    tt = t;
                    continue;
                    }
                    else{
                        break;
                    }
                }


            }
            if(tt>100)
            {
                niuwa = k;
                n3 = matGF2FromUint8_reverse(niuwa);
                Sbox_huifu[m][2][0] = n3;
                for(int k=0 ; k<255 ; k++)
                {
                    Sbox_huifu[m][2][k+1] = um_ceshi2[m][niuwa][k];
                }
                tt = -1;
                break;
            }
        }
    }

     for(int m =0;m<32;m++)
    {
        for(int k=0;k<256;k++)
        {
            for(int t = 0 ;t<255;t++)
            {
                outi13 = Ni_sm4Sbox(um_ceshi3[m][k][t]);
                m1 = t + 1;
                n1 = matGF2FromUint8_reverse(m1);
                n1 = n1+x1;
           
                if(n1 == x[0])
                {
                    outi23 = Ni_sm4Sbox(x[k]);
                }
                else{
                    dop = reverse(n1[0]);
                    VectorCopy(dfp[0],dop,8);
                    m2 = uint8FromMatGF2(transpose(dfp));
                    flag = m2 - 1;
                    outi23 = Ni_sm4Sbox(um_ceshi3[m][k][flag]);

                }
                out3[t] = outi13 + outi23;
    
            
                if(t == 0 || t%2 == 1)
                {
                    continue;
                }
                else
                {
                    if(out3[t] == out3[t-2])
                    {
                    tt = t;
                    continue;
                    }
                    else{
                        break;
                    }
                }


            }
            if(tt>100)
            {
                niuwa = k;
                n3 = matGF2FromUint8_reverse(niuwa);
                Sbox_huifu[m][3][0] = n3;
                for(int k=0 ; k<255 ; k++)
                {
                    Sbox_huifu[m][3][k+1] = um_ceshi3[m][niuwa][k];
                }
     
                tt = -1;
                break;
            }
        }
    }


}

void huifuSboxNi()
{
    mat_GF2 temp0,temp1;
    vec_GF2 ga;
    int dop;

    temp1.SetDims(1,8);

    for(int t = 0 ;t<32;t++)
    {
        for(int k = 0 ;k<256;k++)
        {
            temp0 = Sbox_huifu[t][0][k];
            ga = reverse(temp0[0]);
            VectorCopy(temp1[0],ga,8);
            dop = uint8FromMatGF2(transpose(temp1));
            NiSbox_huifu[t][0][dop] = matGF2FromUint8_reverse(k);
        }

    }

    for(int t = 0 ;t<32;t++)
    {
        for(int k = 0 ;k<256;k++)
        {
            temp0 = Sbox_huifu[t][1][k];
            ga = reverse(temp0[0]);
            VectorCopy(temp1[0],ga,8);
            dop = uint8FromMatGF2(transpose(temp1));
            NiSbox_huifu[t][1][dop] = matGF2FromUint8_reverse(k);
        }
    }

    for(int t = 0 ;t<32;t++)
    {
        for(int k = 0 ;k<256;k++)
        {
            temp0 = Sbox_huifu[t][2][k];
            ga = reverse(temp0[0]);
            VectorCopy(temp1[0],ga,8);
            dop = uint8FromMatGF2(transpose(temp1));
            NiSbox_huifu[t][2][dop] = matGF2FromUint8_reverse(k);
        }
    }

    for(int t = 0 ;t<32;t++)
    {
        for(int k = 0 ;k<256;k++)
        {
            temp0 = Sbox_huifu[t][3][k];
            ga = reverse(temp0[0]);
            VectorCopy(temp1[0],ga,8);
            dop = uint8FromMatGF2(transpose(temp1));
            NiSbox_huifu[t][3][dop] = matGF2FromUint8_reverse(k);
        }
    }

}







static mat_GF2 invEijxing_huifu[31][4];
static mat_GF2 B1;
static mat_GF2 B2;
static mat_GF2 B3;
static mat_GF2 * x;


void huifu_invEij()
{
    B1 = setup_matB1();
    B2 = setup_matB2();
    B3 = setup_matB3();
    x = Produce_x();

    mat_GF2* temp = Produce_temp();//产生1-255的输入
    mat_GF2 x0,x1,op,tmp,out,mat_temp;
    mat_GF2 Ni_B1,Ni_B2,Ni_B3;
    vec_GF2 dop;
    uint8_t flag = -1;
    int n1;
    
    tmp.SetDims(1,8);
    mat_temp.SetDims(8,8);
    x0 = x[0]; //00000000
    x1 = x[1]; //00000001
    //cout << "x[2]"<<x[2] << endl;

    Ni_B1 = inv(B1);
    Ni_B2 = inv(B2);
    Ni_B3 = inv(B3);

    for(int t=0;t<31;t++)
    {
        clear(mat_temp);
        for(int k=0;k<8;k++) //x是从00000001开始
        {   

        if(k==0)
        {
            op = temp[0] * Ni_B1;
        }
        else{
            op = x[1u<<k] * Ni_B1;//1x8
        }
        
        dop = reverse(op[0]);
        VectorCopy(tmp[0],dop,8);//1x8
        flag = uint8FromMatGF2(transpose(tmp));
        out = collision_function(t,NiSbox_huifu[t][0][flag],NiSbox_huifu[t][1][0],NiSbox_huifu[t][2][0],NiSbox_huifu[t][3][0])[0];
        //out3[k] = out;
        VectorCopy(mat_temp[7-k],out[0],8);

        }
        invEijxing_huifu[t][0] = mat_temp;

    }

    for(int t=0;t<31;t++)
    {
        clear(mat_temp);
        for(int k=0;k<8;k++) //x是从00000001开始
        {   

        if(k==0)
        {
            op = temp[0] * Ni_B3;
        }
        else{
            op = x[1u<<k] * Ni_B3;//1x8
        }
        
        dop = reverse(op[0]);
        VectorCopy(tmp[0],dop,8);//1x8
        flag = uint8FromMatGF2(transpose(tmp));
        out = collision_function(t,NiSbox_huifu[t][0][flag],NiSbox_huifu[t][1][0],NiSbox_huifu[t][2][0],NiSbox_huifu[t][3][0])[1];
        //out3[k] = out;
        VectorCopy(mat_temp[7-k],out[0],8);

        }
        invEijxing_huifu[t][1] = mat_temp;

    }

    for(int t=0;t<31;t++)
    {
        clear(mat_temp);
        for(int k=0;k<8;k++) //x是从00000001开始
        {   

        if(k==0)
        {
            op = temp[0] * Ni_B2;
        }
        else{
            op = x[1u<<k] * Ni_B2;//1x8
        }
        
        dop = reverse(op[0]);
        VectorCopy(tmp[0],dop,8);//1x8
        flag = uint8FromMatGF2(transpose(tmp));
        out = collision_function(t,NiSbox_huifu[t][0][flag],NiSbox_huifu[t][1][0],NiSbox_huifu[t][2][0],NiSbox_huifu[t][3][0])[2];
        //out3[k] = out;
        VectorCopy(mat_temp[7-k],out[0],8);

        }
        invEijxing_huifu[t][2] = mat_temp;

    }

    for(int t=0;t<31;t++)
    {
        clear(mat_temp);
        for(int k=0;k<8;k++) //x是从00000001开始
        {   

        if(k==0)
        {
            op = temp[0] * Ni_B2;
        }
        else{
            op = x[1u<<k] * Ni_B2;//1x8
        }
        
        dop = reverse(op[0]);
        VectorCopy(tmp[0],dop,8);//1x8
        flag = uint8FromMatGF2(transpose(tmp));
        out = collision_function(t,NiSbox_huifu[t][0][flag],NiSbox_huifu[t][1][0],NiSbox_huifu[t][2][0],NiSbox_huifu[t][3][0])[3];
        //out3[k] = out;
        VectorCopy(mat_temp[7-k],out[0],8);

        }
        invEijxing_huifu[t][3] = mat_temp;

    }





}




//static mat_GF2** invEij_temp;
static mat_GF2 key_found[5][4];
static mat_GF2 key32_found[5];
static uint32_t rki[5];

void key1_function()
{
    x = Produce_x();//产生0-255的输入
    mat_GF2 x0,x1,n1,ceshi0,ceshi1,temp0,temp1,output0;
    mat_GF2 out3[256];
    vec_GF2 vec_ga;
    int tt=-1;

   // invEij_temp = invEij_ceshi() ;


    x1 = x[1] ;//00000001
    x0 = x[0];//00000000


    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[0][0] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[0][0];
            temp0 = collision_function(1,ceshi0,x0,x0,x0)[0];
            temp1 = collision_function(1,ceshi1,x0,x0,x0)[0];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[1][0] = x[k];
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[0][1] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[0][1] ;
            temp0 = collision_function(1,x0,ceshi0,x0,x0)[1];
            temp1 = collision_function(1,x0,ceshi1,x0,x0)[1];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[1][1] = x[k];
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[0][2] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[0][2] ;
            temp0 = collision_function(1,x0,x0,ceshi0,x0)[2];
            temp1 = collision_function(1,x0,x0,ceshi1,x0)[2];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[1][2] = x[k];
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[0][3] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[0][3] ;
            temp0 = collision_function(1,x0,x0,x0,ceshi0)[3];
            temp1 = collision_function(1,x0,x0,x0,ceshi1)[3];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[1][3] = x[k];
            tt = -1;
            break;
        }
    }
    vec_ga.SetLength(0);
    append(vec_ga,key_found[1][0][0]);
    append(vec_ga,key_found[1][1][0]);
    append(vec_ga,key_found[1][2][0]);
    append(vec_ga,key_found[1][3][0]);
    vec_ga = reverse(vec_ga);
    key32_found[1].SetDims(1,32);
    VectorCopy(key32_found[1][0],vec_ga,32);
    rki[1] = uint32FromMatGF2(transpose(key32_found[1])) ;
    
    //cout << key32_found[1] <<endl;


}

void key2_function()
{
    //x = Produce_x();//产生0-255的输入
    mat_GF2 x0,x1,n1,ceshi0,ceshi1,temp0,temp1,output0;
    mat_GF2 out3[256];
    vec_GF2 vec_ga;
    int tt=-1;

   // invEij_temp = invEij_ceshi() ;


    x1 = x[1] ;//00000001
    x0 = x[0];//00000000


    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[1][0] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[1][0] ;
            temp0 = collision_function(2,ceshi0,x0,x0,x0)[0];
            temp1 = collision_function(2,ceshi1,x0,x0,x0)[0];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            if(k == 0)
            {
                continue;
            }
            else
            {
                key_found[2][0] = x[k];
                //cout << key_found[2][0] <<endl;
                tt = -1;
                break;
            }
            
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[1][1];
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[1][1] ;
            temp0 = collision_function(2,x0,ceshi0,x0,x0)[1];
            temp1 = collision_function(2,x0,ceshi1,x0,x0)[1];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[2][1] = x[k];
            //cout << key_found[2][1] <<endl;
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[1][2] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[1][2] ;
            temp0 = collision_function(2,x0,x0,ceshi0,x0)[2];
            temp1 = collision_function(2,x0,x0,ceshi1,x0)[2];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[2][2] = x[k];
            //cout << key_found[2][2] <<endl;
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[1][3] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) *  invEijxing_huifu[1][3] ;
            temp0 = collision_function(2,x0,x0,x0,ceshi0)[3];
            temp1 = collision_function(2,x0,x0,x0,ceshi1)[3];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[2][3] = x[k];
            //cout << key_found[2][3] <<endl;
            tt = -1;
            break;
        }
    }

    vec_ga.SetLength(0);
    append(vec_ga,key_found[2][0][0]);
    append(vec_ga,key_found[2][1][0]);
    append(vec_ga,key_found[2][2][0]);
    append(vec_ga,key_found[2][3][0]);
    vec_ga = reverse(vec_ga);
    key32_found[2].SetDims(1,32);
    VectorCopy(key32_found[2][0],vec_ga,32);
    rki[2] = uint32FromMatGF2(transpose(key32_found[2])) ;

}

void key3_function()
{
    //x = Produce_x();//产生0-255的输入
    mat_GF2 x0,x1,n1,ceshi0,ceshi1,temp0,temp1,output0;
    mat_GF2 out3[256];
    vec_GF2 vec_ga;
    int tt=-1;

   // invEij_temp = invEij_ceshi() ;


    x1 = x[1] ;//00000001
    x0 = x[0];//00000000


    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[2][0];
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[2][0] ;
            temp0 = collision_function(3,ceshi0,x0,x0,x0)[0];
            temp1 = collision_function(3,ceshi1,x0,x0,x0)[0];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>200)
        {
            if(k == 0)
            {
                continue;
            }
            else
            {
                key_found[3][0] = x[k];
                //cout << key_found[3][0] <<endl;
                tt = -1;
                break;
            }
            
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[2][1] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[2][1] ;
            temp0 = collision_function(3,x0,ceshi0,x0,x0)[1];
            temp1 = collision_function(3,x0,ceshi1,x0,x0)[1];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[3][1] = x[k];
            //cout << key_found[3][1] <<endl;
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[2][2] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[2][2] ;
            temp0 = collision_function(3,x0,x0,ceshi0,x0)[2];
            temp1 = collision_function(3,x0,x0,ceshi1,x0)[2];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[3][2] = x[k];
            //cout << key_found[3][2] <<endl;
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[2][3] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[2][3];
            temp0 = collision_function(3,x0,x0,x0,ceshi0)[3];
            temp1 = collision_function(3,x0,x0,x0,ceshi1)[3];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[3][3] = x[k];
            //cout << key_found[3][3] <<endl;
            tt = -1;
            break;
        }
    }

    vec_ga.SetLength(0);
    append(vec_ga,key_found[3][0][0]);
    append(vec_ga,key_found[3][1][0]);
    append(vec_ga,key_found[3][2][0]);
    append(vec_ga,key_found[3][3][0]);
    vec_ga = reverse(vec_ga);
    key32_found[3].SetDims(1,32);
    VectorCopy(key32_found[3][0],vec_ga,32);
    rki[3] = uint32FromMatGF2(transpose(key32_found[3])) ;

}

void key4_function()
{
    //x = Produce_x();//产生0-255的输入
    mat_GF2 x0,x1,n1,ceshi0,ceshi1,temp0,temp1,output0;
    mat_GF2 out3[256];
    vec_GF2 vec_ga;
    int tt=-1;

   // invEij_temp = invEij_ceshi() ;


    x1 = x[1] ;//00000001
    x0 = x[0];//00000000


    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[3][0];
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[3][0];
            temp0 = collision_function(4,ceshi0,x0,x0,x0)[0];
            temp1 = collision_function(4,ceshi1,x0,x0,x0)[0];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>200)
        {
            if(k == 0)
            {
                continue;
            }
            else
            {
                key_found[4][0] = x[k];
                //cout << key_found[4][0] <<endl;
                tt = -1;
                break;
            }
            
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) *invEijxing_huifu[3][1] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[3][1];
            temp0 = collision_function(4,x0,ceshi0,x0,x0)[1];
            temp1 = collision_function(4,x0,ceshi1,x0,x0)[1];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[4][1] = x[k];
            //cout << key_found[4][1] <<endl;
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[3][2] ;
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[3][2];
            temp0 = collision_function(4,x0,x0,ceshi0,x0)[2];
            temp1 = collision_function(4,x0,x0,ceshi1,x0)[2];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[4][2] = x[k];
            //cout << key_found[4][2] <<endl;
            tt = -1;
            break;
        }
    }

    for(int k = 0;k<256;k++)
    {
        for(int t = 0;t<256;t++)
        {
            ceshi0 = (Ni_sm4Sbox(x[t]) + x[k]) * invEijxing_huifu[3][3];
            n1 = x[t] + x1;//1x8
            ceshi1 = (Ni_sm4Sbox(n1) + x[k]) * invEijxing_huifu[3][3] ;
            temp0 = collision_function(4,x0,x0,x0,ceshi0)[3];
            temp1 = collision_function(4,x0,x0,x0,ceshi1)[3];
            output0 = temp0 + temp1;
            out3[t] = output0;
            if(t == 0 || t%2 == 1)
            {
                continue;
            }
            else
            {
                if(out3[t] == out3[t-2])
                {
                    tt = t;
                continue;
                }
                else
                {
                    break;
                }
            }


        }
        if(tt>100)
        {
            key_found[4][3] = x[k];
            //cout << key_found[4][3] <<endl;
            tt = -1;
            break;
        }
    }

    vec_ga.SetLength(0);
    append(vec_ga,key_found[4][0][0]);
    append(vec_ga,key_found[4][1][0]);
    append(vec_ga,key_found[4][2][0]);
    append(vec_ga,key_found[4][3][0]);
    vec_ga = reverse(vec_ga);
    key32_found[4].SetDims(1,32);
    VectorCopy(key32_found[4][0],vec_ga,32);
    rki[4] = uint32FromMatGF2(transpose(key32_found[4])) ;

}

unsigned int SM4_CK1[32] ={0x00070e15, 0x1c232a31, 0x383f464d, 0x545b6269,
0x70777e85, 0x8c939aa1, 0xa8afb6bd, 0xc4cbd2d9,
0xe0e7eef5, 0xfc030a11, 0x181f262d, 0x343b4249,
0x50575e65, 0x6c737a81, 0x888f969d, 0xa4abb2b9,
0xc0c7ced5, 0xdce3eaf1, 0xf8ff060d, 0x141b2229,
0x30373e45, 0x4c535a61, 0x686f767d, 0x848b9299,
0xa0a7aeb5, 0xbcc3cad1, 0xd8dfe6ed, 0xf4fb0209,
0x10171e25, 0x2c333a41, 0x484f565d, 0x646b7279};

unsigned char SM4_Sbox1[256] =
{0xd6,0x90,0xe9,0xfe,0xcc,0xe1,0x3d,0xb7,0x16,0xb6,0x14,0xc2,0x28,0xfb,0x2c,0x05,
0x2b,0x67,0x9a,0x76,0x2a,0xbe,0x04,0xc3,0xaa,0x44,0x13,0x26,0x49,0x86,0x06,0x99,
0x9c,0x42,0x50,0xf4,0x91,0xef,0x98,0x7a,0x33,0x54,0x0b,0x43,0xed,0xcf,0xac,0x62,
0xe4,0xb3,0x1c,0xa9,0xc9,0x08,0xe8,0x95,0x80,0xdf,0x94,0xfa,0x75,0x8f,0x3f,0xa6,
0x47,0x07,0xa7,0xfc,0xf3,0x73,0x17,0xba,0x83,0x59,0x3c,0x19,0xe6,0x85,0x4f,0xa8,
0x68,0x6b,0x81,0xb2,0x71,0x64,0xda,0x8b,0xf8,0xeb,0x0f,0x4b,0x70,0x56,0x9d,0x35,
0x1e,0x24,0x0e,0x5e,0x63,0x58,0xd1,0xa2,0x25,0x22,0x7c,0x3b,0x01,0x21,0x78,0x87,
0xd4,0x00,0x46,0x57,0x9f,0xd3,0x27,0x52,0x4c,0x36,0x02,0xe7,0xa0,0xc4,0xc8,0x9e,
0xea,0xbf,0x8a,0xd2,0x40,0xc7,0x38,0xb5,0xa3,0xf7,0xf2,0xce,0xf9,0x61,0x15,0xa1,
0xe0,0xae,0x5d,0xa4,0x9b,0x34,0x1a,0x55,0xad,0x93,0x32,0x30,0xf5,0x8c,0xb1,0xe3,
0x1d,0xf6,0xe2,0x2e,0x82,0x66,0xca,0x60,0xc0,0x29,0x23,0xab,0x0d,0x53,0x4e,0x6f,
0xd5,0xdb,0x37,0x45,0xde,0xfd,0x8e,0x2f,0x03,0xff,0x6a,0x72,0x6d,0x6c,0x5b,0x51,
0x8d,0x1b,0xaf,0x92,0xbb,0xdd,0xbc,0x7f,0x11,0xd9,0x5c,0x41,0x1f,0x10,0x5a,0xd8,
0x0a,0xc1,0x31,0x88,0xa5,0xcd,0x7b,0xbd,0x2d,0x74,0xd0,0x12,0xb8,0xe5,0xb4,0xb0,
0x89,0x69,0x97,0x4a,0x0c,0x96,0x77,0x7e,0x65,0xb9,0xf1,0x09,0xc5,0x6e,0xc6,0x84,
0x18,0xf0,0x7d,0xec,0x3a,0xdc,0x4d,0x20,0x79,0xee,0x5f,0x3e,0xd7,0xcb,0x39,0x48};



void key0_function()
{
    unsigned int temp,buf,dop;
    vec_GF2 out_vec;

    temp = rki[1] ^ rki[2] ^ rki[3] ^ SM4_CK1[4];

    buf= (SM4_Sbox1[(temp >> 24) & 0xFF]) << 24
    |(SM4_Sbox1[(temp >> 16) & 0xFF]) << 16
    |(SM4_Sbox1[(temp >> 8) & 0xFF]) << 8
    |(SM4_Sbox1[temp & 0xFF]);

    dop = ((buf)^(SM4_Rotl32((buf),13))^(SM4_Rotl32((buf),23)));

    rki[0] = rki[4] ^ dop;
    out_vec = reverse(transpose(matGF2FromUint32(rki[0]))[0]);
    key32_found[0].SetDims(1,32);
    VectorCopy(key32_found[0][0],out_vec,32);

}

static uint32_t K[4];

void Ki_found()
{

    unsigned int temp0,temp1,temp2,temp3,buf0,buf1,buf2,buf3,dop0,dop1,dop2,dop3;

    temp0 = rki[0] ^ rki[1] ^ rki[2] ^ SM4_CK1[3];

    buf0= (SM4_Sbox1[(temp0 >> 24) & 0xFF]) << 24
    |(SM4_Sbox1[(temp0>> 16) & 0xFF]) << 16
    |(SM4_Sbox1[(temp0 >> 8) & 0xFF]) << 8
    |(SM4_Sbox1[temp0 & 0xFF]);

    dop0 = ((buf0)^(SM4_Rotl32((buf0),13))^(SM4_Rotl32((buf0),23)));

    K[3] = rki[3] ^ dop0;

    temp1 = K[3] ^ rki[0] ^ rki[1] ^ SM4_CK1[2];

    buf1= (SM4_Sbox1[(temp1 >> 24) & 0xFF]) << 24
    |(SM4_Sbox1[(temp1>> 16) & 0xFF]) << 16
    |(SM4_Sbox1[(temp1 >> 8) & 0xFF]) << 8
    |(SM4_Sbox1[temp1 & 0xFF]);

    dop1 = ((buf1)^(SM4_Rotl32((buf1),13))^(SM4_Rotl32((buf1),23)));

    K[2] = rki[2] ^ dop1;

    temp2 = K[2] ^ K[3] ^ rki[0] ^ SM4_CK1[1];

    buf2= (SM4_Sbox1[(temp2 >> 24) & 0xFF]) << 24
    |(SM4_Sbox1[(temp2>> 16) & 0xFF]) << 16
    |(SM4_Sbox1[(temp2 >> 8) & 0xFF]) << 8
    |(SM4_Sbox1[temp2 & 0xFF]);

    dop2 = ((buf2)^(SM4_Rotl32((buf2),13))^(SM4_Rotl32((buf2),23)));

    K[1] = rki[1] ^ dop2;

    temp3 = K[1] ^ K[2] ^ K[3] ^ SM4_CK1[0];

    buf3= (SM4_Sbox1[(temp3 >> 24) & 0xFF]) << 24
    |(SM4_Sbox1[(temp3>> 16) & 0xFF]) << 16
    |(SM4_Sbox1[(temp3 >> 8) & 0xFF]) << 8
    |(SM4_Sbox1[temp3 & 0xFF]);

    dop3 = ((buf3)^(SM4_Rotl32((buf3),13))^(SM4_Rotl32((buf3),23)));

    K[0] = rki[0] ^ dop3;

    
}

unsigned int SM4_FK1[4] = {0xA3B1BAC6, 0x56AA3350, 0x677D9197, 0xB27022DC};

static uint32_t key_output[4];

void key_end()
{
    cout << "碰撞攻击后恢复的密钥:" << endl;

    uint8_t in_ciph[4][4];
    key_output[0] = K[0] ^ SM4_FK1[0];
    key_output[1] = K[1] ^ SM4_FK1[1];
    key_output[2] = K[2] ^ SM4_FK1[2];
    key_output[3] = K[3] ^ SM4_FK1[3];

    *(uint32_t *)&in_ciph[0] = key_output[0];
    *(uint32_t *)&in_ciph[1] = key_output[1];
    *(uint32_t *)&in_ciph[2] = key_output[2];
    *(uint32_t *)&in_ciph[3] = key_output[3];

    for(int k =0 ;k<4;k++)
    {
        printstate(in_ciph[k]);
    }
    cout << endl;

    //cout << transpose(matGF2FromUint32(key_output[0]) ) ;

}

























