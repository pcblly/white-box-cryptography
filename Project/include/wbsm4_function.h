#pragma once
#include "stdint.h"
#include <NTL/mat_GF2E.h>
#include <NTL/mat_GF2.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <cstdint>

using namespace NTL;

mat_GF2 DiagMatrix(mat_GF2 * Eij);

//mat_GF2 DiagMatrix32(mat_GF2 * Fi) ;

//mat_GF2 DiagInvMatrix32(mat_GF2 * Gi);

mat_GF2 Xi123 (const mat_GF2& x1,const mat_GF2& x2,const mat_GF2& x3,int i,const mat_GF2 & Mi,const mat_GF2 &InvDiag_Mat,const mat_GF2* Oi );

mat_GF2 QiMat(const mat_GF2& xi0,int i,const mat_GF2 * Oi) ;

mat_GF2 Tbox_ij(const uint8_t x0,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi,const mat_GF2& aij,const mat_GF2& M_xing,const mat_GF2* Oi);

mat_GF2 Sibox(const mat_GF2 & x, const mat_GF2 & Mi);

mat_GF2 LinearTrans(const mat_GF2 & Mi ,int j ,const mat_GF2& x);

mat_GF2 JiLian(const mat_GF2& Mi,const mat_GF2&x) ;

mat_GF2 aij_part(const mat_GF2& aij,const mat_GF2& M_xing,const mat_GF2& Mi,int i,const mat_GF2* Oi);

mat_GF2 roundF(const mat_GF2& xi0,const mat_GF2& xi1,const mat_GF2& xi2,const mat_GF2& xi3,int i,const mat_GF2& rki,const mat_GF2 * Oi);



void printstate(unsigned char * in);

mat_GF2 * Produce_x() ;

mat_GF2 Tbox_ceshi(const mat_GF2& x,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi,const mat_GF2* Oi);

mat_GF2 Tbox_change(const mat_GF2& x,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi);

void wbsm4_gen(const mat_GF2 * Oi,unsigned char MK[]);

mat_GF2 roundF_ceshi(const mat_GF2& xi0,const mat_GF2& xi1,const mat_GF2& xi2,const mat_GF2& xi3,int i,const mat_GF2 * Oi);

void wbsm4_Encrypt(unsigned char CipherText[],const mat_GF2 * Oi);

mat_GF2 Tbox_change(const mat_GF2& x,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi);

mat_GF2 * collision_function(int i,const mat_GF2& x0,const mat_GF2& x1,const mat_GF2& x2,const mat_GF2& x3);

mat_GF2 * ui0_xing(int i);

mat_GF2 Sbox_change(const mat_GF2& x,int i,int j,mat_GF2 rki , const mat_GF2& Eij,const mat_GF2& Mi);

mat_GF2 * ui1_xing(int i);

mat_GF2 * ui2_xing(int i);

mat_GF2 * ui3_xing(int i);

mat_GF2 ** invEij_ceshi();

//void eij();