#pragma once
#include "stdint.h"
#include <NTL/mat_GF2E.h>
#include <NTL/mat_GF2.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <cstdint>

using namespace NTL;

mat_GF2 * Produce_temp();

void coll();

void  huifuSBOX();

void huifuSboxNi();

//int gauss(int ** a);

void huifu_invEij();

void key1_function();

void key2_function();

void key3_function();

void key4_function();

void key0_function();

void Ki_found();

void key_end();