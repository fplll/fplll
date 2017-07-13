/* Copyright (C) 2017 Marc Stevens

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#ifndef FPLLL_SIMD_AVX256_H
#define FPLLL_SIMD_AVX256_H

#include "simd.h"
#include "simd_generic.h"
#include <immintrin.h>

FPLLL_BEGIN_NAMESPACE

struct SIMD_avx256
{
  typedef __m256d vec_t;
  static const size_t width = 4;

  static vec_t v_zero() { return _mm256_setzero_pd(); }
  static vec_t vv_mul(vec_t x, vec_t y) { return _mm256_mul_pd(x, y); }
  static vec_t vv_add(vec_t x, vec_t y) { return _mm256_add_pd(x, y); }
  static bool detect()
  {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx");
  }
};

template class SIMD_generic_implementation<SIMD_avx256>;
using AVX256 = SIMD_generic_implementation<SIMD_avx256>;

const SIMD_functions SIMD_avx256_functions = {AVX256::detect, AVX256::dot_product};

FPLLL_END_NAMESPACE

#endif
