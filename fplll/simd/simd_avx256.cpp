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

struct SIMD_double_avx256
{
  typedef __m256d vec_t;
  static const size_t width = 4;

  static inline vec_t v_zero() { return _mm256_setzero_pd(); }
  static inline vec_t vv_mul(vec_t x, vec_t y) { return _mm256_mul_pd(x, y); }
  static inline vec_t vv_add(vec_t x, vec_t y) { return _mm256_add_pd(x, y); }
  static bool cpu_supported()
  {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx");
  }
};

using AVX256D = SIMD_double_implementation<SIMD_double_avx256>;
template class SIMD_double_implementation<SIMD_double_avx256>;

const SIMD_double_functions SIMD_double_avx256_functions = {AVX256D::cpu_supported, AVX256D::dot_product};

FPLLL_END_NAMESPACE

#endif
