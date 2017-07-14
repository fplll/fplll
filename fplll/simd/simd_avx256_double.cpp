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

#ifdef FPLLL_HAVE_AVX
#include <cpuid.h>
#include <immintrin.h>

FPLLL_BEGIN_NAMESPACE

struct AVX256D
{
  typedef __m256d vec_t;
  static const size_t width     = 4;
  static const size_t alignment = 32;  // vector byte alignment for aligned load/write

  static inline vec_t v_zero() { return _mm256_setzero_pd(); }
  static inline vec_t vv_mul(vec_t x, vec_t y) { return _mm256_mul_pd(x, y); }
  static inline vec_t vv_add(vec_t x, vec_t y) { return _mm256_add_pd(x, y); }
  static inline vec_t unaligned_load(const double *p) { return _mm256_loadu_pd(p); }
  static inline vec_t aligned_load(const double *p) { return _mm256_load_pd(p); }

  static bool cpu_supported()
  {
    unsigned int a, b, c, d;
    __cpuid_count(0, 0, a, b, c, d);
    unsigned int ID = a;
    if (!(ID >= 0x00000001))
      return false;
    __cpuid_count(0x00000001, 0, a, b, c, d);
    if (0 == (c & (1 << 28)))
      return false;
    // TODO: check if OS has enabled AVX support
    return true;
  }

  static double dot_product(const double *x, const double *y, size_t n)
  {
    if (n < width * 2)
    {
      // for short rows we use the non-simd version
      double r = 0;
      for (size_t i = 0; i < n; ++i)
      {
        r += x[i] * y[i];
      }
      return r;
    }
    if ((size_t(x) - size_t(y)) % alignment == 0)
    {
      // aligned version
      double r = 0;
      size_t i = 0, e = ((size_t(32) - size_t(x)) % alignment) / sizeof(double);
      for (; i < e; ++i)
      {
        r += x[i] * y[i];
      }
      vec_t rv = v_zero();
      for (; i + width <= n; i += width)
      {
        rv = vv_add(rv, vv_mul(aligned_load(x + i), aligned_load(y + i)));
      }
      for (; i < n; ++i)
      {
        r += x[i] * y[i];
      }
      for (i = 0; i < width; ++i)
      {
        r += reinterpret_cast<double *>(&rv)[i];
      }
      return r;
    }
    else
    {
      // one aligned, other unaligned version
      // aligned version
      double r = 0;
      size_t i = 0, e = ((size_t(32) - size_t(x)) % alignment) / sizeof(double);
      for (; i < e; ++i)
      {
        r += x[i] * y[i];
      }
      vec_t rv = v_zero();
      for (; i + width <= n; i += width)
      {
        rv = vv_add(rv, vv_mul(aligned_load(x + i), unaligned_load(y + i)));
      }
      for (; i < n; ++i)
      {
        r += x[i] * y[i];
      }
      for (i = 0; i < width; ++i)
      {
        r += reinterpret_cast<double *>(&rv)[i];
      }
      return r;
    }
  }
};

SIMD_double_functions SIMD_double_avx256_functions({"avx256d", AVX256D::cpu_supported, nonsimd_add,
                                                    nonsimd_sub, nonsimd_addmul, nonsimd_addmul2,
                                                    AVX256D::dot_product, nonsimd_givens_rotation});

FPLLL_END_NAMESPACE

#endif
#endif
