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

#include "simd.h"

FPLLL_BEGIN_NAMESPACE

#ifdef FPLLL_HAVE_AVX
extern const SIMD_double_functions SIMD_double_avx256_functions;
#endif

/* non-SIMD implementation */
static bool always_true() { return true; }
static double dot_product(const double *x, const double *y, size_t n)
{
  double ret = 0;
  for (size_t i = 0; i < n; ++i)
  {
    ret += x[i] * y[i];
  }
  return ret;
}
SIMD_double_functions SIMD_double_nonsimd_functions({"nonsimd", always_true, dot_product});

/* SIMD_operations constructor: decide which SIMD version to use */
SIMD_operations::SIMD_operations(const std::string &disable_versions)
{
#ifdef FPLLL_HAVE_AVX
  if (disable_versions.find("avx256d") == std::string::npos &&
      SIMD_double_avx256_functions.cpu_supported())
  {
    double_functions = SIMD_double_avx256_functions;
    return;
  }
#endif
  double_functions = SIMD_double_nonsimd_functions;
}

FPLLL_END_NAMESPACE
