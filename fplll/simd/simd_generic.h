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

#ifndef FPLLL_SIMD_GENERIC_H
#define FPLLL_SIMD_GENERIC_H

#include <fplll/defs.h>

FPLLL_BEGIN_NAMESPACE

template <typename SIMD_type> class SIMD_double_implementation : public SIMD_type
{
public:
  using typename SIMD_type::vec_t;
  using SIMD_type::width;
  using SIMD_type::v_zero;
  using SIMD_type::vv_mul;
  using SIMD_type::vv_add;
  using SIMD_type::cpu_supported;
  using SIMD_type::load;

  static double &get(vec_t &x, size_t i) { return reinterpret_cast<double *>(&x)[i]; }
  static const double &get(const vec_t &x, size_t i)
  {
    return reinterpret_cast<const double *>(&x)[i];
  }
  static double dot_product(const double *x, const double *y, size_t n)
  {
    vec_t res = v_zero();

    size_t i = 0, e = n & (size_t(0) - size_t(width));
    for (; i != e; i += width)
    {
      res = vv_add(res, vv_mul(load(x + i), load(y + i)));
    }
    double r = 0;
    for (size_t j = 0; j < width; ++j)
    {
      r += get(res, j);
    }
    for (; i != n; ++i)
    {
      r += x[i] * y[i];
    }
    return r;
  }
};

FPLLL_END_NAMESPACE

#endif
