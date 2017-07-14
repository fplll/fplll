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

#ifndef FPLLL_SIMD_H
#define FPLLL_SIMD_H

#include <fplll/defs.h>
#include <string>

FPLLL_BEGIN_NAMESPACE

/* each SIMD version provides a struct of pointers to use */
struct SIMD_double_functions
{
  const char *name;

  typedef bool (*cpu_supported_ptr)();
  cpu_supported_ptr cpu_supported;

  typedef void (*add_ptr)(double *, const double *, const double *, size_t);
  add_ptr add;

  typedef void (*sub_ptr)(double *, const double *, const double *, size_t);
  sub_ptr sub;

  typedef void (*addmul_ptr)(double *, const double *, double, const double *, size_t);
  addmul_ptr addmul;

  typedef void (*addmul2_ptr)(double *, double, const double *, double, const double *, size_t);
  addmul2_ptr addmul2;

  typedef double (*dot_product_ptr)(const double *, const double *, size_t);
  dot_product_ptr dot_product;

  typedef void (*givens_rotation_ptr)(double *, double *, const double *, const double *, double,
                                      double, size_t);
  givens_rotation_ptr givens_rotation;
};

/* these are nonsimd placeholder implementations for unimplemented simd operations */
void nonsimd_add(double *, const double *, const double *, size_t);
void nonsimd_sub(double *, const double *, const double *, size_t);
void nonsimd_addmul(double *, const double *, double, const double *, size_t);
void nonsimd_addmul2(double *, double, const double *, double, const double *, size_t);
double nonsimd_dot_product(const double *, const double *, size_t);
void nonsimd_givens_rotation(double *, double *, const double *, const double *, double, double,
                             size_t);

class SIMD_operations
{
public:
  /* detect which SIMD version to use and initialize functions,
     disable by keywords: sse128d avx256d avx512d */
  SIMD_operations(const std::string &disable_versions = "");

  std::string simd_double_version() const { return std::string(double_functions.name); }

  // add: dstx = x+y (dstx can be equal to x)
  inline void add(double *dstx, const double *x, const double *y, size_t n) const
  {
    (*double_functions.add)(dstx, x, y, n);
  }

  // sub: dstx = x-y (dstx can be equal to x)
  inline void sub(double *dstx, const double *x, const double *y, size_t n) const
  {
    (*double_functions.sub)(dstx, x, y, n);
  }

  // addmul: dstx = x+a*y (dstx can be equal to x) (does not specialize for value a=0,1,-1)
  inline void addmul(double *dstx, const double *x, double a, const double *y, size_t n) const
  {
    (*double_functions.addmul)(dstx, x, a, y, n);
  }

  // addmul: dstx = a*x+b*y (dstx can be equal to x) (does not specialize for value a/b=0,1,-1)
  inline void addmul(double *dstx, double a, const double *x, double b, const double *y,
                     size_t n) const
  {
    (*double_functions.addmul2)(dstx, a, x, b, y, n);
  }

  // dot product: return sum_(i=0)^(n-1) x[i]*y[i]
  inline double dot_product(const double *x, const double *y, size_t n) const
  {
    return (*double_functions.dot_product)(x, y, n);
  }

  // givens rotation: dstx = a*x + b*y, dsty = -b*x + a*y (dstx/dsty can be equal to x/y) (does not
  // specialize for value a/b=0,1,-1)
  inline void givens_rotation(double *dstx, double *dsty, const double *x, const double *y,
                              double a, double b, size_t n) const
  {
    (*double_functions.givens_rotation)(dstx, dsty, x, y, a, b, n);
  }

private:
  SIMD_double_functions double_functions;
};

FPLLL_END_NAMESPACE

#endif
