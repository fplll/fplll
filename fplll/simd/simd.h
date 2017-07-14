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

  typedef double (*dot_product_ptr)(const double *, const double *, size_t n);
  dot_product_ptr dot_product;
};

class SIMD_operations
{
public:
  /* detect which SIMD version to use and initialize functions,
     disable by keywords: sse128d avx256d avx512d */
  SIMD_operations(const std::string &disable_versions = "");

  std::string simd_double_version() const { return std::string(double_functions.name); }

  inline double dot_product(const double *x, const double *y, size_t n) const
  {
    return (*double_functions.dot_product)(x, y, n);
  }

private:
  SIMD_double_functions double_functions;
};

FPLLL_END_NAMESPACE

#endif
