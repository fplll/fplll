/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.

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
/** \file util.h
    Miscellaneous. */

#ifndef FPLLL_UTIL_H
#define FPLLL_UTIL_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

static inline int cputime()
{
#ifdef FPLLL_WITH_GETRUSAGE
  struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
#else
  return time(NULL) * 1000;
#endif
}

/* Computes result = x * m (it is not required that result is initialized,
   x.size() must be equal to m.GetNumRows()) */
template <class ZT>
void vector_matrix_product(vector<ZT> &result, const vector<ZT> &x, const Matrix<ZT> &m)
{
  int nrows = m.get_rows(), ncols = m.get_cols();
  gen_zero_vect(result, ncols);
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
      result[j].addmul(x[i], m(i, j));
}

/* Computes result = x * m (it is not required that result is initialized,
   x.size() must be equal to m.GetNumRows()) */
template <class ZT>
void vector_matrix_product(NumVect<ZT> &result, const NumVect<ZT> &x, const Matrix<ZT> &m)
{
  int nrows = m.get_rows(), ncols = m.get_cols();
  result.gen_zero(ncols);
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
      result[j].addmul(x[i], m(i, j));
}

template <class T>
void scalar_product(T &result, const MatrixRow<T> &v1, const MatrixRow<T> &v2, int n)
{
  FPLLL_DEBUG_CHECK(n <= static_cast<int>(v1.size()) && n <= static_cast<int>(v2.size()));
  T tmp;
  result.mul(v1[0], v2[0]);
  for (int i = 1; i < n; i++)
  {
    tmp.mul(v1[i], v2[i]);
    result.add(result, tmp);
  }
}

template <class T> inline void sqr_norm(T &result, const MatrixRow<T> &v, int n)
{
  scalar_product(result, v, v, n);
}

const double DEF_GSO_PREC_EPSILON = 0.03;

/**
 * Returns the minimum precision required to ensure that error bounds on the
 * GSO are valid. Computes rho such that for all 0 &lt;= i &lt; d and 0 &lt;= j &lt;= i:
 *
 *   |r~_i - r_i| / r_i     <= d * rho ^ (i + 1) * 2 ^ (2 - prec)
 *
 *   |mu~_(i,j) - mu_(i,j)| <= d * rho ^ (i + 1) * 2 ^ (4 - prec)
 */
int gso_min_prec(double &rho, int d, double delta, double eta,
                 double epsilon = DEF_GSO_PREC_EPSILON);

/**
 * Returns the minimum precision for the proved version of LLL.
 */
int l2_min_prec(int d, double delta, double eta, double epsilon);

/**
 * Computes the volume of a d-dimensional hypersphere of radius 1.
 */
void sphere_volume(FP_NR<> &volume, int d);

/**
 * Estimates the cost of the enumeration for SVP.
 */
void cost_estimate(FP_NR<> &cost, const FP_NR<> &bound, const Matrix<FP_NR<>> &r, int dimMax);

template <class ZT> void zeros_first(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv_t);

template <class ZT> void zeros_last(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv_t);

/**
 * Returns the string corresponding to an error code of LLL/BKZ.
 */
const char *get_red_status_str(int status);

FPLLL_END_NAMESPACE

#endif
