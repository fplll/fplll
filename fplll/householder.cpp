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

/* Template source file */

#include "householder.h"

FPLLL_BEGIN_NAMESPACE

// TODO: maybe add a variable available on DEBUG to verify in update_R(i, false) was done
template <class ZT, class FT> void MatHouseholder<ZT, FT>::update_R_last(int i)
{
  // sigma[i] = sign(r[1])
  sigma[i] = (R(i, i).cmp(0.0) < 0) ? -1.0 : 1.0;
  // V(i, i) is used as a temporary variable. In the following, V(i, i) is always modified.
  if (i + 1 == n)
    V(i, i) = 0.0;
  else
  {
    // r^T * r (with first coefficient ignored)
    R[i].dot_product(V(i, i), R[i], i + 1, n);
  }
  ftmp1.mul(R(i, i), R(i, i));
  // ftmp1 = r^T * r
  ftmp1.add(ftmp1, V(i, i));

  if (ftmp1.cmp(0.0) != 0)
  {
    ftmp2.sqrt(ftmp1);
    // s = sigma[i] * ||r|| = sigma[i] * sqrt(r * r^T)
    ftmp0.mul(sigma[i], ftmp2);
    ftmp1.add(R(i, i), ftmp0);
    // vi[1] = (-sum(r[j]^2, j, 2, n-i+1)
    V(i, i).neg(V(i, i));
    // vi[1] = (-sum(r[j]^2, j, 2, n-i+1) / (r[1] + s)
    V(i, i).div(V(i, i), ftmp1);
    if (V(i, i).cmp(0.0) != 0)
    {
      ftmp0.neg(ftmp0);
      ftmp0.mul(ftmp0, V(i, i));
      ftmp0.sqrt(ftmp0);

#ifndef HOUSEHOLDER_PRECOMPUTE_INVERSE
      // Here, ftmp0 = 1 / sqrt(-s * vi[1])
      V(i, i).div(V(i, i), ftmp0);
      R(i, i) = ftmp2;

      for (int k = i + 1; k < n; k++)
      {
        V(i, k).div(R(i, k), ftmp0);
        R(i, k) = 0.0;
      }
#else   // HOUSEHOLDER_PRECOMPUTE_INVERSE
      ftmp0.div(1.0, ftmp0);
      // Here, ftmp0 = 1 / sqrt(-s * vi[1])
      V(i, i).mul(V(i, i), ftmp0);
      R(i, i) = ftmp2;
      for (int k = i + 1; k < n; k++)
      {
        V(i, k).mul(R(i, k), ftmp0);
        R(i, k) = 0.0;
      }

      R_inverse_diag[i].div(1.0, ftmp2);
#endif  // HOUSEHOLDER_PRECOMPUTE_INVERSE
      // Here, vi = vi / ftmp0 and ri[i..n] = (||r||, 0, 0, ..., 0)
    }
    else
    {
      if (R(i, i).cmp(0.0) < 0)
        R(i, i).neg(R(i, i));

      for (int k = i + 1; k < n; k++)
      {
        // if enable_row_expo, R can be not correct at some point of the computation
        FPLLL_DEBUG_CHECK(enable_row_expo ? 1 : R(i, k).is_zero());
        V(i, k) = 0.0;
        R(i, k) = 0.0;
      }

#ifdef HOUSEHOLDER_PRECOMPUTE_INVERSE
      R_inverse_diag[i].div(1.0, R(i, i));
#endif  // HOUSEHOLDER_PRECOMPUTE_INVERSE
    }
  }
  else
  {
    for (int k = i; k < n; k++)
    {
      // if enable_row_expo, R can be not correct at some point of the computation
      FPLLL_DEBUG_CHECK(enable_row_expo ? 1 : R(i, k).is_zero());
      V(i, k) = 0.0;
      R(i, k) = 0.0;
    }

#ifdef HOUSEHOLDER_PRECOMPUTE_INVERSE
    // Result is inf.
    // TODO: set inf instead of doing the computation.
    R_inverse_diag[i].div(1.0, 0.0);
#endif  // HOUSEHOLDER_PRECOMPUTE_INVERSE
  }

  n_known_rows++;
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::update_R(int i, bool last_j)
{
  // To update row i, we need to know n_known_rows rows
  FPLLL_DEBUG_CHECK(i <= n_known_rows);
  if (/*i == n_known_rows && =*/!updated_R)
  {
    int j, k;

    for (j = 0; j < i; j++)
    {
      // vj * ri[j..n]^T
      // TODO: n_known_cols here?
      V[j].dot_product(ftmp0, R[i], j, n);

      //-vj * ri[j..n]^T
      ftmp0.neg(ftmp0);
      for (k = j; k < n; k++)
      {
        // ri[j..n] = ri[j..n] - (vj * ri[j..n]^T) * vj
        R(i, k).addmul(V(j, k), ftmp0);
      }
      // ri[j] = sigma[j] * ri[j]
      R(i, j).mul(sigma[j], R(i, j));

      // Copy R into R_history
      for (int k           = j; k < n; k++)
        R_history[i][j][k] = R(i, k);
    }

    if (last_j)
      update_R_last(i);
  }
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::refresh_R_bf(int i)
{
  int j;

  n_known_cols = max(n_known_cols, init_row_size[i]);

  if (enable_row_expo)
  {
    long max_expo = LONG_MIN;

    for (j = 0; j < n_known_cols; j++)
    {
      b(i, j).get_f_exp(R(i, j), tmp_col_expo[j]);
      max_expo = max(max_expo, tmp_col_expo[j]);
    }

    for (j = 0; j < n_known_cols; j++)
      R(i, j).mul_2si(R(i, j), tmp_col_expo[j] - max_expo);
    for (j = n_known_cols; j < n; j++)
      R(i, j) = 0.0;

    row_expo[i] = max_expo;
    FPLLL_DEBUG_CHECK(row_expo[i] >= 0);
  }
  else
  {
    for (j = 0; j < n_known_cols; j++)
      R(i, j).set_z(b(i, j));
    for (j = n_known_cols; j < n; j++)
      R(i, j) = 0.0;

    b[i].dot_product(ztmp0, b[i], 0, n_known_cols);
    norm_square_b[i].set_z(ztmp0);
  }

  // Copy R[i] in bf[i] (while we copy b[i] in R[i])
  if (enable_bf)
  {
    for (j = 0; j < n_known_cols; j++)
      bf(i, j) = R(i, j);
    for (j = n_known_cols; j < n; j++)
      bf(i, j) = 0.0;

    bf[i].dot_product(norm_square_b[i], bf[i], 0, n_known_cols);
    if (enable_row_expo)
      expo_norm_square_b[i] = 2 * row_expo[i];
  }
  else
  {
    b[i].dot_product(ztmp0, b[i], 0, n_known_cols);
    ztmp0.get_f_exp(norm_square_b[i], expo_norm_square_b[i]);
  }
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::refresh_R(int i)
{
  int j;

  // Copy bf[i] in R[i] (while we have copied b[i] in bf[i])
  if (enable_bf)
  {
    for (j = 0; j < n_known_cols; j++)
      R(i, j) = bf(i, j);
    for (j = n_known_cols; j < n; j++)
      R(i, j) = 0.0;
  }
  else
  {
    if (enable_row_expo)
    {
      for (j = 0; j < n_known_cols; j++)
      {
        b(i, j).get_f_exp(R(i, j), tmp_col_expo[j]);
        R(i, j).mul_2si(R(i, j), tmp_col_expo[j] - row_expo[i]);
      }
      for (j = n_known_cols; j < n; j++)
        R(i, j) = 0.0;
    }
    else
    {
      for (j = 0; j < n_known_cols; j++)
        R(i, j).set_z(b(i, j));
      for (j = n_known_cols; j < n; j++)
        R(i, j) = 0.0;
    }
  }
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::update_R_naively(int i)
{
  FPLLL_DEBUG_CHECK(i <= n_known_rows_naively);

  int j;

  // Set B in R_naively.
  if (enable_row_expo)
  {
    long max_expo = LONG_MIN;

    for (j = 0; j < n; j++)
    {
      b(i, j).get_f_exp(R_naively(i, j), tmp_col_expo[j]);
      max_expo = max(max_expo, tmp_col_expo[j]);
    }

    for (j = 0; j < n; j++)
      R_naively(i, j).mul_2si(R_naively(i, j), tmp_col_expo[j] - max_expo);

    row_expo_naively[i] = max_expo;
    FPLLL_DEBUG_CHECK(row_expo_naively[i] >= 0);
  }
  else
  {
    for (j = 0; j < n; j++)
      R_naively(i, j).set_z(b(i, j));
  }

  for (j = 0; j < i; j++)
  {
    // vj * ri[j..n]^T
    V_naively[j].dot_product(ftmp0, R_naively[i], j, n);

    //-vj * ri[j..n]^T
    ftmp0.neg(ftmp0);
    for (int k = j; k < n; k++)
    {
      // ri[j..n] = ri[j..n] - (vj * ri[j..n]^T) * vj
      R_naively(i, k).addmul(V_naively(j, k), ftmp0);
    }
    // ri[j] = sigma_naively[j] * ri[j]
    R_naively(i, j).mul(sigma_naively[j], R_naively(i, j));
  }

  // Here, ftmp2 is equal to s in [MSV, ISSAC'09].
  //
  // Copy R_naively[i][i..n] in V_naively
  for (j = i; j < n; j++)
  {
    V_naively(i, j) = R_naively(i, j);
  }
  // sigma_naively[i] = sign(r[1])
  sigma_naively[i] = (R_naively(i, i).cmp(0.0) < 0) ? -1.0 : 1.0;
  R_naively[i].dot_product(ftmp2, R_naively[i], i, n);
  ftmp2.sqrt(ftmp2);
  ftmp2.mul(ftmp2, sigma_naively[i]);
  // Here, ftmp2 = sigma_naively[i] * ||r||

  // ftmp0 = (r[1] + ftmp2)
  ftmp0.add(R_naively(i, i), ftmp2);
  if (ftmp0.cmp(0.0) != 0)
  {
    if (i + 1 == n)
      ftmp1 = 0.0;
    else
      R_naively[i].dot_product(ftmp1, R_naively[i], i + 1, n);
    if (ftmp1.cmp(0.0) != 0)
    {
      // ftmp1 = (-sum(r[j]^2, j, 2, n-i+1)
      ftmp1.neg(ftmp1);

      // vi[1] = (-sum(r[j]^2, j, 2, n-i+1) / (r[1] + ftmp2)
      V_naively(i, i).div(ftmp1, ftmp0);

      ftmp2.neg(ftmp2);
      ftmp0.mul(ftmp2, V_naively(i, i));
      // Here, ftmp0 = -ftmp2 * vi[1]

      // ftmp0 = sqrt(-ftmp2 * vi[1])
      ftmp0.sqrt(ftmp0);

      V_naively(i, i).div(V_naively(i, i), ftmp0);
      R_naively(i, i).abs(ftmp2);

      for (j = i + 1; j < n; j++)
      {
        V_naively(i, j).div(V_naively(i, j), ftmp0);
        R_naively(i, j) = 0.0;
      }
    }
    else
    {
      if (R_naively(i, i).cmp(0.0) < 0)
        R_naively(i, i).neg(R_naively(i, i));

      V_naively(i, i) = 0.0;
      for (int k = i + 1; k < n; k++)
      {
        V_naively(i, k) = 0.0;
        R_naively(i, k) = 0.0;
      }
    }
  }
  else
  {
    for (int k = i; k < n; k++)
    {
      V_naively(i, k) = 0.0;
      R_naively(i, k) = 0.0;
    }
  }

  n_known_rows_naively++;
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::swap(int i, int j)
{
  FPLLL_DEBUG_CHECK(0 <= i && i < j && j < d);

  // Invalidate to the min not modified row, that is i
  invalidate_row(i);

  b.swap_rows(i, j);
  if (enable_bf)
    bf.swap_rows(i, j);
  iter_swap(sigma.begin() + i, sigma.begin() + j);
  if (enable_row_expo)
    iter_swap(row_expo.begin() + i, row_expo.begin() + j);
  iter_swap(init_row_size.begin() + i, init_row_size.begin() + j);
  iter_swap(R_history.begin() + i, R_history.begin() + j);
  if (enable_transform)
  {
    u.swap_rows(i, j);
    if (enable_inverse_transform)
      u_inv_t.swap_rows(i, j);
  }
  iter_swap(norm_square_b.begin() + i, norm_square_b.begin() + j);
  iter_swap(expo_norm_square_b.begin() + i, expo_norm_square_b.begin() + j);
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::size_reduce(const FT &xf, int k, int i)
{
  FPLLL_DEBUG_CHECK(k > 0 && k < d);
  FPLLL_DEBUG_CHECK(xf.cmp(0.0) != 0);

  row_addmul_we(k, i, xf, row_expo[k] - row_expo[i]);

  // Invalidate since b has changed, but R is not recomputed. We do operations on R[k], but not the
  // one to get the correct R[k]: the operations are the one mentionned in line 5 of Algorithm 3
  // [MSV, ISSAC'09].
  invalidate_row(k);
}

/* Taken from fplll/gso.cpp (commit 3d0d962)*/
template <class ZT, class FT> void MatHouseholder<ZT, FT>::row_add(int i, int j)
{
  b[i].add(b[j], n_known_cols);

  if (enable_transform)
  {
    u[i].add(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].sub(u_inv_t[i]);
  }
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::row_sub(int i, int j)
{
  b[i].sub(b[j], n_known_cols);

  if (enable_transform)
  {
    u[i].sub(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].add(u_inv_t[i]);
  }
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::row_addmul_si(int i, int j, long x)
{
  b[i].addmul_si(b[j], x, n_known_cols);

  if (enable_transform)
  {
    u[i].addmul_si(u[j], x);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si(u_inv_t[i], -x);
  }
}

template <class ZT, class FT>
void MatHouseholder<ZT, FT>::row_addmul_si_2exp(int i, int j, long x, long expo)
{
  b[i].addmul_si_2exp(b[j], x, expo, n_known_cols, ztmp0);

  if (enable_transform)
  {
    u[i].addmul_si_2exp(u[j], x, expo, ztmp0);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si_2exp(u_inv_t[i], -x, expo, ztmp0);
  }
}

template <class ZT, class FT>
void MatHouseholder<ZT, FT>::row_addmul_2exp(int i, int j, const ZT &x, long expo)
{
  // Cannot use ztmp0 here, since x is ztmp0. Use ztmp1 instead.
  b[i].addmul_2exp(b[j], x, expo, n_known_cols, ztmp1);

  if (enable_transform)
  {
    u[i].addmul_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
    {
      ZT minus_x;
      minus_x.neg(x);
      u_inv_t[j].addmul_2exp(u_inv_t[i], minus_x, expo, ztmp1);
    }
  }
}

template <class ZT, class FT>
void MatHouseholder<ZT, FT>::row_addmul_we(int i, int j, const FT &x, long expo_add)
{
  FPLLL_DEBUG_CHECK(j >= 0 && i == n_known_rows && j < i);
  FPLLL_DEBUG_CHECK(x.cmp(0.0) != 0);

  long expo;
  long lx = x.get_si_exp_we(expo, expo_add);

  if (expo == 0)
  {
    if (lx == 1)
      row_add(i, j);
    else if (lx == -1)
      row_sub(i, j);
    else if (lx != 0)
      row_addmul_si(i, j, lx);
  }
  else if (row_op_force_long)
    row_addmul_si_2exp(i, j, lx, expo);
  else
  {
    x.get_z_exp_we(ztmp0, expo, expo_add);
    row_addmul_2exp(i, j, ztmp0, expo);
  }

  // TODO: is it possible to combine this specialization with the condition on lx?
  // Cannot specialize depending on lx, since x contains the contribution of the 2^row_expo[i] and
  // 2^row_expo[j].
  if (x.cmp(1.0) == 0.0)
    R[i].add(R[j], i);
  else if (x.cmp(-1.0) == 0.0)
    R[i].sub(R[j], i);
  else
    R[i].addmul(R[j], x, i);
}

template class MatHouseholder<Z_NR<long>, FP_NR<double>>;
template class MatHouseholder<Z_NR<double>, FP_NR<double>>;
template class MatHouseholder<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class MatHouseholder<Z_NR<long>, FP_NR<long double>>;
template class MatHouseholder<Z_NR<double>, FP_NR<long double>>;
template class MatHouseholder<Z_NR<mpz_t>, FP_NR<long double>>;
#endif

#ifdef FPLLL_WITH_QD
template class MatHouseholder<Z_NR<long>, FP_NR<dd_real>>;
template class MatHouseholder<Z_NR<double>, FP_NR<dd_real>>;
template class MatHouseholder<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class MatHouseholder<Z_NR<long>, FP_NR<qd_real>>;
template class MatHouseholder<Z_NR<double>, FP_NR<qd_real>>;
template class MatHouseholder<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class MatHouseholder<Z_NR<long>, FP_NR<dpe_t>>;
template class MatHouseholder<Z_NR<double>, FP_NR<dpe_t>>;
template class MatHouseholder<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class MatHouseholder<Z_NR<long>, FP_NR<mpfr_t>>;
template class MatHouseholder<Z_NR<double>, FP_NR<mpfr_t>>;
template class MatHouseholder<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE
