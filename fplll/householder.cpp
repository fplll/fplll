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

template <class ZT, class FT> void MatHouseholder<ZT, FT>::update_R_last(int i)
{
  FT ftmp0, ftmp1, ftmp2;

  // sigma[i] = sign(r[1])
  sigma[i] = (R(i, i).cmp(0) < 0) ? -1.0 : 1.0;
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
      ftmp0.div(1.0, ftmp0);
      // Here, ftmp0 = 1 / sqrt(-s * vi[1])
      V(i, i).mul(V(i, i), ftmp0);
      R(i, i) = ftmp2;
      for (int k = i + 1; k < n; k++)
      {
        V(i, k).mul(R(i, k), ftmp0);
        R(i, k) = 0.0;
      }
      // Here, vi = vi / ftmp0 and ri[i..n] = (||r||, 0, 0, ..., 0)

      R_inverse_diag[i].div(1.0, ftmp2);
    }
    else
    {
      if (R(i, i).cmp(0.0) < 0)
        R(i, i).neg(R(i, i));

      for (int k = i + 1; k < n; k++)
      {
        // if enable_row_expo, R can be not correct at some point of the computation
        FPLLL_DEBUG_CHECK(enable_row_expo ? 1 : R(i, k).is_zero());
        R(i, k) = 0.0;
        V(i, k) = 0.0;
      }

      R_inverse_diag[i].div(1.0, R(i, i));
    }
  }
  else
  {
    for (int k = i; k < n; k++)
    {
      // if enable_row_expo, R can be not correct at some point of the computation
      FPLLL_DEBUG_CHECK(enable_row_expo ? 1 : R(i, k).is_zero());
      R(i, k) = 0.0;
      V(i, k) = 0.0;
    }

    // Result is inf.
    // TODO: set inf instead of doing the computation.
    R_inverse_diag[i].div(1.0, 0.0);
  }

  n_known_rows++;
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::update_R(int i, int last_j)
{
  // To update row i, we need to know n_known_rows rows
  FPLLL_DEBUG_CHECK(i <= n_known_rows);
  if (i == n_known_rows && !updated_R)
  {
    FT ftmp0, ftmp1;
    FPLLL_DEBUG_CHECK(last_j <= i + 1);

    n_known_cols = max(n_known_cols, init_row_size[i]);
    FPLLL_DEBUG_CHECK(last_j <= n_known_cols);

    int j, k;
    int j_stop = last_j == i + 1 ? i : last_j;

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
    }

    // Copy R[i] in bf[i] (while we copy b[i] in R[i])
    if (enable_bf)
    {
      for (j = 0; j < n_known_cols; j++)
        bf(i, j) = R(i, j);
      for (j = n_known_cols; j < n; j++)
        bf(i, j) = 0.0;
    }

    if (j_stop - 1 >= 0)
    {
      for (j = 0; j < j_stop - 1; j++)
      {
        // vj * ri[j..n]^T
        V[j].dot_product(ftmp1, R[i], j, n);

        //-vj * ri[j..n]^T
        ftmp1.neg(ftmp1);
        for (k = j; k < n; k++)
        {
          // ri[j..n] = ri[j..n] - (vj * ri[j..n]^T) * vj
          R(i, k).addmul(V(j, k), ftmp1);
        }
        // ri[j] = sigma[j] * ri[j]
        R(i, j).mul(sigma[j], R(i, j));

        // Copy R into R_history
        for (int k           = j; k < n; k++)
          R_history[i][j][k] = R(i, k);
      }

      V[j_stop - 1].dot_product(ftmp1, R[i], j_stop - 1, n);
      ftmp1.neg(ftmp1);
      R(i, j_stop - 1).addmul(V(j_stop - 1, j_stop - 1), ftmp1);
      R(i, j_stop - 1).mul(sigma[j_stop - 1], R(i, j_stop - 1));

      for (k = j_stop; k < n; k++)
        R(i, k).addmul(V(j_stop - 1, k), ftmp1);

      // Copy R into R_history
      for (k                    = j_stop; k < n; k++)
        R_history[i][j_stop][k] = R(i, k);
    }

    if (last_j == i + 1)
      update_R_last(i);
  }
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::swap(int i, int j)
{
  FPLLL_DEBUG_CHECK(0 <= i && i < j && j < d);

  invalidate_row(i);

  b.swap_rows(i, j);
  if (enable_bf)
    bf.swap_rows(i, j);
  R.swap_rows(i, j);
  V.swap_rows(i, j);
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
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::addmul_b_rows(int k, vector<FT> xf)
{
  FPLLL_DEBUG_CHECK(k > 0 && k < d);

  for (int i = 0; i < k; i++)
  {
    row_addmul_we(k, i, xf[i], row_expo[k] - row_expo[i]);

    if (enable_bf)
      bf[k].addmul(bf[i], xf[i], n_known_cols);
  }

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
  ZT ztmp0;
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
  ZT ztmp0;
  b[i].addmul_2exp(b[j], x, expo, n_known_cols, ztmp0);
  if (enable_transform)
  {
    u[i].addmul_2exp(u[j], x, expo, ztmp0);
    if (enable_inverse_transform)
    {
      ZT minus_x;
      minus_x.neg(x);
      u_inv_t[j].addmul_2exp(u_inv_t[i], minus_x, expo, ztmp0);
    }
  }
}

template <class ZT, class FT>
void MatHouseholder<ZT, FT>::row_addmul_we(int i, int j, const FT &x, long expo_add)
{
  FPLLL_DEBUG_CHECK(j >= 0 && i == n_known_rows && j < i);

  ZT ztmp0;
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
  {
    row_addmul_si_2exp(i, j, lx, expo);
  }
  else
  {
    x.get_z_exp_we(ztmp0, expo, expo_add);
    row_addmul_2exp(i, j, ztmp0, expo);
  }
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
