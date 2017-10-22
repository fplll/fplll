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

template <class ZT, class FT> void MatHouseholder<ZT, FT>::update_R_row(int i, int last_j)
{
  // Restriction on last_j
  FPLLL_DEBUG_CHECK(last_j == i || last_j == i - 1);
  // To update i, we need to know n_known_rows rows
  FPLLL_DEBUG_CHECK(i <= n_known_rows);

  int j, k;
  if (valid_cols[i] == 0)
  {
    for (j = 0; j < n; j++)
    {
      R(i, j) = bf(i, j);
    }
  }
  int j_stop = last_j == i ? i - 1 : last_j;
  for (j = valid_cols[i]; j <= j_stop; j++)
  {
    // vj * ri[j..n]^T
    dot_product(ftmp1, V[j], R[i], j, n);
    //-vj * ri[j..n]^T
    ftmp1.neg(ftmp1);
    for (k = j; k < n; k++)
    {
      // ri[j..n] = ri[j..n] - (vj * ri[j..n]^T) * vj
      R(i, k).addmul(V(j, k), ftmp1);
    }
    // ri[j] = sigma[j] * ri[j]
    R(i, j).mul(sigma[j], R(i, j));
  }
  if (last_j == i)
  {
    // sigma[i] = sign(r[1])
    sigma[i] = (R(i, i).cmp(0) < 0) ? -1.0 : 1.0;
    // r^T * r
    dot_product(ftmp1, R[i], R[i], i, n);
    if (ftmp1.cmp(0) != 0)
    {
      ftmp2.sqrt(ftmp1);
      // s = sigma[i] * ||r|| = sigma[i] * sqrt(r * r^T)
      ftmp0.mul(sigma[i], ftmp2);
      V(i, i).mul(R(i, i), R(i, i));
      V(i, i).sub(V(i, i), ftmp1);
      ftmp1.add(R(i, i), ftmp0);
      V(i, i).div(V(i, i), ftmp1);
      // Here, vi[1] = (-sum(r[j]^2, j, 2, n-i+1) / (r[1] + s)
      if (V(i, i).cmp(0) != 0)
      {
        ftmp0.neg(ftmp0);
        ftmp0.mul(ftmp0, V(i, i));
        ftmp0.sqrt(ftmp0);
        ftmp0.div(1.0, ftmp0);
        // Here, ftmp0 = 1 / sqrt(-s * vi[1])
        V(i, i).mul(V(i, i), ftmp0);
        R(i, i) = ftmp2;
        for (k = i + 1; k < n; k++)
        {
          V(i, k).mul(R(i, k), ftmp0);
#ifdef DEBUG
          R(i, k) = 0.0;
          FPLLL_DEBUG_CHECK(R(i, k).is_zero());
#endif  // DEBUG
        }
        // Here, vi = vi / ftmp0 and ri[i..n] = (||r||, 0, 0, ..., 0)
      }
    }
    n_known_rows++;
  }
  valid_cols[i] = last_j + 1;
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
