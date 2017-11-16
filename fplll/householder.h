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

#ifndef FPLLL_HOUSEHOLDER_H
#define FPLLL_HOUSEHOLDER_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

enum MatHouseholderFlags
{
  HOUSEHOLDER_DEFAULT = 0
};

/**
 * MatHouseholder provides an interface for performing elementary operations on a basis
 * and computing its R matrix using Householder transformations.
 */
template <class ZT, class FT> class MatHouseholder
{
public:
  /**
   * Constructor.
   * The precision of FT must be defined before creating an instance of the
   * class and must remain the same until the object is destroyed (or no longer
   * needed).
   * @param b
   *   The matrix on which row operations are performed. It must not be empty.
   */
  MatHouseholder(Matrix<ZT> &arg_b, int flags) : b(arg_b)
  {
    d            = b.get_rows();
    n            = b.get_cols();
    n_known_rows = 0;
    sigma.resize(d);
    R.resize(d, n);
    V.resize(d, n);

#ifdef DEBUG
    for (int i = 0; i < d; i++)
    {
      for (int j = 0; j < n; j++)
      {
        V(i, j).set_nan();
      }
    }
#endif  // DEBUG
  }

  ~MatHouseholder() {}

  /**
   * Returns f = R(i, j).
   *
   * Returns reference to `f`.
   */
  inline FT &get_R(FT &f, int i, int j);

  inline void set_R(FT &f, int i, int j);

  /**
   * Returns R[i].
   */
  inline MatrixRow<FT> get_R(int i);

  /**
   * Returns the R matrix
   */
  const Matrix<FT> &get_R() { return R; }

  /**
   * Returns b[i].
   */
  MatrixRow<ZT> get_b(int i);

  /**
   * Returns the b matrix
   */
  const Matrix<ZT> &get_b() { return b; }

  /**
   * Apply Householder transformation on row i, from cols 0 to last_j.
   * Restriction: last_j == i - 1 or i.
   */
  void update_R_row(int i, int last_j);

  /**
   * Apply Householder transformation on row i.
   */
  void update_R_row(int i);

  /**
   * Full computation of the matrix R.
   */
  inline void update_R();

  inline int get_d() { return d; }
  inline int get_n() { return n; }

  /**
   * Norm square of b[k].
   */
  inline void norm_square_b_row(FT &f, int k);

  /**
   * Truncated norm square of R[k], with coefficients of R[k][0..end-1].
   */
  inline void norm_square_R_row(FT &f, int k, int end);

  /**
   * b[k] = b[k] - sum_{i = 0}^{k - 1}(x[i] * b[i])
   */
  void add_mul_b_rows(int k, ZT *x);

  /**
   * Swap row i and j of b.
   */
  inline void swap(int i, int j);

  /**
   * Invalidate row k to row n_known_rows - 1.
   * Update n_known_rows to k.
   */
  inline void invalidate_row(int k);

private:
  /**
   * Number of rows of b (dimension of the lattice).
   */
  int d;

  /**
   * Number of columns of b (dimension of the lattice).
   */
  int n;

  /**
   * Basis of the lattice
   */
  Matrix<ZT> &b;

  /**
   * b = R * q_householder.
   * R is lower triangular and the diagonal coefficient are >= 0.
   */
  Matrix<FT> R;

  /**
   * Vector v following [MSV, ISSAC'09].
   */
  Matrix<FT> V;

  /**
   * Sigma values following [MSV, ISSAC'09].
   */
  vector<FT> sigma;

  /**
   * R[i] is invalid for i >= n_known_rows.
   */
  int n_known_rows;
};

template <class ZT, class FT> inline FT &MatHouseholder<ZT, FT>::get_R(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d && j >= 0 && j <= i);
  f = R(i, j);
  return f;
}

template <class ZT, class FT> inline MatrixRow<FT> MatHouseholder<ZT, FT>::get_R(int i)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d);
  return R[i];
}

template <class ZT, class FT> MatrixRow<ZT> MatHouseholder<ZT, FT>::get_b(int i)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d);
  return b[i];
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::set_R(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d && i >= j && j >= 0);
  FPLLL_DEBUG_CHECK(j <= i);
  R(i, j) = f;
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::update_R_row(int i)
{
  update_R_row(i, i);
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::update_R()
{
  for (int i = 0; i < d; i++)
    update_R_row(i);
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::norm_square_b_row(FT &f, int k)
{
  FPLLL_DEBUG_CHECK(k >= 0 && k < d);
  ZT ztmp0;
  b[k].dot_product(ztmp0, b[k], 0, n);
  f.set_z(ztmp0);
}

template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::norm_square_R_row(FT &f, int k, int end)
{
  FPLLL_DEBUG_CHECK(k >= 0 && k < d);
  FPLLL_DEBUG_CHECK(0 <= end && end <= k);
  if (end == 0)
  {
    f = 0.0;
    FPLLL_DEBUG_CHECK(f.is_zero());
  }
  else
  {
    R[k].dot_product(f, R[k], 0, end);
  }
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::invalidate_row(int k)
{
#ifdef DEBUG
  if (k < n_known_rows)
    n_known_rows = k;
#endif  // DEBUG
}

template <class ZT, class FT> void MatHouseholder<ZT, FT>::swap(int i, int j)
{
  FPLLL_DEBUG_CHECK(0 <= i && i < j && j < d);

  invalidate_row(i);

  b.swap_rows(i, j);
}

FPLLL_END_NAMESPACE

#endif  // FPLLL_HOUSEHOLDER_H
