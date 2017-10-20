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
  MatHouseholder(Matrix<ZT> &arg_b) : b(arg_b)
  {
    alloc_dim = 0;
    d         = b.get_rows();
    n         = b.get_cols();
    size_increased();
  }

  ~MatHouseholder() {}

  /**
   * Returns f = householder_r(i, j).
   *
   * Returns reference to `f`.
   */
  inline FT &get_R(FT &f, int i, int j);

  /**
   * Returns the R matrix
   */
  const Matrix<FT> &get_R_matrix() { return R; }

  /**
   * Apply Householder transformation on row i.
   */
  void update_R_row(int i);

  /**
   * Full computation of the matrix R.
   */
  inline void update_R();

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
   * Floating-point representation of the basis.
   */
  Matrix<FT> bf;

  /* Allocates matrices and arrays whose size depends on d.
   When enable_int_gram=false, initializes bf. */
  void size_increased();

  void discover_row();

  /* Upates the i-th row of bf. It does not invalidate anything, so the caller
     must take into account that it might change row_expo. */
  void update_bf(int i);

  /**
   * bf = R * q_householder.
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

  /* Used by update_gso_row (+ update_gso), get_max_mu_exp and row_addmul_we. */
  FT ftmp0, ftmp1, ftmp2;

  int alloc_dim;
};

template <class ZT, class FT> inline FT &MatHouseholder<ZT, FT>::get_R(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && j >= 0 && j <= i);
  f = R(i, j);
  return f;
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::update_R()
{
  for (int i = 0; i < d; i++)
    update_R_row(i);
}

FPLLL_END_NAMESPACE

#endif  // FPLLL_HOUSEHOLDER_H
