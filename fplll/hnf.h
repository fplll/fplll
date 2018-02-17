/* Copyright (C) 2018 Arnaud Sipasseuth.

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

/**
  * WARNING : all HNFs here are considered lower triangular. Not upper.
  */

#ifndef FPLLL_HNF_H
#define FPLLL_HNF_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

/**
   * @brief The class performing HNF reduction.
   *
   * This class seems very similar to BKZReduction, for API consistency
**/

/* I'm bad at templating and can't compile with it so for now, without templates*/

template <class ZT> class HNFReduction
{
  /**
   * @brief Create a BKZObject
   *
   * @param m
   *    GSO object corresponding to the basis to be reduced
   * @param lll_obj
   *    LLL object associated to the same GSO object m
   * @param param
   *    parameter object (see bkz_param.h)
   */
public:
  HNFReduction(ZZ_mat<ZT> &b, HNFMethod method, Z_NR<ZT> &det)
      : basis(b), method(method), det(det){};
  ~HNFReduction() = default;

  /**
   * Status of reduction (see defs.h)
   */
  int status;
  /**
   * Rank of the lattice
   */
  int rank;
  /**
   * Basis of the lattice
   */
  ZZ_mat<ZT> &basis;
  /**
   * Method of reduction
   */
  HNFMethod method;
  /**
   * determinant of the lattice
   */
  Z_NR<ZT> &det;
  /**
     @brief Test the HNF form of basis.

     @return true on success.
  */
  bool is_reduced();
  /**
   * @brief Tests a couple matrix vector
   *    returns true if B is a HNF and if v is in L(B), else 0
   *
   * @param w
   *    vector for the membership test
   */
  bool in_lattice_given_hnf(const vector<Z_NR<ZT>> &w);
  /**
   * @brief Tests a couple matrix matrix
   *    returns 1 if B is a HNF and B is in HNF(A), else 0
   *
   * @param A
   *    basis of the lattice to be tested
   */
  bool in_lattice_given_hnf(const ZZ_mat<ZT> &A);
  /**
   * @brief Runs the main loop of hnf reduction.
   *
   * @return
   *    true if the reduction was successful, false otherwise.
   */
  void hnf();
  /**
   * @brief Compute the row rank of the matrix
   *
   */
  void compute_rank();

private:
  // Temporary data
  double cputime_start;
};

/**
 * @brief Performs hnf reduction using the multiple pgcd algorithm
 *
 * @param B
 *    basis of the lattice to be reduced
 */

template <class ZT> int hnf_xgcd_reduction(ZZ_mat<ZT> &B);

/**
 * @brief Performs hnf reduction using the classical absolute value algorithm
 *
 * @param B
 *    basis of the lattice to be reduced
 */

template <class ZT> int hnf_classical_reduction(ZZ_mat<ZT> &B);

/**
 * @brief Performs hnf reduction using the modulo determinant algorithm (Hafner McCurley)
 *    works only if the biggest square bottom-right matrix has a non-zero determinant
 *    and the number of rows is superior or equal to the number of columns
 *
 * @param B
 *    basis of the lattice to be reduced
 * @param D
 *    multiple of the determinant of the biggest square bottom-right matrix
 */

template <class ZT> int hnf_modular_reduction(ZZ_mat<ZT> &B, const Z_NR<ZT> D);

/**
 * @brief Performs hnf reduction using the minors algorithm (Kannan Bachem)
 *    good algorithm for full rank matrices, errors otherwise
 *
 * @param B
 *    basis of the lattice to be reduced
 */

template <class ZT> int hnf_minors_reduction(ZZ_mat<ZT> &B);

/**
 * @brief add a row to a HNF and reduce by the diagonal
 *
 * @param B
 *    initial HNF basis
 * @param v
 *    vector to be added
 */

template <class ZT> void hnf_addrow(ZZ_mat<ZT> &B, const vector<Z_NR<ZT>> &v);

/**
 * @brief performs HNF reduction, autoselect the best method (in development)
 *
 * @param B
 *    basis of the lattice to be reduced
 */

template <class ZT> int hnf_autoselect(ZZ_mat<ZT> &B);

/**
 * @brief Performs hnf reduction using the selected algorithm
 *
 * @param B
 *    basis of the lattice to be reduced
 * @param method
 *    reduction method
 * @param det
 *    determinant, might be useful
 */

template <class ZT> int hnf_reduction(ZZ_mat<ZT> &B, HNFMethod method, Z_NR<ZT> &det);

FPLLL_END_NAMESPACE

#endif /* FPLLL_HNF_H */
