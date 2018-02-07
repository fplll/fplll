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

class HNFReduction
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
  HNFReduction(ZZ_mat<mpz_t> &b, HNFMethod method, Z_NR<mpz_t> &det):
    b(b),method(method),det(det){};
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
  ZZ_mat<mpz_t> &b;
  /**
   * Method of reduction
   */
  HNFMethod method;
  /**
   * determinant of the lattice
   */
  Z_NR<mpz_t> &det;
  /**
   * modify the status of reduction (see defs.h)
   */
  void is_reduced();
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
  bool set_status(int new_status);

  // Temporary data
  double cputime_start;
};

/**
   @brief Test the HNF form of B.

   @param A
   @return zero on success.
*/

int is_hnf_reduced(const ZZ_mat<mpz_t> &B);

/**
 * @brief Tests a couple matrix vector
 *    returns 1 if B is not a HNF or -1 if v is not in L(B), else 0
 *
 * @param B
 *    HNF basis of the lattice to be tested
 * @param v
 *    vector for the membership test
 */

int in_lattice_given_hnf(const ZZ_mat<mpz_t> &B, const vector<Z_NR<mpz_t>> &v);

/**
 * @brief Tests a couple matrix matrix
 *    returns 1 if B is not a HNF or -1 if B is not HNF(A), else 0
 *
 * @param B
 *    HNF basis of the lattice to be tested
 * @param A
 *    basis of the lattice to be tested
 */

int in_lattice_given_hnf(const ZZ_mat<mpz_t> &B, const ZZ_mat<mpz_t> &A);

/**
 * @brief Performs hnf reduction using the multiple pgcd algorithm
 *
 * @param B
 *    basis of the lattice to be reduced
 */

int hnf_xgcd_reduction(ZZ_mat<mpz_t> &B);

/**
 * @brief Performs hnf reduction using the classical absolute value algorithm
 *
 * @param B
 *    basis of the lattice to be reduced
 */

int hnf_classical_reduction(ZZ_mat<mpz_t> &B);

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

int hnf_modular_reduction(ZZ_mat<mpz_t> &B, const Z_NR<mpz_t> D);

/**
 * @brief Performs hnf reduction using the minors algorithm (Kannan Bachem)
 *    good algorithm for full rank matrices, errors otherwise
 *
 * @param B
 *    basis of the lattice to be reduced
 */

int hnf_minors_reduction(ZZ_mat<mpz_t> &B);

/**
 * @brief Performs hnf reduction using the selected algorithm
 *    This only works for full-rank matrices
 *
 * @param B
 *    basis of the lattice to be reduced
 * @param method
 *    reduction method
 */

int hnf_reduction(ZZ_mat<mpz_t> &B, HNFMethod method, Z_NR<mpz_t> &det);

/**
 * @brief performs HNF reduction, autoselect the best method (in development)
 *
 * @param B
 *    basis of the lattice to be reduced
 */

int hnf_autoselect(ZZ_mat<mpz_t> &B);

FPLLL_END_NAMESPACE

#endif /* FPLLL_HNF_H */
