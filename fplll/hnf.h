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

int in_hnf(const ZZ_mat<mpz_t> &B, const vector<Z_NR<mpz_t>> &v);

/**
 * @brief Tests a couple matrix matrix
 *    returns 1 if B is not a HNF or -1 if B is not HNF(A), else 0
 *
 * @param B
 *    HNF basis of the lattice to be tested
 * @param A
 *    basis of the lattice to be tested
 */

int in_hnf(const ZZ_mat<mpz_t> &B, const ZZ_mat<mpz_t> &A);

/**
 * @brief Performs hnf reduction using the multiple pgcd algorithm
 *
 * @param B
 *    basis of the lattice to be reduced
 */

int hnf_xgcd_reduction(ZZ_mat<mpz_t> &B);

/**
 * @brief Performs hnf reduction using the selected algorithm
 *
 * @param B
 *    basis of the lattice to be reduced
 * @param method
 *    reduction method
 */

int hnf(ZZ_mat<mpz_t> &B, HNFMethod method);

/**
 * @brief performs HNF reduction, autoselect the best method (in development)
 *
 * @param B
 *    basis of the lattice to be reduced*
 */

int hnf_autoselect(ZZ_mat<mpz_t> &B);

FPLLL_END_NAMESPACE

#endif /* FPLLL_HNF_H */
