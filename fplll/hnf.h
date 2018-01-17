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

#ifndef FPLLL_HNF_H
#define FPLLL_HNF_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

/**
 * @brief Tests a couple matrix vector
 *    returns 1 if B is not a HNF or -1 if v is not in L(B), else 0
 *
 * @param B
 *    HNF basis of the lattice to be tested
 * @param v
 *    vector for the membership test
 */

int in_hnf(const ZZ_mat<mpz_t> &B, const vector<mpz_t> &v);

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
 * @param U
 *    transformation matrix (don't pass a parameter to ignore this option)
 */

// template <class ZT> int hnf_xgcd_reduction(ZZ_mat<ZT> &B, ZZ_mat<ZT> &U);
int hnf_xgcd_reduction(ZZ_mat<mpz_t> &B);

FPLLL_END_NAMESPACE

#endif /* FPLLL_HNF_H */
