/* Copyright (C) 2008-2011 Xavier Pujol.

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
/** \file svpcvp.h
    Shortest vector problem. */

#ifndef FPLLL_SVPCVP_H
#define FPLLL_SVPCVP_H

#include "util.h"

FPLLL_BEGIN_NAMESPACE

/**
 * Computes a shortest vector of a lattice.
 * The vectors must be linearly independant and the basis must be LLL-reduced
 * with delta=LLL_DEF_DELTA and eta=LLL_DEF_ETA.
 * The result is guaranteed if method = SVPM_PROVED.
 */
int shortest_vector(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord,
                    SVPMethod method = SVPM_PROVED, int flags = SVP_DEFAULT);

int shortest_vector_pruning(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord,
                            const vector<double> &pruning, int flags = SVP_DEFAULT);

int shortest_vector_pruning(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &subsol_coord, vector<double> &subsol_dist,
                            const vector<double> &pruning, int flags = SVP_DEFAULT);

int shortest_vector_pruning(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &auxsol_coord, vector<double> &auxsol_dist,
                            const int max_aux_sols, const vector<double> &pruning,
                            int flags = SVP_DEFAULT);



int shortest_vector(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso, vector<Z_NR<mpz_t>> &sol_coord,
                    SVPMethod method = SVPM_PROVED, int flags = SVP_DEFAULT);

int shortest_vector_pruning(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso, vector<Z_NR<mpz_t>> &sol_coord,
                            const vector<double> &pruning, int flags = SVP_DEFAULT);

int shortest_vector_pruning(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso, vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &subsol_coord, vector<double> &subsol_dist,
                            const vector<double> &pruning, int flags = SVP_DEFAULT);

int shortest_vector_pruning(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso, vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &auxsol_coord, vector<double> &auxsol_dist,
                            const int max_aux_sols, const vector<double> &pruning,
                            int flags = SVP_DEFAULT);
/**
 * Computes a closest vector of a lattice to a target.
 * The vectors must be linearly independant and the basis must be LLL-reduced
 * with delta=LLL_DEF_DELTA and eta=LLL_DEF_ETA.
 * The result is guaranteed if method = CVPM_PROVED.
 */
int closest_vector(ZZ_mat<mpz_t> &b, const vector<Z_NR<mpz_t>> &int_target,
                   vector<Z_NR<mpz_t>> &sol_coord, int method = CVPM_FAST, int flags = CVP_DEFAULT);

FPLLL_END_NAMESPACE

#endif
