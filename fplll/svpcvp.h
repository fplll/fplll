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
 * @brief SVP solver - computes a shortest vector in a lattice.
 *
 * @param[out] b
 * basis vectors must be linearly independent and LLL-reduced.
 *
 * @param[out] sol_coord
 * coordinates of the shortest vector with respect to the basis b.
 *
 * @param[in] method = SVPM_PROVED
 * The result is guaranteed if method = SVPM_PROVED opposed to SVPM_FAST.
 *
 * @param[in] flags = SVP_DEFAULT
 * one of SVP_DEFAULT, SVP_VERBOSE, SVP_OVERRIDE_BND, and SVP_DUAL with:
 * SVP_DEFAULT the standard behaviour,
 * SVP_VERBOSE writing additional information to console,
 * SVP_OVERRIDE_BND computing and checking initial max error bound, and
 * SVP_DUAL computing the shortest vector of the dual lattice (in the dual basis) instead.
 *
 * @return result Success or failure (due to numerical instability).
 */
int shortest_vector(IntMatrix &b, IntVect &sol_coord, SVPMethod method = SVPM_PROVED,
                    int flags = SVP_DEFAULT);

int shortest_vector_pruning(IntMatrix &b, IntVect &sol_coord, const vector<double> &pruning,
                            int flags = SVP_DEFAULT);

int shortest_vector_pruning(IntMatrix &b, IntVect &sol_coord, vector<IntVect> &subsol_coord,
                            vector<double> &subsol_dist, const vector<double> &pruning,
                            int flags = SVP_DEFAULT);

int shortest_vector_pruning(IntMatrix &b, IntVect &sol_coord, vector<IntVect> &auxsol_coord,
                            vector<double> &auxsol_dist, const int max_aux_sols,
                            const vector<double> &pruning, int flags = SVP_DEFAULT);
/**
 * Computes a closest vector of a lattice to a target.
 * The vectors must be linearly independent and the basis must be LLL-reduced
 * with delta=LLL_DEF_DELTA and eta=LLL_DEF_ETA.
 * The result is guaranteed if method = CVPM_PROVED.
 */
/**
 * @brief CVP solver - computes a closest vector in a lattice with respect to a target.
 *
 * @param[out] b
 * basis vectors must be linearly independent and LLL-reduced.
 *
 * @param[out] sol_coord
 * coordinates of the closest lattice vector with respect to the basis b.
 *
 * @param[in] int_target
 * coordinates of the target vector.
 *
 * @param[in] method = CVPM_FAST
 * or, alternatively CVPM_PROVED which finds a guaranteed solution and resets
 * enumeration below depth with maximal r_i.
 *
 * @param[in] flags = SVP_DEFAULT
 * or, alternatively SVP_VERBOSE(!), both taken from SVPFlags.
 *
 * @return result Success or failure (due to numerical instability).
 */
int closest_vector(IntMatrix &b, const IntVect &int_target, vector<Integer> &sol_coord,
                   int method = CVPM_FAST, int flags = CVP_DEFAULT);

FPLLL_END_NAMESPACE

#endif
