/*
   (C) 2016 Marc Stevens.

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

#ifndef FPLLL_ENUMERATE_EXT_API_H
#define FPLLL_ENUMERATE_EXT_API_H

#include <array>
#include <functional>
#include <memory>

typedef double fplll_extenum_enumf;
#define FPLLL_EXTENUM_MAX_EXTENUM_DIM 256

/* function callback API for external enumeration library (extenum) */

/**
 * Callback function given to external enumeration library.
 *
 * @param[out] mu - this is a pointer to an array 'enumf mu[mudim][mudim]'.
 * 	       Upon return, this array will contain
 * 	       the mu values for the lattice basis.
 * 	       Note that the array pointed to by this argument needs to be contiguous,
 * 	       otherwise there will be write errors.
 * 	       This means that the only acceptable arguments are pointers to
 * 	       variable-length arrays, 1D std::vector of dimension mudim*mudim, or
 * 	       1D raw(resp. std::array) arrays of dimension mudim*mudim.
 * @param[in]  mudim - the number of rows(resp. columns) in the mu array.
 * @param[out] mutranspose - when true, mutranspose is stored in the mu param.
 *                           Otherwise, mu is stored in the mu param.
 * 			Storing mutranspose allows for more efficient memory access, as accesses
 * will be contiguous.
 * @param[out] rdiag - a pointer to an array 'enumf rdiag[mudim]'.
 * 	               Upon return, this will contain the squared norms of the Gram-Schmidt vectors.
 * @param[out] pruning - a pointer to an array 'enumf pruning[mudim]'. Upon return, this will
 * contain the pruning coefficients for enumeration. In rigorous enumeration, this array will
 * consist solely of 1's.
 */
typedef void(extenum_cb_set_config)(fplll_extenum_enumf *mu, std::size_t mudim, bool mutranspose,
                                    fplll_extenum_enumf *rdiag, fplll_extenum_enumf *pruning);

/**
 * Callback function given to external enumeration library.
 * Passes a new solution and its length to Evaluator, returning the new enumeration bound.
 * @param[in] dist - the norm of the new solution.
 * @param[in] sol - a pointer to the new solution.
 * @return The new enumeration bound.
 */
typedef fplll_extenum_enumf(extenum_cb_process_sol)(fplll_extenum_enumf dist,
                                                    fplll_extenum_enumf *sol);

/**
 * Callback function given to external enumeration library.
 *
 * Pass a subsolution and its partial length to Evaluator.
 */
typedef void(extenum_cb_process_subsol)(fplll_extenum_enumf dist, fplll_extenum_enumf *subsol,
                                        int offset);

/**
 * External enumeration function prototype.
 *
 * @param dim         enumeration dimension.
 * @param maxdist     initial enumeration bound.
 * @param cbfunc      given callback function to get mu, rdiag, pruning
 * @param cbsol       given callback function to pass solution and its length to Evaluator,
 *                    it returns new enumeration bound
 * @param cbsubsol    given callback function to pass subsolution and its length to Evaluator
 * @param dual        do dual SVP enumeration
 * @param findsubsols find subsolutions and pass them to Evaluator
 * @return number of nodes visited.
 *         Or ~uint64_t(0) when instance is not supported
 *         in which case fplll falls back to its own enumeration.
 */
typedef std::array<std::uint64_t, FPLLL_EXTENUM_MAX_EXTENUM_DIM>(extenum_fc_enumerate)(
    const int dim, fplll_extenum_enumf maxdist, std::function<extenum_cb_set_config> cbfunc,
    std::function<extenum_cb_process_sol> cbsol, std::function<extenum_cb_process_subsol> cbsubsol,
    bool dual /*=false*/, bool findsubsols /*=false*/
);

#endif
