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

#ifndef FPLLL_ENUMERATE_EXT_H
#define FPLLL_ENUMERATE_EXT_H

#include "enumerate_ext_api.h"

#include <fplll/enum/enumerate_base.h>
#include <fplll/enum/evaluator.h>
#include <fplll/gso_interface.h>

FPLLL_BEGIN_NAMESPACE

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
// typedef void(extenum_cb_set_config)(enumf *mu, size_t mudim, bool mutranspose, enumf *rdiag,
//                                    enumf *pruning);
using ::extenum_cb_set_config;

/**
 * Callback function given to external enumeration library.
 * Passes a new solution and its length to Evaluator, returning the new enumeration bound.
 * @param[in] dist - the norm of the new solution.
 * @param[in] sol - a pointer to the new solution.
 * @return The new enumeration bound.
 */
// typedef enumf(extenum_cb_process_sol)(enumf dist, enumf *sol);
using ::extenum_cb_process_sol;

/**
 * Callback function given to external enumeration library.
 *
 * Pass a subsolution and its partial length to Evaluator.
 */
// typedef void(extenum_cb_process_subsol)(enumf dist, enumf *subsol, int offset);
using ::extenum_cb_process_subsol;

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
// typedef array<uint64_t, FPLLL_EXTENUM_MAX_EXTENUM_DIM>(extenum_fc_enumerate)(
//    const int dim, enumf maxdist, std::function<extenum_cb_set_config> cbfunc,
//    std::function<extenum_cb_process_sol> cbsol, std::function<extenum_cb_process_subsol>
//    cbsubsol, bool dual /*=false*/, bool findsubsols /*=false*/
//);
using ::extenum_fc_enumerate;

/* set & get external enumerator. If extenum = nullptr then this interface is disabled,
                                  and fplll will use the internal enumerator.
                                  Otherwise, fplll will use the enumeration function pointed to
                                  by extenum.
*/
void set_external_enumerator(std::function<extenum_fc_enumerate> extenum = nullptr);
std::function<extenum_fc_enumerate> get_external_enumerator();

template <typename ZT, typename FT> class ExternalEnumeration
{
public:
  ExternalEnumeration(MatGSOInterface<ZT, FT> &gso, Evaluator<FT> &evaluator)
      : _gso(gso), _evaluator(evaluator)
  {
  }

  bool enumerate(int first, int last, FT &fmaxdist, long fmaxdistexpo,
                 const vector<enumf> &pruning = vector<enumf>(), bool dual = false);

  // get_nodes. This returns the number of nodes visited by the external enumeration process.
  // If this returns 0, then fplll will fall back to the internal enumerator.
  inline uint64_t get_nodes(const int level = -1) const
  {
    if (level == -1)
    {
      return std::accumulate(_nodes.begin(), _nodes.end(), 0);
    }
    return _nodes[level];
  }

  inline std::array<uint64_t, FPLLL_EXTENUM_MAX_EXTENUM_DIM> get_nodes_array() const { return _nodes; }

private:
  void callback_set_config(enumf *mu, size_t mudim, bool mutranspose, enumf *rdiag, enumf *pruning);

  enumf callback_process_sol(enumf dist, enumf *sol);

  void callback_process_subsol(enumf dist, enumf *subsol, int offset);

  MatGSOInterface<ZT, FT> &_gso;
  Evaluator<FT> &_evaluator;
  vector<enumf> _pruning;
  long _normexp;

  array<uint64_t, FPLLL_EXTENUM_MAX_EXTENUM_DIM> _nodes;
  bool _dual;
  int _d, _first;
  enumf _maxdist;
  vector<FT> _fx;
};

FPLLL_END_NAMESPACE

#endif
