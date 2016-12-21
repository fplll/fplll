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

#include <array>
#include <fplll/enum/enumerate_base.h>
#include <fplll/enum/evaluator.h>
#include <fplll/gso.h>
#include <functional>
#include <memory>

FPLLL_BEGIN_NAMESPACE

/* function callback API for external enumeration library (extenum) */

/**
 * Callback function given to external enumeration library.
 *
 * You have to pass a pointer to an array 'enumf mu[mudim][mudim]'.
 * When flag mutranspose is true then mutransposed is actually stored there.
 * You have to pass a pointer to an array 'enumf rdiag[mudim]'
 * and an array 'enumf pruning[mudim]'.
 * Note: for dual SVP you also get mu and rdiag as is,
 * so the external library must make the respective changes to Mu and Rdiag itself.
 */
typedef void(extenum_cb_set_config)(enumf *mu, size_t mudim, bool mutranspose, enumf *rdiag,
                                    enumf *pruning);

/**
 * Callback function given to external enumeration library.
 *
 * Pass a new solution and its length to Evaluator, it returns the new enumeration bound.
 */
typedef enumf(extenum_cb_process_sol)(enumf dist, enumf *sol);

/**
 * Callback function given to external enumeration library.
 *
 * Pass a subsolution and its partial length to Evaluator.
 */
typedef void(extenum_cb_process_subsol)(enumf dist, enumf *subsol, int offset);

/**
 * External enumeration function prototype.
 *
 * @param dim         enumeration dimension
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
typedef uint64_t(extenum_fc_enumerate)(int dim, enumf maxdist,
                                       std::function<extenum_cb_set_config> cbfunc,
                                       std::function<extenum_cb_process_sol> cbsol,
                                       std::function<extenum_cb_process_subsol> cbsubsol,
                                       bool dual /*=false*/, bool findsubsols /*=false*/
                                       );

// set & get external enumerator (nullptr => disabled)
void set_external_enumerator(std::function<extenum_fc_enumerate> extenum = nullptr);
std::function<extenum_fc_enumerate> get_external_enumerator();

template <typename FT> class ExternalEnumeration
{
public:
  ExternalEnumeration(MatGSO<Integer, FT> &gso, Evaluator<FT> &evaluator)
      : _gso(gso), _evaluator(evaluator)
  {
  }

  bool enumerate(int first, int last, FT &fmaxdist, long fmaxdistexpo,
                 const vector<enumf> &pruning = vector<enumf>(), bool dual = false);

  inline uint64_t get_nodes() const { return _nodes; }

private:
  void callback_set_config(enumf *mu, size_t mudim, bool mutranspose, enumf *rdiag, enumf *pruning);

  enumf callback_process_sol(enumf dist, enumf *sol);

  void callback_process_subsol(enumf dist, enumf *subsol, int offset);

  MatGSO<Integer, FT> &_gso;
  Evaluator<FT> &_evaluator;
  vector<enumf> _pruning;
  long _normexp;
  uint64_t _nodes;

  bool _dual;
  int _d, _first;
  enumf _maxdist;
  vector<FT> _fx;
};

FPLLL_END_NAMESPACE

#endif
