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

#include "enumerate_ext.h"
#include "../fplll_config.h"

#ifdef FPLLL_WITH_PARALLEL_ENUM
#include "../enum-parallel/enumlib.h"
#endif

FPLLL_BEGIN_NAMESPACE

// set & get external enumerator (nullptr => disabled)
#ifdef FPLLL_WITH_PARALLEL_ENUM
std::function<extenum_fc_enumerate> fplll_extenum = enumlib::enumlib_enumerate;
#else
std::function<extenum_fc_enumerate> fplll_extenum = nullptr;
#endif

void set_external_enumerator(std::function<extenum_fc_enumerate> extenum)
{
  fplll_extenum = extenum;
}

std::function<extenum_fc_enumerate> get_external_enumerator()
{
	return fplll_extenum;
}

template <typename ZT, typename FT>
bool ExternalEnumeration<ZT, FT>::enumerate(int first, int last, FT &fmaxdist, long fmaxdistexpo,
                                            const vector<enumf> &pruning, bool dual)
{
  using namespace std::placeholders;
  if (fplll_extenum == nullptr)
    return false;
  if (last == -1)
    last = _gso.d;

  _first   = first;
  _dual    = dual;
  _pruning = pruning;
  _d       = last - _first;
  _fx.resize(_d);

  FPLLL_CHECK(_pruning.empty() || int(_pruning.size()) == _d,
              "ExternalEnumeration: non-empty pruning vector dimension does not match");

  FT fr, fmu, fmaxdistnorm;
  long rexpo;
  _normexp = -1;
  for (int i = 0; i < _d; ++i)
  {
    fr       = _gso.get_r_exp(i + first, i + first, rexpo);
    _normexp = max(_normexp, rexpo + fr.exponent());
  }
  fmaxdistnorm.mul_2si(fmaxdist, dual ? _normexp - fmaxdistexpo : fmaxdistexpo - _normexp);

  _maxdist = fmaxdistnorm.get_d(GMP_RNDU);
  _evaluator.set_normexp(_normexp);

  // clang-format off
  _nodes = fplll_extenum(_d, _maxdist,
                         std::bind(&ExternalEnumeration<ZT,FT>::callback_set_config, this, _1, _2, _3, _4, _5),
                         std::bind(&ExternalEnumeration<ZT,FT>::callback_process_sol, this, _1, _2),
                         std::bind(&ExternalEnumeration<ZT,FT>::callback_process_subsol, this, _1, _2, _3),
               _dual, _evaluator.findsubsols
               );
  // clang-format on
  return _nodes != ~uint64_t(0);
}

template <typename ZT, typename FT>
void ExternalEnumeration<ZT, FT>::callback_set_config(enumf *mu, size_t mudim, bool mutranspose,
                                                      enumf *rdiag, enumf *pruning)
{
  FT fr, fmu;
  long rexpo;

  for (int i = 0; i < _d; ++i)
  {
    fr = _gso.get_r_exp(i + _first, i + _first, rexpo);
    fr.mul_2si(fr, rexpo - _normexp);
    rdiag[i] = fr.get_d();
  }

  if (mutranspose)
  {
    size_t offs = 0;
    for (int i = 0; i < _d; ++i, offs += mudim)
    {
      for (int j = 0; j < _d; ++j)
      {
        _gso.get_mu(fmu, j + _first, i + _first);
        /* mu[i][j]= */
        mu[offs + j] = fmu.get_d();
      }
    }
  }
  else
  {
    size_t offs = 0;
    for (int j = 0; j < _d; ++j, offs += mudim)
    {
      for (int i = 0; i < _d; ++i)
      {
        _gso.get_mu(fmu, j + _first, i + _first);
        /* mu[j][i] = */
        mu[offs + i] = fmu.get_d();
      }
    }
  }

  if (_pruning.empty())
  {
    for (int i = 0; i < _d; ++i)
      pruning[i] = 1.0;
  }
  else
  {
    for (int i = 0; i < _d; ++i)
      pruning[i] = _pruning[i];
  }
}

template <typename ZT, typename FT>
enumf ExternalEnumeration<ZT, FT>::callback_process_sol(enumf dist, enumf *sol)
{
  for (int i = 0; i < _d; ++i)
    _fx[i] = sol[i];
  _evaluator.eval_sol(_fx, dist, _maxdist);
  return _maxdist;
}

template <typename ZT, typename FT>
void ExternalEnumeration<ZT, FT>::callback_process_subsol(enumf dist, enumf *subsol, int offset)
{
  for (int i = 0; i < offset; ++i)
    _fx[i] = 0.0;
  for (int i = offset; i < _d; ++i)
    _fx[i] = subsol[i];
  _evaluator.eval_sub_sol(offset, _fx, dist);
}

template class ExternalEnumeration<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class ExternalEnumeration<Z_NR<mpz_t>, FP_NR<long double>>;
#endif

#ifdef FPLLL_WITH_QD
template class ExternalEnumeration<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class ExternalEnumeration<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class ExternalEnumeration<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class ExternalEnumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

template class ExternalEnumeration<Z_NR<long>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class ExternalEnumeration<Z_NR<long>, FP_NR<long double>>;
#endif

#ifdef FPLLL_WITH_QD
template class ExternalEnumeration<Z_NR<long>, FP_NR<dd_real>>;

template class ExternalEnumeration<Z_NR<long>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class ExternalEnumeration<Z_NR<long>, FP_NR<dpe_t>>;
#endif

template class ExternalEnumeration<Z_NR<long>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE
