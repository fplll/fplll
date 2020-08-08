/* Copyright (C) 2008-2011 Xavier Pujol
   (C) 2015 Michael Walter.
   (C) 2016 Marc Stevens. (generic improvements, auxiliary solutions, subsolutions)
   (C) 2016 Guillaume Bonnoron. (CVP improvements)

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

#include "enumerate.h"

FPLLL_BEGIN_NAMESPACE

template <typename ZT, typename FT, typename CounterClass>
void EnumerationDyn<ZT, FT, CounterClass>::reset(enumf cur_dist, int cur_depth)
{
  // FPLLL_TRACE("Reset level " << cur_depth);
  int new_dim = cur_depth + 1;
  auto &d     = this->d;

  vector<enumxt> partial_sol(d - cur_depth - 1);
  for (int i = cur_depth + 1; i < d; ++i)
    partial_sol[i - cur_depth - 1] = this->x[i];

  FT new_dist = 0.0;
  for (int i = 0; i < new_dim; i++)
    new_dist.add(new_dist, _gso.get_r_exp(i, i));

  FastEvaluator<FT> new_evaluator;
  Enumeration<ZT, FT, CounterClass> enumobj(_gso, new_evaluator, this->_max_indices);
  enumobj.enumerate(0, d, new_dist, 0, target, partial_sol, pruning_bounds, false, true);

  if (!new_evaluator.empty())
  {
    FT sol_dist2 = new_evaluator.begin()->first;
    sol_dist2.mul_2si(sol_dist2, -new_evaluator.normExp);
    enumf sol_dist = sol_dist2.get_d();
    // FPLLL_TRACE("Recovering sub-solution at level: " << cur_depth <<" soldist: " << sol_dist);

    if (sol_dist + cur_dist < this->partdistbounds[0])
    {
      // FPLLL_TRACE("Saving it.");
      for (int i = 0; i < new_dim; ++i)
        this->x[i] = new_evaluator.begin()->second[i].get_d();
      process_solution(sol_dist + cur_dist);
    }
  }
}

template <typename ZT, typename FT, typename CounterClass>
void EnumerationDyn<ZT, FT, CounterClass>::enumerate(
    int first, int last, FT &fmaxdist, long fmaxdistexpo, const vector<FT> &target_coord,
    const vector<enumxt> &subtree, const vector<enumf> &pruning, bool _dual, bool subtree_reset)
{
  auto &d         = this->d;
  bool solvingsvp = target_coord.empty();
  this->dual      = _dual;
  pruning_bounds  = pruning;
  target          = target_coord;
  if (last == -1)
    last = _gso.d;
  d = last - first;
  fx.resize(d);
  FPLLL_CHECK(d < this->maxdim, "enumerate: dimension is too high");
  FPLLL_CHECK((solvingsvp || !this->dual),
              "CVP for dual not implemented! What does that even mean? ");
  FPLLL_CHECK((subtree.empty() || !this->dual), "Subtree enumeration for dual not implemented!");

  this->resetflag = !this->_max_indices.empty();
  if (this->resetflag)
    this->reset_depth = this->_max_indices[last - subtree.size() - 1];

  if (solvingsvp)
  {
    for (int i = 0; i < d; ++i)
      this->center_partsum[i] = 0.0;
  }
  else
  {
    for (int i = 0; i < d; ++i)
      this->center_partsum[i] = target_coord[i + first].get_d();
  }

  FT fr, fmu, fmaxdistnorm;
  long rexpo, normexp = -1;
  for (int i = 0; i < d; ++i)
  {
    fr      = _gso.get_r_exp(i + first, i + first, rexpo);
    normexp = max(normexp, rexpo + fr.exponent());
  }
  // normalization is multiplication by 2^(-normexp). in case of dual, we are normalizing
  // the inverse of r, so we negate the normalizing exponent
  if (this->dual)
  {
    normexp *= -1;
  }
  fmaxdistnorm.mul_2si(fmaxdist, fmaxdistexpo - normexp);
  maxdist = fmaxdistnorm.get_d(GMP_RNDU);
  _evaluator.set_normexp(normexp);

  if (this->dual)
  {
    for (int i = 0; i < d; ++i)
    {
      fr = _gso.get_r_exp(i + first, i + first, rexpo);
      fr.mul_2si(fr, rexpo + normexp);
      this->rdiag[d - i - 1] = enumf(1.0) / fr.get_d();
    }
    for (int i = 0; i < d; ++i)
    {
      for (int j = i + 1; j < d; ++j)
      {
        _gso.get_mu(fmu, j + first, i + first);
        this->mut[d - j - 1][d - i - 1] = -fmu.get_d();
      }
    }
  }
  else
  {
    for (int i = 0; i < d; ++i)
    {
      fr = _gso.get_r_exp(i + first, i + first, rexpo);
      fr.mul_2si(fr, rexpo - normexp);
      this->rdiag[i] = fr.get_d();
    }
    for (int i = 0; i < d; ++i)
    {
      for (int j = i + 1; j < d; ++j)
      {
        _gso.get_mu(fmu, j + first, i + first);
        this->mut[i][j] = fmu.get_d();
      }
    }
  }

  this->subsoldists = this->rdiag;

  this->save_rounding();
  prepare_enumeration(subtree, solvingsvp, subtree_reset);
  do_enumerate();
  this->restore_rounding();

  fmaxdistnorm = maxdist;  // Exact

  fmaxdist.mul_2si(fmaxdistnorm, normexp - fmaxdistexpo);

  if (this->dual && !_evaluator.empty())
  {
    for (auto it = _evaluator.begin(), itend = _evaluator.end(); it != itend; ++it)
      reverse_by_swap(it->second, 0, d - 1);
  }
}

template <typename ZT, typename FT, typename CounterClass>
void EnumerationDyn<ZT, FT, CounterClass>::prepare_enumeration(const vector<enumxt> &subtree,
                                                               bool solvingsvp, bool subtree_reset)
{
  this->is_svp = solvingsvp;
  auto &k      = this->k;
  auto &d      = this->d;
  auto &k_end  = this->k_end;

  enumf newdist = 0.0;
  this->k_end   = d - subtree.size();
  for (k = d - 1; k >= 0 && newdist <= maxdist; --k)
  {
    enumf newcenter = this->center_partsum[k];
    if (k >= k_end)
    {
      // use subtree
      this->x[k] = subtree[k - k_end];

      if (this->x[k] != 0)
        this->is_svp = false;

      for (int j = 0; j < k; ++j)
        this->center_partsum[j] -= this->x[k] * this->mut[j][k];
    }
    else
    {
      if (this->dual)
      {
        for (int j = k + 1; j < k_end; ++j)
          newcenter -= this->alpha[j] * this->mut[k][j];
      }
      else
      {
        for (int j = k + 1; j < k_end; ++j)
          newcenter -= this->x[j] * this->mut[k][j];
      }
      roundto(this->x[k], newcenter);  // newX = rint(newCenter) / lrint(newCenter)
      this->center[k]   = newcenter;
      this->partdist[k] = newdist;
      this->dx[k] = this->ddx[k] = (((int)(newcenter >= this->x[k]) & 1) << 1) - 1;
    }
    if (!subtree_reset || k < k_end)
    {
      this->alpha[k] = this->x[k] - newcenter;
      newdist += this->alpha[k] * this->alpha[k] * this->rdiag[k];
    }
  }
  if (!this->is_svp)
  {
    this->k_max = k_end;  // The last non-zero coordinate of x will not stay positive
  }
  else
  {
    this->k_max = 0;
    this->x[0]  = 1;  // Excludes (0,...,0) from the enumeration
  }
  ++k;
}

template <typename ZT, typename FT, typename CounterClass>
void EnumerationDyn<ZT, FT, CounterClass>::set_bounds()
{
  auto &d = this->d;
  if (pruning_bounds.empty())
  {
    fill(&this->partdistbounds[0] + 0, &this->partdistbounds[0] + d, maxdist);
  }
  else
  {
    for (int i = 0; i < d; ++i)
      this->partdistbounds[i] = this->pruning_bounds[i] * maxdist;
  }
}

template <typename ZT, typename FT, typename CounterClass>
void EnumerationDyn<ZT, FT, CounterClass>::process_solution(enumf newmaxdist)
{
  FPLLL_TRACE("Sol dist: " << newmaxdist << " (nodes:" << nodes << ")");
  auto &d = this->d;
  for (int j = 0; j < d; ++j)
    fx[j] = this->x[j];
  _evaluator.eval_sol(fx, newmaxdist, maxdist);

  set_bounds();
}

template <typename ZT, typename FT, typename CounterClass>
void EnumerationDyn<ZT, FT, CounterClass>::process_subsolution(int offset, enumf newdist)
{
  auto &d = this->d;
  for (int j = 0; j < offset; ++j)
    fx[j] = 0.0;
  for (int j = offset; j < d; ++j)
    fx[j] = this->x[j];
  _evaluator.eval_sub_sol(offset, fx, newdist);
}

template <typename ZT, typename FT, typename CounterClass>
void EnumerationDyn<ZT, FT, CounterClass>::do_enumerate()
{

  this->nodes_counter.reset();
  set_bounds();
  auto dual      = this->dual;
  auto resetflag = this->resetflag;

  if (dual && _evaluator.findsubsols && !resetflag)
    (*this).template enumerate_loop<true, true, false>();
  else if (!dual && _evaluator.findsubsols && !resetflag)
    (*this).template enumerate_loop<false, true, false>();
  else if (dual && !_evaluator.findsubsols && !resetflag)
    (*this).template enumerate_loop<true, false, false>();
  else if (!dual && !_evaluator.findsubsols && !resetflag)
    (*this).template enumerate_loop<false, false, false>();
  else if (!dual && _evaluator.findsubsols && resetflag)
    (*this).template enumerate_loop<false, true, true>();
  else if (!dual && !_evaluator.findsubsols && resetflag)
    (*this).template enumerate_loop<false, false, true>();
}

template class Enumeration<Z_NR<mpz_t>, FP_NR<double>, WholeTreeCounter>;
template class Enumeration<Z_NR<mpz_t>, FP_NR<double>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<double>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<double>, LevelTreeCounter>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class Enumeration<Z_NR<mpz_t>, FP_NR<long double>, WholeTreeCounter>;
template class Enumeration<Z_NR<mpz_t>, FP_NR<long double>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<long double>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<long double>, LevelTreeCounter>;
#endif

#ifdef FPLLL_WITH_QD
template class Enumeration<Z_NR<mpz_t>, FP_NR<dd_real>, WholeTreeCounter>;
template class Enumeration<Z_NR<mpz_t>, FP_NR<dd_real>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<dd_real>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<dd_real>, LevelTreeCounter>;

template class Enumeration<Z_NR<mpz_t>, FP_NR<qd_real>, WholeTreeCounter>;
template class Enumeration<Z_NR<mpz_t>, FP_NR<qd_real>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<qd_real>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<qd_real>, LevelTreeCounter>;
#endif

#ifdef FPLLL_WITH_DPE
template class Enumeration<Z_NR<mpz_t>, FP_NR<dpe_t>, WholeTreeCounter>;
template class Enumeration<Z_NR<mpz_t>, FP_NR<dpe_t>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<dpe_t>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<dpe_t>, LevelTreeCounter>;
#endif

template class Enumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>, WholeTreeCounter>;
template class Enumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<mpfr_t>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<mpz_t>, FP_NR<mpfr_t>, LevelTreeCounter>;

template class Enumeration<Z_NR<long>, FP_NR<double>, WholeTreeCounter>;
template class Enumeration<Z_NR<long>, FP_NR<double>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<double>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<double>, LevelTreeCounter>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class Enumeration<Z_NR<long>, FP_NR<long double>, WholeTreeCounter>;
template class Enumeration<Z_NR<long>, FP_NR<long double>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<long double>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<long double>, LevelTreeCounter>;
#endif

#ifdef FPLLL_WITH_QD
template class Enumeration<Z_NR<long>, FP_NR<dd_real>, WholeTreeCounter>;
template class Enumeration<Z_NR<long>, FP_NR<dd_real>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<dd_real>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<dd_real>, LevelTreeCounter>;

template class Enumeration<Z_NR<long>, FP_NR<qd_real>, WholeTreeCounter>;
template class Enumeration<Z_NR<long>, FP_NR<qd_real>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<qd_real>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<qd_real>, LevelTreeCounter>;
#endif

#ifdef FPLLL_WITH_DPE
template class Enumeration<Z_NR<long>, FP_NR<dpe_t>, WholeTreeCounter>;
template class Enumeration<Z_NR<long>, FP_NR<dpe_t>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<dpe_t>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<dpe_t>, LevelTreeCounter>;
#endif

template class Enumeration<Z_NR<long>, FP_NR<mpfr_t>, WholeTreeCounter>;
template class Enumeration<Z_NR<long>, FP_NR<mpfr_t>, LevelTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<mpfr_t>, WholeTreeCounter>;
template class EnumerationDyn<Z_NR<long>, FP_NR<mpfr_t>, LevelTreeCounter>;

FPLLL_END_NAMESPACE
