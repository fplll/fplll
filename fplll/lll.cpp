/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.

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

/* Template source file */

#include "lll.h"
#include "util.h"

FPLLL_BEGIN_NAMESPACE

static inline bool is_power_of_2(int i) { return (i & (i - 1)) == 0; }

template <class ZT, class FT>
LLLReduction<ZT, FT>::LLLReduction(MatGSOInterface<ZT, FT> &m, double delta, double eta, int flags)
    : status(RED_SUCCESS), final_kappa(0), last_early_red(0), n_swaps(0), m(m)
{
  /* No early reduction in proved mode (i.e. enable_int_gram=true).
     NOTE: To make this possible, the hypothesis "g(i, j) is valid if
     0 <= i < n_known_rows and j <= i" in gso.h should be changed and
     MatGSOInterface<ZT, FT>::discover_row() should be rewritten. */
  enable_early_red = (flags & LLL_EARLY_RED) && !m.enable_int_gram;
  siegel           = flags & LLL_SIEGEL;
  verbose          = flags & LLL_VERBOSE;
  this->delta      = delta;
  this->eta        = eta;
  swap_threshold   = siegel ? delta - eta * eta : delta;
  zeros            = 0;
}

template <class ZT, class FT>
bool LLLReduction<ZT, FT>::lll(int kappa_min, int kappa_start, int kappa_end,
                               int size_reduction_start)
{
  if (kappa_end == -1)
    kappa_end = m.d;

  FPLLL_DEBUG_CHECK(kappa_min <= kappa_start && kappa_start < kappa_end && kappa_end <= m.d);
  int start_time = cputime();
  int kappa      = kappa_start + 1;
  int kappa_max  = 0;
  int d          = kappa_end - kappa_min;

  zeros       = 0;
  n_swaps     = 0;
  final_kappa = 0;
  if (verbose)
    print_params();
  extend_vect(lovasz_tests, kappa_end);
  extend_vect(babai_mu, kappa_end);
  extend_vect(babai_expo, kappa_end);

  for (; zeros < d && m.b_row_is_zero(0); zeros++)
  {
    m.move_row(kappa_min, kappa_end - 1 - zeros);
  }

  if (zeros < d && ((kappa_start > 0 && !babai(kappa_start, kappa_start, size_reduction_start)) ||
                    !m.update_gso_row(kappa_start)))
  {
    final_kappa = kappa_start;
    return false;
  }

  long long iter, max_iter;
  max_iter = static_cast<long long>(d - 2 * d * (d + 1) *
                                            ((m.get_max_exp_of_b() + 3) / std::log(delta.get_d())));

  for (iter = 0; iter < max_iter && kappa < kappa_end - zeros; iter++)
  {
    if (kappa > kappa_max)
    {
      if (verbose)
      {
        cerr << "Discovering vector " << kappa - kappa_min + 1 + zeros << "/" << d
             << " cputime=" << cputime() - start_time << endl;
      }
      kappa_max = kappa;
      if (enable_early_red && is_power_of_2(kappa) && kappa > last_early_red)
      {
        if (!early_reduction(kappa, size_reduction_start))
        {
          final_kappa = kappa;
          return false;
        }
      }
    }

    // Lazy size reduction
    if (!babai(kappa, kappa, size_reduction_start))
    {
      final_kappa = kappa;
      return false;
    }

    // Tests Lovasz's condition
    m.get_gram(lovasz_tests[0], kappa, kappa);
    for (int i = 1; i <= kappa; i++)
    {
      ftmp1.mul(m.get_mu_exp(kappa, i - 1), m.get_r_exp(kappa, i - 1));
      lovasz_tests[i].sub(lovasz_tests[i - 1], ftmp1);
    }
    ftmp1.mul(m.get_r_exp(kappa - 1, kappa - 1), swap_threshold);
    if (m.enable_row_expo)
    {
      ftmp1.mul_2si(ftmp1, 2 * (m.row_expo[kappa - 1] - m.row_expo[kappa]));
    }

    if (ftmp1 > lovasz_tests[siegel ? kappa : kappa - 1])
    {
      n_swaps++;
      // Failure, computes the insertion index
      int old_k = kappa;
      for (kappa--; kappa > kappa_min; kappa--)
      {
        ftmp1.mul(m.get_r_exp(kappa - 1, kappa - 1), swap_threshold);
        if (m.enable_row_expo)
        {
          ftmp1.mul_2si(ftmp1, 2 * (m.row_expo[kappa - 1] - m.row_expo[old_k]));
        }
        if (ftmp1 < lovasz_tests[siegel ? kappa : kappa - 1])
          break;
      }
      // FPLLL_TRACE("Lovasz's condition is not satisfied, kappa=" << kappa << " old_kappa=" <<
      // old_k);
      // Moves the vector
      if (lovasz_tests[kappa] > 0)
      {
        m.move_row(old_k, kappa);
      }
      else
      {
        zeros++;
        m.move_row(old_k, kappa_end - zeros);
        kappa = old_k;
        continue;
      }
    }

    m.set_r(kappa, kappa, lovasz_tests[kappa]);
    kappa++;
  }

  if (kappa < kappa_end - zeros)
    return set_status(RED_LLL_FAILURE);
  else
    return set_status(RED_SUCCESS);
}

template <class ZT, class FT>
bool LLLReduction<ZT, FT>::babai(int kappa, int size_reduction_end, int size_reduction_start)
{
  // FPLLL_TRACE_IN("kappa=" << kappa);
  long max_expo = LONG_MAX;

  for (int iter = 0;; iter++)
  {
    if (!m.update_gso_row(kappa, size_reduction_end - 1))
      return set_status(RED_GSO_FAILURE);

    bool loop_needed = false;
    for (int j = size_reduction_end - 1; j >= size_reduction_start && !loop_needed; j--)
    {
      m.get_mu(ftmp1, kappa, j);
      ftmp1.abs(ftmp1);
      loop_needed |= (ftmp1 > eta);
    }
    if (!loop_needed)
      break;

    if (iter >= 2)
    {
      long new_max_expo = m.get_max_mu_exp(kappa, size_reduction_end);
      if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
      {
        return set_status(RED_BABAI_FAILURE);
      }
      max_expo = new_max_expo;
    }

    for (int j = size_reduction_start; j < size_reduction_end; j++)
    {
      babai_mu[j] = m.get_mu_exp(kappa, j, babai_expo[j]);
    }
    m.row_op_begin(kappa, kappa + 1);
    for (int j = size_reduction_end - 1; j >= size_reduction_start; j--)
    {
      mu_m_ant.rnd_we(babai_mu[j], babai_expo[j]);
      if (mu_m_ant.zero_p())
        continue;
      // Approximate update of the mu_(kappa,k)'s
      for (int k = size_reduction_start; k < j; k++)
      {
        ftmp1.mul(mu_m_ant, m.get_mu_exp(j, k));
        /* When enable_row_expo=true, the following line relies on the fact that
           get_mu_exp(a, b, expo) returns expo = row_expo[a] - row_expo[b]. */
        babai_mu[k].sub(babai_mu[k], ftmp1);
      }
      // Operation on the basis
      // FPLLL_TRACE("Babai : row[" << kappa << "] += " << mu_m_ant << " * 2^" << babai_expo[j] << "
      // * row[" << j << "]");
      mu_m_ant.neg(mu_m_ant);
      m.row_addmul_we(kappa, j, mu_m_ant, babai_expo[j]);
    }
    m.row_op_end(kappa, kappa + 1);
  }
  return true;
}

template <class ZT, class FT>
bool is_lll_reduced(MatGSOInterface<ZT, FT> &m, double delta, double eta)
{
  FT ftmp1;
  FT ftmp2;
  FT delta_ = delta;
  m.update_gso();
  for (int i = 0; i < m.d; i++)
  {
    for (int j = 0; j < i; j++)
    {
      m.get_mu(ftmp1, i, j);
      ftmp1.abs(ftmp1);
      if (ftmp1 > eta)
        return false;
    }
  }
  for (int i = 1; i < m.d; i++)
  {
    m.get_mu(ftmp2, i, i - 1);
    ftmp2.mul(ftmp2, ftmp2);  // μ^2

    ftmp2.sub(delta_, ftmp2);  // δ - μ^2
    m.get_r(ftmp1, i - 1, i - 1);
    ftmp2.mul(ftmp1, ftmp2);  // (δ - μ^2) ⋅ r_{i-1,i-1}

    m.get_r(ftmp1, i, i);  // r_{i,i}

    if (ftmp1 < ftmp2)
      return false;
  }
  return true;
}

/** instantiate functions **/

template class LLLReduction<Z_NR<long>, FP_NR<double>>;
template class LLLReduction<Z_NR<double>, FP_NR<double>>;
template class LLLReduction<Z_NR<mpz_t>, FP_NR<double>>;
template bool
is_lll_reduced<Z_NR<mpz_t>, FP_NR<double>>(MatGSOInterface<Z_NR<mpz_t>, FP_NR<double>> &m,
                                           double delta, double eta);
template bool
is_lll_reduced<Z_NR<long>, FP_NR<double>>(MatGSOInterface<Z_NR<long>, FP_NR<double>> &m,
                                          double delta, double eta);
template bool
is_lll_reduced<Z_NR<double>, FP_NR<double>>(MatGSOInterface<Z_NR<double>, FP_NR<double>> &m,
                                            double delta, double eta);

#ifdef FPLLL_WITH_LONG_DOUBLE
template class LLLReduction<Z_NR<long>, FP_NR<long double>>;
template class LLLReduction<Z_NR<double>, FP_NR<long double>>;
template class LLLReduction<Z_NR<mpz_t>, FP_NR<long double>>;

template bool
is_lll_reduced<Z_NR<mpz_t>, FP_NR<long double>>(MatGSOInterface<Z_NR<mpz_t>, FP_NR<long double>> &m,
                                                double delta, double eta);
template bool
is_lll_reduced<Z_NR<long>, FP_NR<long double>>(MatGSOInterface<Z_NR<long>, FP_NR<long double>> &m,
                                               double delta, double eta);
template bool is_lll_reduced<Z_NR<double>, FP_NR<long double>>(
    MatGSOInterface<Z_NR<double>, FP_NR<long double>> &m, double delta, double eta);
#endif

#ifdef FPLLL_WITH_QD
template class LLLReduction<Z_NR<long>, FP_NR<dd_real>>;
template class LLLReduction<Z_NR<double>, FP_NR<dd_real>>;
template class LLLReduction<Z_NR<mpz_t>, FP_NR<dd_real>>;

template bool
is_lll_reduced<Z_NR<mpz_t>, FP_NR<dd_real>>(MatGSOInterface<Z_NR<mpz_t>, FP_NR<dd_real>> &m,
                                            double delta, double eta);
template bool
is_lll_reduced<Z_NR<long>, FP_NR<dd_real>>(MatGSOInterface<Z_NR<long>, FP_NR<dd_real>> &m,
                                           double delta, double eta);
template bool
is_lll_reduced<Z_NR<double>, FP_NR<dd_real>>(MatGSOInterface<Z_NR<double>, FP_NR<dd_real>> &m,
                                             double delta, double eta);

template class LLLReduction<Z_NR<long>, FP_NR<qd_real>>;
template class LLLReduction<Z_NR<double>, FP_NR<qd_real>>;
template class LLLReduction<Z_NR<mpz_t>, FP_NR<qd_real>>;

template bool
is_lll_reduced<Z_NR<mpz_t>, FP_NR<qd_real>>(MatGSOInterface<Z_NR<mpz_t>, FP_NR<qd_real>> &m,
                                            double delta, double eta);
template bool
is_lll_reduced<Z_NR<long>, FP_NR<qd_real>>(MatGSOInterface<Z_NR<long>, FP_NR<qd_real>> &m,
                                           double delta, double eta);
template bool
is_lll_reduced<Z_NR<double>, FP_NR<qd_real>>(MatGSOInterface<Z_NR<double>, FP_NR<qd_real>> &m,
                                             double delta, double eta);
#endif

#ifdef FPLLL_WITH_DPE
template class LLLReduction<Z_NR<long>, FP_NR<dpe_t>>;
template class LLLReduction<Z_NR<double>, FP_NR<dpe_t>>;
template class LLLReduction<Z_NR<mpz_t>, FP_NR<dpe_t>>;
template bool
is_lll_reduced<Z_NR<mpz_t>, FP_NR<dpe_t>>(MatGSOInterface<Z_NR<mpz_t>, FP_NR<dpe_t>> &m,
                                          double delta, double eta);
template bool is_lll_reduced<Z_NR<long>, FP_NR<dpe_t>>(MatGSOInterface<Z_NR<long>, FP_NR<dpe_t>> &m,
                                                       double delta, double eta);
template bool
is_lll_reduced<Z_NR<double>, FP_NR<dpe_t>>(MatGSOInterface<Z_NR<double>, FP_NR<dpe_t>> &m,
                                           double delta, double eta);
#endif

template class LLLReduction<Z_NR<long>, FP_NR<mpfr_t>>;
template class LLLReduction<Z_NR<double>, FP_NR<mpfr_t>>;
template class LLLReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>>;
template bool
is_lll_reduced<Z_NR<mpz_t>, FP_NR<mpfr_t>>(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &m,
                                           double delta, double eta);
template bool
is_lll_reduced<Z_NR<long>, FP_NR<mpfr_t>>(MatGSOInterface<Z_NR<long>, FP_NR<mpfr_t>> &m,
                                          double delta, double eta);
template bool
is_lll_reduced<Z_NR<double>, FP_NR<mpfr_t>>(MatGSOInterface<Z_NR<double>, FP_NR<mpfr_t>> &m,
                                            double delta, double eta);

FPLLL_END_NAMESPACE
