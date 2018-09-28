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

#include "hlll.h"
#include "util.h"

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT> void HLLLReduction<ZT, FT>::lll()
{
  int k     = 1;
  int k_max = 0;
  FT s;
  FT tmp;
  FT sum;

#ifdef HOUSEHOLDER_NAIVELY
  FT dR;
  long expo_dR;
#endif  // HOUSEHOLDER_NAIVELY

  /* TODO: not exactly the good value
   * delta_ in (delta + 2^(-p + p0), 1 - 2^(-p + p0))
   */
  FT delta_      = delta;
  int start_time = cputime();
  long expo_k1_k1, expo_k_k1, expo_k_k;

#ifndef HOUSEHOLDER_NAIVELY
  m.refresh_R_bf(0);
  m.update_R(0);
  compute_dR(0, delta_);
#else   // HOUSEHOLDER_NAIVELY
  m.update_R_naively(0);
#endif  // HOUSEHOLDER_NAIVELY

  if (verbose)
  {
    print_params();
    // Discover vector 1
    cerr << "Discovering vector " << k + 1 << "/" << m.get_d()
         << " cputime=" << cputime() - start_time << endl;
  }

  m.refresh_R_bf(1);

  while (true)
  {
    size_reduction(k);

// s = R(k, k-1)
#ifndef HOUSEHOLDER_NAIVELY
    m.get_R(s, k, k - 1, expo_k_k1);
#else   //  HOUSEHOLDER_NAIVELY
    m.get_R_naively(s, k, k - 1, expo_k_k1);
#endif  // HOUSEHOLDER_NAIVELY

    s.mul(s, s);  // s = R(k, k - 1)^2

// tmp = R(k, k)
#ifndef HOUSEHOLDER_NAIVELY
    m.get_R(tmp, k, k, expo_k_k);
#else   //  HOUSEHOLDER_NAIVELY
    m.get_R_naively(tmp, k, k, expo_k_k);
#endif  // HOUSEHOLDER_NAIVELY

    tmp.mul(tmp, tmp);  // tmp = R(k, k)^2
    s.add(tmp, s);      // s = R(k, k - 1)^2 + R(k, k)^2
// Here, s = R(k, k - 1)^2 + R(k, k)^2 = ||b_k||^2 - sum_{i in [0, k-2)} R(k, i)^2

// Get expo of row k - 1
#ifndef HOUSEHOLDER_NAIVELY
    expo_k1_k1 = m.get_row_expo(k - 1);
#else   //  HOUSEHOLDER_NAIVELY
    expo_k1_k1 = m.get_row_expo_naively(k - 1);
#endif  // HOUSEHOLDER_NAIVELY

    if (expo_k1_k1 > -1)
      s.mul_2si(s, 2 * (expo_k_k - expo_k1_k1));

// Test if delta_ * R(k, k)^2 <= s
#ifndef HOUSEHOLDER_NAIVELY
    if (dR[k - 1].cmp(s) <= 0)
#else   //  HOUSEHOLDER_NAIVELY
    m.get_R_naively(dR, k - 1, k - 1, expo_dR);
    dR.mul(dR, dR);
    dR.mul(delta_, dR);  // dR[k] = delta_ * R(k, k)^2

    if (dR.cmp(s) <= 0)
#endif  // HOUSEHOLDER_NAIVELY

    {
#ifndef HOUSEHOLDER_NAIVELY
      compute_dR(k, delta_);
#endif  // HOUSEHOLDER_NAIVELY
      k++;

      if (k < m.get_d())
      {
        if (k > k_max)
        {
          if (verbose)
          {
            cerr << "Discovering vector " << k + 1 << "/" << m.get_d()
                 << " cputime=" << cputime() - start_time << endl;
          }
          k_max = k;
        }
        // TODO: do not need to refresh bf if B[k] was already seen?
        m.refresh_R_bf(k);
      }
      else
        return;
    }
    else
    {
      m.swap(k - 1, k);

      if (k - 1 == 0)
      {
// Update row 0
#ifndef HOUSEHOLDER_NAIVELY
        m.refresh_R_bf(0);
        m.update_R(0);
        compute_dR(0, delta_);
#else   // HOUSEHOLDER_NAIVELY
        m.update_R_naively(0);
#endif  // HOUSEHOLDER_NAIVELY

        // TODO: do not need to refresh bf?
        m.refresh_R_bf(1);
        k = 1;
      }
      else
      {
        k--;
        m.recover_R(k);
      }
    }
  }
}

template <class ZT, class FT> void HLLLReduction<ZT, FT>::size_reduction(int kappa)
{
  vector<FT> xf(kappa);
  FT ftmp0, ftmp1;
  long expo0 = -1;
  long expo1 = -1;
  // for all i > max_index, xf[i] == 0.
  int max_index = -1;

#ifndef HOUSEHOLDER_NAIVELY
  m.update_R(kappa, false);

  /* Most likely, at this step, the next update_R(kappa, false) must modify some coefficients since
   * b will most likely
   * be changed. If b is not modified during the size reduction, there will be only a call to
   * update_R_last(kappa),
   * which automatically must set updated_R to false.
   * TODO: find a best place to use this function.
   */
  m.set_updated_R_false();
#else   // HOUSEHOLDER_NAIVELY
  m.update_R_naively(kappa, false);
#endif  // HOUSEHOLDER_NAIVELY

  do
  {
    max_index = -1;

    for (int i = kappa - 1; i >= 0; i--)
    {
#ifndef HOUSEHOLDER_NAIVELY
      m.get_R(ftmp1, kappa, i, expo0);  // expo0 = row_expo[kappa]
      m.get_R(ftmp0, i, i, expo1);      // expo1 = row_expo[i]

#ifndef HOUSEHOLDER_PRECOMPUTE_INVERSE
      ftmp1.div(ftmp1, ftmp0);  // x[i] = R(kappa, i) / R(i, i)
#else                           // HOUSEHOLDER_PRECOMPUTE_INVERSE
      ftmp1.mul(ftmp1, m.get_R_inverse_diag(i));  // x[i] = R(kappa, i) / R(i, i)
#endif                          // HOUSEHOLDER_PRECOMPUTE_INVERSE
#else                           // HOUSEHOLDER_NAIVELY
      m.get_R_naively(ftmp1, kappa, i, expo0);  // expo0 = row_expo[kappa]
      m.get_R_naively(ftmp0, i, i, expo1);      // expo1 = row_expo[i]

      ftmp1.div(ftmp1, ftmp0);                      // x[i] = R(kappa, i) / R(i, i)
#endif                          // HOUSEHOLDER_NAIVELY

      /* If T = mpfr or dpe, enable_row_expo must be false and then, expo0 - expo1 == 0 (required by
       * rnd_we with this types) */
      ftmp1.rnd_we(ftmp1, expo0 - expo1);
      xf[i].neg(ftmp1);

      if (ftmp1.cmp(0.0) != 0)
      {
        for (int j = 0; j < i; j++)
        {
#ifndef HOUSEHOLDER_NAIVELY
          m.get_R(ftmp1, i, j, expo0);  // expo0 = row_expo[i]
#else  // HOUSEHOLDER_NAIVELY
          m.get_R_naively(ftmp1, i, j, expo0);      // expo0 = row_expo[i]
#endif  // HOUSEHOLDER_NAIVELY

          ftmp1.mul(xf[i], ftmp1);  // ftmp1 = x[i] * R(i, j)

#ifndef HOUSEHOLDER_NAIVELY
          m.get_R(ftmp0, kappa, j, expo0);  // expo0 = row_expo[kappa]
#else                                       // HOUSEHOLDER_NAIVELY
          m.get_R_naively(ftmp0, kappa, j, expo0);  // expo0 = row_expo[kappa]
#endif                                      // HOUSEHOLDER_NAIVELY

          ftmp0.add(ftmp0, ftmp1);  // ftmp0 = R(kappa, j) + x[i] * R(i, j)

#ifndef HOUSEHOLDER_NAIVELY
          m.set_R(ftmp0, kappa, j);
#else   // HOUSEHOLDER_NAIVELY
          m.set_R_naively(ftmp0, kappa, j);
#endif  // HOUSEHOLDER_NAIVELY
        }
        max_index = max(max_index, i);
      }
    }

#ifndef HOUSEHOLDER_PRECOMPUTE_INVERSE
    if (max_index == -1)
    {
// If max_index == -1, b(kappa) has not changed. Computing ||b[kappa]||^2 is not necessary.
// 1 > 2^(-cd)=sr since cd > 0. Then, compute the last coefficient of R and stop.
#ifndef HOUSEHOLDER_NAIVELY
      m.update_R_last(kappa);
#else   // HOUSEHOLDER_NAIVELY
      m.update_R_last_naively(kappa);
#endif  // HOUSEHOLDER_NAIVELY

      return;
    }
    else
#endif  // HOUSEHOLDER_PRECOMPUTE_INVERSE
    {
#ifndef HOUSEHOLDER_NAIVELY
      m.norm_square_b_row(ftmp1, kappa, expo0);  // ftmp1 = ||b[kappa]||^2
      m.addmul_b_rows(kappa, xf);
      m.norm_square_b_row(ftmp0, kappa, expo1);  // ftmp0 = ||b[kappa]||^2
      m.refresh_R_bf(kappa);
#else   // HOUSEHOLDER_NAIVELY
      m.norm_square_b_row_naively(ftmp1, kappa, expo0);  // ftmp1 = ||b[kappa]||^2
      m.addmul_b_rows_naively(kappa, xf);
      m.norm_square_b_row_naively(ftmp0, kappa, expo1);  // ftmp0 = ||b[kappa]||^2
#endif  // HOUSEHOLDER_NAIVELY

      ftmp1.mul(sr, ftmp1);  // ftmp1 = 2^(-cd) * ftmp1 = sr * ftmp1

      if (expo1 > -1)
        ftmp0.mul_2si(ftmp0, expo1 - expo0);

      if (ftmp0.cmp(ftmp1) <= 0)
      {
// Continue to try to reduce b(kappa).
#ifndef HOUSEHOLDER_NAIVELY
        m.update_R(kappa, false);
#else   // HOUSEHOLDER_NAIVELY
        m.update_R_naively(kappa, false);
#endif  // HOUSEHOLDER_NAIVELY
      }
      else
      {
// Compute the last coefficients of R and stop.
#ifndef HOUSEHOLDER_NAIVELY
        m.update_R(kappa);
#else   // HOUSEHOLDER_NAIVELY
        m.update_R_naively(kappa);
#endif  // HOUSEHOLDER_NAIVELY

        return;
      }
    }
  } while (true);
}

/* d is to perform the check only for row [0, d[. If d = -1, d = m.get_d(). */
/* If compute is true and d > -1, apply Householder transpose on row [0, d[. */
template <class ZT, class FT>
bool is_hlll_reduced(MatHouseholder<ZT, FT> &m, double delta, double eta, int d, bool compute)
{
  FT ftmp0;
  FT ftmp1;
  FT delta_ = delta;
  FT eta_   = eta;

  if (d == -1)
  {
    d = m.get_d();
    m.update_R_naively();
  }
  else if (compute)
    for (int i = 0; i < d; i++)
      m.update_R_naively(i);

  long expo0 = -1;
  long expo1 = -1;

  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < i; j++)
    {
      m.get_R_naively(ftmp0, i, j, expo0);
      m.get_R_naively(ftmp1, j, j, expo1);
      ftmp1.div(ftmp0, ftmp1);
      ftmp1.abs(ftmp1);

      if (expo0 > -1)
        ftmp1.mul_2si(ftmp1, expo0 - expo1);

      if (ftmp1.cmp(eta_) > 0)
        return false;

      // Since eta_ is not involved in the test of the algorithm, this test is probably the one we
      // want
      if (ftmp1.cmp(0.5) > 0)
        return false;
    }
  }

  for (int i = 1; i < d; i++)
  {
    m.norm_square_b_row_naively(ftmp0, i, expo0);  // ftmp0 = ||b[i]||^2
    m.norm_square_R_row_naively(ftmp1, i, i - 1,
                                expo1);  // ftmp1 = sum_{i = 0}^{i < i - 1}R[i][i]^2

    if (expo0 > -1)
      ftmp0.mul_2si(ftmp0, expo0 - expo1);

    ftmp1.sub(ftmp0, ftmp1);  // ftmp1 = ||b[i]||^2 - sum_{i = 0}^{i < i - 1}R[i][i]^2
    m.get_R_naively(ftmp0, i - 1, i - 1, expo0);
    ftmp0.mul(ftmp0, ftmp0);
    ftmp0.mul(delta_, ftmp0);  // ftmp0 = delta_ * R(i - 1, i - 1)^2

    if (expo0 > -1)
      ftmp1.mul_2si(ftmp1, expo1 - expo0);

    if (ftmp0.cmp(ftmp1) > 0)
      return false;
  }
  return true;
}

/** instantiate functions **/

template class HLLLReduction<Z_NR<long>, FP_NR<double>>;
template class HLLLReduction<Z_NR<double>, FP_NR<double>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<double>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<double>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<double>> &m,
                                            double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<double>>(MatHouseholder<Z_NR<long>, FP_NR<double>> &m,
                                           double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<double>>(MatHouseholder<Z_NR<double>, FP_NR<double>> &m,
                                             double delta, double eta, int d, bool compute);

#ifdef FPLLL_WITH_LONG_DOUBLE
template class HLLLReduction<Z_NR<long>, FP_NR<long double>>;
template class HLLLReduction<Z_NR<double>, FP_NR<long double>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<long double>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<long double>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<long double>> &m,
                                                 double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<long double>>(MatHouseholder<Z_NR<long>, FP_NR<long double>> &m,
                                                double delta, double eta, int d, bool compute);
template bool is_hlll_reduced<Z_NR<double>, FP_NR<long double>>(
    MatHouseholder<Z_NR<double>, FP_NR<long double>> &m, double delta, double eta, int d,
    bool compute);
#endif

#ifdef FPLLL_WITH_QD
template class HLLLReduction<Z_NR<long>, FP_NR<dd_real>>;
template class HLLLReduction<Z_NR<double>, FP_NR<dd_real>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class HLLLReduction<Z_NR<long>, FP_NR<qd_real>>;
template class HLLLReduction<Z_NR<double>, FP_NR<qd_real>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<qd_real>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<qd_real>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<qd_real>> &m,
                                             double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<qd_real>>(MatHouseholder<Z_NR<long>, FP_NR<qd_real>> &m,
                                            double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<qd_real>>(MatHouseholder<Z_NR<double>, FP_NR<qd_real>> &m,
                                              double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<dd_real>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<dd_real>> &m,
                                             double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<dd_real>>(MatHouseholder<Z_NR<long>, FP_NR<dd_real>> &m,
                                            double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<dd_real>>(MatHouseholder<Z_NR<double>, FP_NR<dd_real>> &m,
                                              double delta, double eta, int d, bool compute);
#endif

#ifdef FPLLL_WITH_DPE
template class HLLLReduction<Z_NR<long>, FP_NR<dpe_t>>;
template class HLLLReduction<Z_NR<double>, FP_NR<dpe_t>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<dpe_t>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<dpe_t>> &m,
                                           double delta, double eta, int d, bool compute);
template bool is_hlll_reduced<Z_NR<long>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<long>, FP_NR<dpe_t>> &m,
                                                        double delta, double eta, int d,
                                                        bool compute);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<double>, FP_NR<dpe_t>> &m,
                                            double delta, double eta, int d, bool compute);
#endif

template class HLLLReduction<Z_NR<long>, FP_NR<mpfr_t>>;
template class HLLLReduction<Z_NR<double>, FP_NR<mpfr_t>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<mpfr_t>> &m,
                                            double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<long>, FP_NR<mpfr_t>> &m,
                                           double delta, double eta, int d, bool compute);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<double>, FP_NR<mpfr_t>> &m,
                                             double delta, double eta, int d, bool compute);

FPLLL_END_NAMESPACE
