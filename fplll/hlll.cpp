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
  /* TODO: not exactly the good value
   * delta_ in (delta + 2^(-p + p0), 1 - 2^(-p + p0))
   */
  FT delta_      = delta;
  int start_time = cputime();
  long expo_k1_k1, expo_k_k1, expo_k_k;

  m.update_R(0);
  compute_dR(0, delta_);

  if (verbose)
    print_params();

  while (k < m.get_d())
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

    size_reduction(k);

    m.get_R(s, k, k - 1, expo_k_k1);
    s.mul(s, s);  // s = R(k, k - 1)^2
    m.get_R(tmp, k, k, expo_k_k);
    tmp.mul(tmp, tmp);  // tmp = R(k, k)^2
    s.add(tmp, s);      // s = R(k, k - 1)^2 + R(k, k)^2
    // Here, s = R(k, k - 1)^2 + R(k, k)^2 = ||b_k||^2 - sum_{i in [0, k-2)} R(k, i)^2
    expo_k1_k1 = m.get_row_expo(k - 1);

    if (expo_k1_k1 > -1)
      s.mul_2si(s, 2 * (expo_k_k - expo_k1_k1));

    if (dR[k - 1] <= s)
    {
      compute_dR(k, delta_);
      k++;
    }
    else
    {
      m.swap(k - 1, k);

      if (k - 1 == 0)
      {
        m.update_R(0);
        compute_dR(0, delta_);
      }
      else
        m.recover_R(k - 1);

      k = max(k - 1, 1);
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

  m.update_R(kappa, kappa);

  /* Most likely, at this step, the next update_R(kappa, kappa) must modify some coefficients since
   * b will most likely
   * be changed. If b is not modified during the size reduction, there will be only a call to
   * update_R_last(kappa),
   * which automatically must set updated_R to false.
   * TODO: find a best place to use this function.
   */
  m.set_updated_R_false();

  do
  {
    max_index = -1;

    for (int i = kappa - 1; i >= 0; i--)
    {
      m.get_R(ftmp1, kappa, i, expo0);  // expo0 = row_expo[kappa]
      m.get_R(ftmp0, i, i, expo1);      // expo1 = row_expo[i]

      ftmp1.mul(ftmp1, m.get_R_inverse_diag(i));  // x[i] = R(kappa, i) / R(i, i)
      /* If T = mpfr or dpe, enable_row_expo must be false and then, expo0 - expo1 == 0 (required by
       * rnd_we with this types) */
      ftmp1.rnd_we(ftmp1, expo0 - expo1);
      xf[i].neg(ftmp1);

      if (ftmp1.cmp(0.0) != 0)
      {
        for (int j = 0; j < i; j++)
        {
          m.get_R(ftmp1, i, j, expo0);      // expo0 = row_expo[i]
          ftmp1.mul(xf[i], ftmp1);          // ftmp1 = x[i] * R(i, j)
          m.get_R(ftmp0, kappa, j, expo0);  // expo0 = row_expo[kappa]
          ftmp0.add(ftmp0, ftmp1);          // ftmp0 = R(kappa, j) + x[i] * R(i, j)
          m.set_R(ftmp0, kappa, j);
        }
        max_index = max(max_index, i);
      }
    }

    if (max_index == -1)
    {
      // If max_index == -1, b(kappa) has not changed. Computing ||b[kappa]||^2 is not necessary.
      // 1 > 2^(-cd)=sr since cd > 0. Then, compute the last coefficient of R and stop.
      m.update_R_last(kappa);
      return;
    }
    else
    {
      m.norm_square_b_row(ftmp1, kappa, expo0);  // ftmp1 = ||b[kappa]||^2
      m.addmul_b_rows(kappa, xf);
      m.norm_square_b_row(ftmp0, kappa, expo1);  // ftmp0 = ||b[kappa]||^2
      ftmp1.mul(sr, ftmp1);                      // ftmp1 = 2^(-cd) * ftmp1 = sr * ftmp1

      if (expo1 > -1)
        ftmp0.mul_2si(ftmp0, expo1 - expo0);

      if (ftmp0.cmp(ftmp1) <= 0)
      {
        // Continue to try to reduce b(kappa).
        m.update_R(kappa, kappa);
      }
      else
      {
        // Compute the last coefficients of R and stop.
        m.update_R(kappa);
        return;
      }
    }
  } while (true);
}

template <class ZT, class FT>
bool is_hlll_reduced(MatHouseholder<ZT, FT> &m, double delta, double eta)
{
  FT ftmp0;
  FT ftmp1;
  FT delta_ = delta;
  FT eta_   = eta;

  m.compute_R_naively();

  long expo0 = -1;
  long expo1 = -1;

  for (int i = 0; i < m.get_d(); i++)
  {
    for (int j = 0; j < i; j++)
    {
      m.get_R(ftmp0, i, j, expo0);
      m.get_R(ftmp1, j, j, expo1);
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

  for (int i = 1; i < m.get_d(); i++)
  {
    m.norm_square_b_row(ftmp0, i, expo0);  // ftmp0 = ||b[i]||^2
    m.norm_square_R_row(ftmp1, i, i - 1, expo1);

    if (expo0 > -1)
      ftmp0.mul_2si(ftmp0, expo0 - expo1);

    ftmp1.sub(ftmp0, ftmp1);  // ftmp1 = ||b[i]||^2 - sum_{i = 0}^{i < i - 1}R[i][i]^2
    m.get_R(ftmp0, i - 1, i - 1, expo0);
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
                                            double delta, double eta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<double>>(MatHouseholder<Z_NR<long>, FP_NR<double>> &m,
                                           double delta, double eta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<double>>(MatHouseholder<Z_NR<double>, FP_NR<double>> &m,
                                             double delta, double eta);

#ifdef FPLLL_WITH_LONG_DOUBLE
template class HLLLReduction<Z_NR<long>, FP_NR<long double>>;
template class HLLLReduction<Z_NR<double>, FP_NR<long double>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<long double>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<long double>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<long double>> &m,
                                                 double delta, double eta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<long double>>(MatHouseholder<Z_NR<long>, FP_NR<long double>> &m,
                                                double delta, double eta);
template bool is_hlll_reduced<Z_NR<double>, FP_NR<long double>>(
    MatHouseholder<Z_NR<double>, FP_NR<long double>> &m, double delta, double eta);
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
                                             double delta, double eta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<qd_real>>(MatHouseholder<Z_NR<long>, FP_NR<qd_real>> &m,
                                            double delta, double eta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<qd_real>>(MatHouseholder<Z_NR<double>, FP_NR<qd_real>> &m,
                                              double delta, double eta);
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<dd_real>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<dd_real>> &m,
                                             double delta, double eta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<dd_real>>(MatHouseholder<Z_NR<long>, FP_NR<dd_real>> &m,
                                            double delta, double eta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<dd_real>>(MatHouseholder<Z_NR<double>, FP_NR<dd_real>> &m,
                                              double delta, double eta);
#endif

#ifdef FPLLL_WITH_DPE
template class HLLLReduction<Z_NR<long>, FP_NR<dpe_t>>;
template class HLLLReduction<Z_NR<double>, FP_NR<dpe_t>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<dpe_t>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<dpe_t>> &m,
                                           double delta, double eta);
template bool is_hlll_reduced<Z_NR<long>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<long>, FP_NR<dpe_t>> &m,
                                                        double delta, double eta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<double>, FP_NR<dpe_t>> &m,
                                            double delta, double eta);
#endif

template class HLLLReduction<Z_NR<long>, FP_NR<mpfr_t>>;
template class HLLLReduction<Z_NR<double>, FP_NR<mpfr_t>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<mpfr_t>> &m,
                                            double delta, double eta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<long>, FP_NR<mpfr_t>> &m,
                                           double delta, double eta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<double>, FP_NR<mpfr_t>> &m,
                                             double delta, double eta);

FPLLL_END_NAMESPACE
