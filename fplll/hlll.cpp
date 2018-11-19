/*
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

template <class ZT, class FT> bool HLLLReduction<ZT, FT>::hlll()
{
  /* TODO: not exactly the good value
   * delta_ in (delta + 2^(-p + p0), 1 - 2^(-p + p0))
   */
  FT delta_      = delta;
  int start_time = cputime();
  long expo0 = 0, expo1 = 0;

  if (verbose)
  {
    // Print the parameters of the computation
    print_params();
    // Discover b[0]
    cerr << "Discovering vector 1/" << m.get_d() << " cputime=" << cputime() - start_time << endl;
  }

  // Set R[0] and bf[0] to b[0], precompute ||b[0]||^2
  m.refresh_R_bf(0);
  // Compute R[0]
  m.update_R_last(0);  // In this case, update_R(0) is exactly equal to update_R_last(0)
  // Precompute dR[0]: R[0]^2 * delta_ = dR[0] * 2^(2*row_expo[0])
  compute_dR(0, delta_);

  int k = 1;
  // Remember which was the largest b[k_max] that is tried to be size-reduced
  int k_max = 1;

  int prev_k = -1;
  vector<FT> prev_R;
  vector<long> prev_expo;
  prev_R.resize(m.get_d());
  prev_expo.resize(m.get_d());
  // prev_R[0] is never used: m.get_R(prev_R[0], 0, 0, prev_expo[0]);

  if (verbose)
  {
    // Discover b[1]
    cerr << "Discovering vector 2/" << m.get_d() << " cputime=" << cputime() - start_time << endl;
  }

  // Set R[1] and bf[1] to b[1], precompute ||b[1]||^2
  m.refresh_R_bf(1);

  while (true)
  {
    // Size reduce b[k] thanks to b[0] to b[k - 1]
    size_reduction(k, k, 0);

#ifndef MODIFIED_LOVASZ_TEST
    // This Lovasz test is the one proposed in [MSV, ISSAC'09]
    //
    // Prior, we used another test, which is dR[k-1].cmp(R(k, k - 1)^2 + R(k, k)^2) <= 0, if R(k, k)
    // is known (which is not the case here at this step of the computation, but can be computed
    // thanks to sqrt(sum_{i=k}^{i<n}R(k, i)^2) (indices must be checked)). In the prior version,
    // since R(k, k)^2 was known, we directly used the test. However, the test was probably not as
    // accurate as what we hope. However, since this formula is used in hplll
    // (https://github.com/gilvillard/hplll/releases), this can maybe be retested. An example of
    // matrices that was not HLLL reduced because of a fail in this test can be generated thanks to
    // latticegen -randseed 122 r 300 30000.
    // Such a test can be activate by compiling with -DMODIFIED_LOVASZ_TEST
    //
    // TODO: this section must be improved and investigated, i.e.:
    //   * probably other improvement to be done
    m.get_norm_square_b(ftmp0, k, expo0);  // ||b[k]||^2 = ftmp0 * 2^expo0
    m.norm_square_R_row(ftmp1, k, 0, k - 1,
                        expo1);  // sum_{i = 0}^{i < k - 1}R[k][i]^2 = ftmp1 * 2^expo1

    // If this check is false, we need to reenable
    // ftmp0.mul_2si(ftmp0, expo0 - expo1);
    FPLLL_DEBUG_CHECK(expo0 == expo1);

    ftmp1.sub(ftmp0, ftmp1);  // ||b[k]||^2 - sum_{i = 0}^{i < k - 1}R[k][i]^2 = ftmp1 * 2^expo1

    expo0 = m.get_row_expo(k - 1);

    // Here, delta * R(k - 1, k - 1)^2 = dR[k-1] * 2^(2*expo0). We want to compare
    //   delta * R(k - 1, k - 1)^2 <= ||b[k]||^2 - sum_{i = 0}^{i < k - 1}R[k][i]^2
    //   dR[k-1] * 2^(2*expo0) <= ftmp1 * 2^expo1
    //   dR[k-1] <= ftmp1 * 2^(expo1 - 2*expo0)
    ftmp1.mul_2si(ftmp1, expo1 - 2 * expo0);
#else  // MODIFIED_LOVASZ_TEST
    // Modified Lovasz test, following the comment above.
    // FIXME: probably not maintened.

    m.norm_square_R_row(ftmp1, k, k, m.get_n(),
                        expo1);       // sum_{i = k}^{i < n}R[k][i]^2 = ftmp1 * 2^expo1
#ifdef DEBUG
    m.get_R(ftmp0, k, k - 1, expo0);  // R(k, k - 1) = ftmp0 * 2^expo0
#else   // DEBUG
    m.get_R(ftmp0, k, k - 1);  // R(k, k - 1) = ftmp0 * 2^expo0
#endif  // DEBUG

    ftmp0.mul(ftmp0, ftmp0);  // R(k, k - 1)^2 = ftmp0 * 2^(2 * expo 0)

    // If this check is false, we need to reenable
    // ftmp0.mul_2si(ftmp0, 2 * expo0 - expo1);  // 2 * expo0 since R(k, k-1)^2 = ftmp0 *
    // (2^expo0)^2
    FPLLL_DEBUG_CHECK(2 * expo0 == expo1);

    ftmp1.add(ftmp0, ftmp1);  // sum_{i = k}^{i < n}R[k][i]^2 + R(k, k-1)^2 = ftmp1 * 2^expo1

    expo0 = m.get_row_expo(k - 1);

    // Here, delta * R(k - 1, k - 1)^2 = dR[k-1] * 2^(2*expo0). We want to compare
    //   delta * R(k - 1, k - 1)^2 <= sum_{i = k}^{i < n}R[k][i]^2 + R(k, k-1)^2
    //   dR[k-1] * 2^(2*expo0) <= ftmp1 * 2^expo1
    //   dR[k-1] <= ftmp1 * 2^(expo1 - 2*expo0)
    ftmp1.mul_2si(ftmp1, expo1 - 2 * expo0);
#endif  // MODIFIED_LOVASZ_TEST

    // Test if delta * R(k - 1, k - 1)^2 <= ||b[k]||^2 - sum_{i = 0}^{i < k - 1}R[k][i]^2 (depending
    // on the way ftmp1 is computed, this test can be slightly different, but the purpose keeps the
    // same)
    if (dR[k - 1].cmp(ftmp1) <= 0)
    {
      // Fully compute R[k], since all the coefficient except one of R[k] were computed during
      // size_reduction.
      m.update_R_last(k);
      // Compute delta_ * R(k, k)^2 = dR[k] * 2^(2*row_expo[k])
      compute_dR(k, delta_);

      // Heuristic precision check : when R(kappa-1,kappa-1) increases in a 2x2 up and down (see
      // hplll)
      if (prev_k == k + 1)
      {
        // R(k, k) = ftmp0 * 2^expo0
        m.get_R(ftmp0, k, k, expo0);
        ftmp1.mul_2si(prev_R[k], prev_expo[k] - expo0);

        if (ftmp0.cmp(ftmp1) > 0)
        {
          // cerr << "Anomaly: the norm increases for kappa = " << k << endl;
          return set_status(RED_NORM_HLLL_FAILURE);
        }
      }

      prev_k = k;
      // R(k, k) = prev_R[k] * 2^prev_expo[k]
      m.get_R(prev_R[k], k, k, prev_expo[k]);

      // b[k] is size reduced, now, size reduce b[k + 1]
      k++;

      if (k < m.get_d())
      {
        if (k > k_max)
        {
          // First time b[k] is discovered
          if (verbose)
          {
            cerr << "Discovering vector " << k + 1 << "/" << m.get_d()
                 << " cputime=" << cputime() - start_time << endl;
          }
          k_max = k;
          // Set R[k] and bf[k] to b[k], precompute ||b[k]||^2
          m.refresh_R_bf(k);
        }
        else
          // Set R[k] to b[k]. Indeed, it is not necessary to refresh bf[k], since b[k] has not
          // changed. However, it is mandatory to refresh R[k], since b[0] to b[k - 1] may have
          // changed, and then, it is necessary to recompute R[k].
          m.refresh_R(k);
      }
      else
        // if k == m.get_d(), then b[k] does not exist and the computation is ended
        return set_status(RED_SUCCESS);
    }
    else
    {
      // Swap b[k-1] and b[k] and other usefull variables
      m.swap(k - 1, k);
      prev_k = k;

      if (k - 1 == 0)
      {
        // Set R[0] to b[0] (bf[0] and other usefull variables were swaped previously)
        m.refresh_R(0);
        // Compute R[0]
        m.update_R_last(0);  // In this case, update_R(0) is exactly equal to update_R_last(0)
        // Precompute dR[0]: R[0]^2 * delta_ = dR[0] * 2^(2*row_expo[0])
        compute_dR(0, delta_);

        // Set R[1] to b[1] (bf[1] and other usefull variables were swaped previously)
        m.refresh_R(1);
        k = 1;
      }
      else
      {
        // Size reduce b[k - 1]
        k--;
        // Since b[k] was not changed, a previous computation of R[k][0..k-1] can be used instead of
        // recomputing R[k][0..k-1], which is the only interesting part to begin the size reduction
        // of b[k]
        m.recover_R(k);
      }
    }
  }
}

template <class ZT, class FT>
void HLLLReduction<ZT, FT>::size_reduction(int kappa, int size_reduction_end,
                                           int size_reduction_start)
{
  FPLLL_DEBUG_CHECK(kappa >= size_reduction_end);
  FPLLL_DEBUG_CHECK(size_reduction_start < size_reduction_end);
  FPLLL_DEBUG_CHECK(0 <= size_reduction_start);

  long expo0 = 0;
  long expo1 = 0;
  // If b[kappa] is reduced by at least one b[i], then reduced will be set to true.
  bool reduced = false;

  /*
   * Variables introduced in hplll (https://github.com/gilvillard/hplll)
   * See commit a6b29d1a23ca34000264e22608ef23a64e3cac9d
   * It seems that testing only one time the condition on line 7 of Algorithm 3 (incomplete size
   * reduction) is not sufficient to ensure that b[kappa] is size-reduced. We then test it two
   * times. If two times in a row, the condition is reached, then we consider that b[kappa] is
   * size-reduced.
   */
  bool not_stop      = true;
  bool prev_not_stop = true;

#ifndef HOUSEHOLDER_USE_SIZE_REDUCTION_TEST
  /*
   * Variables introduced in hplll (https://github.com/gilvillard/hplll)
   * See commit a6b29d1a23ca34000264e22608ef23a64e3cac9d
   * The loop is done while (||bk||^2 <= a * t, where t is the squared norm of bk before size
   * reduction and a = 2^(-cd).
   *   The condition is equivalent to -1/2 * log(a) <= 1/2*log(t) - log(||bk||), where 1/2*log(t) is
   *   the log of the squared norm of bk before size reduction. For c=0.1 and d=300,
   *   -1/2 * log_2(a) = 15, then the loop continues when the length of bk decreases from more than
   *   15 bits. For a=0.1, the same thing occurs with less than 2 bits. It experimentally allows to
   *   size reduce the vectors for dimension >= 260.
   * TODO: this hard coded value must be theoretically analyzed.
   */
  FT approx = 0.1;
#else   // HOUSEHOLDER_USE_SIZE_REDUCTION_TEST
  FT approx = sr;
#endif  // HOUSEHOLDER_USE_SIZE_REDUCTION_TEST

  m.update_R(kappa, false);

  /* Most likely, at this step, the next update_R(kappa, false) must modify some coefficients since
   * b will most likely be changed. If b is not modified during the size reduction, there will be
   * only a call to update_R_last(kappa), which automatically must set updated_R to false. We set
   * updated_R to false since it can be equal to true if recover_R was called that allows to avoid
   * an unuseful recomputation of R[kappa] with update_R.
   * TODO: maybe find a best place to use this function.
   */
  m.set_updated_R_false();

  do
  {
    // No b[i] reduced b[kappa]
    reduced = false;

    for (int i = size_reduction_end - 1; i >= size_reduction_start; i--)
    {
      m.get_R(ftmp1, kappa, i, expo1);  // R(kappa, i) = ftmp1 * 2^expo1
      m.get_R(ftmp0, i, i, expo0);      // R(i, i) = ftmp0 * 2^expo0

#ifndef HOUSEHOLDER_PRECOMPUTE_INVERSE
      ftmp1.div(ftmp1, ftmp0);  // R(kappa, i) / R(i, i) = ftmp1 * 2^(expo1 - expo0)
#else                           // HOUSEHOLDER_PRECOMPUTE_INVERSE
      ftmp1.mul(ftmp1,
                m.get_R_inverse_diag(i));  // R(kappa, i) / R(i, i) = ftmp1 * 2^(expo1 - expo0)
#endif                          // HOUSEHOLDER_PRECOMPUTE_INVERSE

      /* If T = mpfr or dpe, enable_row_expo must be false and then, expo1 - expo0 == 0 (required by
       * rnd_we with these types) */
      ftmp1.rnd_we(ftmp1, expo1 - expo0);  // rnd(R(kappa, i) / R(i, i)) = ftmp1 * 2^(expo1 - expo0)

      // ftmp1 * 2^(expo1 - expo0) is equal to -X[i] in Algorithm 3 of [MSV, ISSAC'09]
      ftmp1.neg(ftmp1);

      if (ftmp1.sgn() != 0)  // Equivalent to test if ftmp1 == 0
      {
        // Reduce b[kappa] and R[kappa] accordingly (Step 5 and Step 6 of Algorithm 3 of [MSV,
        // ISSAC'09])
        m.size_reduce(ftmp1, kappa, i);
        // b[kappa] was reduced by -ftmp1 * b[i]
        reduced = true;
      }
    }

    // If not reduced, b(kappa) has not changed. Computing ||b[kappa]||^2 is not necessary.
    // 1 > 2^(-cd)=sr since cd > 0.
    if (!reduced)
      return;
    else
    {
      // At this point, even if b has changed, the precomputed squared norm of b was for b before
      // the reduction
      m.get_norm_square_b(ftmp0, kappa, expo0);  // ||b[kappa]||^2 = t = ftmp0 * 2^expo
      // Since b has changed, R must be recomputed (latter in the implementation) and then R[kappa]
      // and bf[kappa] are set to b[kappa]. The squared norm of b is updated, then, the next call to
      // get_norm_square_b(..., kappa, ...) will get the squared norm of the current b.
      m.refresh_R_bf(kappa);
      m.get_norm_square_b(ftmp1, kappa, expo1);  // ||b[kappa]||^2 = ftmp1 * 2^expo1

      ftmp0.mul(approx, ftmp0);  // approx * t = ftmp0 * 2^expo0

      // TODO: Why not doing ftmp1.mul_2si(ftmp1, expo1 - expo0); ?
      // We want to compare
      //   ||b[k]||^2 > approx * t
      //   ftmp1 * 2^expo1 > ftmp0 * 2^expo0
      //   ftmp1 > ftmp0 * 2^(expo0-expo1)
      ftmp0.mul_2si(ftmp0, expo0 - expo1);

      // If (||b(kappa)||^2 > approx * t => ftmp1 > ftmp0), stop the loop.
      not_stop = (ftmp1.cmp(ftmp0) <= 0);

      // Update R(kappa, 0..kappa-1).
      m.update_R(kappa, false);

      if (prev_not_stop || not_stop)
        prev_not_stop = not_stop;  // Continue to try to reduce b(kappa).
      else
        return;  // b[kappa] should be size_reduced.
    }
  } while (true);
}

/*
 * Verify if the basis b inside m is (delta, eta)-hlll reduced.
 * Use a different implementation of the Householder transformation to compute R in this test than
 * the one used to reduced the basis.
 */
template <class ZT, class FT>
bool is_hlll_reduced(MatHouseholder<ZT, FT> &m, double delta, double eta, double theta)
{
  int i, j;
  // Temporary variables
  FT ftmp0, ftmp1, ftmp2;
  // FT version of delta and eta
  FT delta_ = delta;
  FT eta_   = eta;
  FT theta_ = theta;

  // Compute the R coefficients of b
  m.update_R_naively();

  // Exponent associated to ftmp0 and ftmp1 (respectively)
  long expo0 = 0;
  long expo1 = 0;
  long expo2 = 0;

  // Verify if |R(j, i)| <= eta * R(i, i) + theta * R(j, j) (weak size-reduction of Definition 2
  // [MSV'09].
  for (j = 0; j < m.get_d(); j++)
  {
    for (i = 0; i < j; i++)
    {
      m.get_R_naively(ftmp0, j, i, expo0);  // R(j, i) = ftmp0 * 2^expo0
      ftmp0.abs(ftmp0);                     // |R(j, i)| = |ftmp0| * 2^expo0
      m.get_R_naively(ftmp1, j, j, expo1);  // R(j, j) = ftmp1 * 2^expo1
      m.get_R_naively(ftmp2, i, i, expo2);  // R(i, i) = ftmp2 * 2^expo2

      FPLLL_DEBUG_CHECK(expo0 == expo1);

      ftmp1.mul(ftmp1, theta_);  // eta_ * R(i, i) = ftmp2 * 2^expo2
      ftmp2.mul(ftmp2, eta_);    // eta_ * R(i, i) = ftmp2 * 2^expo2

      // We want to test if
      //   |R(j, i)| <= eta * R(i, i) + theta * R(j, j)
      //   ftmp0 * 2^expo0 <= ftmp1 * 2^expo0 + ftmp2 * 2^expo2, since expo0 == expo1
      //   ftmp0 <= ftmp1 + ftmp2 * 2^(expo2 - expo0)
      ftmp2.mul_2si(ftmp2, expo2 - expo0);
      ftmp1.add(ftmp1, ftmp2);

      if (ftmp0.cmp(ftmp1) > 0)
        return false;
    }
  }

  // At this step, we verify if two consecutive vectors must be swapped during the hlll-reduction or
  // not (Lovasz's condition)
  for (i = 1; i < m.get_d(); i++)
  {
    m.get_R_naively(ftmp0, i - 1, i - 1, expo0);  // R(i - 1, i - 1) = ftmp0 * 2^expo0
    m.get_R_naively(ftmp1, i, i - 1, expo1);      // R(i, i - 1) = ftmp1 * 2^expo1
    m.get_R_naively(ftmp2, i, i, expo2);          // R(i, i) = ftmp2 * 2^expo2
    FPLLL_DEBUG_CHECK(expo0 == expo1);

    ftmp0.mul(ftmp0, ftmp0);
    ftmp1.mul(ftmp1, ftmp1);
    ftmp2.mul(ftmp2, ftmp2);
    expo0 = 2 * expo0;
    // expo1 = 2 * expo1; : not necessary, since expo0 == expo1
    expo2 = 2 * expo2;
    // Here: R(i - 1, i - 1)^2 = ftmp0 * 2^expo0
    // Here: R(i, i - 1)^2 = ftmp1 * 2^expo1
    // Here: R(i, i)^2 = ftmp2 * 2^expo2

    ftmp0.mul(ftmp0, delta);  // delta * R(i - 1, i - 1)^2 = delta * ftmp0 * 2^expo0

    // We want to test if
    //   delta * R(i - 1, i - 1)^2 <= R(i, i - 1)^2 + R(i, i)^2
    //   ftmp0 * 2^expo0 <= ftmp1 * 2^expo0 + ftmp2 * 2^expo2
    //   ftmp0 <= ftmp1 + ftmp2 * 2^(expo2 - expo0)
    ftmp2.mul_2si(ftmp2, expo2 - expo0);
    ftmp1.add(ftmp1, ftmp2);

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
                                            double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<double>>(MatHouseholder<Z_NR<long>, FP_NR<double>> &m,
                                           double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<double>>(MatHouseholder<Z_NR<double>, FP_NR<double>> &m,
                                             double delta, double eta, double theta);

#ifdef FPLLL_WITH_LONG_DOUBLE
template class HLLLReduction<Z_NR<long>, FP_NR<long double>>;
template class HLLLReduction<Z_NR<double>, FP_NR<long double>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<long double>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<long double>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<long double>> &m,
                                                 double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<long double>>(MatHouseholder<Z_NR<long>, FP_NR<long double>> &m,
                                                double delta, double eta, double theta);
template bool is_hlll_reduced<Z_NR<double>, FP_NR<long double>>(
    MatHouseholder<Z_NR<double>, FP_NR<long double>> &m, double delta, double eta, double theta);
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
                                             double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<qd_real>>(MatHouseholder<Z_NR<long>, FP_NR<qd_real>> &m,
                                            double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<qd_real>>(MatHouseholder<Z_NR<double>, FP_NR<qd_real>> &m,
                                              double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<dd_real>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<dd_real>> &m,
                                             double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<dd_real>>(MatHouseholder<Z_NR<long>, FP_NR<dd_real>> &m,
                                            double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<dd_real>>(MatHouseholder<Z_NR<double>, FP_NR<dd_real>> &m,
                                              double delta, double eta, double theta);
#endif

#ifdef FPLLL_WITH_DPE
template class HLLLReduction<Z_NR<long>, FP_NR<dpe_t>>;
template class HLLLReduction<Z_NR<double>, FP_NR<dpe_t>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<dpe_t>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<dpe_t>> &m,
                                           double delta, double eta, double theta);
template bool is_hlll_reduced<Z_NR<long>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<long>, FP_NR<dpe_t>> &m,
                                                        double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<dpe_t>>(MatHouseholder<Z_NR<double>, FP_NR<dpe_t>> &m,
                                            double delta, double eta, double theta);
#endif

template class HLLLReduction<Z_NR<long>, FP_NR<mpfr_t>>;
template class HLLLReduction<Z_NR<double>, FP_NR<mpfr_t>>;
template class HLLLReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>>;
template bool
is_hlll_reduced<Z_NR<mpz_t>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<mpz_t>, FP_NR<mpfr_t>> &m,
                                            double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<long>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<long>, FP_NR<mpfr_t>> &m,
                                           double delta, double eta, double theta);
template bool
is_hlll_reduced<Z_NR<double>, FP_NR<mpfr_t>>(MatHouseholder<Z_NR<double>, FP_NR<mpfr_t>> &m,
                                             double delta, double eta, double theta);

FPLLL_END_NAMESPACE
