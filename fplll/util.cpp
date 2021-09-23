/* Copyright (C) 2011 Xavier Pujol.

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

#include "util.h"

#ifdef DEBUG
int debug_depth = 0;
#endif

FPLLL_BEGIN_NAMESPACE

enum MinPrecAlgo
{
  MINPREC_GSO,
  MINPREC_L2
};

/* State of the random generator (declared in nr.h, must be defined in exactly
   one source file) */
bool RandGen::initialized = false;
gmp_randstate_t RandGen::gmp_state;

static int compute_min_prec(double &rho, int d, double delta, double eta, double epsilon,
                            MinPrecAlgo algo)
{
  int old_prec = FP_NR<mpfr_t>::set_prec(53);
  FP_NR<mpfr_t> f_minprec, f_rho, f_d, f_eta, f_delta, f_epsilon, tmp1, tmp2;

  // These four conversions are exact
  f_d       = static_cast<double>(d);
  f_eta     = eta;
  f_delta   = delta;
  f_epsilon = epsilon;
  if (algo == MINPREC_L2)
  {
    // eta - 0.5 is an exact fp operation
    if (f_epsilon > eta - 0.5)
      f_epsilon = eta - 0.5;
    tmp1 = 1.0;
    tmp1.sub(tmp1, f_delta, GMP_RNDD);
    if (f_epsilon > tmp1)
      f_epsilon = tmp1;
    // now fEpsilon <= min(epsilon, eta - 0.5, 1 - delta);
  }
  // Computes tmp1 >= (1 + eta) ^ 2 + epsilon
  tmp1 = 1.0;                       // exact
  tmp1.add(f_eta, tmp1, GMP_RNDU);  // >= 1 + eta
  tmp1.mul(tmp1, tmp1, GMP_RNDU);   // >= (1 + eta) ^ 2
  tmp1.add(tmp1, f_epsilon, GMP_RNDU);
  // Computes tmp2 <= delta - eta ^ 2
  tmp2.mul(f_eta, f_eta, GMP_RNDU);
  tmp2.sub(f_delta, tmp2, GMP_RNDD);
  FPLLL_CHECK(tmp2 > 0, "invalid LLL parameters, eta must be < sqrt(delta)");
  // Computes rho >= ((1 + eta) ^ 2 + epsilon) / (delta - eta ^ 2)
  f_rho.div(tmp1, tmp2, GMP_RNDU);
  rho = f_rho.get_d(GMP_RNDU);

  /* Computes minprec >= constant + 2 * log2(d) - log2(epsilon) + d * log2(rho)
     (constant = 5 for GSO, 10 for LLL) */
  tmp1.log(f_d, GMP_RNDU);         // >= log(d)
  tmp1.mul_2si(tmp1, 1);           // >= 2 * log(d)
  tmp2.log(f_epsilon, GMP_RNDD);   // <= log(epsilon) (<= 0)
  tmp1.sub(tmp1, tmp2, GMP_RNDU);  // >= 2 * log(d) - log(epsilon)
  tmp2.log(f_rho, GMP_RNDU);       // >= log(rho)
  tmp2.mul(f_d, tmp2, GMP_RNDU);   // >= d * log(rho)
  tmp1.add(tmp1, tmp2, GMP_RNDU);  // >= 2*log(d)-log(epsilon)+d*log(rho)
  tmp2 = 2.0;                      // exact
  tmp2.log(tmp2, GMP_RNDD);        // <= log(2)
  tmp1.div(tmp1, tmp2, GMP_RNDU);  // >= 2*log2(d)-log2(epsilon)+d*log2(rho)
  tmp2 = (algo == MINPREC_L2) ? 10.0 : 5.0;
  f_minprec.add(tmp1, tmp2, GMP_RNDU);
  int minprec = static_cast<int>(ceil(f_minprec.get_d(GMP_RNDU)));
  mpfr_free_cache();
  FP_NR<mpfr_t>::set_prec(old_prec);
  return minprec;
}

int gso_min_prec(double &rho, int d, double delta, double eta, double epsilon)
{
  return compute_min_prec(rho, d, delta, eta, epsilon, MINPREC_GSO);
}

int l2_min_prec(int d, double delta, double eta, double epsilon)
{
  double rho;
  return compute_min_prec(rho, d, delta, eta, epsilon, MINPREC_L2);
}

int hlll_min_prec(int d_i, int n_i, double delta_d, double eta_d, double theta_d, double c_d)
{
  FPLLL_CHECK(delta_d < 1.0 && delta_d >= 0.25, "delta must be in [1/4, 1).");
  FPLLL_CHECK(theta_d >= 0.0, "theta must be positive.");
  FPLLL_CHECK(eta_d >= 0.5, "theta must be larger than or equal to 0.5.");
  FPLLL_CHECK(eta_d - theta_d > 0.5, "eta - theta must be larger than 0.5.");

  int old_prec = FP_NR<mpfr_t>::set_prec(53);
  FP_NR<mpfr_t> d, n, delta, eta, theta, c, alpha, c0, c1, rho, phi, p0, p;
  FP_NR<mpfr_t> ftmp0, ftmp1, ftmp2, ftmp3, ftmp4;

  d     = d_i;
  n     = n_i;
  delta = delta_d;
  eta   = eta_d;
  theta = theta_d;
  c     = c_d;

  // ftmp0 = (1 + theta^2) * delta - eta^2
  ftmp0 = (1.0 + theta * theta) * delta - eta * eta;
  // ftmp0 = sqrt((1 + theta^2) * delta - eta^2)
  ftmp0.sqrt(ftmp0);

  // alpha = theta * eta + sqrt((1 + theta^2) * delta - eta^2) / (delta - eta^2)
  alpha = (theta * eta + ftmp0) / (delta - eta * eta);

  // ftmp0 = 3 / 2
  ftmp0 = 3.0 / 2.0;
  // ftmp0 = sqrt(3 / 2)
  ftmp0.sqrt(ftmp0);
  // ftmp1 = 1 - eta - theta
  ftmp1 = 1.0 - eta - theta;
  // ftmp1 = |1 - eta - theta|
  ftmp1.abs(ftmp1);
  ftmp2 = 6.0;
  // ftmp2 = sqrt(6)
  ftmp2.sqrt(ftmp2);
  // ftmp3 = 1 + d * eta^2
  ftmp3 = 1.0 + d * eta * eta;
  // ftmp3 = sqrt(1 + d * eta^2)
  ftmp3.sqrt(ftmp3);
  // ftmp4 = sqrt(d)
  ftmp4.sqrt(d);
  // ftmp0 = 1 + |1 - eta - theta| * alpha / ((eta + theta) * (-1 + sqrt(3/2)))
  ftmp0 = (1.0 + ftmp1 * alpha) / ((eta + theta) * (-1.0 + ftmp0));
  // ftmp1 = 4 * sqrt(6) / (1 + eta) * sqrt(1 + d * eta^2)
  ftmp1 = 4.0 * ftmp2 / (1.0 + eta) * ftmp3;
  // ftmp0 = max(1 + |1 - eta - theta| * alpha / ((eta + theta) * (-1 + sqrt(3/2))),
  //             4 * sqrt(6) / (1 + eta) * sqrt(1 + d * eta^2))
  ftmp0.max_f(ftmp1);
  // c0 = max(...) * n * sqrt(d)
  c0 = ftmp0 * n * ftmp4;

  // c1 = 8 * d * (n + 9) * c0
  c1 = 8.0 * d * (n + 9.0) * c0;

  // rho = (1 + eta + theta) * alpha
  rho = (1.0 + eta + theta) * alpha;

  // ftmp0 = rho^(d + 1) (since we want to compute phi(d))
  ftmp0.pow_si(rho, d_i + 1);
  // phi(d) = c1 * (1 + 1 / theta) * ftmp0
  phi = c1 * (1.0 + 1.0 / theta) * ftmp0;

  // ftmp0 = alpha^d
  ftmp0.pow_si(alpha, d_i);
  // ftmp0 = log(d^3 * phi(d) * alpha^d / theta)
  ftmp0.log(d * d * d * phi * ftmp0 / theta);
  // ftmp1 = log(2)
  ftmp1.log(2);
  // ftmp0 = log(d^3 * phi(d) * alpha^d / theta) / log(2)
  ftmp0 = ftmp0 / ftmp1;
  // p0 = log2(d^3 * phi(d) * alpha^d / theta) + 16 + c * d / 2
  p0 = ftmp0 + 16.0 + c * d / 2.0;

  // ftmp0 = log(1 - delta)
  ftmp0.log(1.0 - delta);
  // ftmp0 = log(1 - delta) / log(2)
  ftmp0 = ftmp0 / ftmp1;
  // ftmp2 = log(eta - theta - 1/2)
  ftmp2.log(eta - theta - 0.5);
  // ftmp2 = log(eta - theta - 1/2) / log(2)
  ftmp2 = ftmp2 / ftmp1;
  // p = p0 + 1 - log2(1 - delta) - log2(eta - theta - 1 / 2)
  p = p0 + 1.0 - ftmp0 - ftmp2;

  // Convert p in int
  int p_i = static_cast<int>(ceil(p.get_d(GMP_RNDU)));

  FP_NR<mpfr_t>::set_prec(old_prec);

  return p_i;
}

/**
 * Computes the volume of a d-dimensional hypersphere of radius 1.
 */
void sphere_volume(FP_NR<mpfr_t> &volume, int d)
{
  FP_NR<mpfr_t> rtmp1;
  volume = pow(M_PI, (double)(d / 2));

  if (d % 2 == 0)
    for (int i = 1; i <= d / 2; i++)
    {
      rtmp1 = (double)i;
      volume.div(volume, rtmp1);
    }
  else
    for (int i = 0; i <= d / 2; i++)
    {
      rtmp1 = 2.0 / (double)(2 * i + 1);
      volume.mul(volume, rtmp1);
    }
}

/**
 * Estimates the cost of the enumeration for SVP.
 */
void cost_estimate(FP_NR<mpfr_t> &cost, const FP_NR<mpfr_t> &bound, const Matrix<FP_NR<mpfr_t>> &r,
                   int dimMax)
{
  FP_NR<mpfr_t> det, level_cost, tmp1;
  det  = 1.0;
  cost = 0.0;

  for (int i = dimMax - 1; i >= 0; i--)
  {
    tmp1.div(bound, r(i, i));
    det.mul(det, tmp1);

    level_cost.sqrt(det);
    sphere_volume(tmp1, dimMax - i);
    level_cost.mul(level_cost, tmp1);

    cost.add(cost, level_cost);
  }
}

const char *get_red_status_str(int status)
{
  if (status >= 0 && status < RED_STATUS_MAX)
    return RED_STATUS_STR[status];
  else
    return "unknown error";
}

template <class ZT> void zeros_first(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv_t)
{
  int i, d = b.get_rows();
  for (i = d; i > 0 && b[i - 1].is_zero(); i--)
  {
  }
  if (i > 0 && i < d)
  {
    b.rotate(0, i, d - 1);
    if (!u.empty())
      u.rotate(0, i, d - 1);
    if (!u_inv_t.empty())
      u_inv_t.rotate(0, i, d - 1);
  }
}

template <class ZT> void zeros_last(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv_t)
{
  int i, d = b.get_rows();
  for (i = 0; i < d && b[i].is_zero(); i++)
  {
  }
  if (i > 0 && i < d)
  {
    b.rotate(0, i, d - 1);
    if (!u.empty())
      u.rotate(0, i, d - 1);
    if (!u_inv_t.empty())
      u_inv_t.rotate(0, i, d - 1);
  }
}

template void zeros_first<mpz_t>(ZZ_mat<mpz_t> &, ZZ_mat<mpz_t> &, ZZ_mat<mpz_t> &);
template void zeros_last<mpz_t>(ZZ_mat<mpz_t> &, ZZ_mat<mpz_t> &, ZZ_mat<mpz_t> &);

#ifdef FPLLL_WITH_ZLONG
template void zeros_first<long>(ZZ_mat<long> &, ZZ_mat<long> &, ZZ_mat<long> &);
template void zeros_last<long>(ZZ_mat<long> &, ZZ_mat<long> &, ZZ_mat<long> &);
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template void zeros_first<double>(ZZ_mat<double> &, ZZ_mat<double> &, ZZ_mat<double> &);
template void zeros_last<double>(ZZ_mat<double> &, ZZ_mat<double> &, ZZ_mat<double> &);
#endif

FPLLL_END_NAMESPACE
