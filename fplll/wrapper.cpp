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

#include "wrapper.h"
#include "lll.h"
#include "util.h"

FPLLL_BEGIN_NAMESPACE

/* prec=53, eta=0.501, dim < dim_double_max [ (delta / 100.0) + 25 ] */
const double dim_double_max[75] = {
    0,     26,    29.6,  28.1,  31.1,  32.6,  34.6,  34,    37.7,  38.8,  39.6,  41.8,  40.9,
    43.6,  44.2,  47,    46.8,  50.6,  49.1,  51.5,  52.5,  54.8,  54.6,  57.4,  57.6,  59.9,
    61.8,  62.3,  64.5,  67.1,  68.8,  68.3,  69.9,  73.1,  74,    76.1,  76.8,  80.9,  81.8,
    83,    85.3,  87.9,  89,    90.1,  89,    94.6,  94.8,  98.7,  99,    101.6, 104.9, 106.8,
    108.2, 107.4, 110,   112.7, 114.6, 118.1, 119.7, 121.8, 122.9, 126.6, 128.6, 129,   133.6,
    126.9, 135.9, 139.5, 135.2, 137.2, 139.3, 142.8, 142.4, 142.5, 145.4};

const double eta_dep[10] = {1.,       // 0.5
                            1.,       // 0.55
                            1.0521,   // 0.6
                            1.1254,   // 0.65
                            1.2535,   // 0.7
                            1.3957,   // 0.75
                            1.6231,   // 0.8
                            1.8189,   // 0.85
                            2.1025,   // 0.9
                            2.5117};  // 0.95

Wrapper::Wrapper(IntMatrix &b, IntMatrix &u, IntMatrix &u_inv, double delta, double eta, int flags)
    : status(RED_SUCCESS), b(b), u(u), u_inv(u_inv), delta(delta), eta(eta), use_long(false),
      last_early_red(0)
{
  n            = b.get_cols();
  d            = b.get_rows();
  this->flags  = flags;
  max_exponent = b.get_max_exp() + (int)ceil(0.5 * log2((double)d * n));

  // Computes the parameters required for the proved version
  good_prec = l2_min_prec(d, delta, eta, LLL_DEF_EPSILON);
}

bool Wrapper::little(int kappa, int precision)
{
  /*one may add here dimension arguments with respect to eta and delta */
  int dm = (int)(delta * 100. - 25.);
  if (dm < 0)
    dm = 0;
  if (dm > 74)
    dm = 74;

  int em = (int)((eta - 0.5) * 20);
  if (em < 0)
    em = 0;
  if (em > 9)
    em = 9;

  double p = max(1.0, precision / 53.0);

  p *= eta_dep[em]; /* eta dependance */
  p *= dim_double_max[dm];
  // cerr << kappa << " compared to " << p << endl;
  return kappa < p;
}

/**
 * main function. Method determines whether heuristic, fast or proved
 */
template <class Z, class F>
int Wrapper::call_lll(ZZ_mat<Z> &bz, ZZ_mat<Z> &uz, ZZ_mat<Z> &u_invZ, LLLMethod method,
                      int precision, double delta, double eta)
{
  typedef Z_NR<Z> ZT;
  typedef FP_NR<F> FT;

  if (flags & LLL_VERBOSE)
  {
    cerr << "====== Wrapper: calling " << LLL_METHOD_STR[method] << "<" << num_type_str<Z>() << ","
         << num_type_str<F>() << "> method";
    if (precision > 0)
    {
      cerr << " (precision=" << precision << ")";
    }
    cerr << " ======" << endl;
  }

  int gso_flags = 0;
  if (method == LM_PROVED)
    gso_flags |= GSO_INT_GRAM;
  if (method == LM_FAST)
    gso_flags |= GSO_ROW_EXPO;
  if (method != LM_PROVED && precision == 0)
    gso_flags |= GSO_OP_FORCE_LONG;

  int oldprec = Float::get_prec();
  if (precision > 0)
  {
    Float::set_prec(precision);
  }
  MatGSO<ZT, FT> m_gso(bz, uz, u_invZ, gso_flags);
  LLLReduction<ZT, FT> lll_obj(m_gso, delta, eta, flags);
  lll_obj.last_early_red = last_early_red;
  lll_obj.lll();
  status         = lll_obj.status;
  last_early_red = max(last_early_red, lll_obj.last_early_red);
  if (precision > 0)
  {
    Float::set_prec(oldprec);
  }

  if (flags & LLL_VERBOSE)
  {
    cerr << "====== Wrapper: end of " << LLL_METHOD_STR[method] << " method ======\n" << endl;
  }

  if (lll_obj.status == RED_SUCCESS)
    return 0;
  else if (lll_obj.status == RED_GSO_FAILURE || lll_obj.status == RED_BABAI_FAILURE)
    return lll_obj.final_kappa;
  else
    return -1;
}

/**
 * pass the method to call_lll()
 */
template <class F> int Wrapper::fast_lll(double delta, double eta)
{
  return call_lll<mpz_t, F>(b, u, u_inv, LM_FAST, 0, delta, eta);
}

template <class Z, class F>
int Wrapper::heuristic_lll(ZZ_mat<Z> &bz, ZZ_mat<Z> &uz, ZZ_mat<Z> &u_invZ, int precision,
                           double delta, double eta)
{
  return call_lll<Z, F>(bz, uz, u_invZ, LM_HEURISTIC, precision, delta, eta);
}

template <class Z, class F>
int Wrapper::proved_lll(ZZ_mat<Z> &bz, ZZ_mat<Z> &uz, ZZ_mat<Z> &u_invZ, int precision, double delta,
                        double eta)
{
  return call_lll<Z, F>(bz, uz, u_invZ, LM_PROVED, precision, delta, eta);
}

/**
 * In heuristic_loop(), we only use double or dpe_t or mpfr_t.
 */
int Wrapper::heuristic_loop(int precision)
{
  int kappa;

  if (precision > numeric_limits<double>::digits)
    kappa = heuristic_lll<mpz_t, mpfr_t>(b, u, u_inv, precision, delta, eta);
  else
  {
#ifdef FPLLL_WITH_DPE
    kappa = heuristic_lll<mpz_t, dpe_t>(b, u, u_inv, 0, delta, eta);
#else
    kappa = heuristic_lll<mpz_t, mpfr_t>(b, u, u_inv, precision, delta, eta);
#endif
  }

  if (kappa == 0)
    return 0;  // Success
  else if (precision < good_prec && !little(kappa, precision))
    return heuristic_loop(increase_prec(precision));
  else
    return proved_loop(precision);
}

int Wrapper::proved_loop(int precision)
{
  int kappa;
#ifdef FPLLL_WITH_QD
  if (precision > PREC_DD)
#else
  if (precision > numeric_limits<double>::digits)
#endif
    kappa = proved_lll<mpz_t, mpfr_t>(b, u, u_inv, precision, delta, eta);
  else if (max_exponent * 2 > MAX_EXP_DOUBLE)
  {
#ifdef FPLLL_WITH_DPE
    kappa = proved_lll<mpz_t, dpe_t>(b, u, u_inv, 0, delta, eta);
#else
    kappa                = proved_lll<mpz_t, mpfr_t>(b, u, u_inv, precision, delta, eta);
#endif
  }
#ifdef FPLLL_WITH_QD
  else if (precision > numeric_limits<double>::digits)
    kappa = proved_lll<mpz_t, dd_real>(b, u, u_inv, precision, delta, eta);
#endif
  else
    kappa = proved_lll<mpz_t, double>(b, u, u_inv, 0, delta, eta);

  if (kappa == 0)
    return 0;  // Success
  else if (precision < good_prec)
    return proved_loop(increase_prec(precision));
  else
    return -1;  // This point should never be reached
}

/**
 * last call to LLL. Need to be proved_lll.
 */
int Wrapper::last_lll()
{

/* <long, FT> */
#ifdef FPLLL_WITH_ZLONG
  if (use_long)
  {
    int kappa;
    if (good_prec <= numeric_limits<double>::digits)
      kappa = proved_lll<long, double>(b_long, u_long, u_inv_long, good_prec, delta, eta);
#ifdef FPLLL_WITH_QD
    else if (good_prec <= PREC_DD)
      kappa = proved_lll<long, dd_real>(b_long, u_long, u_inv_long, good_prec, delta, eta);
#endif
    else
      kappa = proved_lll<long, mpfr_t>(b_long, u_long, u_inv_long, good_prec, delta, eta);
    return kappa;
  }
#endif

/* <mpfr, FT> */
#ifdef FPLLL_WITH_DPE
  if (good_prec <= numeric_limits<double>::digits)
    return proved_lll<mpz_t, dpe_t>(b, u, u_inv, good_prec, delta, eta);
#ifdef FPLLL_WITH_QD
  else if (good_prec <= PREC_DD)
    return proved_lll<mpz_t, dd_real>(b, u, u_inv, good_prec, delta, eta);
#endif
#endif
  return proved_lll<mpz_t, mpfr_t>(b, u, u_inv, good_prec, delta, eta);
}

/**
 * Wrapper.lll() calls
 *  - heuristic_lll()
 *  - fast_lll()
 *  - proved_lll()
 */
bool Wrapper::lll()
{
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;

#ifdef FPLLL_WITH_ZLONG
  bool heuristic_with_long =
      max_exponent < numeric_limits<long>::digits - 2 && u.empty() && u_inv.empty();
  bool proved_with_long =
      2 * max_exponent < numeric_limits<long>::digits - 2 && u.empty() && u_inv.empty();
#else
  bool heuristic_with_long = false, proved_with_long = false;
#endif

  int kappa;

  /* small matrix */
  if (heuristic_with_long)
  {
#ifdef FPLLL_WITH_ZLONG
    set_use_long(true);
    /* try heuristic_lll <long, double> */
    heuristic_lll<long, double>(b_long, u_long, u_inv_long, 0, delta, eta);
#endif
  }
  /* large matrix */
  else
  {

    /* try fast_lll<mpz_t, double> */
    kappa           = fast_lll<double>(delta, eta);
    bool lll_failure = (kappa != 0);
    int last_prec;

/* try fast_lll<mpz_t, long double> */
#ifdef FPLLL_WITH_LONG_DOUBLE
    if (lll_failure)
    {
      kappa      = fast_lll<long double>(delta, eta);
      lll_failure = kappa != 0;
    }
    last_prec = numeric_limits<long double>::digits;
#else
    last_prec = numeric_limits<double>::digits;
#endif

/* try fast_lll<mpz_t, dd_real> */
#ifdef FPLLL_WITH_QD
    if (lll_failure)
    {
      kappa      = fast_lll<dd_real>(delta, eta);
      lll_failure = kappa != 0;
    }
    last_prec = PREC_DD;
#else
#ifdef FPLLL_WITH_LONG_DOUBLE
    last_prec = numeric_limits<long double>::digits;
#else
    last_prec = numeric_limits<double>::digits;
#endif
#endif

    /* loop */
    if (lll_failure)
    {
      int prec_d = numeric_limits<double>::digits;
      if (little(kappa, last_prec))
        kappa = proved_loop(prec_d);
      else
        kappa = heuristic_loop(increase_prec(prec_d));
    }
  }

  set_use_long(proved_with_long);
  /* final LLL */
  kappa = last_lll();
  set_use_long(false);
  return kappa == 0;
}

/**
 * set blong <-- b
 */
void Wrapper::set_use_long(bool value)
{
#ifdef FPLLL_WITH_ZLONG
  if (!use_long && value)
  {
    if (b_long.empty())
    {
      b_long.resize(d, n);
    }
    for (int i = 0; i < d; i++)
      for (int j = 0; j < n; j++)
        b_long(i, j) = b(i, j).get_si();
  }
  else if (use_long && !value)
  {
    for (int i = 0; i < d; i++)
      for (int j = 0; j < n; j++)
        b(i, j) = b_long(i, j).get_si();
  }
  use_long = value;
#endif
}

int Wrapper::increase_prec(int precision) { return min(precision * 2, good_prec); }

FPLLL_END_NAMESPACE
