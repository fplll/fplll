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
#include "hlll.h"
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

Wrapper::Wrapper(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, ZZ_mat<mpz_t> &u_inv, double delta, double eta,
                 int flags)
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

// Constructor for HLLL
Wrapper::Wrapper(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, ZZ_mat<mpz_t> &u_inv, double delta, double eta,
                 double theta, double c, int flags)
    : status(RED_SUCCESS), b(b), u(u), u_inv(u_inv), delta(delta), eta(eta), use_long(false),
      last_early_red(-1), theta(theta), c(c)
{
  n           = b.get_cols();
  d           = b.get_rows();
  this->flags = flags;

  // Computes the parameters required for the proved version
  good_prec = hlll_min_prec(d, n, delta, eta, theta, c);
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

  int old_prec = FP_NR<mpfr_t>::get_prec();
  if (precision > 0)
  {
    FP_NR<mpfr_t>::set_prec(precision);
  }
  MatGSO<ZT, FT> m_gso(bz, uz, u_invZ, gso_flags);
  LLLReduction<ZT, FT> lll_obj(m_gso, delta, eta, flags);
  lll_obj.last_early_red = last_early_red;
  lll_obj.lll();
  status         = lll_obj.status;
  last_early_red = max(last_early_red, lll_obj.last_early_red);
  if (precision > 0)
  {
    FP_NR<mpfr_t>::set_prec(old_prec);
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
int Wrapper::proved_lll(ZZ_mat<Z> &bz, ZZ_mat<Z> &uz, ZZ_mat<Z> &u_invZ, int precision,
                        double delta, double eta)
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
    kappa = proved_lll<mpz_t, mpfr_t>(b, u, u_inv, precision, delta, eta);
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
  {
    max_exponent = b.get_max_exp() + (int)ceil(0.5 * log2((double)d * n));
    if (max_exponent * 2 < MAX_EXP_DOUBLE)
    {
      return proved_lll<mpz_t, dd_real>(b, u, u_inv, good_prec, delta, eta);
    }
  }
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
    kappa            = fast_lll<double>(delta, eta);
    bool lll_failure = (kappa != 0);
    int last_prec;

    /* try fast_lll<mpz_t, long double> */
#ifdef FPLLL_WITH_LONG_DOUBLE
    if (lll_failure)
    {
      kappa       = fast_lll<long double>(delta, eta);
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
      kappa       = fast_lll<dd_real>(delta, eta);
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

/* Set of methods for the HLLL wrapper. */

/**
 * main function to call hlll
 * Return true if success, false otherwise
 */
template <class F> bool Wrapper::call_hlll(LLLMethod method, int precision)
{
  typedef FP_NR<F> FT;

  if (flags & LLL_VERBOSE)
  {
    cerr << "====== Wrapper: calling " << HLLL_METHOD_STR[method] << "<mpz_t," << num_type_str<F>()
         << "> method";
    if (precision > 0)
    {
      cerr << " (precision=" << precision << ")";
    }
    cerr << " ======" << endl;
  }

  int householder_flags = 0;
  if (method == LM_PROVED)
    householder_flags |= HOUSEHOLDER_DEFAULT;
  if (method == LM_FAST)
    householder_flags |= HOUSEHOLDER_ROW_EXPO | HOUSEHOLDER_OP_FORCE_LONG;

  int old_prec = FP_NR<mpfr_t>::get_prec();

  if (precision > 0)
    FP_NR<mpfr_t>::set_prec(precision);

  MatHouseholder<Z_NR<mpz_t>, FT> m(b, u, u_inv, householder_flags);
  HLLLReduction<Z_NR<mpz_t>, FT> hlll_obj(m, delta, eta, theta, c, flags);
  hlll_obj.hlll();
  int status = hlll_obj.get_status();

  if (precision > 0)
    FP_NR<mpfr_t>::set_prec(old_prec);

  if (flags & LLL_VERBOSE)
  {
    cerr << "====== Wrapper: end of " << HLLL_METHOD_STR[method] << " method ======\n" << endl;
  }

  if (status == RED_SUCCESS)
    return true;
  else
    return false;
}

/**
 * last call to LLL. Need to be proved_lll.
 */
bool Wrapper::last_hlll()
{
/* <mpfr, FT> */
#ifdef FPLLL_WITH_DPE
  if (good_prec <= numeric_limits<double>::digits)
    return proved_hlll<dpe_t>(good_prec);
#ifdef FPLLL_WITH_QD
  else if (good_prec <= PREC_DD)
    return proved_hlll<dd_real>(good_prec);
#endif  // FPLLL_WITH_QD
#endif  // FPLLL_WITH_DPE
  return proved_hlll<mpfr_t>(good_prec);
}
int Wrapper::hlll_proved_loop(int precision)
{
  bool status = proved_hlll<mpfr_t>(precision);

  if (status)
    return 0;  // Success
  else if (precision < good_prec)
    return hlll_proved_loop(increase_prec(precision));
  else
    return -1;  // This point should never be reached
}

/**
 * pass the method to call_hlll()
 */
template <class F> bool Wrapper::fast_hlll() { return call_hlll<F>(LM_FAST, 0); }

template <class F> bool Wrapper::proved_hlll(int precision)
{
  return call_hlll<F>(LM_PROVED, precision);
}

bool Wrapper::hlll()
{
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;

  int last_prec      = numeric_limits<double>::digits;
  bool hlll_complete = false;

// TODO: since classical lll is faster than hlll for dim <~ 160, maybe we can
// use fast_lll<double>() at the beginning of hlll, before calling
// fast_hlll<double>()
// Something like the one used in #if 0 should work
#if 0
  /* try fast_lll<mpz_t, double> */
  int kappa        = fast_lll<double>(delta, eta);
  bool lll_failure = (kappa != 0);

  /* try fast_lll<mpz_t, double> */
  if (lll_failure)
    hlll_complete = fast_hlll<double>();
  else
    hlll_complete = true;
#else   // 0
  hlll_complete = fast_hlll<double>();
#endif  // 0

  /* try fast_hlll<mpz_t, long double> */
#ifdef FPLLL_WITH_LONG_DOUBLE
  if (!hlll_complete)
  {
    hlll_complete = fast_hlll<long double>();
    last_prec     = numeric_limits<long double>::digits;
  }
#endif  // FPLLL_WITH_LONG_DOUBLE

  /* try fast_hlll<mpz_t, dd_real> */
#ifdef FPLLL_WITH_QD
  if (!hlll_complete)
  {
    hlll_complete = fast_hlll<dd_real>();
    last_prec     = PREC_DD;
  }
#endif  // FPLLL_WITH_QD

  /* loop */
  if (!hlll_complete)
    hlll_complete = hlll_proved_loop(last_prec);

  hlll_complete = last_hlll();

  return hlll_complete == RED_SUCCESS;
}

/**
 * LLL with a typical method "proved or heuristic or fast".
 * @proved:     exact gram +   exact rowexp +   exact rowaddmul
 * @heuristic:  approx. gram +   exact rowexp +   exact rowaddmul
 * @fast:       approx. gram + approx. rowexp + approx. rowaddmul
 *    (double, long double, dd_real, qd_real)
 */
template <class ZT, class FT>
int lll_reduction_zf(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv, double delta, double eta,
                     LLLMethod method, int flags)
{
  int gso_flags = 0;
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;
  if (method == LM_PROVED)
    gso_flags |= GSO_INT_GRAM;
  if (method == LM_FAST)
    gso_flags |= GSO_ROW_EXPO | GSO_OP_FORCE_LONG;
  MatGSO<Z_NR<ZT>, FP_NR<FT>> m_gso(b, u, u_inv, gso_flags);
  LLLReduction<Z_NR<ZT>, FP_NR<FT>> lll_obj(m_gso, delta, eta, flags);
  lll_obj.lll();
  return lll_obj.status;
}

template <class ZT>
int lll_reduction_wrapper(ZZ_mat<ZT> &, ZZ_mat<ZT> &, ZZ_mat<ZT> &, double, double, FloatType, int,
                          int)
{
  FPLLL_ABORT("The wrapper method works only with integer type mpz");
  return RED_LLL_FAILURE;
}

template <>
int lll_reduction_wrapper(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, ZZ_mat<mpz_t> &u_inv, double delta,
                          double eta, FloatType float_type, int precision, int flags)
{
  FPLLL_CHECK(float_type == FT_DEFAULT,
              "The floating point type cannot be specified with the wrapper method");
  FPLLL_CHECK(precision == 0, "The precision cannot be specified with the wrapper method");
  Wrapper wrapper(b, u, u_inv, delta, eta, flags);
  wrapper.lll();
  zeros_first(b, u, u_inv);
  return wrapper.status;
}

/**
 * Main function called from call_lll().
 */
template <class ZT>
int lll_reduction_z(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv, double delta, double eta,
                    LLLMethod method, IntType int_type, FloatType float_type, int precision,
                    int flags)
{

  /* switch to wrapper */
  if (method == LM_WRAPPER)
    return lll_reduction_wrapper(b, u, u_inv, delta, eta, float_type, precision, flags);

  FPLLL_CHECK(!(method == LM_PROVED && (flags & LLL_EARLY_RED)),
              "LLL method 'proved' with early reduction is not implemented");

  /* computes the parameters required for the proved version */
  int good_prec = l2_min_prec(b.get_rows(), delta, eta, LLL_DEF_EPSILON);

  /* sets the parameters and checks the consistency */
  int sel_prec = 0;
  if (method == LM_PROVED)
  {
    sel_prec = (precision != 0) ? precision : good_prec;
  }
  else
  {
    sel_prec = (precision != 0) ? precision : PREC_DOUBLE;
  }

  FloatType sel_ft = float_type;

  /* if manually input precision */
  if (precision != 0)
  {
    if (sel_ft == FT_DEFAULT)
    {
      sel_ft = FT_MPFR;
    }
    FPLLL_CHECK(sel_ft == FT_MPFR,
                "The floating type must be mpfr when the precision is specified");
  }

  if (sel_ft == FT_DEFAULT)
  {
    if (method == LM_FAST)
      sel_ft = FT_DOUBLE;
#ifdef FPLLL_WITH_DPE
    else if (sel_prec <= static_cast<int>(FP_NR<dpe_t>::get_prec()))
      sel_ft = FT_DPE;
#endif
#ifdef FPLLL_WITH_QD
    else if (sel_prec <= static_cast<int>(FP_NR<dd_real>::get_prec()))
      sel_ft = FT_DD;
    else if (sel_prec <= static_cast<int>(FP_NR<qd_real>::get_prec()))
      sel_ft = FT_QD;
#endif
    else
      sel_ft = FT_MPFR;
  }
  else if (method == LM_FAST &&
           (sel_ft != FT_DOUBLE && sel_ft != FT_LONG_DOUBLE && sel_ft != FT_DD && sel_ft != FT_QD))
  {
    FPLLL_ABORT("'double' or 'long double' or 'dd' or 'qd' required for "
                << LLL_METHOD_STR[method]);
  }

  if (sel_ft == FT_DOUBLE)
    sel_prec = FP_NR<double>::get_prec();
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (sel_ft == FT_LONG_DOUBLE)
    sel_prec = FP_NR<long double>::get_prec();
#endif
#ifdef FPLLL_WITH_DPE
  else if (sel_ft == FT_DPE)
    sel_prec = FP_NR<dpe_t>::get_prec();
#endif
#ifdef FPLLL_WITH_QD
  else if (sel_ft == FT_DD)
    sel_prec = FP_NR<dd_real>::get_prec();
  else if (sel_ft == FT_QD)
    sel_prec = FP_NR<qd_real>::get_prec();
#endif

  if (flags & LLL_VERBOSE)
  {
    cerr << "Starting LLL method '" << LLL_METHOD_STR[method] << "'" << endl
         << "  integer type '" << INT_TYPE_STR[int_type] << "'" << endl
         << "  floating point type '" << FLOAT_TYPE_STR[sel_ft] << "'" << endl;
    if (method != LM_PROVED || int_type != ZT_MPZ || sel_ft == FT_DOUBLE)
    {
      cerr << "  The reduction is not guaranteed";
    }
    else if (sel_prec < good_prec)
    {
      cerr << "  prec < " << good_prec << ", the reduction is not guaranteed";
    }
    else
    {
      cerr << "  prec >= " << good_prec << ", the reduction is guaranteed";
    }
    cerr << endl;
  }

  // Applies the selected method
  int status;
  if (sel_ft == FT_DOUBLE)
  {
    status = lll_reduction_zf<ZT, double>(b, u, u_inv, delta, eta, method, flags);
  }
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (sel_ft == FT_LONG_DOUBLE)
  {
    status = lll_reduction_zf<ZT, long double>(b, u, u_inv, delta, eta, method, flags);
  }
#endif
#ifdef FPLLL_WITH_DPE
  else if (sel_ft == FT_DPE)
  {
    status = lll_reduction_zf<ZT, dpe_t>(b, u, u_inv, delta, eta, method, flags);
  }
#endif
#ifdef FPLLL_WITH_QD
  else if (sel_ft == FT_DD)
  {
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    status = lll_reduction_zf<ZT, dd_real>(b, u, u_inv, delta, eta, method, flags);
    fpu_fix_end(&old_cw);
  }
  else if (sel_ft == FT_QD)
  {
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    status = lll_reduction_zf<ZT, qd_real>(b, u, u_inv, delta, eta, method, flags);
    fpu_fix_end(&old_cw);
  }
#endif
  else if (sel_ft == FT_MPFR)
  {
    int old_prec = FP_NR<mpfr_t>::set_prec(sel_prec);
    status       = lll_reduction_zf<ZT, mpfr_t>(b, u, u_inv, delta, eta, method, flags);
    FP_NR<mpfr_t>::set_prec(old_prec);
  }
  else
  {
    if (0 <= sel_ft && sel_ft <= FT_MPFR)
    {
      // it's a valid choice but we don't have support for it
      FPLLL_ABORT("Compiled without support for LLL reduction with " << FLOAT_TYPE_STR[sel_ft]);
    }
    else
    {
      FPLLL_ABORT("Floating point type " << sel_ft << "not supported in LLL");
    }
  }
  zeros_first(b, u, u_inv);
  return status;
}

// Verify if b is hlll reduced according to delta and eta
// For FT != dpe and FT != mpfr
// This function is not used, but can be used during a testing step.
template <class ZT, class FT>
int is_hlll_reduced_zf(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv, double delta, double eta,
                       double theta)
{
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;

  int householder_flags = HOUSEHOLDER_DEFAULT | HOUSEHOLDER_ROW_EXPO;
  MatHouseholder<Z_NR<ZT>, FP_NR<FT>> m(b, u, u_inv, householder_flags);

  return is_hlll_reduced<Z_NR<ZT>, FP_NR<FT>>(m, delta, eta, theta);
}

// Verify if b is hlll reduced according to delta and eta
// For FT == dpe or FT == mpfr
template <class ZT, class FT>
int is_hlll_reduced_pr(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv, double delta, double eta,
                       double theta)
{
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;

  int householder_flags = HOUSEHOLDER_DEFAULT;
  MatHouseholder<Z_NR<ZT>, FP_NR<FT>> m(b, u, u_inv, householder_flags);

  return is_hlll_reduced<Z_NR<ZT>, FP_NR<FT>>(m, delta, eta, theta);
}

template <class ZT>
int hlll_reduction_wrapper(ZZ_mat<ZT> &, ZZ_mat<ZT> &, ZZ_mat<ZT> &, double, double, double, double,
                           FloatType, int, int)
{
  FPLLL_ABORT("The wrapper method works only with integer type mpz");
  return RED_LLL_FAILURE;
}

template <>
int hlll_reduction_wrapper(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, ZZ_mat<mpz_t> &u_inv, double delta,
                           double eta, double theta, double c, FloatType float_type, int precision,
                           int flags)
{
  FPLLL_CHECK(float_type == FT_DEFAULT,
              "The floating point type cannot be specified with the wrapper method");
  FPLLL_CHECK(precision == 0, "The precision cannot be specified with the wrapper method");
  Wrapper wrapper(b, u, u_inv, delta, eta, theta, c, flags);
  wrapper.hlll();
  zeros_first(b, u, u_inv);
  return wrapper.status;
}

template <class ZT, class FT>
int hlll_reduction_zf(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv, double delta, double eta,
                      double theta, double c, LLLMethod method, int flags)
{
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;
  int householder_flags = HOUSEHOLDER_DEFAULT;
  if (method == LM_FAST)
  {
    householder_flags |= HOUSEHOLDER_ROW_EXPO | HOUSEHOLDER_OP_FORCE_LONG;
    // householder_flags |= HOUSEHOLDER_ROW_EXPO;
  }
  MatHouseholder<Z_NR<ZT>, FP_NR<FT>> m(b, u, u_inv, householder_flags);
  HLLLReduction<Z_NR<ZT>, FP_NR<FT>> hlll_obj(m, delta, eta, theta, c, flags);
  hlll_obj.hlll();

  return hlll_obj.get_status();
}

template <class ZT>
int hlll_reduction_z(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv, double delta, double eta,
                     double theta, double c, LLLMethod method, IntType int_type,
                     FloatType float_type, int precision, int flags, bool nolll)
{
  FPLLL_CHECK(method != LM_HEURISTIC, "HLLL heuristic is not implementated.");

  int sel_prec = 0;
  int status   = -1;
  /* computes the parameters required for the proved version */
  int good_prec = hlll_min_prec(b.get_rows(), b.get_cols(), delta, eta, theta, c);

  // If nolll, just verify if the basis is reduced or not
  /*
   * FIXME: for an unknow reason, this test can use a lot of RAM. Example:
   *   latticegen n 613 2048 q | fplll -a hlll | fplll -a hlll -nolll
   * uses more than 16 GB of RAM. It is advised to use for now the binary
   * isreduced provided inside hplll
   *   latticegen n 613 2048 q | fplll -a hlll | isreduced
   */
  if (nolll)
  {
    sel_prec = (precision != 0) ? precision : good_prec;

    if (flags & LLL_VERBOSE)
    {
      cerr << "Starting HLLL method 'verification'" << endl
           << "  integer type '" << INT_TYPE_STR[int_type] << "'" << endl
           << "  floating point type 'mpfr_t'" << endl;

      if (sel_prec < good_prec)
        cerr << "  prec < " << good_prec << ", the verification is not guaranteed";
      else
        cerr << "  prec >= " << good_prec << ", the verification is guaranteed";

      cerr << endl;
    }

    int old_prec = FP_NR<mpfr_t>::set_prec(sel_prec);

    status = is_hlll_reduced_pr<ZT, mpfr_t>(b, u, u_inv, delta, eta, theta);

    if (flags & LLL_VERBOSE)
    {
      if (status == RED_SUCCESS)
        cerr << "Basis is reduced";
      else
        cerr << "Basis is not reduced";
      cerr << endl;
    }

    FP_NR<mpfr_t>::set_prec(old_prec);

    return status;
  }

  /* switch to wrapper */
  if (method == LM_WRAPPER)
    return hlll_reduction_wrapper(b, u, u_inv, delta, eta, theta, c, float_type, precision, flags);

  /* sets the parameters and checks the consistency */
  if (method == LM_PROVED)
  {
    sel_prec = (precision != 0) ? precision : good_prec;
  }
  else
  {
    sel_prec = (precision != 0) ? precision : PREC_DOUBLE;
  }

  FloatType sel_ft = float_type;

  /* if manually input precision */
  if (precision != 0)
  {
    if (sel_ft == FT_DEFAULT)
    {
      sel_ft = FT_MPFR;
    }
    FPLLL_CHECK(sel_ft == FT_MPFR,
                "The floating type must be mpfr when the precision is specified");
  }

  if (sel_ft == FT_DEFAULT)
  {
    if (method == LM_FAST)
      sel_ft = FT_DOUBLE;
#ifdef FPLLL_WITH_DPE
    else if (sel_prec <= static_cast<int>(FP_NR<dpe_t>::get_prec()))
      sel_ft = FT_DPE;
#endif
#ifdef FPLLL_WITH_QD
    else if (sel_prec <= static_cast<int>(FP_NR<dd_real>::get_prec()))
      sel_ft = FT_DD;
    else if (sel_prec <= static_cast<int>(FP_NR<qd_real>::get_prec()))
      sel_ft = FT_QD;
#endif
    else
      sel_ft = FT_MPFR;
  }
  else if (method == LM_FAST &&
           (sel_ft != FT_DOUBLE && sel_ft != FT_LONG_DOUBLE && sel_ft != FT_DD && sel_ft != FT_QD))
  {
    FPLLL_ABORT("'double' or 'long double' or 'dd' or 'qd' required for "
                << LLL_METHOD_STR[method]);
  }

  if (sel_ft == FT_DOUBLE)
    sel_prec = FP_NR<double>::get_prec();
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (sel_ft == FT_LONG_DOUBLE)
    sel_prec = FP_NR<long double>::get_prec();
#endif
#ifdef FPLLL_WITH_DPE
  else if (sel_ft == FT_DPE)
    sel_prec = FP_NR<dpe_t>::get_prec();
#endif
#ifdef FPLLL_WITH_QD
  else if (sel_ft == FT_DD)
    sel_prec = FP_NR<dd_real>::get_prec();
  else if (sel_ft == FT_QD)
    sel_prec = FP_NR<qd_real>::get_prec();
#endif

  if (flags & LLL_VERBOSE)
  {
    cerr << "Starting HLLL method '" << LLL_METHOD_STR[method] << "'" << endl
         << "  integer type '" << INT_TYPE_STR[int_type] << "'" << endl
         << "  floating point type '" << FLOAT_TYPE_STR[sel_ft] << "'" << endl;
    if (method != LM_PROVED || int_type != ZT_MPZ || sel_ft == FT_DOUBLE)
    {
      cerr << "  The reduction is not guaranteed";
    }
    else if (sel_prec < good_prec)
    {
      cerr << "  prec < " << good_prec << ", the reduction is not guaranteed";
    }
    else
    {
      cerr << "  prec >= " << good_prec << ", the reduction is guaranteed";
    }
    cerr << endl;
  }

  // Applies the selected method
  if (sel_ft == FT_DOUBLE)
    status = hlll_reduction_zf<ZT, double>(b, u, u_inv, delta, eta, theta, c, method, flags);
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (sel_ft == FT_LONG_DOUBLE)
    status = hlll_reduction_zf<ZT, long double>(b, u, u_inv, delta, eta, theta, c, method, flags);
#endif
#ifdef FPLLL_WITH_DPE
  else if (sel_ft == FT_DPE)
    status = hlll_reduction_zf<ZT, dpe_t>(b, u, u_inv, delta, eta, theta, c, method, flags);
#endif
#ifdef FPLLL_WITH_QD
  else if (sel_ft == FT_DD)
  {
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    status = hlll_reduction_zf<ZT, dd_real>(b, u, u_inv, delta, eta, theta, c, method, flags);
    fpu_fix_end(&old_cw);
  }
  else if (sel_ft == FT_QD)
  {
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    status = hlll_reduction_zf<ZT, qd_real>(b, u, u_inv, delta, eta, theta, c, method, flags);
    fpu_fix_end(&old_cw);
  }
#endif
  else if (sel_ft == FT_MPFR)
  {
    int old_prec = FP_NR<mpfr_t>::set_prec(sel_prec);
    status       = hlll_reduction_zf<ZT, mpfr_t>(b, u, u_inv, delta, eta, theta, c, method, flags);
    FP_NR<mpfr_t>::set_prec(old_prec);
  }
  else
  {
    if (0 <= sel_ft && sel_ft <= FT_MPFR)
    {
      // it's a valid choice but we don't have support for it
      FPLLL_ABORT("Compiled without support for LLL reduction with " << FLOAT_TYPE_STR[sel_ft]);
    }
    else
    {
      FPLLL_ABORT("Floating point type " << sel_ft << "not supported in LLL");
    }
  }
  zeros_first(b, u, u_inv);

  return status;
}

/**
 * We define LLL for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
#define FPLLL_DEFINE_LLL(T, id_t)                                                                  \
  int lll_reduction(ZZ_mat<T> &b, double delta, double eta, LLLMethod method,                      \
                    FloatType float_type, int precision, int flags)                                \
  {                                                                                                \
    ZZ_mat<T> empty_mat; /* Empty u -> transform disabled */                                       \
    return lll_reduction_z<T>(b, empty_mat, empty_mat, delta, eta, method, id_t, float_type,       \
                              precision, flags);                                                   \
  }                                                                                                \
                                                                                                   \
  int lll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, double delta, double eta, LLLMethod method,        \
                    FloatType float_type, int precision, int flags)                                \
  {                                                                                                \
    ZZ_mat<T> empty_mat;                                                                           \
    if (!u.empty())                                                                                \
      u.gen_identity(b.get_rows());                                                                \
    return lll_reduction_z<T>(b, u, empty_mat, delta, eta, method, id_t, float_type, precision,    \
                              flags);                                                              \
  }                                                                                                \
                                                                                                   \
  int lll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, ZZ_mat<T> &u_inv, double delta, double eta,        \
                    LLLMethod method, FloatType float_type, int precision, int flags)              \
  {                                                                                                \
    if (!u.empty())                                                                                \
      u.gen_identity(b.get_rows());                                                                \
    if (!u_inv.empty())                                                                            \
      u_inv.gen_identity(b.get_rows());                                                            \
    u_inv.transpose();                                                                             \
    int status =                                                                                   \
        lll_reduction_z<T>(b, u, u_inv, delta, eta, method, id_t, float_type, precision, flags);   \
    u_inv.transpose();                                                                             \
    return status;                                                                                 \
  }

FPLLL_DEFINE_LLL(mpz_t, ZT_MPZ)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DEFINE_LLL(long, ZT_LONG)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DEFINE_LLL(double, ZT_DOUBLE)
#endif

// HLLL

/**
 * We define HLLL for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
#define FPLLL_DEFINE_HLLL(T, id_t)                                                                 \
  int hlll_reduction(ZZ_mat<T> &b, double delta, double eta, double theta, double c,               \
                     LLLMethod method, FloatType float_type, int precision, int flags, bool nolll) \
  {                                                                                                \
    ZZ_mat<T> empty_mat; /* Empty u -> transform disabled */                                       \
    return hlll_reduction_z<T>(b, empty_mat, empty_mat, delta, eta, theta, c, method, id_t,        \
                               float_type, precision, flags, nolll);                               \
  }                                                                                                \
                                                                                                   \
  int hlll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, double delta, double eta, double theta, double c, \
                     LLLMethod method, FloatType float_type, int precision, int flags, bool nolll) \
  {                                                                                                \
    ZZ_mat<T> empty_mat;                                                                           \
    if (!u.empty())                                                                                \
      u.gen_identity(b.get_rows());                                                                \
    return hlll_reduction_z<T>(b, u, empty_mat, delta, eta, theta, c, method, id_t, float_type,    \
                               precision, flags, nolll);                                           \
  }                                                                                                \
                                                                                                   \
  int hlll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, ZZ_mat<T> &u_inv, double delta, double eta,       \
                     double theta, double c, LLLMethod method, FloatType float_type,               \
                     int precision, int flags, bool nolll)                                         \
  {                                                                                                \
    if (!u.empty())                                                                                \
      u.gen_identity(b.get_rows());                                                                \
    if (!u_inv.empty())                                                                            \
      u_inv.gen_identity(b.get_rows());                                                            \
    u_inv.transpose();                                                                             \
    int status = hlll_reduction_z<T>(b, u, u_inv, delta, eta, theta, c, method, id_t, float_type,  \
                                     precision, flags, nolll);                                     \
    u_inv.transpose();                                                                             \
    return status;                                                                                 \
  }

FPLLL_DEFINE_HLLL(mpz_t, ZT_MPZ)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DEFINE_HLLL(long, ZT_LONG)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DEFINE_HLLL(double, ZT_DOUBLE)
#endif

FPLLL_END_NAMESPACE
