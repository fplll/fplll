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

#include "fplll.h"
#include "wrapper.h"
#include "svpcvp.cpp"

FPLLL_BEGIN_NAMESPACE

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
int lll_reduction_wrapper(ZZ_mat<ZT> &b, ZZ_mat<ZT> &u, ZZ_mat<ZT> &u_inv, double delta, double eta,
                          FloatType float_type, int precision, int flags)
{
  FPLLL_ABORT("The wrapper method works only with integer type mpz");
  return RED_LLL_FAILURE;
}

template <>
int lll_reduction_wrapper(IntMatrix &b, IntMatrix &u, IntMatrix &u_inv, double delta, double eta,
                          FloatType float_type, int precision, int flags)
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
                  LLLMethod method, IntType int_type, FloatType float_type, int precision, int flags)
{

  /* switch to wrapper */
  if (method == LM_WRAPPER)
    return lll_reduction_wrapper(b, u, u_inv, delta, eta, float_type, precision, flags);

  FPLLL_CHECK(!(method == LM_PROVED && (flags & LLL_EARLY_RED)),
              "LLL method 'proved' with early reduction is not implemented");

  /* computes the parameters required for the proved version */
  int good_prec = l2MinPrec(b.get_rows(), delta, eta, LLL_DEF_EPSILON);

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
    FPLLL_CHECK(sel_ft == FT_MPFR, "The floating type must be mpfr when the precision is specified");
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
    status = lll_reduction_zf<ZT, mpfr_t>(b, u, u_inv, delta, eta, method, flags);
    FP_NR<mpfr_t>::set_prec(old_prec);
  }
  else
  {
    FPLLL_ABORT("Compiled without support for LLL reduction with " << FLOAT_TYPE_STR[sel_ft]);
  }
  zeros_first(b, u, u_inv);
  return status;
}

/**
 * We define LLL for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
#define FPLLL_DEFINE_LLL(T, id_t)                                                                   \
  int lll_reduction(ZZ_mat<T> &b, double delta, double eta, LLLMethod method, FloatType float_type,  \
                   int precision, int flags)                                                       \
  {                                                                                                \
    ZZ_mat<T> empty_mat; /* Empty u -> transform disabled */                                        \
    return lll_reduction_z<T>(b, empty_mat, empty_mat, delta, eta, method, id_t, float_type, precision,  \
                              flags);                                   \
  }                                                                                                \
                                                                                                   \
  int lll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, double delta, double eta, LLLMethod method,         \
                   FloatType float_type, int precision, int flags)                                  \
  {                                                                                                \
    ZZ_mat<T> empty_mat;                                                                            \
    if (!u.empty())                                                                                \
      u.gen_identity(b.get_rows());                                                                 \
    return lll_reduction_z<T>(b, u, empty_mat, delta, eta, method, id_t, float_type, precision, flags); \
  }                                                                                                \
                                                                                                   \
  int lll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, ZZ_mat<T> &u_inv, double delta, double eta,          \
                   LLLMethod method, FloatType float_type, int precision, int flags)                \
  {                                                                                                \
    if (!u.empty())                                                                                \
      u.gen_identity(b.get_rows());                                                                 \
    if (!u_inv.empty())                                                                             \
      u_inv.gen_identity(b.get_rows());                                                              \
    u_inv.transpose();                                                                              \
    int status =                                                                                   \
        lll_reduction_z<T>(b, u, u_inv, delta, eta, method, id_t, float_type, precision, flags);        \
    u_inv.transpose();                                                                              \
    return status;                                                                                 \
  }

FPLLL_DEFINE_LLL(mpz_t, ZT_MPZ)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DEFINE_LLL(long, ZT_LONG)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DEFINE_LLL(double, ZT_DOUBLE)
#endif

/***************************************************************/

/**
 * call LLLReduction() and then BKZReduction.
 */
template <class FT>
int bkz_reduction_f(IntMatrix &b, const BKZParam &param, int sel_ft, double lll_delta, IntMatrix &u,
                    IntMatrix &u_inv)
{
  int gso_flags = 0;
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;
  if (sel_ft == FT_DOUBLE || sel_ft == FT_LONG_DOUBLE)
    gso_flags |= GSO_ROW_EXPO;
  MatGSO<Integer, FT> m_gso(b, u, u_inv, gso_flags);
  LLLReduction<Integer, FT> lll_obj(m_gso, lll_delta, LLL_DEF_ETA, LLL_DEFAULT);
  BKZReduction<FT> bkz_obj(m_gso, lll_obj, param);
  bkz_obj.bkz();
  return bkz_obj.status;
}

/**
 * interface called from call_bkz() from main.cpp.
 */
int bkz_reduction(IntMatrix *B, IntMatrix *U, const BKZParam &param, FloatType float_type,
                 int precision)
{
  IntMatrix empty_mat;
  IntMatrix &u    = U ? *U : empty_mat;
  IntMatrix &u_inv = empty_mat;
  FPLLL_CHECK(B, "B == NULL in bkzReduction");

  if (U && (!u.empty()))
  {
    u.gen_identity(B->get_rows());
  }

  double lll_delta = param.delta < 1 ? param.delta : LLL_DEF_DELTA;

  FloatType sel_ft = (float_type != FT_DEFAULT) ? float_type : FT_DOUBLE;
  FPLLL_CHECK(!(sel_ft == FT_MPFR && precision == 0),
              "Missing precision for BKZ with floating point type mpfr");

  /* lllwrapper (no FloatType needed, -m ignored) */
  if (param.flags & BKZ_NO_LLL)
    zeros_last(*B, u, u_inv);
  else
  {
    Wrapper wrapper(*B, u, u_inv, lll_delta, LLL_DEF_ETA, LLL_DEFAULT);
    if (!wrapper.lll())
      return wrapper.status;
  }

  /* bkz (with float_type) */
  int status;
  if (sel_ft == FT_DOUBLE)
  {
    status = bkz_reduction_f<FP_NR<double>>(*B, param, sel_ft, lll_delta, u, u_inv);
  }
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (sel_ft == FT_LONG_DOUBLE)
  {
    status = bkz_reduction_f<FP_NR<long double>>(*B, param, sel_ft, lll_delta, u, u_inv);
  }
#endif
#ifdef FPLLL_WITH_DPE
  else if (sel_ft == FT_DPE)
  {
    status = bkz_reduction_f<FP_NR<dpe_t>>(*B, param, sel_ft, lll_delta, u, u_inv);
  }
#endif
#ifdef FPLLL_WITH_QD
  else if (sel_ft == FT_DD)
  {
    status = bkz_reduction_f<FP_NR<dd_real>>(*B, param, sel_ft, lll_delta, u, u_inv);
  }
  else if (sel_ft == FT_QD)
  {
    status = bkz_reduction_f<FP_NR<qd_real>>(*B, param, sel_ft, lll_delta, u, u_inv);
  }
#endif
  else if (sel_ft == FT_MPFR)
  {
    int old_prec = FP_NR<mpfr_t>::set_prec(precision);
    status = bkz_reduction_f<FP_NR<mpfr_t>>(*B, param, sel_ft, lll_delta, u, u_inv);
    FP_NR<mpfr_t>::set_prec(old_prec);
  }
  else
  {
    FPLLL_ABORT("Compiled without support for BKZ reduction with " << FLOAT_TYPE_STR[sel_ft]);
  }
  zeros_first(*B, u, u_inv);
  return status;
}

/**
 * We define BKZ/HKZ for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
int bkz_reduction(IntMatrix &b, int block_size, int flags, FloatType float_type, int precision)
{
  vector<Strategy> strategies;
  BKZParam param(block_size, strategies);
  param.flags = flags;
  return bkz_reduction(&b, NULL, param, float_type, precision);
}

int bkz_reduction(IntMatrix &b, IntMatrix &u, int block_size, int flags, FloatType float_type,
                  int precision)
{
  vector<Strategy> strategies;
  BKZParam param(block_size, strategies);
  param.flags = flags;
  return bkz_reduction(&b, &u, param, float_type, precision);
}

int hkz_reduction(IntMatrix &b, int flags, FloatType float_type, int precision)
{
  vector<Strategy> strategies;
  BKZParam param(b.get_rows(), strategies);
  param.block_size = b.get_rows();
  param.delta      = 1;
  if (flags & HKZ_VERBOSE)
    param.flags |= BKZ_VERBOSE;
  return bkz_reduction(&b, NULL, param, float_type, precision);
}

const char *get_red_status_str(int status)
{
  if (status >= 0 && status < RED_STATUS_MAX)
    return RED_STATUS_STR[status];
  else
    return "unknown error";
}

FPLLL_END_NAMESPACE
