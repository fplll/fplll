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
#include "fplll.h"

#include "gso.cpp"
#include "lll.cpp"
#include "wrapper.cpp"
#include "svpcvp.cpp"
#include "bkz.cpp"
#include "enum/enumerate.cpp"
#include "pruner.cpp"

FPLLL_BEGIN_NAMESPACE

template<class ZT>
void zerosFirst(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInvT) {
  int i, d = b.getRows();
  for (i = d; i > 0 && b[i - 1].is_zero(); i--) {}
  if (i > 0 && i < d) {
    b.rotate(0, i, d - 1);
    if (!u.empty())
      u.rotate(0, i, d - 1);
    if (!uInvT.empty())
      uInvT.rotate(0, i, d - 1);
  }
}

template<class ZT>
void zerosLast(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInvT) {
  int i, d = b.getRows();
  for (i = 0; i < d && b[i].is_zero(); i++) {}
  if (i > 0 && i < d) {
    b.rotate(0, i, d - 1);
    if (!u.empty())
      u.rotate(0, i, d - 1);
    if (!uInvT.empty())
      uInvT.rotate(0, i, d - 1);
  }
}


/**
 * LLL with a typical method "proved or heuristic or fast".
 * @ proved:       exact gram +   exact rowexp +   exact rowaddmul
 * @ heuristic:  approx. gram +   exact rowexp +   exact rowaddmul
 * @ fast:       approx. gram + approx. rowexp + approx. rowaddmul
 *    (double, long double, dd_real, qd_real)
 */
template<class ZT, class FT>
int lllReductionZF(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInv,
                   double delta, double eta, LLLMethod method, int flags) {
  int gsoFlags = 0;
  if (b.getRows() == 0 || b.getCols() == 0) return RED_SUCCESS;
  if (method == LM_PROVED) gsoFlags |= GSO_INT_GRAM;
  if (method == LM_FAST)   gsoFlags |= GSO_ROW_EXPO | GSO_OP_FORCE_LONG;
  MatGSO<Z_NR<ZT>, FP_NR<FT> > mGSO(b, u, uInv, gsoFlags);
  LLLReduction<Z_NR<ZT>, FP_NR<FT> > lllObj(mGSO, delta, eta, flags);
  lllObj.lll();
  return lllObj.status;
}

template<class ZT>
int lllReductionWrapper(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInv,
        double delta, double eta, FloatType floatType,
        int precision, int flags) {
  FPLLL_ABORT("The wrapper method works only with integer type mpz");
  return RED_LLL_FAILURE;
}

template<>
int lllReductionWrapper(IntMatrix& b, IntMatrix& u, IntMatrix& uInv,
        double delta, double eta, FloatType floatType,
        int precision, int flags) {
  FPLLL_CHECK(floatType == FT_DEFAULT,
      "The floating point type cannot be specified with the wrapper method");
  FPLLL_CHECK(precision == 0,
      "The precision cannot be specified with the wrapper method");
  Wrapper wrapper(b, u, uInv, delta, eta, flags);
  wrapper.lll();
  zerosFirst(b, u, uInv);
  return wrapper.status;
}


/**
 * Main function called from call_lll().
 */
template<class ZT>
int lllReductionZ(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInv,
                  double delta, double eta, LLLMethod method,
                  IntType intType, FloatType floatType,
                  int precision, int flags)
{

  /* switch to wrapper */
  if (method == LM_WRAPPER)
    return lllReductionWrapper(b, u, uInv, delta, eta,
            floatType, precision, flags);

  FPLLL_CHECK(!(method == LM_PROVED && (flags & LLL_EARLY_RED)),
    "LLL method 'proved' with early reduction is not implemented");

  /* computes the parameters required for the proved version */
  int goodPrec = l2MinPrec(b.getRows(), delta, eta, LLL_DEF_EPSILON);

  /* sets the parameters and checks the consistency */
  int selPrec = 0;
  if (method == LM_PROVED) {
    selPrec = (precision != 0) ? precision : goodPrec;
  }
  else {
    selPrec = (precision != 0) ? precision : PREC_DOUBLE;
  }

  FloatType selFT = floatType;

  /* if manually input precision */
  if (precision != 0) {
    if (selFT == FT_DEFAULT) {
      selFT = FT_MPFR;
    }
    FPLLL_CHECK(selFT == FT_MPFR,
            "The floating type must be mpfr when the precision is specified");
  }

  if (selFT == FT_DEFAULT) {
    if (method == LM_FAST)
      selFT = FT_DOUBLE;
#ifdef FPLLL_WITH_DPE
    else if (selPrec <= static_cast<int>(FP_NR<dpe_t>::getprec()))
      selFT = FT_DPE;
#endif
#ifdef FPLLL_WITH_QD
    else if (selPrec <= static_cast<int>(FP_NR<dd_real>::getprec()))
      selFT = FT_DD;
    else if (selPrec <= static_cast<int>(FP_NR<qd_real>::getprec()))
      selFT = FT_QD;
#endif
    else
      selFT = FT_MPFR;
  }
  else if (method == LM_FAST
           && (selFT != FT_DOUBLE && selFT != FT_LONG_DOUBLE
               && selFT != FT_DD && selFT != FT_QD)) {
    FPLLL_ABORT("'double' or 'long double' or 'dd' or 'qd' required for "
                  << LLL_METHOD_STR[method]);
  }

  if (selFT == FT_DOUBLE)
    selPrec = FP_NR<double>::getprec();
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (selFT == FT_LONG_DOUBLE)
    selPrec = FP_NR<long double>::getprec();
#endif
#ifdef FPLLL_WITH_DPE
  else if (selFT == FT_DPE)
    selPrec = FP_NR<dpe_t>::getprec();
#endif
#ifdef FPLLL_WITH_QD
  else if (selFT == FT_DD)
    selPrec = FP_NR<dd_real>::getprec();
  else if (selFT == FT_QD)
    selPrec = FP_NR<qd_real>::getprec();
#endif


  if (flags & LLL_VERBOSE) {
    cerr << "Starting LLL method '" << LLL_METHOD_STR[method] << "'" << endl
         << "  integer type '" << INT_TYPE_STR[intType] << "'" << endl
         << "  floating point type '" << FLOAT_TYPE_STR[selFT] << "'" << endl;
    if (method != LM_PROVED || intType != ZT_MPZ || selFT == FT_DOUBLE)
      cerr << "  The reduction is not guaranteed";
    else if (selPrec < goodPrec)
      cerr << "  prec < " << goodPrec << ", the reduction is not guaranteed";
    else {
      cerr << "  prec >= " << goodPrec << ", the reduction is guaranteed";
    }
    cerr << endl;
  }

  // Applies the selected method
  int status;
  if (selFT == FT_DOUBLE) {
    status = lllReductionZF<ZT, double>(b, u, uInv,
            delta, eta, method, flags);
  }
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (selFT == FT_LONG_DOUBLE) {
    status = lllReductionZF<ZT, long double>(b, u, uInv,
            delta, eta, method, flags);
  }
#endif
#ifdef FPLLL_WITH_DPE
  else if (selFT == FT_DPE) {
    status = lllReductionZF<ZT, dpe_t>(b, u, uInv,
            delta, eta, method, flags);
  }
#endif
#ifdef FPLLL_WITH_QD
  else if (selFT == FT_DD) {
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    status = lllReductionZF<ZT, dd_real>(b, u, uInv,
                                         delta, eta, method, flags);
    fpu_fix_end(&old_cw);
  }
  else if (selFT == FT_QD) {
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    status = lllReductionZF<ZT, qd_real>(b, u, uInv,
                                         delta, eta, method, flags);
    fpu_fix_end(&old_cw);
  }
#endif
  else if (selFT == FT_MPFR) {
    int oldPrec = FP_NR<mpfr_t>::setprec(selPrec);
    status = lllReductionZF<ZT, mpfr_t>(b, u, uInv,
            delta, eta, method, flags);
    FP_NR<mpfr_t>::setprec(oldPrec);
  }
  else {
    FPLLL_ABORT("Compiled without support for LLL reduction with " 
                 << FLOAT_TYPE_STR[selFT]);
  }
  zerosFirst(b, u, uInv);
  return status;
}


/**
 * We define LLL for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
#define FPLLL_DEFINE_LLL(T, idT)                                        \
int lllReduction(ZZ_mat<T>& b, double delta, double eta,                \
                 LLLMethod method, FloatType floatType,                 \
                 int precision, int flags) {                            \
  ZZ_mat<T> emptyMat; /* Empty u -> transform disabled */               \
  return lllReductionZ<T>(b, emptyMat, emptyMat, delta, eta, method,    \
            idT, floatType, precision, flags);                          \
}                                                                       \
                                                                        \
int lllReduction(ZZ_mat<T>& b, ZZ_mat<T>& u, double delta, double eta,  \
                 LLLMethod method, FloatType floatType,                 \
                 int precision, int flags) {                            \
  ZZ_mat<T> emptyMat;                                                   \
  if (!u.empty()) u.gen_identity(b.getRows());                           \
  return lllReductionZ<T>(b, u, emptyMat, delta, eta, method,           \
            idT, floatType, precision, flags);                          \
}                                                                       \
                                                                        \
int lllReduction(ZZ_mat<T>& b, ZZ_mat<T>& u, ZZ_mat<T>& uInv,           \
                 double delta, double eta,                              \
                 LLLMethod method, FloatType floatType,                 \
                 int precision, int flags) {                            \
  if (!u.empty()) u.gen_identity(b.getRows());                           \
  if (!uInv.empty()) uInv.gen_identity(b.getRows());                     \
  uInv.transpose();                                                     \
  int status = lllReductionZ<T>(b, u, uInv, delta, eta, method,         \
            idT, floatType, precision, flags);                          \
  uInv.transpose();                                                     \
  return status;                                                        \
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
template<class FT>
int bkzReductionF ( IntMatrix &b,
                    const BKZParam& param,
                    int selFT,
                    double lllDelta,
                    IntMatrix& u,
                    IntMatrix& uInv ) 
{
  int gsoFlags = 0;
  if (b.getRows() == 0 || b.getCols() == 0)
    return RED_SUCCESS;
  if (selFT == FT_DOUBLE || selFT == FT_LONG_DOUBLE)
    gsoFlags |= GSO_ROW_EXPO;
  MatGSO<Integer, FT> mGSO(b, u, uInv, gsoFlags);
  LLLReduction<Integer, FT> lllObj(mGSO, lllDelta, LLL_DEF_ETA, LLL_DEFAULT);
  BKZReduction<FT> bkzObj(mGSO, lllObj, param);
  bkzObj.bkz();
  return bkzObj.status;
}


/**
 * interface called from call_bkz() from main.cpp.
 */
int bkzReduction(IntMatrix* B, IntMatrix *U, const BKZParam& param, FloatType floatType, int precision) {
  IntMatrix emptyMat;
  IntMatrix& u = U ? *U : emptyMat;
  IntMatrix& uInv = emptyMat;
  FPLLL_CHECK(B, "B == NULL in bkzReduction");

  if (U && (!u.empty())) {
    u.gen_identity(B->getRows());
  }

  double lllDelta = param.delta < 1 ? param.delta : LLL_DEF_DELTA;

  FloatType selFT = (floatType != FT_DEFAULT) ? floatType : FT_DOUBLE;
  FPLLL_CHECK(!(selFT == FT_MPFR && precision == 0),
               "Missing precision for BKZ with floating point type mpfr");

  /* lllwrapper (no FloatType needed, -m ignored) */
  if (param.flags & BKZ_NO_LLL)
    zerosLast(*B, u, uInv);
  else {
    Wrapper wrapper(*B, u, uInv, lllDelta, LLL_DEF_ETA, LLL_DEFAULT);
    if (!wrapper.lll())
      return wrapper.status;
  }

  
  /* bkz (with FloatType) */
  int status;
  if (selFT == FT_DOUBLE) {
    status = bkzReductionF< FP_NR<double> >(*B, param, selFT,
                                            lllDelta, u, uInv);
  }
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (selFT == FT_LONG_DOUBLE) {
    status = bkzReductionF< FP_NR<long double> >(*B, param, selFT,
                                                 lllDelta, u, uInv);
  }
#endif
#ifdef FPLLL_WITH_DPE
  else if (selFT == FT_DPE) {
    status = bkzReductionF< FP_NR<dpe_t> >(*B, param, selFT,
                                           lllDelta, u, uInv);
  }
#endif
#ifdef FPLLL_WITH_QD
  else if (selFT == FT_DD) {
    status = bkzReductionF< FP_NR<dd_real> >(*B, param, selFT,
                                             lllDelta, u, uInv);
  }
  else if (selFT == FT_QD) {
    status = bkzReductionF< FP_NR<qd_real> >(*B, param, selFT,
                                             lllDelta, u, uInv);
  }
#endif
  else if (selFT == FT_MPFR) {
    int oldPrec = FP_NR<mpfr_t>::setprec(precision);
    status = bkzReductionF< FP_NR<mpfr_t> >(*B, param, selFT,
                                            lllDelta, u, uInv);
    FP_NR<mpfr_t>::setprec(oldPrec);
  }
  else {
    FPLLL_ABORT("Compiled without support for BKZ reduction with " 
                 << FLOAT_TYPE_STR[selFT]);
  }
  zerosFirst(*B, u, uInv);
  return status;
}


/**
 * We define BKZ/HKZ for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
int bkzReduction(IntMatrix& b, int blockSize, int flags, FloatType floatType, int precision) {
  BKZParam param;
  param.blockSize = blockSize;
  param.flags = flags;
  return bkzReduction(&b, NULL, param, floatType, precision);
}

int bkzReduction(IntMatrix& b, IntMatrix& u, int blockSize, int flags, FloatType floatType, int precision) {
  BKZParam param;
  param.blockSize = blockSize;
  param.flags = flags;
  return bkzReduction(&b, &u, param, floatType, precision);
}

int hkzReduction(IntMatrix& b, int flags) {
  BKZParam param;
  param.blockSize = b.getRows();
  param.delta = 1;
  if (flags & HKZ_VERBOSE)
    param.flags |= BKZ_VERBOSE;
  return bkzReduction(&b, NULL, param);
}

const char* getRedStatusStr(int status) {
  if (status >= 0 && status < RED_STATUS_MAX)
    return RED_STATUS_STR[status];
  else
    return "unknown error";
}

/** instantiate functions **/

template bool isLLLReduced<Z_NR<mpz_t>, FP_NR<double> >(MatGSO<Z_NR<mpz_t>, FP_NR<double> >& m, double delta, double eta);
#ifdef FPLLL_WITH_LONG_DOUBLE
template bool isLLLReduced<Z_NR<mpz_t>, FP_NR<long double> >(MatGSO<Z_NR<mpz_t>, FP_NR<long double> >& m, double delta, double eta);
#endif
template bool isLLLReduced<Z_NR<mpz_t>, FP_NR<mpfr_t> >(MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t> >& m, double delta, double eta);

#ifdef FPLLL_WITH_QD
template bool isLLLReduced<Z_NR<mpz_t>, FP_NR<dd_real> >(MatGSO<Z_NR<mpz_t>, FP_NR<dd_real> >& m, double delta, double eta);
template bool isLLLReduced<Z_NR<mpz_t>, FP_NR<qd_real> >(MatGSO<Z_NR<mpz_t>, FP_NR<qd_real> >& m, double delta, double eta);
#endif

#ifdef FPLLL_WITH_DPE
template bool isLLLReduced<Z_NR<mpz_t>, FP_NR<dpe_t> >(MatGSO<Z_NR<mpz_t>, FP_NR<dpe_t> >& m, double delta, double eta);
#endif

/*
template void Enumeration::enumerate<FP_NR<double> >(MatGSO<Integer, FP_NR<double> >& gso,
                                                     FP_NR<double>& fMaxDist, long maxDistExpo,
                                                     Evaluator<FP_NR<double> >& evaluator,
                                                     const vector<FP_NR<double> >& targetCoord,
                                                     const vector<FP_NR<double> >& subTree,
                                                     int first, int last, const vector<double>& pruning, bool dual);
*/
FPLLL_END_NAMESPACE
