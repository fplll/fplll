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

#include "util.h"
#include "lll.h"
#include "wrapper.h"

FPLLL_BEGIN_NAMESPACE

/* prec=53, eta=0.501, dim < dim_double_max [ (delta / 100.0) + 25 ] */
const double dim_double_max[75]=
  {0,26,29.6,28.1,31.1,32.6,34.6,34,37.7,38.8,39.6,41.8,40.9,43.6,44.2,47,46.8,
   50.6,49.1,51.5,52.5,54.8,54.6,57.4,57.6,59.9,61.8,62.3,64.5,67.1,68.8,68.3,
   69.9,73.1,74,76.1,76.8,80.9,81.8,83,85.3,87.9,89,90.1,89,94.6,94.8,98.7,99,
   101.6,104.9,106.8,108.2,107.4,110,112.7,114.6,118.1,119.7,121.8,122.9,126.6,
   128.6,129,133.6,126.9,135.9,139.5,135.2,137.2,139.3,142.8,142.4,142.5,145.4};

const double eta_dep[10]=
  {1.,//0.5
   1.,//0.55
   1.0521,//0.6
   1.1254,//0.65
   1.2535,//0.7
   1.3957,//0.75
   1.6231,//0.8
   1.8189,//0.85
   2.1025,//0.9
   2.5117};//0.95

Wrapper::Wrapper(IntMatrix& b, IntMatrix& u, IntMatrix& uInv,
                 double delta, double eta, int flags) :
  status(RED_SUCCESS), b(b), u(u), uInv(uInv), delta(delta), eta(eta),
  useLong(false), lastEarlyRed(0)
{
  n = b.get_cols();
  d = b.get_rows();
  this->flags = flags;
  maxExponent = b.getMaxExp() + (int) ceil(0.5 * log2((double) d * n));

  // Computes the parameters required for the proved version
  goodPrec = l2MinPrec(d, delta, eta, LLL_DEF_EPSILON);
}

bool Wrapper::little(int kappa, int precision) {
  /*one may add here dimension arguments with respect to eta and delta */
  int dm=(int)(delta*100.-25.);
  if (dm<0) dm=0;
  if (dm>74) dm=74;

  int em=(int) ((eta-0.5)*20);
  if (em<0) em=0;
  if (em>9) em=9;

  double p = max(1.0, precision / 53.0);

  p *= eta_dep[em]; /* eta dependance */
  p *= dim_double_max[dm];
  //cerr << kappa << " compared to " << p << endl;
  return kappa < p;
}


/**
 * main function. Method determines whether heuristic, fast or proved
 */
template<class Z, class F>
int Wrapper::callLLL(ZZ_mat<Z>& bz, ZZ_mat<Z>& uz, ZZ_mat<Z>& uInvZ,
        LLLMethod method, int precision, double delta, double eta) {
  typedef Z_NR<Z> ZT;
  typedef FP_NR<F> FT;

  if (flags & LLL_VERBOSE) {
    cerr << "====== Wrapper: calling " << LLL_METHOD_STR[method] << "<"
         << numTypeStr<Z>() << "," << numTypeStr<F>() << "> method";
    if (precision > 0) {
      cerr << " (precision=" << precision << ")";
    }
    cerr << " ======" << endl;
  }

  int gsoFlags = 0;
  if (method == LM_PROVED) gsoFlags |= GSO_INT_GRAM;
  if (method == LM_FAST) gsoFlags |= GSO_ROW_EXPO;
  if (method != LM_PROVED && precision == 0) gsoFlags |= GSO_OP_FORCE_LONG;

  int oldprec = Float::getprec();
  if (precision > 0) {
    Float::setprec(precision);
  }
  MatGSO<ZT, FT> mGSO(bz, uz, uInvZ, gsoFlags);
  LLLReduction<ZT, FT> lllObj(mGSO, delta, eta, flags);
  lllObj.lastEarlyRed = lastEarlyRed;
  lllObj.lll();
  status = lllObj.status;
  lastEarlyRed = max(lastEarlyRed, lllObj.lastEarlyRed);
  if (precision > 0) {
    Float::setprec(oldprec);
  }

  if (flags & LLL_VERBOSE) {
    cerr << "====== Wrapper: end of " << LLL_METHOD_STR[method]
        << " method ======\n" << endl;
  }

  if (lllObj.status == RED_SUCCESS)
    return 0;
  else if (lllObj.status == RED_GSO_FAILURE
            || lllObj.status == RED_BABAI_FAILURE)
    return lllObj.finalKappa;
  else
    return -1;
}


/**
 * pass the method to callLLL()
 */
template<class F>
int Wrapper::fastLLL(double delta, double eta) {
  return callLLL<mpz_t, F>(b, u, uInv, LM_FAST, 0, delta, eta);
}

template<class Z, class F>
int Wrapper::heuristicLLL(ZZ_mat<Z>& bz, ZZ_mat<Z>& uz,
                          ZZ_mat<Z>& uInvZ, int precision,
                          double delta, double eta) {
  return callLLL<Z, F>(bz, uz, uInvZ, LM_HEURISTIC, precision, delta, eta);
}

template<class Z, class F>
int Wrapper::provedLLL(ZZ_mat<Z>& bz, ZZ_mat<Z>& uz,
                       ZZ_mat<Z>& uInvZ, int precision,
                       double delta, double eta) {
  return callLLL<Z, F>(bz, uz, uInvZ, LM_PROVED, precision, delta, eta);
}


/**
 * In heuristicLoop(), we only use double or dpe_t or mpfr_t.
 */
int Wrapper::heuristicLoop(int precision) {
  int kappa;

  if (precision > numeric_limits<double>::digits)
    kappa = heuristicLLL<mpz_t, mpfr_t>(b, u, uInv, precision, delta, eta);
  else {
#ifdef FPLLL_WITH_DPE
    kappa = heuristicLLL<mpz_t, dpe_t>(b, u, uInv, 0, delta, eta);
#else
    kappa = heuristicLLL<mpz_t, mpfr_t>(b, u, uInv, precision, delta, eta);
#endif
  }

  if (kappa == 0)
    return 0; // Success
  else if (precision < goodPrec && !little(kappa, precision))
    return heuristicLoop(increasePrec(precision));
  else
    return provedLoop(precision);
}


int Wrapper::provedLoop(int precision) {
  int kappa;
#ifdef FPLLL_WITH_QD
  if (precision > PREC_DD)
#else
  if (precision > numeric_limits<double>::digits)
#endif
    kappa = provedLLL<mpz_t, mpfr_t>(b, u, uInv, precision, delta, eta);
  else if (maxExponent * 2 > MAX_EXP_DOUBLE) {
#ifdef FPLLL_WITH_DPE
    kappa = provedLLL<mpz_t, dpe_t>(b, u, uInv, 0, delta, eta);
#else
    kappa = provedLLL<mpz_t, mpfr_t>(b, u, uInv, precision, delta, eta);
#endif
  }
#ifdef FPLLL_WITH_QD
  else if (precision > numeric_limits<double>::digits)
    kappa = provedLLL<mpz_t, dd_real>(b, u, uInv, precision, delta, eta);  
#endif
  else
    kappa = provedLLL<mpz_t, double>(b, u, uInv, 0, delta, eta);

  if (kappa == 0)
    return 0; // Success
  else if (precision < goodPrec)
    return provedLoop(increasePrec(precision));
  else
    return -1; // This point should never be reached
}


/**
 * last call to LLL. Need to be provedLLL.
 */
int Wrapper::lastLLL() {

  /* <long, FT> */
#ifdef FPLLL_WITH_ZLONG
  if (useLong) {
    int kappa;
    if (goodPrec <= numeric_limits<double>::digits)
      kappa = provedLLL<long, double>(bLong, uLong, uInvLong,
                                      goodPrec, delta, eta);
#ifdef FPLLL_WITH_QD
    else if (goodPrec <= PREC_DD)
      kappa = provedLLL<long, dd_real>(bLong, uLong, uInvLong,
                                       goodPrec, delta, eta);
#endif
    else
      kappa = provedLLL<long, mpfr_t>(bLong, uLong, uInvLong,
                                      goodPrec, delta, eta);
    return kappa;
  }
#endif

  /* <mpfr, FT> */
#ifdef FPLLL_WITH_DPE
  if (goodPrec <= numeric_limits<double>::digits)
    return provedLLL<mpz_t, dpe_t>(b, u, uInv, goodPrec, delta, eta);
#ifdef FPLLL_WITH_QD
  else if (goodPrec <= PREC_DD)
    return provedLLL<mpz_t, dd_real>(b, u, uInv, goodPrec, delta, eta);
#endif
#endif
  return provedLLL<mpz_t, mpfr_t>(b, u, uInv, goodPrec, delta, eta);
}


/**
 * Wrapper.lll() calls
 *  - heuristicLLL()
 *  - fastLLL()
 *  - provedLLL()
 */
bool Wrapper::lll() {
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;

#ifdef FPLLL_WITH_ZLONG
  bool heuristicWithLong = maxExponent < numeric_limits<long>::digits - 2
            && u.empty() && uInv.empty();
  bool provedWithLong = 2 * maxExponent < numeric_limits<long>::digits - 2
            && u.empty() && uInv.empty();
#else
  bool heuristicWithLong = false, provedWithLong = false;
#endif

  int kappa;

  /* small matrix */
  if (heuristicWithLong) {
#ifdef FPLLL_WITH_ZLONG
    setUseLong(true);
    /* try heuristicLLL <long, double> */
    heuristicLLL<long, double>(bLong, uLong, uInvLong, 0, delta, eta);
#endif
  }
  /* large matrix */
  else {

    /* try fastLLL<mpz_t, double> */
    kappa = fastLLL<double>(delta, eta);
    bool lllFailure = (kappa != 0);
    int lastPrec;

    /* try fastLLL<mpz_t, long double> */
#ifdef FPLLL_WITH_LONG_DOUBLE
    if (lllFailure) {
      kappa = fastLLL<long double>(delta, eta);
      lllFailure = kappa != 0;
    }
    lastPrec = numeric_limits<long double>::digits;
#else
    lastPrec = numeric_limits<double>::digits;
#endif

    /* try fastLLL<mpz_t, dd_real> */
#ifdef FPLLL_WITH_QD
    if (lllFailure) {
      kappa = fastLLL<dd_real>(delta, eta);
      lllFailure = kappa != 0;
    }
    lastPrec = PREC_DD;
#else
#ifdef FPLLL_WITH_LONG_DOUBLE
    lastPrec = numeric_limits<long double>::digits;
#else
    lastPrec = numeric_limits<double>::digits;
#endif
#endif

    /* loop */
    if (lllFailure) {
      int precD = numeric_limits<double>::digits;
      if (little(kappa, lastPrec))
        kappa = provedLoop(precD);
      else
        kappa = heuristicLoop(increasePrec(precD));
    }
  }

  setUseLong(provedWithLong);
  /* final LLL */
  kappa = lastLLL();
  setUseLong(false);
  return kappa == 0;
}


/**
 * set blong <-- b
 */
void Wrapper::setUseLong(bool value) {
#ifdef FPLLL_WITH_ZLONG
  if (!useLong && value) {
    if (bLong.empty()) {
      bLong.resize(d, n);
    }
    for (int i = 0; i < d; i++)
      for (int j = 0; j < n; j++)
        bLong(i, j) = b(i, j).get_si();
  }
  else if (useLong && !value) {
    for (int i = 0; i < d; i++)
      for (int j = 0; j < n; j++)
        b(i, j) = bLong(i, j).get_si();
  }
  useLong = value;
#endif
}

int Wrapper::increasePrec(int precision) {
  return min(precision * 2, goodPrec);
}

FPLLL_END_NAMESPACE
