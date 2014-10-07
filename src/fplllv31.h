/* Copyright (C) 2011 Xavier Pujol.
   Copyright (C) 2013 Damien Stehle.

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

/* Compatibility layer for version 3.1 */

#ifndef FPLLL_V31_H
#define FPLLL_V31_H

FPLLL_BEGIN_NAMESPACE

template<class T> inline FloatType getFloatType()         {return FT_DEFAULT;}
template<>        inline FloatType getFloatType<double>() {return FT_DOUBLE;}
template<>        inline FloatType getFloatType<dpe_t>()  {return FT_DPE;}
template<>        inline FloatType getFloatType<mpfr_t>() {return FT_MPFR;}

template<class ZT, class FT>
struct lll31method {
  lll31method(LLLMethod method, int flags, ZZ_mat<ZT>* B, int precision,
           double eta, double delta, int siegel) :
    B(B), method(method), flags(flags), precision(precision),
    eta(eta), delta(delta) 
  {
    if (siegel) flags = flags | LLL_SIEGEL;
  }
  int LLL() {
    return lllReduction(*B, delta, eta, method, getFloatType<FT>(),
            precision, flags) == RED_SUCCESS ? 0 : -1;
  }
  ZZ_mat<ZT>* GetBase() {return B;}
  ZZ_mat<ZT>* B;
  LLLMethod method;
  int flags, precision;
  double eta, delta;
};

template<class ZT, class FT>
struct fast : public lll31method<ZT, FT> {
  fast(ZZ_mat<ZT>* B, int precision = 0, double eta = LLL_DEF_ETA,
       double delta = LLL_DEF_DELTA) :
    lll31method<ZT, FT>(LM_FAST, 0, B, precision, eta, delta, 0) {}
};

template<class ZT, class FT>
struct fast_early_red : public lll31method<ZT, FT> {
  fast_early_red(ZZ_mat<ZT>* B, int precision = 0, double eta = LLL_DEF_ETA,
                 double delta = LLL_DEF_DELTA) :
    lll31method<ZT, FT>(LM_FAST, LLL_EARLY_RED, B, precision, eta, delta, 0) {}
};

template<class ZT, class FT>
struct heuristic : public lll31method<ZT, FT> {
  heuristic(ZZ_mat<ZT>*B, int precision = 0, double eta = LLL_DEF_ETA,
            double delta = LLL_DEF_DELTA, int siegel = 0) :
    lll31method<ZT, FT>(LM_HEURISTIC, 0, B, precision,
                      eta, delta, siegel) {}
};

template<class ZT, class FT>
struct heuristic_early_red : public lll31method<ZT, FT> {
  heuristic_early_red(ZZ_mat<ZT>*B, int precision = 0, double eta = LLL_DEF_ETA,
                      double delta = LLL_DEF_DELTA, int siegel = 0) :
    lll31method<ZT, FT>(LM_HEURISTIC, LLL_EARLY_RED, B, precision,
                      eta, delta, siegel) {}
};

template<class ZT, class FT>
struct proved : public lll31method<ZT, FT> {
  proved(ZZ_mat<ZT>*B, int precision = 0, double eta = LLL_DEF_ETA,
            double delta = LLL_DEF_DELTA, int siegel = 0) :
    lll31method<ZT, FT>(LM_PROVED, 0, B, precision,
                      eta, delta, siegel) {}
};

struct wrapper {
  wrapper(ZZ_mat<mpz_t>* B, int /*precision*/ = 0, double eta = LLL_DEF_ETA,
          double delta = LLL_DEF_DELTA) : B(B), eta(eta), delta(delta) {}
  int LLL() {
    int s = lllReduction(*B, delta, eta, LM_WRAPPER);
    return s != RED_SUCCESS ? -1 : 0;
  }
  ZZ_mat<mpz_t>* GetBase() {return B;}
  ZZ_mat<mpz_t>* B;
  double eta, delta;
};

FPLLL_END_NAMESPACE

#endif
