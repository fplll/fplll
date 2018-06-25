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

#ifndef FPLLL_WRAPPER_H
#define FPLLL_WRAPPER_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

/* The matrix b must not be modified before calling lll().
   lll() must be called only once. */

class Wrapper
{
public:
  /* u must be either empty or the identity matrix */
  Wrapper(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, ZZ_mat<mpz_t> &u_inv, double delta, double eta,
          int flags);

  bool lll();

  int status;

private:
  ZZ_mat<mpz_t> &b;
  ZZ_mat<mpz_t> &u;
  ZZ_mat<mpz_t> &u_inv;

#ifdef FPLLL_WITH_ZLONG
  ZZ_mat<long> b_long;
  ZZ_mat<long> u_long;      // Always empty
  ZZ_mat<long> u_inv_long;  // Always empty
#endif

  double delta;
  double eta;
  int good_prec;
  bool use_long;
  int flags;

  bool little(int kappa, int precision);

  template <class Z, class F>
  int call_lll(ZZ_mat<Z> &bz, ZZ_mat<Z> &uz, ZZ_mat<Z> &u_inv_z, LLLMethod method, int precision,
               double delta, double eta);

  template <class F> int fast_lll(double delta, double eta);

  template <class Z, class F>
  int heuristic_lll(ZZ_mat<Z> &bz, ZZ_mat<Z> &uz, ZZ_mat<Z> &u_inv_z, int precision, double delta,
                    double eta);

  template <class Z, class F>
  int proved_lll(ZZ_mat<Z> &bz, ZZ_mat<Z> &uz, ZZ_mat<Z> &u_inv_z, int precision, double delta,
                 double eta);

  int heuristic_loop(int precision);
  int proved_loop(int precision);
  int last_lll();

  void set_use_long(bool value);
  int increase_prec(int precision);

  int max_exponent;
  int n;
  int d;
  int last_early_red;
};

#define FPLLL_DECLARE_LLL(T)                                                                       \
  int lll_reduction(ZZ_mat<T> &b, double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,          \
                    LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,               \
                    int precision = 0, int flags = LLL_DEFAULT);                                   \
                                                                                                   \
  int lll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, double delta = LLL_DEF_DELTA,                      \
                    double eta = LLL_DEF_ETA, LLLMethod method = LM_WRAPPER,                       \
                    FloatType floatType = FT_DEFAULT, int precision = 0, int flags = LLL_DEFAULT); \
                                                                                                   \
  int lll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, ZZ_mat<T> &u_inv, double delta = LLL_DEF_DELTA,    \
                    double eta = LLL_DEF_ETA, LLLMethod method = LM_WRAPPER,                       \
                    FloatType floatType = FT_DEFAULT, int precision = 0, int flags = LLL_DEFAULT);

FPLLL_DECLARE_LLL(mpz_t)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DECLARE_LLL(long)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DECLARE_LLL(double)
#endif

// H-LLL

/**
 * We define H-LLL for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
#define FPLLL_DECLARE_HLLL(T)                                                                      \
  int hlll_reduction(ZZ_mat<T> &b, double delta = LLL_DEF_DELTA, double eta = HLLL_DEF_ETA,        \
                     double theta = HLLL_DEF_THETA, double c = HLLL_DEF_C,                         \
                     LLLMethod method = LM_WRAPPER, FloatType float_type = FT_DEFAULT,             \
                     int precision = 0, int flags = LLL_DEFAULT, bool is_reduced = false);         \
  int hlll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, double delta = LLL_DEF_DELTA,                     \
                     double eta = HLLL_DEF_ETA, double theta = HLLL_DEF_THETA,                     \
                     double c = HLLL_DEF_C, LLLMethod method = LM_WRAPPER,                         \
                     FloatType float_type = FT_DEFAULT, int precision = 0,                         \
                     int flags = LLL_DEFAULT, bool is_reduced = false);                            \
  int hlll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, ZZ_mat<T> &u_inv, double delta = LLL_DEF_DELTA,   \
                     double eta = HLLL_DEF_ETA, double theta = HLLL_DEF_THETA,                     \
                     double c = HLLL_DEF_C, LLLMethod method = LM_WRAPPER,                         \
                     FloatType float_type = FT_DEFAULT, int precision = 0,                         \
                     int flags = LLL_DEFAULT, bool is_reduced = false);

FPLLL_DECLARE_HLLL(mpz_t)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DECLARE_HLLL(long)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DECLARE_HLLL(double)
#endif

FPLLL_END_NAMESPACE

#endif
