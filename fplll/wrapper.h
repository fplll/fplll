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

/**
 * @brief Wrapper. This class provides an externally callable API for LLL reducing some basis ``b``.
 * This class forcibly instantiates some template declarations (see FPLLL_DECLARE_LLL(T) for more
 * information at the bottom of this class) and provides an interface for calling provable,
 * heuristic and fast variants of both LLL and HLLL.
 *
 * User note: all of the parameters passed into this class are non-const references. Thus, this
 * class may overwrite these variables.
 *
 * In particular, the methods in this class typically
 * over-write the parameter ``u`` with the transformation matrix applied to ``b`` to produce LLL(B).
 * In other words, let C = LLL(B). Then we can write C = U B for some unimodular transformation
 * matrix U. For this reason, ``u`` should either be an empty matrix or the identity matrix when it
 * is passed in to this function, so that no information is lost.
 *
 * Similarly, the parameter ``u_inv`` will contain the inverse matrix ``u ^ -1``.
 * This allows you to recover the original basis ``b`` by multiplying ``C`` by ``u_inv``.
 * Note that all operations on this object are carried out on the transpose of this matrix
 * (see gso_interface.h) for speed: hence the discrepancy between the names (see the lll_reduction
 * method in wraper.cpp for this). This and other such behaviour on these parameters is described in
 * some detail in gso_interface.h.
 */
class Wrapper
{
public:
  /* u must be either empty or the identity matrix */
  Wrapper(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, ZZ_mat<mpz_t> &u_inv, double delta, double eta,
          int flags);

  // Used for HLLL
  Wrapper(ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, ZZ_mat<mpz_t> &u_inv, double delta, double eta,
          double theta, double c, int flags);

  bool lll();

  // Call HLLL on the wrapper object
  // TODO: this wrapper does not offers as many options as the one of LLL (for example, use long
  // instead of mpz_t, modify delta, ...)
  bool hlll();

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

  // For HLLL
  double theta;
  double c;

  // High level function, select depending on method and precision the best way to call
  // HLLL (types, flags enabled, ...)
  template <class F> bool call_hlll(LLLMethod method, int precision);

  // = call_hlll(LM_FAST, 0);
  template <class F> bool fast_hlll();

  // = call_hlll(LM_PROVED, precision);
  template <class F> bool proved_hlll(int precision);

  // Perform proved_hlll until grood_prec is reached or the lattice is reduced
  int hlll_proved_loop(int precision);

  // Perform proved version of HLLL with good_prec
  bool last_hlll();
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

// HLLL

/**
 * We define HLLL for each input type instead of using a template,
 * in order to force the compiler to instantiate the functions.
 */
#define FPLLL_DECLARE_HLLL(T)                                                                      \
  int hlll_reduction(ZZ_mat<T> &b, double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,         \
                     double theta = HLLL_DEF_THETA, double c = HLLL_DEF_C,                         \
                     LLLMethod method = LM_WRAPPER, FloatType float_type = FT_DEFAULT,             \
                     int precision = 0, int flags = LLL_DEFAULT, bool nolll = false);              \
  int hlll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, double delta = LLL_DEF_DELTA,                     \
                     double eta = LLL_DEF_ETA, double theta = HLLL_DEF_THETA,                      \
                     double c = HLLL_DEF_C, LLLMethod method = LM_WRAPPER,                         \
                     FloatType float_type = FT_DEFAULT, int precision = 0,                         \
                     int flags = LLL_DEFAULT, bool nolll = false);                                 \
  int hlll_reduction(ZZ_mat<T> &b, ZZ_mat<T> &u, ZZ_mat<T> &u_inv, double delta = LLL_DEF_DELTA,   \
                     double eta = LLL_DEF_ETA, double theta = HLLL_DEF_THETA,                      \
                     double c = HLLL_DEF_C, LLLMethod method = LM_WRAPPER,                         \
                     FloatType float_type = FT_DEFAULT, int precision = 0,                         \
                     int flags = LLL_DEFAULT, bool nolll = false);

FPLLL_DECLARE_HLLL(mpz_t)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DECLARE_HLLL(long)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DECLARE_HLLL(double)
#endif

FPLLL_END_NAMESPACE

#endif
