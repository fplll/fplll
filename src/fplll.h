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

#ifndef FPLLL_H
#define FPLLL_H

#include "util.h"
#include "svpcvp.h"

FPLLL_BEGIN_NAMESPACE

#define FPLLL_DECLARE_LLL(T)                                                \
int lllReduction(ZZ_mat<T>& b,                                              \
        double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,             \
        LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,    \
        int precision = 0, int flags = LLL_DEFAULT);                        \
                                                                            \
int lllReduction(ZZ_mat<T>& b, ZZ_mat<T>& u,                                \
        double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,             \
        LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,    \
        int precision = 0, int flags = LLL_DEFAULT);                        \
                                                                            \
int lllReduction(ZZ_mat<T>& b, ZZ_mat<T>& u, ZZ_mat<T>& uInv,               \
        double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,             \
        LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,    \
        int precision = 0, int flags = LLL_DEFAULT);

FPLLL_DECLARE_LLL(mpz_t)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DECLARE_LLL(long)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DECLARE_LLL(double)
#endif

struct BKZParam {
  BKZParam() : b(NULL), u(NULL), blockSize(0), delta(LLL_DEF_DELTA),
    floatType(FT_DEFAULT), precision(0), flags(BKZ_DEFAULT),
    maxLoops(0), maxTime(0) {
  }
  IntMatrix* b;
  IntMatrix* u;
  int blockSize;
  double delta;
  FloatType floatType;
  int precision;
  int flags;
  int maxLoops;
  double maxTime;
  vector<double> pruning;
};

int bkzReduction(const BKZParam& param);
int bkzReduction(IntMatrix& b, int blockSize, int flags = BKZ_DEFAULT);
int bkzReduction(IntMatrix& b, IntMatrix& u, int blockSize, int flags = BKZ_DEFAULT);

int hkzReduction(IntMatrix& b, int flags = HKZ_DEFAULT);

/**
 * Returns the string corresponding to an error code of LLL/BKZ.
 */
const char* getRedStatusStr(int status);

FPLLL_END_NAMESPACE

#ifdef FPLLL_V3_COMPAT
#include "fplllv31.h"
#endif

#endif
