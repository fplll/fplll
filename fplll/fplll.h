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

#if __cplusplus < 201103L
#error fplll needs at least a C++11 compliant compiler
#endif

#include "util.h"
#include "svpcvp.h"
#include "bkz_param.h"
#include "bkz.h"
#include "wrapper.h"
#include "pruner.h"

FPLLL_BEGIN_NAMESPACE

#define FPLLL_DECLARE_LLL(T)                                                \
int lll_reduction(ZZ_mat<T>& b,                                              \
        double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,             \
        LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,    \
        int precision = 0, int flags = LLL_DEFAULT);                        \
                                                                            \
int lll_reduction(ZZ_mat<T>& b, ZZ_mat<T>& u,                                \
        double delta = LLL_DEF_DELTA, double eta = LLL_DEF_ETA,             \
        LLLMethod method = LM_WRAPPER, FloatType floatType = FT_DEFAULT,    \
        int precision = 0, int flags = LLL_DEFAULT);                        \
                                                                            \
int lll_reduction(ZZ_mat<T>& b, ZZ_mat<T>& u, ZZ_mat<T>& uInv,               \
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


int bkz_reduction(IntMatrix* B, IntMatrix* U, const BKZParam& param, FloatType float_type=FT_DEFAULT, int precision=0);
int bkz_reduction(IntMatrix& b, int block_size, int flags = BKZ_DEFAULT, FloatType float_type=FT_DEFAULT, int precision=0);
int bkz_reduction(IntMatrix& b, IntMatrix& u, int block_size, int flags = BKZ_DEFAULT, FloatType float_type=FT_DEFAULT, int precision=0);

int hkz_reduction(IntMatrix& b, int flags = HKZ_DEFAULT, FloatType float_type=FT_DEFAULT, int precision=0);

/**
 * Returns the string corresponding to an error code of LLL/BKZ.
 */
const char* get_red_status_str(int status);

FPLLL_END_NAMESPACE

#ifdef FPLLL_V3_COMPAT
#include "fplllv31.h"
#endif

#endif
