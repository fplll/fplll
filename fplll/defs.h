/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.
   Copyright (C) 2013 Damien Stehle

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

#ifndef FPLLL_DEFS_H
#define FPLLL_DEFS_H

#ifndef __CYGWIN__
#define FPLLL_WITH_LONG_DOUBLE
#endif

#define FPLLL_WITH_DPE
#define FPLLL_WITH_ZDOUBLE
#define FPLLL_WITH_ZLONG
#define FPLLL_WITH_GETRUSAGE

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef FPLLL_WITH_GETRUSAGE
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include "fplll_config.h"
#ifdef HAVE_LIBMPIR
#include <mpir.h>
#endif
#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif
#include <mpfr.h>
#ifdef FPLLL_WITH_DPE
#include "nr/dpe.h"
#endif

#if defined(__sun) || defined(__CYGWIN__)
#include <ieeefp.h>
extern "C" long double ldexpl(long double x, int exp);
#ifndef NAN
#define NAN __builtin_nanf("")
#endif
#endif

#define FPLLL_INFO(x)                                                                              \
  {                                                                                                \
    cerr << x << endl;                                                                             \
  }
#define FPLLL_ABORT(x)                                                                             \
  {                                                                                                \
    cerr << "fplll: " << x << endl;                                                                \
    abort();                                                                                       \
  }
#define FPLLL_CHECK(x, y)                                                                          \
  {                                                                                                \
    if (!(x))                                                                                      \
      FPLLL_ABORT(y);                                                                              \
  }

#ifdef DEBUG
#include <cassert>
extern int debug_depth;
#define FPLLL_TRACE(x) std::cerr << "TRACE: " << std::string(debug_depth * 2, ' ') << x << std::endl
struct DebugTracer
{
  DebugTracer(const char *f) : f(f) { debug_depth++; }
  ~DebugTracer()
  {
    debug_depth--;
    FPLLL_TRACE("</" << f << ">");
  }
  std::string f;
};
#define FPLLL_DEBUG_ABORT(x) FPLLL_ABORT(x)
#define FPLLL_DEBUG_CHECK(x) assert(x);
#define FPLLL_TRACE_IN(x)                                                                          \
  FPLLL_TRACE("<" << __func__ << " " << x << ">");                                                 \
  DebugTracer debugTracer(__func__);
#define FPLLL_DEBUG_SAFEVECT
#else
#define FPLLL_DEBUG_ABORT(x)
#define FPLLL_DEBUG_CHECK(x)
#define FPLLL_TRACE(x)
#define FPLLL_TRACE_IN(x)
#endif

#define FPLLL_BEGIN_NAMESPACE                                                                      \
  namespace fplll                                                                                  \
  {
#define FPLLL_END_NAMESPACE }

/** \namespace fplll
    The fplll namespace */
FPLLL_BEGIN_NAMESPACE

using namespace std;

/* this trick will not work on 16-bit machines*/
#if (LONG_MAX == 2147483647L)
const int CPU_SIZE         = 32;
const int CPU_SIZE_1       = 30;
const double MAX_LONG_FAST = 0x1p30;
const long int EXPO_MAX    = 30;
#else
const int CPU_SIZE         = 64;
const int CPU_SIZE_1       = 53;
const double MAX_LONG_FAST = 0x1p53;
const long int EXPO_MAX    = 53;
#endif

const int MAX_EXP_DOUBLE = 1000;
const int PREC_DOUBLE    = 53;
const int PREC_DD        = 106;
const int PREC_QD        = 212;

const double LLL_DEF_DELTA        = 0.99;
const double LLL_DEF_ETA          = 0.51;
const double LLL_DEF_EPSILON      = 0.01;
const int SIZE_RED_FAILURE_THRESH = 5;

enum RedStatus
{
  RED_SUCCESS = 0,
  // Skips value 1
  RED_GSO_FAILURE     = 2,
  RED_BABAI_FAILURE   = 3,
  RED_LLL_FAILURE     = 4,
  RED_ENUM_FAILURE    = 5,
  RED_BKZ_FAILURE     = 6,
  RED_BKZ_TIME_LIMIT  = 7,
  RED_BKZ_LOOPS_LIMIT = 8,
  RED_STATUS_MAX      = 9
};

const char *const RED_STATUS_STR[RED_STATUS_MAX] = {"success",
                                                    "",
                                                    "infinite number in GSO",
                                                    "infinite loop in babai",
                                                    "infinite loop in LLL",
                                                    "error in SVP solver",
                                                    "error in BKZ",
                                                    "time limit exceeded in BKZ",
                                                    "loops limit exceeded in BKZ"};

enum LLLMethod
{
  LM_WRAPPER   = 0,
  LM_PROVED    = 1,
  LM_HEURISTIC = 2,
  LM_FAST      = 3
};

const char *const LLL_METHOD_STR[6] = {"wrapper", "proved", "heuristic", "fast"};

enum IntType
{
  ZT_MPZ    = 0,
  ZT_LONG   = 1,
  ZT_DOUBLE = 2
};

const char *const INT_TYPE_STR[5] = {"mpz", "long", "double"};

enum FloatType
{
  FT_DEFAULT     = 0,
  FT_DOUBLE      = 1,
  FT_LONG_DOUBLE = 2,
  FT_DPE         = 3,
  FT_DD          = 4,
  FT_QD          = 5,
  FT_MPFR        = 6
};

const char *const FLOAT_TYPE_STR[7] = {"", "double", "long double", "dpe", "dd", "qd", "mpfr"};

enum LLLFlags
{
  LLL_VERBOSE   = 1,
  LLL_EARLY_RED = 2,
  LLL_SIEGEL    = 4,
  LLL_DEFAULT   = 0
};

enum SVPMethod
{
  SVPM_FAST   = 0,
  SVPM_PROVED = 2
};

enum CVPMethod
{
  CVPM_FAST   = 0,
  CVPM_PROVED = 2
};

enum SVPFlags
{
  SVP_DEFAULT      = 0,
  SVP_VERBOSE      = 1,
  SVP_OVERRIDE_BND = 2,
  SVP_DUAL         = 4
};

enum CVPFlags
{
  CVP_DEFAULT = SVP_DEFAULT,
  CVP_VERBOSE = SVP_VERBOSE
};

const double BKZ_DEF_AUTO_ABORT_SCALE        = 1.0;
const int BKZ_DEF_AUTO_ABORT_MAX_NO_DEC      = 5;
const double BKZ_DEF_GH_FACTOR               = 1.1;
const double BKZ_DEF_MIN_SUCCESS_PROBABILITY = 0.5;
const int BKZ_DEF_RERANDOMIZATION_DENSITY    = 3;

enum BKZFlags
{
  BKZ_DEFAULT     = 0,
  BKZ_VERBOSE     = 1,
  BKZ_NO_LLL      = 2,
  BKZ_MAX_LOOPS   = 4,
  BKZ_MAX_TIME    = 8,
  BKZ_BOUNDED_LLL = 0x10,
  BKZ_AUTO_ABORT  = 0x20,
  BKZ_DUMP_GSO    = 0x40,
  BKZ_GH_BND      = 0x80,
  BKZ_SD_VARIANT  = 0x100,
  BKZ_SLD_RED     = 0x200
};

enum HKZFlags
{
  HKZ_DEFAULT = 0,
  HKZ_VERBOSE = 1
};

#ifndef FPLLL_DEFAULT_STRATEGY_PATH
#define FPLLL_DEFAULT_STRATEGY_PATH ""
#endif

#ifndef FPLLL_DEFAULT_STRATEGY
#define FPLLL_DEFAULT_STRATEGY ""
#endif

enum PrunerMetric
{
  PRUNER_METRIC_PROBABILITY_OF_SHORTEST = 0,
  PRUNER_METRIC_EXPECTED_SOLUTIONS      = 1
};

enum PrunerFlags
{
  PRUNER_CVP =
      0x1,  // Do not Halve the count of nodes, according to enumeration optimization from symmetry
  // Descent methods. If several activated, pruner will execute them in the order below.
  PRUNER_START_FROM_INPUT = 0x2,
  PRUNER_GRADIENT         = 0x4,  // Activate the gradient descent
  PRUNER_NELDER_MEAD      = 0x8,
  // Verbosity
  PRUNER_VERBOSE = 0x10,
  // Optimize w.r.t to half of the coefficients (those of even indices)
  // (by default this is not enabled)
  PRUNER_HALF = 0x20,
  // Optimize goal set to single enumeration cost while fixing the probability ~ target. Note that flags PRUNER_HALF and PRUNER_SINGLE are mutually exclusive.
  PRUNER_SINGLE = 0x40
};

#define PRUNER_ZEALOUS (PRUNER_GRADIENT | PRUNER_NELDER_MEAD)

FPLLL_END_NAMESPACE

#endif
