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
#define FPLLL_V3_COMPAT
#define FPLLL_WITH_GETRUSAGE

#include <iostream>
#include <cstring>
#include <cstdio>
#include <climits>
#include <ctime>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#ifdef FPLLL_WITH_GETRUSAGE
  #include <sys/time.h>
  #include <sys/resource.h>
  #include <unistd.h>
#endif

#include <gmp.h>
#include <mpfr.h>
#ifdef FPLLL_WITH_DPE
  #include "dpe.h"
#endif

#if defined(__sun) || defined(__CYGWIN__)
  #include <ieeefp.h>
  extern "C" long double ldexpl(long double x, int exp);
  #ifndef NAN
    #define NAN __builtin_nanf("")
  #endif
#endif

#define FPLLL_INFO(x) {cerr << x << endl;}
#define FPLLL_ABORT(x) {cerr << "fplll: " << x << endl; abort();}
#define FPLLL_CHECK(x, y) {if (!(x)) FPLLL_ABORT(y);}

#ifdef DEBUG
  #include <cassert>
  extern int debugDepth;
  #define FPLLL_TRACE(x) std::cerr << "TRACE: " << std::string(debugDepth * 2, ' ') << x << std::endl
  struct DebugTracer {
    DebugTracer(const char* f) : f(f) {
      debugDepth++;
    }
    ~DebugTracer() {
      debugDepth--;
      FPLLL_TRACE("</" << f << ">");
    }
    std::string f;
  };
  #define FPLLL_DEBUG_ABORT(x) FPLLL_ABORT(x)
  #define FPLLL_DEBUG_CHECK(x) assert(x)
  #define FPLLL_TRACE_IN(x) FPLLL_TRACE("<" << __func__ << " " << x << ">"); DebugTracer debugTracer(__func__);
  #define FPLLL_DEBUG_SAFEVECT
#else
  #define FPLLL_DEBUG_ABORT(x)
  #define FPLLL_DEBUG_CHECK(x)
  #define FPLLL_TRACE(x)
  #define FPLLL_TRACE_IN(x)
#endif

#define FPLLL_BEGIN_NAMESPACE namespace fplll {
#define FPLLL_END_NAMESPACE }

/** \namespace fplll
    The fplll namespace */
FPLLL_BEGIN_NAMESPACE

using namespace std;

/* this trick will not work on 16-bit machines*/
#if (LONG_MAX==2147483647L)
  const int CPU_SIZE = 32;
  const int CPU_SIZE_1 = 30;
  const double MAX_LONG_FAST = 0x1p30;
  const long int EXPO_MAX = 30;
#else
  const int CPU_SIZE = 64;
  const int CPU_SIZE_1 = 53;
  const double MAX_LONG_FAST = 0x1p53;
  const long int EXPO_MAX = 53;
#endif

const int MAX_EXP_DOUBLE = 1000;
const int PREC_DOUBLE = 53;

const double LLL_DEF_DELTA = 0.99;
const double LLL_DEF_ETA = 0.51;
const double LLL_DEF_EPSILON = 0.01;
const int SIZE_RED_FAILURE_THRESH = 5;


enum RedStatus {
  RED_SUCCESS = 0,
  // Skips value 1
  RED_GSO_FAILURE = 2,
  RED_BABAI_FAILURE,
  RED_LLL_FAILURE,
  RED_ENUM_FAILURE,
  RED_BKZ_FAILURE,
  RED_BKZ_TIME_LIMIT,
  RED_BKZ_LOOPS_LIMIT,
  RED_STATUS_MAX
};

const char* const RED_STATUS_STR[RED_STATUS_MAX] =
  {"success",
   "",
   "infinite number in GSO",
   "infinite loop in babai",
   "infinite loop in LLL",
   "error in SVP solver",
   "error in BKZ",
   "time limit exceeded in BKZ",
   "loops limit exceeded in BKZ"};

enum LLLMethod {
  LM_WRAPPER,
  LM_PROVED,
  LM_HEURISTIC,
  LM_FAST
};

const char* const LLL_METHOD_STR[6] =
  {"wrapper", "proved", "heuristic", "fast"};

enum IntType {
  ZT_MPZ,
  ZT_LONG,
  ZT_DOUBLE
};

const char* const INT_TYPE_STR[5] =
  {"mpz", "long", "double"};

enum FloatType {
  FT_DEFAULT,
  FT_DOUBLE,
  FT_LONG_DOUBLE,
  FT_DPE,
  FT_MPFR
};

const char* const FLOAT_TYPE_STR[5] =
  {"", "double", "long double", "dpe", "mpfr"};

enum LLLFlags {
  LLL_VERBOSE = 1,
  LLL_EARLY_RED = 2,
  LLL_SIEGEL = 4,
  LLL_DEFAULT = 0
};

enum SVPMethod {
  SVPM_FAST = 0,
  SVPM_PROVED = 2
};

enum SVPFlags {
  SVP_DEFAULT = 0,
  SVP_VERBOSE = 1,
  SVP_OVERRIDE_BND = 2
};

enum CVPFlags {
  CVP_DEFAULT = SVP_DEFAULT,
  CVP_VERBOSE = SVP_VERBOSE
};

enum BKZFlags {
  BKZ_DEFAULT = 0,
  BKZ_VERBOSE = 1,
  BKZ_NO_LLL = 2,
  BKZ_MAX_LOOPS = 4,
  BKZ_MAX_TIME = 8,
  BKZ_BOUNDED_LLL = 0x10,
  BKZ_AUTO_ABORT = 0x20,
  BKZ_DUMP_GSO = 0x40,
  BKZ_GH_BND = 0x80
};

enum HKZFlags {
  HKZ_DEFAULT = 0,
  HKZ_VERBOSE = 1
};

FPLLL_END_NAMESPACE

#endif
