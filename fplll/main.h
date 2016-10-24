/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2008-2011 Xavier Pujol.

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

#ifndef FPLLL_MAIN_H
#define FPLLL_MAIN_H

#include "fplll.h"
#include <cstring>

#define ABORT_MSG(y)                                                                               \
  {                                                                                                \
    cerr << "fplll: " << y << endl;                                                                \
    exit(1);                                                                                       \
  }
#define CHECK(x, y)                                                                                \
  if (!(x))                                                                                        \
  ABORT_MSG(y)

using namespace std;
using namespace fplll;

enum Action
{
  ACTION_LLL,
  ACTION_HKZ,
  ACTION_BKZ,
  ACTION_SVP,
  ACTION_CVP
};

struct Options
{
  Options()
      : action(ACTION_LLL), method(LM_WRAPPER), int_type(ZT_MPZ), float_type(FT_DEFAULT),
        delta(LLL_DEF_DELTA), eta(LLL_DEF_ETA), precision(0), early_red(false), siegel(false),
        no_lll(false), block_size(0), bkz_gh_factor(1.1), verbose(false), input_file(NULL),
        output_format(NULL)
  {
    bkz_flags     = 0;
    bkz_max_loops = 0;
    bkz_max_time  = 0;
  }
  Action action;
  LLLMethod method;
  IntType int_type;
  FloatType float_type;
  double delta;
  double eta;
  int precision;
  bool early_red;
  bool siegel;

  bool no_lll;
  int block_size;
  int bkz_flags;
  int bkz_max_loops;
  double bkz_max_time;
  string bkz_dump_gso_filename;
  double bkz_gh_factor;
  string bkz_strategy_file;

  bool verbose;
  const char *input_file;
  const char *output_format;
};

#endif
