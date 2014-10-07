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

#include <cstring>
#include "fplll.h"

#define ABORT_MSG(y) {cerr << "fplll: " << y << endl; exit(1);}
#define CHECK(x, y) if (!(x)) ABORT_MSG(y)

using namespace std;
using namespace fplll;


enum Action {
  ACTION_LLL,
  ACTION_HKZ,
  ACTION_BKZ,
  ACTION_SVP,
  ACTION_CVP
};

struct Options {
  Options() : action(ACTION_LLL), method(LM_WRAPPER), intType(ZT_MPZ),
              floatType(FT_DEFAULT), delta(LLL_DEF_DELTA), eta(LLL_DEF_ETA),
              precision(0), earlyRed(false), siegel(false), noLLL(false),
              blockSize(0), verbose(false), inputFile(NULL),
              outputFormat(NULL), pruningFile(NULL) {
    bkzFlags = 0;
    bkzMaxLoops = 0;
    bkzMaxTime = 0;
  }
  Action action;
  LLLMethod method;
  IntType intType;
  FloatType floatType;
  double delta;
  double eta;
  int precision;
  bool earlyRed;
  bool siegel;

  bool noLLL;
  int blockSize;
  int bkzFlags;
  int bkzMaxLoops;
  double bkzMaxTime;
  string bkzDumpGSOFilename;

  bool verbose;
  const char* inputFile;
  const char* outputFormat;
  const char* pruningFile;
};

#endif
