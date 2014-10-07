/* Copyright (C) 2011 Xavier Pujol.

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

#ifndef FPLLL_BKZ_H
#define FPLLL_BKZ_H

#include "lll.h"

FPLLL_BEGIN_NAMESPACE

template<class FT>
class BKZAutoAbort {
public:
  BKZAutoAbort(MatGSO<Integer, FT>& m, int numRows) : m(m),
      oldSlope(numeric_limits<double>::max()), noDec(-1), numRows(numRows) {}
  bool testAbort();

private:
  MatGSO<Integer, FT>& m;
  vector<double> x;
  double oldSlope;
  int noDec;
  int numRows;
};

/* The matrix must be LLL-reduced */
template<class FT>
class BKZReduction {
public:
  BKZReduction(MatGSO<Integer, FT>& m, LLLReduction<Integer, FT>& lllObj,
               const BKZParam& param);
  ~BKZReduction();

  bool svpReduction(int kappa, int blockSize, bool& clean);
  bool bkzLoop(int& kappaMax, int& nIter, bool& clean);
  bool bkz();

  int status;

private:
  void printParams();
  bool setStatus(int newStatus);

  const BKZParam& param;
  int numRows;
  MatGSO<Integer, FT>& m;
  LLLReduction<Integer, FT>& lllObj;
  FastEvaluator<FT> evaluator;
  FT delta;

  // Temporary data
  const vector<FT> emptyTarget, emptySubTree;
  FT maxDist, deltaMaxDist;
  double cputimeStart;
};

FPLLL_END_NAMESPACE

#endif
