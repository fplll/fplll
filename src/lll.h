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

#ifndef FPLLL_LLL_H
#define FPLLL_LLL_H

#include "gso.h"

FPLLL_BEGIN_NAMESPACE

/* The precision of FT must be defined before creating an instance of
   LLLReduction. The matrix b can be modified between each call to lll. */

template<class ZT, class FT>
class LLLReduction {
public:
  /**
   * Constructor.
   * The precision of FT must be defined before creating an instance of the
   * class and must remain the same until the object is destroyed (or no longer
   * needed).
   */
  LLLReduction(MatGSO<ZT, FT>& m, double delta, double eta,
               int flags);

  bool lll(int kappaMin = 0, int kappaStart = 0, int kappaEnd = -1);
  bool sizeReduction(int kappaMin = 0, int kappaEnd = -1);

  int status;
  int finalKappa;
  int lastEarlyRed;
  int zeros;
  int nSwaps;

private:
  bool babai(int kappa, int nCols);
  bool earlyReduction(int start);
  void printParams();
  bool setStatus(int newStatus);

  MatGSO<ZT, FT>& m;
  FT delta, eta, swapThreshold;

  bool enableEarlyRed;
  bool siegel;
  bool verbose;

  vector<FT> lovaszTests;
  vector<FT> babaiMu;
  vector<long> babaiExpo;
  ZT ztmp1;
  FT muMant, ftmp1;
};

template<class ZT, class FT>
bool isLLLReduced(MatGSO<ZT, FT>& m, double delta, double eta);

FPLLL_END_NAMESPACE

#endif
