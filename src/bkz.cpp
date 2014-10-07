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

/* Template source file */
#include "bkz.h"
#include "enumerate.h"

FPLLL_BEGIN_NAMESPACE

static double getSlope(const vector<double>& x) {
  int n = x.size();
  double iMean = (n - 1) * 0.5, xMean = 0, v1 = 0, v2 = 0;
  for (int i = 0; i < n; i++) {
    xMean += x[i];
  }
  xMean /= n;
  for (int i = 0; i < n; i++) {
    v1 += (i - iMean) * (x[i] - xMean);
    v2 += (i - iMean) * (i - iMean);
  }
  return v1 / v2;
}

template<class FT>
bool BKZAutoAbort<FT>::testAbort() {
  const int MAX_NO_DEC = 5;
  FT f, logF; 
  long expo;
  if (x.empty()) x.resize(numRows);
  for (int i = 0; i < numRows; i++) {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    x[i] = logF.get_d() + expo * log(2.0);
  }
  double newSlope = -getSlope(x);
  if (noDec == -1 || newSlope < oldSlope)
    noDec = 0;
  else
    noDec++;
  oldSlope = min(oldSlope, newSlope);
  return noDec >= MAX_NO_DEC;
}

template<class FT>
BKZReduction<FT>::BKZReduction(MatGSO<Integer, FT>& m,
        LLLReduction<Integer, FT>& lllObj, const BKZParam& param) :
  status(RED_SUCCESS), param(param), m(m), lllObj(lllObj)
{
  for (numRows = m.d; numRows > 0 && m.b[numRows - 1].is_zero(); numRows--) {}
  this->delta = param.delta;
}

template<class FT>
BKZReduction<FT>::~BKZReduction() {
}

template<class FT>
bool BKZReduction<FT>::svpReduction(int kappa, int blockSize, bool& clean) {
  long maxDistExpo;

  int lllStart = (param.flags & BKZ_BOUNDED_LLL) ? kappa : 0;

  if (!lllObj.lll(lllStart, kappa, kappa + blockSize)) {
    return setStatus(lllObj.status);
  }
  if (lllObj.nSwaps > 0) {
    clean = false;
  }
  maxDist = m.getRExp(kappa, kappa, maxDistExpo);
  deltaMaxDist.mul(delta, maxDist);
  vector<FT>& solCoord = evaluator.solCoord;
  solCoord.clear();
  Enumeration::enumerate(m, maxDist, maxDistExpo, evaluator, emptySubTree,
            emptySubTree, kappa, kappa + blockSize, param.pruning);
  if (solCoord.empty()) {
    return setStatus(RED_ENUM_FAILURE);
  }

  // Is it already in the basis ?
  int nzVectors = 0, iVector = -1;
  for (int i = 0; i < blockSize; i++) {
    if (!solCoord[i].is_zero()) {
      nzVectors++;
      if (iVector == -1 && (solCoord[i].get_d() == 1 || solCoord[i].get_d() == -1))
        iVector = i;
    }
  }
  FPLLL_DEBUG_CHECK(nzVectors > 0);

  if (maxDist >= deltaMaxDist) {
    return true; // Do nothing
  }

  if (nzVectors == 1) {
    // Yes, it is another vector
    FPLLL_DEBUG_CHECK(iVector != -1 && iVector != 0);
    m.moveRow(kappa + iVector, kappa);
    if (!lllObj.sizeReduction(kappa, kappa + 1))
      return setStatus(lllObj.status);
  }
  else {
    // No, general case
    int d = m.d;
    m.createRow();
    m.rowOpBegin(d, d + 1);
    for (int i = 0; i < blockSize; i++) {
      m.row_addmul(d, kappa + i, solCoord[i]);
    }
    m.rowOpEnd(d, d + 1);
    m.moveRow(d, kappa);
    if (!lllObj.lll(kappa, kappa, kappa + blockSize + 1))
      return setStatus(lllObj.status);
    FPLLL_DEBUG_CHECK(m.b[kappa + blockSize].is_zero());
    m.moveRow(kappa + blockSize, d);
    m.removeLastRow();
  }
  clean = false;
  return true;
}

template<class FT>
bool BKZReduction<FT>::bkzLoop(int& kappaMax, bool& clean) {
  int flags = param.flags;
  for (int kappa = 0; kappa < numRows-1; kappa++) {
    // SVP-reduces a block
    int blockSize = min(param.blockSize, numRows - kappa);
    if (!svpReduction(kappa, blockSize, clean)) return false;
    if ((flags & BKZ_VERBOSE) && kappaMax < kappa && clean) {
      cerr << "Block [1-" << kappa + 1
           << "] BKZ-reduced for the first time" << endl;
      kappaMax = kappa;
    }
  }
  return true;
}

template<class FT>
bool BKZReduction<FT>::bkz() {
  int flags = param.flags;
  int finalStatus = RED_SUCCESS;

  if (param.blockSize < 2)
    return setStatus(RED_SUCCESS);

  int kappaMax = 0;
  int iLoop =0;
  BKZAutoAbort<FT> autoAbort(m, numRows);

  if (flags & BKZ_VERBOSE) printParams();
  cputimeStart = cputime();

  m.discoverAllRows();

  for (iLoop = 0;; iLoop++) {
    if ((flags & BKZ_MAX_LOOPS) && iLoop >= param.maxLoops) {
      finalStatus = RED_BKZ_LOOPS_LIMIT;
      break;
    }
    if ((flags & BKZ_MAX_TIME) && (cputime() - cputimeStart) * 0.001 >= param.maxTime) { 
      finalStatus = RED_BKZ_TIME_LIMIT;
      break;
    }
    if ((flags & BKZ_AUTO_ABORT) && autoAbort.testAbort()) break;
    bool clean = true;
    if (!bkzLoop(kappaMax, clean)) return false;
    if (clean || param.blockSize >= numRows) break;
    if (flags & BKZ_VERBOSE) {
      FT r0;
      Float fr0;
      long expo;
      r0 = m.getRExp(0, 0, expo);
      fr0 = r0.get_d();
      fr0.mul_2si(fr0, expo);
      cerr << "End of BKZ loop, time=" << cputime() * 0.001 << ", r_0 = " << fr0 << endl;
    }
  }
  return setStatus(finalStatus);
}

template<class FT>
void BKZReduction<FT>::printParams() {
  cerr << "Entering BKZ"
       << "\nblocksize = " << param.blockSize << endl;
}

template<class FT>
bool BKZReduction<FT>::setStatus(int newStatus) {
  status = newStatus;
  if (param.flags & BKZ_VERBOSE) {
    if (status == RED_SUCCESS)
      cerr << "End of BKZ: success" << endl;
    else
      cerr << "End of BKZ: failure: " << RED_STATUS_STR[status] << endl;
  }
  return status == RED_SUCCESS;
}

FPLLL_END_NAMESPACE
