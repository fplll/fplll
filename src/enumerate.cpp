/* Copyright (C) 2008-2011 Xavier Pujol.

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

#include "enumerate.h"

FPLLL_BEGIN_NAMESPACE

enumf Enumeration::mut[DMAX][DMAX];
enumf Enumeration::centerPartSums[DMAX][DMAX + 1];
EnumfVect Enumeration::rdiag;
EnumfVect Enumeration::x;
EnumfVect Enumeration::dx;
EnumfVect Enumeration::ddx;
EnumfVect Enumeration::dist;
EnumfVect Enumeration::center;
EnumfVect Enumeration::centerPartSum;
EnumfVect Enumeration::maxDists;
EnumfVect Enumeration::centerLoopBg;
int Enumeration::d;
int Enumeration::k;
int Enumeration::kEnd;
int Enumeration::kMax;

static const vector<FP_NR<double> > EMPTY_DOUBLE_VECT;

template<class FT>
void Enumeration::prepareEnumeration(enumf maxDist, const vector<FT>& subTree, bool solvingSVP) {
  //FPLLL_TRACE_IN("maxDist=" << maxDist);
  // true->SVP and all coordinates in subTree are null
  bool svpBeginning = solvingSVP;
  enumf newX, newDist = enumf(0.0);

  kEnd = d - subTree.size();
  // Prepares the loop (goes to the first vector)
  for (k = d - 1; k >= 0 && newDist <= maxDist; k--) {
    enumf newCenter = centerPartSum[k];
    for (int j = k + 1; j < kEnd; j++) {
      newCenter -= x[j] * mut[k][j];
    }

    if (k >= kEnd) {
      newX = subTree[k - kEnd].get_d();
      if (newX != 0.0) svpBeginning = false;
      for (int j = 0; j < k; j++) {
        centerPartSum[j] -= newX * mut[j][k];
      }
    }
    else {
      newX = rint(newCenter);
      center[k] = newCenter;
      dist[k] = newDist;
      dx[k] = enumf(0.0);
      ddx[k] = newCenter < newX ? enumf(1.0) : enumf(-1.0);
      /*FPLLL_TRACE("k=" << k << " x_k=" << newX << " center_k=" << center[k]
              << " dist_k=" << dist[k] << " dx_k=" << dx[k]
              << " ddx_k=" << ddx[k]);*/
    }
    x[k] = newX;
    enumf y = newCenter - newX;
    newDist += y * y * rdiag[k];
  }
  if (!svpBeginning) {
    kMax = kEnd; // The last non-zero coordinate of x will not stay positive
  }
  else {
    kMax = 0;
    x[0] = enumf(1.0); // Excludes (0,...,0) from the enumeration
  }
  k++;
  // now, 0 <= k <= kEnd - 1
}

/* Input: rdiag, center, dist, centerPartSum, x, dx, ddx, maxDists, k, kEnd, kMax
   Output: center, dist, centerPartSum, x, dx, ddx, k, kMax, newMaxDist, newKMax */
bool Enumeration::enumerateLoop(enumf& newMaxDist, int& newKMax) {
  //FPLLL_TRACE_IN("k=" << k);
  if (k >= kEnd) return false;

  for (int i = 0; i < kEnd; i++) {
    centerLoopBg[i] = kEnd - 1;
    centerPartSums[i][kEnd] = centerPartSum[i];
  }

  while (true) {
    enumf y = center[k] - x[k];
    enumf newDist = dist[k] + y * y * rdiag[k];
    /*FPLLL_TRACE("k=" << k << " x_k=" << x[k] << " center_k=" << center[k]
            << " dist_k=" << dist[k] << " r_k=" << rdiag[k]
            << " y=" << y << " newDist=" << newDist);*/
    if (newDist <= maxDists[k]) {
      k--;
      if (k < 0) {
        newMaxDist = newDist;
        newKMax = kMax;
        return true; // New solution found
      }

      for (int j = centerLoopBg[k]; j > k; j--) {
        centerPartSums[k][j] = centerPartSums[k][j + 1] - x[j] * mut[k][j];
      }
      enumf newCenter = centerPartSums[k][k + 1];
      if (k > 0) centerLoopBg[k - 1] = max(centerLoopBg[k - 1], centerLoopBg[k]);
      centerLoopBg[k] = k + 1;

      center[k] = newCenter;
      dist[k] = newDist;
      x[k] = rint(newCenter);
      dx[k] = enumf(0.0);
      ddx[k] = newCenter < x[k] ? enumf(1.0) : enumf(-1.0);
    }
    else {
      if (!nextPosUp()) {
        // End of the enumeration
        return false;
      }
    }
  }
}

/* Input: d, rdiag, center, dist, centerPartSum, x, dx, ddx, k, kEnd
   Output: center, dist, centerPartSum, x, dx, ddx, k
   Internal use: kMax */
template<class FT>
void Enumeration::enumerate(enumf& maxDist, long normExp, Evaluator<FT>& evaluator, const vector<double>& pruning) {
  vector<FT> fX(d);
  enumf newMaxDist;
  while (true) {
    if (pruning.empty()) {
      fill(maxDists, maxDists + d, maxDist);
    }
    else {
      for (int i = 0; i < d; i++) {
        maxDists[i] = pruning[i] * maxDist;
      }
    }

    if (!enumerateLoop(newMaxDist, kMax)) {
      break;
    }

    // We have found a solution
    for (int j = 0; j < d; j++) {
      fX[j] = x[j];
    }
    evaluator.evalSol(fX, newMaxDist, maxDist, normExp);
    k = -1;
    // Goes to the next step and continues the loop
    nextPosUp();
  }
}

template<class FT>
void Enumeration::enumerate(MatGSO<Integer, FT>& gso, FT& fMaxDist, long maxDistExpo,
               Evaluator<FT>& evaluator, const vector<FT>& targetCoord,
               const vector<FT>& subTree, int first, int last,
               const vector<double>& pruning) {
  bool solvingSVP;    // true->SVP, false->CVP
  enumf maxDist;
  FT fR, fMu, fMaxDistNorm;
  long rExpo, normExp = LONG_MIN;

  if (last == -1) last = gso.d;
  d = last - first;

  solvingSVP = targetCoord.empty();
  FPLLL_CHECK(d <= DMAX, "enumerate: dimension is too high");

  // FT->enumf conversion and transposition of mu
  for (int i = 0; i < d; i++) {
    fR = gso.getRExp(i + first, i + first, rExpo);
    normExp = max(normExp, rExpo + fR.exponent());
  }

  fMaxDistNorm.mul_2si(fMaxDist, maxDistExpo - normExp);
  maxDist = fMaxDistNorm.get_d(GMP_RNDU);

  for (int i = 0; i < d; i++) {
    fR = gso.getRExp(i + first, i + first, rExpo);
    fR.mul_2si(fR, rExpo - normExp);
    rdiag[i] = fR.get_d();

    if (solvingSVP)
      centerPartSum[i] = 0.0;
    else
      centerPartSum[i] = targetCoord[i + first].get_d();

    for (int j = i + 1; j < d; j++) {
      gso.getMu(fMu, j + first, i + first);
      mut[i][j] = fMu.get_d();
    }
  }

  prepareEnumeration(maxDist, subTree, solvingSVP);
  enumerate(maxDist, normExp, evaluator, pruning);

  fMaxDistNorm = maxDist; // Exact
  fMaxDist.mul_2si(fMaxDistNorm, normExp - maxDistExpo);
}


void Enumeration::enumerateDouble(MatGSO<Z_NR<double>, FP_NR<double> >& gso,
            FP_NR<double>& fMaxDist, Evaluator<FP_NR<double> >& evaluator, int first, int last,
            const vector<double>& pruning) {
  enumf maxDist = fMaxDist.get_d();
  FP_NR<double> fR, fMu;

  if (last == -1) last = gso.d;
  d = last - first;
  FPLLL_CHECK(d <= DMAX, "enumerate: dimension is too high");

  // FT->enumf conversion and transposition of mu
  for (int i = 0; i < d; i++) {
    gso.getR(fR, i + first, i + first);
    rdiag[i] = fR.get_d();
    centerPartSum[i] = 0.0;
    for (int j = i + 1; j < d; j++) {
      gso.getMu(fMu, j + first, i + first);
      mut[i][j] = fMu.get_d();
    }
  }

  prepareEnumeration(maxDist, EMPTY_DOUBLE_VECT, true);
  enumerate(maxDist, 0, evaluator, pruning);
  fMaxDist = maxDist;
}

FPLLL_END_NAMESPACE
