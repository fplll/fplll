/* Copyright (C) 2008-2011 Xavier Pujol
   (C) 2015 Michael Walter.

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

template<typename FT, int MaxDimension>
void Enumeration2<FT,MaxDimension>::enumerate(int first, int last, FT& fMaxDist, long maxDistExpo,
                                const vector<FT>& targetCoord,
                                const vector<enumxt>& subTree,
                                const vector<enumf>& pruning,
                                bool _dual)
{
    bool solvingSVP = targetCoord.empty();
    dual = _dual;
    if (last == -1)
        last = _gso.d;
    d = last - first;
    FPLLL_CHECK(d <= DMAX, "enumerate: dimension is too high");
    FPLLL_CHECK((solvingSVP || !dual), "CVP for dual not implemented! What does that even mean? ");
    FPLLL_CHECK((subTree.empty() || !dual), "Subtree enumeration for dual not implemented!");

    if (solvingSVP)
    {
        for (int i = 0; i < d; ++i)
            centerPartSum[i] = 0.0;
    }
    else
    {
        for (int i = 0; i < d; ++i)
            centerPartSum[i] = targetCoord[i + first].get_d();
    }
    
    FT fR, fMu, fMaxDistNorm;
    long rExpo, normExp = -1;
    for (int i = 0; i < d; ++i)
    {
        fR = _gso.getRExp(i + first, i + first, rExpo);
        normExp = max(normExp, rExpo + fR.exponent());
    }
    fMaxDistNorm.mul_2si(fMaxDist, dual ? normExp-maxDistExpo : maxDistExpo-normExp);
    enumf maxDist = fMaxDistNorm.get_d(GMP_RNDU);

    if (dual)
    {
        for (int i = 0; i < d; ++i)
        {
            fR = _gso.getRExp(i + first, i + first, rExpo);
            fR.mul_2si(fR, rExpo - normExp);
            rdiag[d-i-1] = enumf(1.0)/fR.get_d();
        }
        for (int i = 0; i < d; ++i)
        {
            for (int j = i + 1; j < d; ++j)
            {
                _gso.getMu(fMu, j + first, i + first);
                mut[d-j-1][d-i-1] = -fMu.get_d();
            }
        }
    }
    else
    {
        for (int i = 0; i < d; ++i)
        {
            fR = _gso.getRExp(i + first, i + first, rExpo);
            fR.mul_2si(fR, rExpo - normExp);
            rdiag[i] = fR.get_d();
        }
        for (int i = 0; i < d; ++i)
        {
            for (int j = i + 1; j < d; ++j)
            {
                _gso.getMu(fMu, j + first, i + first);
                mut[i][j] = fMu.get_d();
            }
        }
    }
    
    prepareEnumeration(maxDist, subTree, solvingSVP);
    enumerate(maxDist, normExp, pruning);
  
    fMaxDistNorm = maxDist; // Exact
  
    fMaxDist.mul_2si(fMaxDistNorm, dual ? maxDistExpo-normExp : normExp-maxDistExpo);
  
    if (dual && !_evaluator.solCoord.empty())
        reverseBySwap(_evaluator.solCoord, 0, d-1);
}

template<typename FT, int MaxDimension>
void Enumeration2<FT,MaxDimension>::prepareEnumeration(enumf maxDist, const vector<enumxt>& subTree, bool solvingSVP)
{
    bool svpBeginning = solvingSVP;
    
    enumxt newX;
    enumf newDist = 0.0;
    kEnd = d - subTree.size();
    for (k = d - 1; k >= 0 && newDist <= maxDist; --k)
    {
        enumf newCenter = centerPartSum[k];
        for (int j = k + 1; j < kEnd; ++j)
        {
            newCenter -= (dual ? alpha[j] : x[j]) * mut[k][j];
        }
        if (k >= kEnd)
        {
            newX = subTree[k - kEnd];
            if (newX != 0)
            {
                svpBeginning = false;
            }
            for (int j = 0; j < k; ++j)
            {
                centerPartSum[j] -= newX * mut[j][k];
            }
        }
        else
        {
            roundto(newX, newCenter); // newX = rint(newCenter) / lrint(newCenter)
            center[k] = newCenter;
            dist[k] = newDist;
            dx[k] = ddx[k] = (((int)(newCenter >= newX) & 1) << 1) - 1;
//            dx[k] = (enumxt)(0);
//            ddx[k] = newCenter < newX ? enumxt(1) : enumxt(-1);
        }
        x[k] = newX;
        alpha[k] = newX - newCenter;
        newDist += alpha[k] * alpha[k] * rdiag[k];
    }
    if (!svpBeginning) 
    {
        kMax = kEnd; // The last non-zero coordinate of x will not stay positive
    }
    else 
    {
        kMax = 0;
        x[0] = enumxt(1); // Excludes (0,...,0) from the enumeration
    }
    ++k;
}

template<typename FT, int MaxDimension>
bool Enumeration2<FT,MaxDimension>::enumerateLoop(enumf& newMaxDist, int& newKMax) 
{
    if (k >= kEnd)
        return false;
    
    for (int i = 0; i < kEnd; ++i)
    {
        centerLoopBg[i] = kEnd - 1;
        centerPartSums[i][kEnd] = centerPartSum[i];
    }
    
    while (true)
    {
        alpha[k] = x[k] - center[k];
        enumf newDist = dist[k] + alpha[k]*alpha[k]*rdiag[k];
        if (newDist <= maxDists[k])
        {
            --k;
            ++nodes;
            if (k < 0)
            {
                newMaxDist = newDist;
                newKMax = kMax;
                return true;
            }
            
            if (dual)
            {
                for (int j = centerLoopBg[k]; j > k; --j)
                    centerPartSums[k][j] = centerPartSums[k][j+1] - alpha[j] * mut[k][j];
            }
            else
            {
                for (int j = centerLoopBg[k]; j > k; --j)
                    centerPartSums[k][j] = centerPartSums[k][j+1] - x[j] * mut[k][j];
            }
            
            enumf newCenter = centerPartSums[k][k+1];
            if (k > 0)
            {
                centerLoopBg[k-1] = max(centerLoopBg[k-1], centerLoopBg[k]);
            }
            centerLoopBg[k] = k+1;
            
            center[k] = newCenter;
            dist[k] = newDist;
            roundto(x[k], newCenter);
            dx[k] = ddx[k] = (((int)(newCenter >= x[k]) & 1) << 1) - 1;
//            dx[k] = 0;
//            ddx[k] = newCenter < x[k] ? enumxt(1) : enumxt(-1);
        }
        else
        {
            if (!nextPosUp())
                return false;
        }
    }
}

template<typename FT, int MaxDimension>
void Enumeration2<FT,MaxDimension>::enumerate(enumf& maxDist, long normExp, const vector<double>& pruning) 
{
    vector<FT> fX(d);
    enumf newMaxDist;
    nodes = 0;
    while (true)
    {
        if (pruning.empty())
        {
            fill(&maxDists[0] + 0, &maxDists[0] + d, maxDist);
        }
        else
        {
            for (int i = 0; i < d; ++i)
                maxDists[i] = pruning[i] * maxDist;
        }
        if (!enumerateLoop(newMaxDist, kMax))
            break;
            
        for (int j = 0; j < d; ++j)
            fX[j] = x[j];
            
        _evaluator.evalSol(fX, newMaxDist, maxDist, normExp);
        k = -1;
        
        nextPosUp();
    }
}


#if 0
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
EnumfVect Enumeration::alpha;
int Enumeration::d;
int Enumeration::k;
int Enumeration::kEnd;
int Enumeration::kMax;
bool Enumeration::dual;
long Enumeration::nodes; // uint64_t is only a thing starting with C++11

static const vector<FP_NR<double> > EMPTY_DOUBLE_VECT;

template<class FT>
void Enumeration::prepareEnumeration(enumf maxDist, const vector<FT>& subTree, bool solvingSVP) {
  //FPLLL_TRACE_IN("maxDist=" << maxDist);
  // true->SVP and all coordinates in subTree are null
  bool svpBeginning = solvingSVP;
  enumf newX, newDist = enumf(0.0);
  enumf *co;
  
  kEnd = d - subTree.size();
  // Prepares the loop (goes to the first vector)
  for (k = d - 1; k >= 0 && newDist <= maxDist; k--) {
    enumf newCenter = centerPartSum[k];
    for (int j = k + 1; j < kEnd; j++) {
      co = dual ? &alpha[j] : &x[j];
      newCenter -= (*co) * mut[k][j];
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
    alpha[k] = newX - newCenter;
    newDist += alpha[k] * alpha[k] * rdiag[k];
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
template< class EnumType >
bool Enumeration::enumerateLoop(enumf& newMaxDist, int& newKMax) {
  //FPLLL_TRACE_IN("k=" << k);
  if (k >= kEnd) return false;

  for (int i = 0; i < kEnd; i++) {
    centerLoopBg[i] = kEnd - 1;
    centerPartSums[i][kEnd] = centerPartSum[i];
  }

  while (true) {
    alpha[k] = x[k] - center[k];
    enumf newDist = dist[k] + alpha[k] * alpha[k] * rdiag[k];
    /*FPLLL_TRACE("k=" << k << " x_k=" << x[k] << " center_k=" << center[k]
            << " dist_k=" << dist[k] << " r_k=" << rdiag[k]
            << " y=" << y << " newDist=" << newDist);*/
    if (newDist <= maxDists[k]) {
      k--;
      nodes++;
      if (k < 0) {
        newMaxDist = newDist;
        newKMax = kMax;
        return true; // New solution found
      }
      
      for (int j = centerLoopBg[k]; j > k; j--) {
        centerPartSums[k][j] = centerPartSums[k][j + 1] - EnumType::center_summand_coefficient(alpha[j], x[j]) * mut[k][j];
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
  nodes = 0;
  while (true) {
    if (pruning.empty()) {
      fill(maxDists, maxDists + d, maxDist);
    }
    else {
      for (int i = 0; i < d; i++) {
        maxDists[i] = pruning[i] * maxDist;
      }
    }
    
    if (dual) {
      if (!enumerateLoop<DualEnum>(newMaxDist, kMax)) {
        break;
      }
    } else {
      if (!enumerateLoop<PrimalEnum>(newMaxDist, kMax)) {
        break;
      }
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
               const vector<double>& pruning, bool dual) {
  bool solvingSVP;    // true->SVP, false->CVP
  Enumeration::dual = dual;
  enumf maxDist;
  FT fR, fMu, fMaxDistNorm;
  long rExpo, normExp = LONG_MIN;

  if (last == -1) last = gso.d;
  d = last - first;

  solvingSVP = targetCoord.empty();
  FPLLL_CHECK(d <= DMAX, "enumerate: dimension is too high");

  FPLLL_CHECK((solvingSVP || !dual), "CVP for dual not implemented! What does that even mean? ");
  FPLLL_CHECK((subTree.empty() || !dual), "Subtree enumeration for dual not implemented!");

  // FT->enumf conversion and transposition of mu
  for (int i = 0; i < d; i++) {
    fR = gso.getRExp(i + first, i + first, rExpo);
    normExp = max(normExp, rExpo + fR.exponent());
  }
  
  if (dual) {
    fMaxDistNorm.mul_2si(fMaxDist, normExp - maxDistExpo);
  } else {
    fMaxDistNorm.mul_2si(fMaxDist, maxDistExpo - normExp);
  }
  maxDist = fMaxDistNorm.get_d(GMP_RNDU);

  for (int i = 0; i < d; i++) {
    fR = gso.getRExp(i + first, i + first, rExpo);
    fR.mul_2si(fR, rExpo - normExp);
    if (dual) {
      rdiag[d-i-1] = enumf(1.0)/fR.get_d();
    } else {
      rdiag[i] = fR.get_d();
    }

    if (solvingSVP)
      centerPartSum[i] = 0.0;
    else
      centerPartSum[i] = targetCoord[i + first].get_d();

    for (int j = i + 1; j < d; j++) {
      gso.getMu(fMu, j + first, i + first);
      if (dual) {
        mut[d-j-1][d-i-1] = -fMu.get_d();
      } else {
        mut[i][j] = fMu.get_d();
      }
    }
  }

  prepareEnumeration(maxDist, subTree, solvingSVP);
  enumerate(maxDist, normExp, evaluator, pruning);
  
  fMaxDistNorm = maxDist; // Exact
  
  if (dual) {
    fMaxDist.mul_2si(fMaxDistNorm, maxDistExpo - normExp);
  } else {
    fMaxDist.mul_2si(fMaxDistNorm, normExp - maxDistExpo);
  }
  
  if (dual && !evaluator.solCoord.empty()) reverseBySwap(evaluator.solCoord, 0, d-1);
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

#endif

FPLLL_END_NAMESPACE
