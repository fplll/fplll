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

#ifndef FPLLL_ENUMERATE_H
#define FPLLL_ENUMERATE_H

#include <array>
#include "evaluator.h"

FPLLL_BEGIN_NAMESPACE

inline void roundto(int& dest, const double& src) { dest = lrint(src); }
inline void roundto(double& dest, const double& src) { dest = rint(src); }

template<typename FT, int MaxDimension = 160>
class Enumeration2
{
public:
    static const int DMAX = MaxDimension;
    Enumeration2(MatGSO<Integer, FT>& gso, Evaluator<FT>& evaluator)
        : _gso(gso), _evaluator(evaluator)
    {
    }
    
    void enumerate(int first, int last,
                FT& fMaxDist, long maxDistExpo, 
                const vector<FT>& targetCoord,
                const vector<enumxt>& subTree,
                const vector<enumf>& pruning,
                bool dual = false);

    inline uint64_t getNodes() const { return nodes; }
    
private:
    MatGSO<Integer, FT>& _gso; 
    Evaluator<FT>& _evaluator;
    
    enumf mut[DMAX][DMAX];
    enumf centerPartSums[DMAX][DMAX];
    array<enumf, DMAX> rdiag, dist, center, centerPartSum, maxDists, alpha;
    array<enumxt, DMAX> x, dx, ddx, centerLoopBg;
    
    int d, k, kEnd, kMax;
    bool dual;
    uint64_t nodes;
    
    inline bool nextPosUp()
    {
        ++k;
        if (k < kMax)
        {
            x[k] += dx[k];
            ddx[k] = -ddx[k];
            dx[k] = ddx[k] - dx[k];
//            x[k] += dx[k];
            return true;
        }
        else if (k < kEnd)
        {
            kMax = k;
            ++x[k];
            return true;
        }
        return false;
    }

    void prepareEnumeration(enumf maxDist, const vector<enumxt>& subTree, bool solvingSVP);
  
    bool enumerateLoop(enumf& newMaxDist, int& newKMax);

    void enumerate(enumf& maxDist, long normExp, const vector<double>& pruning);
};

#if 0
static const int DMAX = 150;
typedef enumf EnumfVect[DMAX];

class Enumeration {
public:
  static void enumerateDouble(MatGSO<Z_NR<double>, FP_NR<double> >& gso,
            FP_NR<double>& fMaxDist, Evaluator<FP_NR<double> >& evaluator, int first, int last,
            const vector<double>& pruning);

  template<class FT>
  static void enumerate(MatGSO<Integer, FT>& gso, FT& fMaxDist, long maxDistExpo,
               Evaluator<FT>& evaluator, const vector<FT>& targetCoord,
               const vector<FT>& subTree, int first, int last,
               const vector<double>& pruning, bool dual = false);

  static inline long getNodes() { return nodes; }

private:
  static const int DMAX = 150;
  static enumf mut[DMAX][DMAX];
  static enumf centerPartSums[DMAX][DMAX + 1];
  static EnumfVect rdiag, x, dx, ddx, dist, center, centerPartSum, maxDists, centerLoopBg, alpha;
  static int d;
  static int k;              // Current level in the enumeration
  static int kEnd;           // The algorithm stops when k = kEnd
  static int kMax;           // Index of the last non-zero value of x (<= kEnd)
  static bool dual;
  static long nodes;

  // Input: x, dx, ddx, k, kMax, kEnd
  // Output: k, kMax
  static inline bool nextPosUp() {
    //FPLLL_TRACE_IN("k=" << k << " kMax=" << kMax << " kEnd=" << kEnd);
    bool result = true;
    k++;
    if (k < kMax) {
      ddx[k] = -ddx[k];
      dx[k] = ddx[k] - dx[k];
      x[k] += dx[k];
      /*FPLLL_TRACE("x_" << k << "=" << x[k] << ", dx_" << k << "="
              << dx[k] << " ddx_" << k << "=" << ddx[k]);*/
    }
    else if (k < kEnd) {
      kMax = k;
      x[k]++;
    }
    else {
      result = false;
    }
    return result;
  }

  template<class FT>
  static void prepareEnumeration(enumf maxDist, const vector<FT>& subTree, bool solvingSVP);
  
  template<class EnumType>
  static bool enumerateLoop(enumf& newMaxDist, int& newKMax);

  template<class FT>
  static void enumerate(enumf& maxDist, long normExp, Evaluator<FT>& evaluator, const vector<double>& pruning);
  
};

class PrimalEnum {
  public:
    static enumf& center_summand_coefficient(enumf& alpha, enumf& x) {
      return x;
    }
};

class DualEnum {
  public:
    static enumf& center_summand_coefficient(enumf& alpha, enumf& x) {
      return alpha;
    }
};
#endif

FPLLL_END_NAMESPACE

#endif
