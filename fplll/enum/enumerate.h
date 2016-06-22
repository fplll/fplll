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
#include <fplll/gso.h>
#include <fplll/enum/evaluator.h>

FPLLL_BEGIN_NAMESPACE

inline void roundto(int& dest, const double& src) { dest = lrint(src); }
inline void roundto(double& dest, const double& src) { dest = rint(src); }

template<typename FT, int MaxDimension = 256>
class Enumeration
{
public:
    static const int DMAX = MaxDimension;

    Enumeration(MatGSO<Integer, FT>& gso, Evaluator<FT>& evaluator)
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
    array<enumf, DMAX> subsolDists;
    bool findbettersubsols;
    
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
  
    template<bool dualenum, bool dosubsols>
    bool enumerateLoop(enumf& newMaxDist, int& newKMax);

    void enumerate(enumf& maxDist, long normExp, const vector<double>& pruning);
    
};


FPLLL_END_NAMESPACE

#endif
