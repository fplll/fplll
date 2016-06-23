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

#include "../nr/nr.h"

FPLLL_BEGIN_NAMESPACE

inline void roundto(int& dest, const double& src) { dest = lrint(src); }
inline void roundto(double& dest, const double& src) { dest = rint(src); }

#define MAXDIMENSION           256
#define MAXTEMPLATEDDIMENSION  80

class enumeration_base
{
public:
    static const int maxdim = MAXDIMENSION;
    static const int maxdim_opt = MAXTEMPLATEDDIMENSION;

    inline uint64_t getNodes() const { return nodes; }

protected:
    enumf mut[maxdim][maxdim];
    enumf centerPartSums[maxdim][maxdim];
    array<enumf, maxdim> rdiag, dist, center, centerPartSum, maxDists, alpha;
    array<enumxt, maxdim> x, dx, ddx, centerLoopBg;
    
    int d, k, kEnd, kMax;
    bool dual;
    uint64_t nodes;
    
    template<bool dualenum>
    bool enumerateLoop(enumf& newMaxDist, int& newKMax);
    
};


FPLLL_END_NAMESPACE

#endif
