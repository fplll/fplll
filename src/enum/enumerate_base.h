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
    /* configuration */
    bool dual, subsols;

    /* enumeration input */
    enumf mut[maxdim][maxdim];
    array<enumf, maxdim> rdiag, partdistbounds;
    int d, kEnd; // dimension, subtreelevel

    /* partial sum cache */
    enumf center_partsums[maxdim][maxdim];
    array<enumf, maxdim> center_partsum;
    array<int,   maxdim> center_partsum_begin;

    /* enumeration data for each level */
    array<enumf, maxdim> partdist, center, alpha;
    array<enumxt,maxdim> x, dx, ddx;

    int k, kMax;
   
    /* nodes count */
    uint64_t nodes;
    
    template<bool dualenum, bool findsubsols>
    bool enumerateLoop(enumf& newMaxDist, int& newKMax);

    virtual process_solution() = 0;
    virtual process_subsolution(int offset) = 0;
};


FPLLL_END_NAMESPACE

#endif
