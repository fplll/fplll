/* Copyright (C) 2008-2011 Xavier Pujol
   (C) 2015 Michael Walter.
   (C) 2016 Marc Stevens. (generic improvements, auxiliary solutions, subsolutions)
   
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

#ifndef FPLLL_ENUMERATE_BASE_H
#define FPLLL_ENUMERATE_BASE_H

#include <array>
#include "../nr/nr.h"

FPLLL_BEGIN_NAMESPACE

inline void roundto(int& dest, const double& src) { dest = lrint(src); }
inline void roundto(double& dest, const double& src) { dest = rint(src); }

#define MAXDIMENSION           256
#define MAXTEMPLATEDDIMENSION  80

class EnumerationBase
{
public:
    static const int maxdim = MAXDIMENSION;
    static const int maxdim_opt = MAXTEMPLATEDDIMENSION;
    
    inline uint64_t getNodes() const { return nodes; }
    
protected:
    /* configuration */
    bool dual;

    /* enumeration input */
    enumf mut[maxdim][maxdim];
    array<enumf, maxdim> rdiag, partdistbounds;
    int d, k_end; // dimension, subtreelevel

    /* partial sum cache */
    enumf center_partsums[maxdim][maxdim];
    array<enumf, maxdim> center_partsum;
    array<int,   maxdim> center_partsum_begin;

    /* enumeration data for each level */
    array<enumf, maxdim> partdist, center, alpha;
    array<enumxt,maxdim> x, dx, ddx;
    array<enumf, maxdim> subsoldists;
    
    int k, k_max;
   
    /* nodes count */
    uint64_t nodes;
    
    template<bool dualenum, bool findsubsols>
    bool enumerate_loop();

    virtual void process_solution(enumf newmaxdist) = 0; 
    virtual void process_subsolution(int offset, enumf newdist) = 0;
    
    inline bool next_pos_up()
    {
        ++k;
        if (k >= k_end)
            return false;
        if (k < k_max)
        {
            x[k] += dx[k];
            ddx[k] = -ddx[k];
            dx[k] = ddx[k] - dx[k];
        }
        else
        {
            k_max = k;
            ++x[k];
        }
        return true;
    }
};


FPLLL_END_NAMESPACE

#endif
