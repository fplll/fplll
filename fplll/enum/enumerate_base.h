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
#include <cfenv>
#include <cmath>
#include "../fplll_config.h"
#include "../nr/nr.h"

FPLLL_BEGIN_NAMESPACE

inline void roundto(int& dest, const double& src) { dest = std::lrint(src); }
inline void roundto(double& dest, const double& src) { dest = std::rint(src); }

/* config */
#define FPLLL_WITH_RECURSIVE_ENUM 1
#define MAXTEMPLATEDDIMENSION  80 // unused
//#define FORCE_ENUM_INLINE // not recommended
/* end config */


#ifndef __has_attribute
#define __has_attribute(x) 0  // Compatibility with non - GCC/clang compilers.
#endif
#if __has_attribute(always_inline)
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define ALWAYS_INLINE
#endif

#ifndef FORCE_ENUM_INLINE
#define ENUM_ALWAYS_INLINE
#else
#define ENUM_ALWAYS_INLINE ALWAYS_INLINE
#endif

class EnumerationBase
{
public:
    static const int maxdim = FPLLL_MAX_ENUM_DIMENSION;
    
    inline uint64_t get_nodes() const { return nodes; }
    virtual ~EnumerationBase() {}

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

    template<int kk, int kk_start, bool dualenum, bool findsubsols>
    struct opts
    {};
    
    /* need templated function argument for support of integer specialization for kk==-1 */
    template<int kk, int kk_start, bool dualenum, bool findsubsols>
    inline void enumerate_recursive( opts<kk, kk_start, dualenum, findsubsols> ) ENUM_ALWAYS_INLINE;
    template<int kk_start, bool dualenum, bool findsubsols>
    inline void enumerate_recursive( opts<-1, kk_start, dualenum, findsubsols> ) 
    {
    }

    /* simple wrapper with no function argument as helper for dispatcher */
    template<int kk, bool dualenum, bool findsubsols>
    void enumerate_recursive_wrapper()
    {
        enumerate_recursive( opts<(kk<maxdim ? kk : -1),0,dualenum,findsubsols>() );
    }
        
    template<bool dualenum, bool findsubsols>
    inline void enumerate_recursive_dispatch(int kk);
        
    template<bool dualenum, bool findsubsols>
    void enumerate_loop();

    virtual void process_solution(enumf newmaxdist) = 0; 
    virtual void process_subsolution(int offset, enumf newdist) = 0;
    
    int rounding_backup;
    void save_rounding()
    {
        rounding_backup = std::fegetround();
        std::fesetround(FE_TONEAREST);
    }
    void restore_rounding()
    {
        std::fesetround(rounding_backup);
    }
    
    inline bool next_pos_up()
    {
        ++k;
        if (partdist[k] != 0.0)
        {
            x[k] += dx[k];
            ddx[k] = -ddx[k];
            dx[k] = ddx[k] - dx[k];
        }
        else
        {
            if (k >= k_end)
                return false;
            k_max = k;
            ++x[k];
        }
        return true;
    }
};


FPLLL_END_NAMESPACE

#endif
