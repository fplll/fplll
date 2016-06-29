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

#include "enumerate_base.h"

FPLLL_BEGIN_NAMESPACE

#ifdef USE_RECURSIVE_ENUM
template<int kk, int kk_start, bool dualenum, bool findsubsols>
inline void EnumerationBase::enumerate_recursive( EnumerationBase::opts<kk, kk_start, dualenum, findsubsols> )
{
    enumf alphak = x[kk] - center[kk];
    enumf newdist = partdist[kk] + alphak*alphak*rdiag[kk];

    if (!(newdist <= partdistbounds[kk]))
        return;
    ++nodes;
    
    alpha[kk] = alphak;
    if (findsubsols && newdist < subsoldists[kk])
    {
        subsoldists[kk] = newdist;
        process_subsolution(kk, newdist);
    }
    
    if (kk == 0)
    {
        if (newdist > 0.0 /* || is_cvp */)
            process_solution(newdist);
    }
    else
    {
        partdist[kk-1] = newdist;
        if (dualenum)
        {
            for (int j = center_partsum_begin[kk]; j > kk-1; --j)
                center_partsums[kk-1][j] = center_partsums[kk-1][j+1] - alpha[j] * mut[kk-1][j];
        }
        else
        {
            for (int j = center_partsum_begin[kk]; j > kk-1; --j)
                center_partsums[kk-1][j] = center_partsums[kk-1][j+1] - x[j] * mut[kk-1][j];
        }
        if (center_partsum_begin[kk] > center_partsum_begin[kk-1])
            center_partsum_begin[kk-1] = center_partsum_begin[kk];
        center_partsum_begin[kk] = kk;
        center[kk-1] = center_partsums[kk-1][kk];    
        roundto(x[kk-1], center[kk-1]);
        dx[kk-1] = ddx[kk-1] = (((int)(center[kk-1] >= x[kk-1]) & 1) << 1) - 1;
    }

    
    while (true)
    {
        enumerate_recursive( opts<kk-1,kk_start,dualenum,findsubsols>() );

        if (partdist[kk] != 0.0)
        {
            x[kk] += dx[kk];
            ddx[kk] = -ddx[kk];
            dx[kk] = ddx[kk] - dx[kk];

            enumf alphak2 = x[kk] - center[kk];
            enumf newdist2 = partdist[kk] + alphak2*alphak2*rdiag[kk];
            if (!(newdist2 <= partdistbounds[kk]))
                return;
            ++nodes;
            alpha[kk] = alphak2;
            if (kk == 0)
            {
                if (newdist2 > 0.0 /* || is_cvp */)
                    process_solution(newdist2);
            }
            else
            {
                partdist[kk-1] = newdist2;
                if (dualenum)
                    center_partsums[kk-1][kk-1+1] = center_partsums[kk-1][kk-1+1+1] - alpha[kk-1+1] * mut[kk-1][kk-1+1];
                else
                    center_partsums[kk-1][kk-1+1] = center_partsums[kk-1][kk-1+1+1] - x[kk-1+1] * mut[kk-1][kk-1+1];
                if (kk > center_partsum_begin[kk-1])
                    center_partsum_begin[kk-1] = kk;
                center[kk-1] = center_partsums[kk-1][kk-1+1];    
                roundto(x[kk-1], center[kk-1]);
                dx[kk-1] = ddx[kk-1] = (((int)(center[kk-1] >= x[kk-1]) & 1) << 1) - 1;
            }
        }
        else
        {
            ++x[kk];

            enumf alphak2 = x[kk] - center[kk];
            enumf newdist2 = partdist[kk] + alphak2*alphak2*rdiag[kk];
            if (!(newdist2 <= partdistbounds[kk]))
                return;
            ++nodes;
            alpha[kk] = alphak2;
            if (kk == 0)
            {
                if (newdist2 > 0.0 /* || is_cvp */)
                    process_solution(newdist2);
            }
            else
            {
                partdist[kk-1] = newdist2;
                if (dualenum)
                    center_partsums[kk-1][kk-1+1] = center_partsums[kk-1][kk-1+1+1] - alpha[kk-1+1] * mut[kk-1][kk-1+1];
                else
                    center_partsums[kk-1][kk-1+1] = center_partsums[kk-1][kk-1+1+1] - x[kk-1+1] * mut[kk-1][kk-1+1];
                if (kk > center_partsum_begin[kk-1])
                    center_partsum_begin[kk-1] = kk;
                center[kk-1] = center_partsums[kk-1][kk-1+1];    
                roundto(x[kk-1], center[kk-1]);
                dx[kk-1] = ddx[kk-1] = (((int)(center[kk-1] >= x[kk-1]) & 1) << 1) - 1;
            }
        }
    }
}

template<int kk_start, bool dualenum, bool findsubsols>
inline void EnumerationBase::enumerate_recursive( EnumerationBase::opts<-1, kk_start, dualenum, findsubsols> )
{
}

template<bool dualenum, bool findsubsols>
inline void EnumerationBase::enumerate_recursive_dispatch(int kk)
{
    typedef void (EnumerationBase::*enum_recur_type)();
    static const enum_recur_type lookup[] = 
        { 
            & EnumerationBase::enumerate_recursive_wrapper<  7,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 15,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 23,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 31,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 39,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 47,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 55,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 63,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 71,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 79,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 87,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper< 95,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper<103,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper<111,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper<119,dualenum,findsubsols>,
            & EnumerationBase::enumerate_recursive_wrapper<127,dualenum,findsubsols>
        };
    (this->*lookup[kk/8])();
}
#endif

template<bool dualenum, bool findsubsols>
void EnumerationBase::enumerate_loop()
{
    if (k >= k_end)
        return;
    
    for (int i = 0; i < k_end; ++i)
    {
        center_partsum_begin[i+1] = k_end - 1;
        center_partsums[i][k_end] = center_partsum[i];
    }
    nodes -= (k_end - 1 - k) + 1;
    k = k_end - 1;

#ifdef USE_RECURSIVE_ENUM
    if ((k & 7) == 7 && k <= 127)
    {
        enumerate_recursive_dispatch<dualenum, findsubsols>(k);
        if (!next_pos_up())
        {
//            std::cout << "Nodes: " << nodes << std::endl;
            return;
        }
    }
#endif
    
    while (true)
    {
        enumf alphak = x[k] - center[k];
        enumf newdist = partdist[k] + alphak * alphak * rdiag[k];
        if (newdist <= partdistbounds[k])
        {
            ++nodes;
            alpha[k] = alphak;
            if (findsubsols && newdist < subsoldists[k])
            {
                subsoldists[k] = newdist;
                process_subsolution(k, newdist);
            }
            --k;
            if (k < 0)
            {
                if (newdist > 0.0)
                    process_solution(newdist);
                if (!next_pos_up())
                    break;
                continue;
            }
            if (dualenum)
            {
                for (int j = center_partsum_begin[k+1]; j > k; --j)
                    center_partsums[k][j] = center_partsums[k][j+1] - alpha[j] * mut[k][j];
            }
            else
            {
                for (int j = center_partsum_begin[k+1]; j > k; --j)
                    center_partsums[k][j] = center_partsums[k][j+1] - x[j] * mut[k][j];
            }
            center_partsum_begin[k] = max(center_partsum_begin[k], center_partsum_begin[k+1]);
            center_partsum_begin[k+1] = k+1;
            
            enumf newcenter = center_partsums[k][k+1];
            center[k] = newcenter;
            partdist[k] = newdist;
            roundto(x[k], newcenter);
            dx[k] = ddx[k] = (((int)(newcenter >= x[k]) & 1) << 1) - 1;

#ifdef USE_RECURSIVE_ENUM
            if ((k & 7) == 7 && k <= 127)
            {
                enumerate_recursive_dispatch<dualenum, findsubsols>(k);
                if (!next_pos_up())
                    break;
            }
#endif

        }
        else
        {
            if (!next_pos_up())
                break;
        }
    }
//    std::cout << "Nodes: " << nodes << std::endl;
}

template void EnumerationBase::enumerate_loop<false,false>();
template void EnumerationBase::enumerate_loop<false,true>();
template void EnumerationBase::enumerate_loop<true,false>();
template void EnumerationBase::enumerate_loop<true,true>();

FPLLL_END_NAMESPACE
