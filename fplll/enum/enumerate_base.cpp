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

template<int kk, bool dualenum, bool findsubsols>
inline void EnumerationBase::enumerate_recursive( EnumerationBase::opts<kk, dualenum, findsubsols> )
{

    enumf newcenter = center_partsums[kk][kk+1];
    enumf newx;
    roundto(newx, newcenter);
    enumf alphak = newx - newcenter;
    enumf newdist = partdist[kk] + alphak*alphak*rdiag[kk];

    if (newdist > partdistbounds[kk])
        return;

    center[kk] = newcenter;
    x[kk] = newx;
    alpha[kk] = alphak;
    dx[kk] = ddx[kk] = (((int)(newcenter >= newx) & 1) << 1) - 1;
    
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
            for (int j = center_partsum_begin[kk-1+1]; j > kk-1; --j)
                center_partsums[kk-1][j] = center_partsums[kk-1][j+1] - alpha[j] * mut[kk-1][j];
        }
        else
        {
            for (int j = center_partsum_begin[kk-1+1]; j > kk-1; --j)
                center_partsums[kk-1][j] = center_partsums[kk-1][j+1] - x[j] * mut[kk-1][j];
        }
        if (center_partsum_begin[kk-1+1] > center_partsum_begin[kk-1])
            center_partsum_begin[kk-1] = center_partsum_begin[kk-1+1];
        center_partsum_begin[kk-1+1] = kk-1+1;
    }

    if (findsubsols && newdist < subsoldists[kk])
    {
        subsoldists[kk] = newdist;
        process_subsolution(kk, newdist);
    }
    
    while (true)
    {
        enumerate_recursive( opts<kk-1,dualenum,findsubsols>() );

        if (partdist[kk] != 0.0)
        {
            x[kk] += dx[kk];
            ddx[kk] = -ddx[kk];
            dx[kk] = ddx[kk] - dx[kk];
            enumf alphak2 = x[kk] - center[kk];
            enumf newdist2 = partdist[kk] + alphak2*alphak2*rdiag[kk];
            if (newdist2 > partdistbounds[kk])
                return;
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
                center_partsum_begin[kk-1] = kk-1+1;
            }
        }
        else
        {
            ++x[kk];
            enumf alphak2 = x[kk] - center[kk];
            enumf newdist2 = partdist[kk] + alphak2*alphak2*rdiag[kk];
            if (newdist2 > partdistbounds[kk])
                return;
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
                center_partsum_begin[kk-1] = kk-1+1;
            }
        }
    }
}

template<bool dualenum, bool findsubsols>
inline void EnumerationBase::enumerate_recursive( EnumerationBase::opts<-1, dualenum, findsubsols> )
{
}

template<bool dualenum, bool findsubsols>
inline void EnumerationBase::enumerate_recursive_dispatch(int kk)
{
//    typedef int (Foo::*Member)(int, int)
//    static const 
}


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
    
    if (k_end==52 && 0)
    {
        cout << k_end << endl;
        enumerate_recursive( opts<51, dualenum, findsubsols>() );
        k=52;
        return;
    }
    
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
            enumf newcenter = center_partsums[k][k+1];
            center_partsum_begin[k] = max(center_partsum_begin[k], center_partsum_begin[k+1]);
            center_partsum_begin[k+1] = k+1;
            
            center[k] = newcenter;
            partdist[k] = newdist;
            roundto(x[k], newcenter);
            dx[k] = ddx[k] = (((int)(newcenter >= x[k]) & 1) << 1) - 1;
        }
        else
        {
            if (!next_pos_up())
                break;
        }
    }
}

template void EnumerationBase::enumerate_loop<false,false>();
template void EnumerationBase::enumerate_loop<false,true>();
template void EnumerationBase::enumerate_loop<true,false>();
template void EnumerationBase::enumerate_loop<true,true>();

FPLLL_END_NAMESPACE
