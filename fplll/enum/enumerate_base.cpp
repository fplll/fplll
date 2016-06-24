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

template<bool dualenum, bool findsubsols>
bool EnumerationBase::enumerate_loop()
{
    if (k >= k_end)
        return false;
    
    for (int i = 0; i < k_end; ++i)
    {
        center_partsum_begin[i+1] = k_end - 1;
        center_partsums[i][k_end] = center_partsum[i];
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
                    return false;
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
                return false;
        }
    }
}

template bool EnumerationBase::enumerate_loop<false,false>();
template bool EnumerationBase::enumerate_loop<false,true>();
template bool EnumerationBase::enumerate_loop<true,false>();
template bool EnumerationBase::enumerate_loop<true,true>();

FPLLL_END_NAMESPACE
