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

#include "enumerate.h"

FPLLL_BEGIN_NAMESPACE

template<typename FT>
void Enumeration<FT>::enumerate(int first, int last, FT& fmaxdist, long fmaxdistexpo,
                                const vector<FT>& targetcoord,
                                const vector<enumxt>& subtree,
                                const vector<enumf>& pruning,
                                bool _dual)
{
    bool solvingsvp = targetcoord.empty();
    dual = _dual;
    pruning_bounds = pruning;
    if (last == -1)
        last = _gso.d;
    d = last - first;
    fx.resize(d);
    FPLLL_CHECK(d < maxdim, "enumerate: dimension is too high");
    FPLLL_CHECK((solvingsvp || !dual), "CVP for dual not implemented! What does that even mean? ");
    FPLLL_CHECK((subtree.empty() || !dual), "Subtree enumeration for dual not implemented!");

    if (solvingsvp)
    {
        for (int i = 0; i < d; ++i)
            center_partsum[i] = 0.0;
    }
    else
    {
        for (int i = 0; i < d; ++i)
            center_partsum[i] = targetcoord[i + first].get_d();
    }
    
    FT fr, fmu, fmaxdistnorm;
    long rexpo, normexp = -1;
    for (int i = 0; i < d; ++i)
    {
        fr = _gso.getRExp(i + first, i + first, rexpo);
        normexp = max(normexp, rexpo + fr.exponent());
    }
    fmaxdistnorm.mul_2si(fmaxdist, dual ? normexp-fmaxdistexpo : fmaxdistexpo-normexp);
    maxdist = fmaxdistnorm.get_d(GMP_RNDU);

    _evaluator.set_normexp(normexp);

    if (dual)
    {
        for (int i = 0; i < d; ++i)
        {
            fr = _gso.getRExp(i + first, i + first, rexpo);
            fr.mul_2si(fr, rexpo - normexp);
            rdiag[d-i-1] = enumf(1.0)/fr.get_d();
        }
        for (int i = 0; i < d; ++i)
        {
            for (int j = i + 1; j < d; ++j)
            {
                _gso.getMu(fmu, j + first, i + first);
                mut[d-j-1][d-i-1] = -fmu.get_d();
            }
        }
    }
    else
    {
        for (int i = 0; i < d; ++i)
        {
            fr = _gso.getRExp(i + first, i + first, rexpo);
            fr.mul_2si(fr, rexpo - normexp);
            rdiag[i] = fr.get_d();
        }
        for (int i = 0; i < d; ++i)
        {
            for (int j = i + 1; j < d; ++j)
            {
                _gso.getMu(fmu, j + first, i + first);
                mut[i][j] = fmu.get_d();
            }
        }
    }
    subsoldists = rdiag;
    
    prepare_enumeration(subtree, solvingsvp);
    do_enumerate();
  
    fmaxdistnorm = maxdist; // Exact
  
    fmaxdist.mul_2si(fmaxdistnorm, dual ? fmaxdistexpo-normexp : normexp-fmaxdistexpo);
  
    if (dual && !_evaluator.solCoord.empty())
        reverseBySwap(_evaluator.solCoord, 0, d-1);
}

template<typename FT>
void Enumeration<FT>::prepare_enumeration(const vector<enumxt>& subtree, bool solvingsvp)
{
    bool svpbeginning = solvingsvp;
    
    enumf newdist = 0.0;
    k_end = d - subtree.size();
    for (k = d - 1; k >= 0 && newdist <= maxdist; --k)
    {
        enumf newcenter = center_partsum[k];
        if (k >= k_end)
        {
            // use subtree
            x[k] = subtree[k - k_end];
            if (x[k] != 0)
                svpbeginning = false;

            for (int j = 0; j < k; ++j)
                center_partsum[j] -= x[k] * mut[j][k];
        }
        else
        {
            if (dual)
            {
                for (int j = k + 1; j < k_end; ++j)
                    newcenter -= alpha[j] * mut[k][j];
            }
            else
            {
                for (int j = k + 1; j < k_end; ++j)
                    newcenter -= x[j] * mut[k][j];
            }
            roundto(x[k], newcenter); // newX = rint(newCenter) / lrint(newCenter)
            center[k] = newcenter;
            partdist[k] = newdist;
            dx[k] = ddx[k] = (((int)(newcenter >= x[k]) & 1) << 1) - 1;
        }
        alpha[k] = x[k] - newcenter;
        newdist += alpha[k] * alpha[k] * rdiag[k];
    }
    if (!svpbeginning) 
    {
        k_max = k_end; // The last non-zero coordinate of x will not stay positive
    }
    else 
    {
        k_max = 0;
        x[0] = 1; // Excludes (0,...,0) from the enumeration
    }
    ++k;
}

template<typename FT>
void Enumeration<FT>::set_bounds()
{
    if (pruning_bounds.empty())
    {
        fill(&partdistbounds[0] + 0, &partdistbounds[0] + d, maxdist);
    }
    else
    {
        for (int i = 0; i < d; ++i)
            partdistbounds[i] = pruning_bounds[i] * maxdist;
    }
}

template<typename FT>
void Enumeration<FT>::process_solution(enumf newmaxdist)
{
    for (int j = 0; j < d; ++j)
        fx[j] = x[j];
    _evaluator.evalSol(fx, newmaxdist, maxdist);
    
    set_bounds();
}

template<typename FT>
void Enumeration<FT>::process_subsolution(int offset, enumf newdist)
{
    for (int j = 0; j < offset; ++j)
        fx[j] = 0;
    for (int j = offset; j < d; ++j)
        fx[j] = x[j];
    _evaluator.evalSubSol(k, fx, newdist);
}

template<typename FT>
void Enumeration<FT>::do_enumerate()
{
    nodes = 0;

    set_bounds();
    
    if      ( dual &&  _evaluator.findsubsols) 
        enumerate_loop<true,true>();
    else if (!dual &&  _evaluator.findsubsols)
        enumerate_loop<false,true>();
    else if ( dual && !_evaluator.findsubsols)
        enumerate_loop<true,false>();
    else if (!dual && !_evaluator.findsubsols)
        enumerate_loop<false,false>();
            
}

template class Enumeration<FP_NR<double> >;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class Enumeration<FP_NR<long double> >;
#endif

#ifdef FPLLL_WITH_QD
template class Enumeration<FP_NR<dd_real> >;

template class Enumeration<FP_NR<qd_real> >;
#endif

#ifdef FPLLL_WITH_DPE
template class Enumeration<FP_NR<dpe_t> >;
#endif

template class Enumeration<FP_NR<mpfr_t> >;

FPLLL_END_NAMESPACE
