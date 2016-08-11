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
                                const vector<FT>& target_coord,
                                const vector<enumxt>& subtree,
                                const vector<enumf>& pruning,
                                bool _dual)
{
    bool solvingsvp = target_coord.empty();
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
            center_partsum[i] = target_coord[i + first].get_d();
    }
    
    FT fr, fmu, fmaxdistnorm;
    long rexpo, normexp = -1;
    for (int i = 0; i < d; ++i)
    {
        fr = _gso.get_r_exp(i + first, i + first, rexpo);
        normexp = max(normexp, rexpo + fr.exponent());
    }
    fmaxdistnorm.mul_2si(fmaxdist, dual ? normexp-fmaxdistexpo : fmaxdistexpo-normexp);
    maxdist = fmaxdistnorm.get_d(GMP_RNDU);

    if (!solvingsvp)
    {
        // For a proper CVP, we need to reset enum below depth with maximal r_i
        
        //find the indices of the maxs
        _max_indices = vector<int>(d);
        int cur, max_index=d, previous_max_index=d;
        auto max_val = _gso.getRExp(d - 1, d - 1);

        for (cur = 0 ; cur < d ; ++cur)
            FPLLL_TRACE("gso[" << cur << "] = "<< _gso.getRExp(cur, cur));

        while (max_index > 0)
        {
            previous_max_index = max_index;
            --max_index;
            for (cur = max_index ; cur >= 0  ; --cur)
            {
                if (max_val <= _gso.getRExp(cur, cur))
                {
                    max_val = _gso.getRExp(cur, cur);
                    max_index = cur;
                }
            }
            for (cur = max_index ; cur < previous_max_index ; ++cur)
                _max_indices[cur] = max_index;
        }
        FPLLL_TRACE("max_indices " << _max_indices);

        //set the correct subbounds
        //pruning_bounds = vector<enumf>(d);  //TODO FIX change between ok or nok to test3
        int i = 0;
        FT cumul_dist = 0;
        enumf cur_dist = 0;
        while (i < d)
        {
            int max = i;
            while (_max_indices[i] == _max_indices[max])
            {
                cumul_dist += _gso.getRExp(i, i);
                ++i;
            }
            for (int j = max ; j < i ; ++j)
            {
                cur_dist = (cumul_dist/fmaxdist).get_d(GMP_RNDU);
            //    pruning_bounds[j] = cur_dist;
            }
        }
        //FPLLL_TRACE("pruning_bounds " << pruning_bounds);
    }

    _evaluator.set_normexp(normexp);

    if (dual)
    {
        for (int i = 0; i < d; ++i)
        {
            fr = _gso.get_r_exp(i + first, i + first, rexpo);
            fr.mul_2si(fr, rexpo - normexp);
            rdiag[d-i-1] = enumf(1.0)/fr.get_d();
        }
        for (int i = 0; i < d; ++i)
        {
            for (int j = i + 1; j < d; ++j)
            {
                _gso.get_mu(fmu, j + first, i + first);
                mut[d-j-1][d-i-1] = -fmu.get_d();
            }
        }
    }
    else
    {
        for (int i = 0; i < d; ++i)
        {
            fr = _gso.get_r_exp(i + first, i + first, rexpo);
            fr.mul_2si(fr, rexpo - normexp);
            rdiag[i] = fr.get_d();
        }
        for (int i = 0; i < d; ++i)
        {
            for (int j = i + 1; j < d; ++j)
            {
                _gso.get_mu(fmu, j + first, i + first);
                mut[i][j] = fmu.get_d();
            }
        }
    }
    subsoldists = rdiag;
    
    prepare_enumeration(subtree, solvingsvp);
    do_enumerate();
  
    fmaxdistnorm = maxdist; // Exact
  
    fmaxdist.mul_2si(fmaxdistnorm, dual ? fmaxdistexpo-normexp : normexp-fmaxdistexpo);
  
    if (dual && !_evaluator.sol_coord.empty())
        reverse_by_swap(_evaluator.sol_coord, 0, d-1);
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
                    newcenter -= x[j] * mut[k][j];  //not (x[j] - center_partsum[j] ???)
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
//    std::cout << "Sol dist: " << newmaxdist << " (nodes:" << nodes << ")" << endl;
    for (int j = 0; j < d; ++j)
        fx[j] = x[j];
    _evaluator.eval_sol(fx, newmaxdist, maxdist);
    
    set_bounds();
}

template<typename FT>
void Enumeration<FT>::process_subsolution(int offset, enumf newdist)
{
    for (int j = 0; j < offset; ++j)
        fx[j] = 0.0;
    for (int j = offset; j < d; ++j)
        fx[j] = x[j];
    _evaluator.eval_sub_sol(k, fx, newdist);
}

template<typename FT>
void Enumeration<FT>::do_enumerate()
{
    nodes = 0;

    save_rounding();
    
    set_bounds();
    
    if      ( dual &&  _evaluator.findsubsols) 
        enumerate_loop<true,true>();
    else if (!dual &&  _evaluator.findsubsols)
        enumerate_loop<false,true>();
    else if ( dual && !_evaluator.findsubsols)
        enumerate_loop<true,false>();
    else if (!dual && !_evaluator.findsubsols)
        enumerate_loop<false,false>();
    
    restore_rounding();        
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
