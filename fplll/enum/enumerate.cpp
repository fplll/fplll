/* Copyright (C) 2008-2011 Xavier Pujol
   (C) 2015 Michael Walter.
   (C) 2016 Marc Stevens. (generic improvements, auxiliary solutions, subsolutions)
   (C) 2016 Guillaume Bonnoron. (CVP improvements)
   
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

template<>
void EnumerationDyn<Float>::reset_rec(enumf cur_dist, int kk)
{
    //FPLLL_TRACE("Reset level " << kk);
    int new_dim = kk+1;

    vector<enumxt> partial_sol(d-kk-1);
    for (int i=kk+1 ; i < d ; ++i)
        partial_sol[i-kk-1] = x[i];

    Float new_dist = 0.0;
    for (int i = 0; i < new_dim; i++)
        new_dist.add(new_dist, _gso.get_r_exp(i, i));

    FastEvaluator<Float> new_evaluator(new_dim, _gso.get_mu_matrix(), _gso.get_r_matrix(), EVALMODE_CV);
    Enumeration<Float> enumobj(_gso, new_evaluator);
    enumobj.enumerate(0, d, new_dist, 0, target, partial_sol, pruning_bounds, false, true);

    if (!new_evaluator.sol_coord.empty())
    {
        Float fsoldistnorm;
        fsoldistnorm.mul_2si(new_evaluator.sol_dist, new_evaluator.normExp - _evaluator.normExp);
        enumf sol_dist = fsoldistnorm.get_d(GMP_RNDU);
        //FPLLL_TRACE("Recovering sub-solution at level: " << kk <<" soldist: " << sol_dist);

        if (sol_dist+cur_dist < partdistbounds[0])
        {
            //FPLLL_TRACE("Saving it.");
            for (int i = 0 ; i < new_dim ; ++i)
                x[i] = new_evaluator.sol_coord[i].get_d();
            process_solution(sol_dist + cur_dist);
        }
    }
}
template<>
void EnumerationDyn<Float>::reset(enumf cur_dist)
{
    //FPLLL_TRACE("Reset level " << k);
    int new_dim = k+1;

    vector<enumxt> partial_sol(d-k-1);
    for (int i=k+1 ; i < d ; ++i)
        partial_sol[i-k-1] = x[i];

    Float new_dist = 0.0;
    for (int i = 0; i < new_dim; i++)
        new_dist.add(new_dist, _gso.get_r_exp(i, i));

    FastEvaluator<Float> new_evaluator(new_dim, _gso.get_mu_matrix(), _gso.get_r_matrix(), EVALMODE_CV);
    Enumeration<Float> enumobj(_gso, new_evaluator);
    enumobj.enumerate(0, d, new_dist, 0, target, partial_sol, pruning_bounds, false, true);

    if (!new_evaluator.sol_coord.empty())
    {
        Float fsoldistnorm;
        fsoldistnorm.mul_2si(new_evaluator.sol_dist, new_evaluator.normExp - _evaluator.normExp);
        enumf sol_dist = fsoldistnorm.get_d(GMP_RNDU);
        //FPLLL_TRACE("Recovering sub-solution at level: " << k <<" soldist: " << sol_dist);

        if (sol_dist+cur_dist < partdistbounds[0])
        {
            //FPLLL_TRACE("Saving it.");
            for (int i = 0 ; i < new_dim ; ++i)
                x[i] = new_evaluator.sol_coord[i].get_d();
            process_solution(sol_dist + cur_dist);
        }
    }
}

template<typename FT>
void EnumerationDyn<FT>::enumerate(int first, int last, FT& fmaxdist, long fmaxdistexpo,
                                const vector<FT>& target_coord,
                                const vector<enumxt>& subtree,
                                const vector<enumf>& pruning,
                                bool _dual,
                                bool subtree_reset)
{
    bool solvingsvp = target_coord.empty();
    dual = _dual;
    pruning_bounds = pruning;
    target = target_coord;
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

#ifdef FPLLL_BIG_ENUM
    memcpy(subsoldists, rdiag, (_gso.d+1)*sizeof(enumf));
#else
    subsoldists = rdiag;
#endif
    
    save_rounding();
    prepare_enumeration(subtree, solvingsvp, subtree_reset);
    do_enumerate();
    restore_rounding();
  
    fmaxdistnorm = maxdist; // Exact
  
    fmaxdist.mul_2si(fmaxdistnorm, dual ? fmaxdistexpo-normexp : normexp-fmaxdistexpo);
  
    if (dual && !_evaluator.sol_coord.empty())
        reverse_by_swap(_evaluator.sol_coord, 0, d-1);
}

template<typename FT>
void EnumerationDyn<FT>::prepare_enumeration(const vector<enumxt>& subtree, bool solvingsvp, bool subtree_reset)
{
    is_svp = solvingsvp;
    
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
                is_svp = false;

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
        if (!subtree_reset || k < k_end)
        {
            alpha[k] = x[k] - newcenter;
            newdist += alpha[k] * alpha[k] * rdiag[k];
        }
    }
    if (!is_svp)
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
void EnumerationDyn<FT>::set_bounds()
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
void EnumerationDyn<FT>::process_solution(enumf newmaxdist)
{
    FPLLL_TRACE("Sol dist: " << newmaxdist << " (nodes:" << nodes << ")");
    for (int j = 0; j < d; ++j)
        fx[j] = x[j];
    _evaluator.eval_sol(fx, newmaxdist, maxdist);
    
    set_bounds();
}

template<typename FT>
void EnumerationDyn<FT>::process_subsolution(int offset, enumf newdist)
{
    for (int j = 0; j < offset; ++j)
        fx[j] = 0.0;
    for (int j = offset; j < d; ++j)
        fx[j] = x[j];
    _evaluator.eval_sub_sol(offset, fx, newdist);
}

template<typename FT>
void EnumerationDyn<FT>::do_enumerate()
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
