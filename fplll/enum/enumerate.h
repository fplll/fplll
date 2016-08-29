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

#ifndef FPLLL_ENUMERATE_H
#define FPLLL_ENUMERATE_H

#include <array>
#include <fplll/gso.h>
#include <fplll/enum/evaluator.h>
#include <fplll/enum/enumerate_base.h>

FPLLL_BEGIN_NAMESPACE

template<typename FT>
class Enumeration : public EnumerationBase
{
public:
#ifdef FPLLL_BIG_ENUM
    ~Enumeration()
    {
        for (int i=0 ; i < _gso.d+1 ; ++i)
        {
            delete[] mut[i];
            delete[] center_partsums[i];
        }
        delete[] mut;
        delete[] rdiag;
        delete[] partdistbounds;
        delete[] center_partsums;
        delete[] center_partsum;
        delete[] center_partsum_begin;
        delete[] partdist;
        delete[] center;
        delete[] alpha;
        delete[] x;
        delete[] dx;
        delete[] ddx;
        delete[] subsoldists;
    }
#endif
    Enumeration(MatGSO<Integer, FT>& gso, Evaluator<FT>& evaluator, vector<int> max_indices=vector<int>())
        : _gso(gso), _evaluator(evaluator)
    {
        _max_indices = max_indices;
#ifdef FPLLL_BIG_ENUM
        int dim = gso.d+1;
        mut = new enumf*[dim];
        for (int i=0 ; i<dim ; ++i)
            mut[i] = new enumf[dim];

        rdiag = new enumf[dim];
        partdistbounds = new enumf[dim];

        center_partsums = new enumf*[dim];
        for (int i=0 ; i<dim ; ++i)
            center_partsums[i] = new enumf[dim];

        center_partsum = new enumf[dim];
        center_partsum_begin = new int[dim];

        partdist = new enumf[dim];
        center= new enumf[dim];
        alpha = new enumf[dim];

        x = new enumxt[dim];
        dx = new enumxt[dim];
        ddx = new enumxt[dim];

        subsoldists = new enumf[dim];
#endif
    }

    void enumerate(int first, int last,
                FT& fmaxdist, long fmaxdistexpo, 
                const vector<FT>& target_coord = vector<FT>(),
                const vector<enumxt>& subtree = vector<enumxt>(),
                const vector<enumf>& pruning = vector<enumf>(),
                bool dual = false,
                bool subtree_reset = false);

    inline uint64_t get_nodes() const { return nodes; }
    
private:
    MatGSO<Integer, FT>& _gso; 
    Evaluator<FT>& _evaluator;
    vector<FT> target;

    vector<enumf> pruning_bounds;
    enumf maxdist;
    vector<FT> fx;
    
    void prepare_enumeration(const vector<enumxt>& subtree, bool solvingsvp, bool subtree_reset);
  
    void do_enumerate();

    void set_bounds();    
    void reset(enumf cur_dist) {};
    void reset_rec(enumf cur_dist, int kk) {};
    virtual void process_solution(enumf newmaxdist);
    virtual void process_subsolution(int offset, enumf newdist);
    
};


FPLLL_END_NAMESPACE

#endif
