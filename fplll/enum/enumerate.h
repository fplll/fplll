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
    Enumeration(MatGSO<Integer, FT>& gso, Evaluator<FT>& evaluator, vector<int> max_indices=vector<int>())
        : _gso(gso), _evaluator(evaluator)
    {
        _max_indices = max_indices;
    }

    void reset(enumf cur_dist) {};

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
    virtual void process_solution(enumf newmaxdist);
    virtual void process_subsolution(int offset, enumf newdist);
    
};


FPLLL_END_NAMESPACE

#endif
