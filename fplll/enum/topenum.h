/* Copyright (C) 2008-2011 Xavier Pujol.

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

#ifndef FPLLL_TOP_ENUM_H
#define FPLLL_TOP_ENUM_H

#include "../util.h"

FPLLL_BEGIN_NAMESPACE

const double ENUM_MAX_VOLUME = 20000000;
const int ENUM_MIN_LEVEL     = 20;

class Enumerator
{
public:
  Enumerator(int d, const Matrix<Float> &mu, const Matrix<Float> &r,
             double max_volume = ENUM_MAX_VOLUME, int min_level = ENUM_MIN_LEVEL);
  bool enum_next(const Float &max_sqr_length);
  inline const vector<enumxt> &get_sub_tree() { return sub_tree; }
private:
  const Matrix<Float> &mu;
  const Matrix<Float> &r;
  int k, kmin, kmax, d;
  FloatVect center, dist;
  FloatVect x, dx, ddx;
  //  FloatVect sub_tree;
  vector<enumxt> sub_tree;
  Float max_volume;
  bool svp_init_needed;
};

FPLLL_END_NAMESPACE

#endif
