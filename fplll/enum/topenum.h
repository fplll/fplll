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
const int ENUM_MIN_LEVEL = 20;

class Enumerator {
public:
  Enumerator(int d, const Matrix<Float>& mu, const Matrix<Float>& r,
             double maxVolume = ENUM_MAX_VOLUME,
             int minLevel = ENUM_MIN_LEVEL);
  bool enum_next(const Float& maxsqrlength);
  inline const vector<enumxt>& get_sub_tree() {
    return subTree;
  }
private:
  const Matrix<Float>& mu;
  const Matrix<Float>& r;
  int k, kmin, kmax, d;
  FloatVect center, dist;
  FloatVect x, dx, ddx;
//  FloatVect subTree;
  vector<enumxt> subTree;
  Float maxVolume;
  bool svpInitNeeded;
};

FPLLL_END_NAMESPACE

#endif
