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

#include "topenum.h"

FPLLL_BEGIN_NAMESPACE

Enumerator::Enumerator(int d, const Matrix<Float>& mu, const Matrix<Float>& r,
                       double argMaxVolume, int minLevel) :
  mu(mu), r(r), kmin(minLevel), d(d)
{
  maxVolume = argMaxVolume > 0 ? argMaxVolume : ENUM_MAX_VOLUME;
  genZeroVect(center, d);
  genZeroVect(dist, d);
  genZeroVect(x, d);
  genZeroVect(dx, d);
  genZeroVect(ddx, d);
  svpInitNeeded = true;
}

bool Enumerator::enumNext(const Float& maxsqrlength) {
  Float newdist, newcenter, y, volume, rtmp1;
  bool notFound = true;

  if (svpInitNeeded) {
    for (k = d - 1; k > kmin; k--) {
      costEstimate(volume, maxsqrlength, r, k - 1);
      if (volume <= maxVolume) break;
    }
    kmax = k;
    svpInitNeeded = false;
  }
  if (k >= d) return false;

  while (notFound) {
    //FPLLL_TRACE("Level k=" << k << " dist_k=" << dist[k] << " x_k=" << x[k]);
    y.sub(center[k], x[k]);
    newdist.mul(y, y);
    newdist.mul(newdist, r(k, k));
    newdist.add(newdist, dist[k]);

    if (newdist <= maxsqrlength) {
      rtmp1.sub(maxsqrlength, newdist);
      costEstimate(volume, rtmp1, r, k - 1);
      if (k > kmin && volume >= maxVolume) {
        k--;
        //FPLLL_TRACE("  Go down, newdist=" << newdist);

        newcenter = 0.0;
        for (int j = d - 1; j > k; j--)
          newcenter.submul(x[j], mu(j, k));

        center[k] = newcenter;
        dist[k] = newdist;
        x[k].rnd(newcenter);
        dx[k] = 0.0;
        ddx[k] = newcenter >= x[k] ? -1.0 : 1.0;
        continue;
      }
      subTree.resize(d - k);
      for (size_t j = 0; j < subTree.size(); j++)
        subTree[j] = x[j + k];
      //FPLLL_TRACE("  SubTree approx_size=" << volume << " coord=" << subTree);
      notFound = false;
    }
    else {
      //FPLLL_TRACE("  Go up");
      k++;
    }
    if (k < kmax) {
      ddx[k].neg(ddx[k]);
      dx[k].sub(ddx[k], dx[k]);
      x[k].add(x[k], dx[k]);
    }
    else {
      if (k >= d) break;
      kmax = k;
      rtmp1 = 1.0;
      x[k].add(x[k], rtmp1);
    }
    //FPLLL_TRACE("  x[" << k << "]=" << x[k]);
  }
  return !notFound;
}

FPLLL_END_NAMESPACE
