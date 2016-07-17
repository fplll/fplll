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

Enumerator::Enumerator(int d, const Matrix<Float> &mu, const Matrix<Float> &r, double argMaxVolume,
                       int min_level)
    : mu(mu), r(r), kmin(min_level), d(d)
{
  max_volume = argMaxVolume > 0 ? argMaxVolume : ENUM_MAX_VOLUME;
  gen_zero_vect(center, d);
  gen_zero_vect(dist, d);
  gen_zero_vect(x, d);
  gen_zero_vect(dx, d);
  gen_zero_vect(ddx, d);
  svp_init_needed = true;
}

bool Enumerator::enum_next(const Float &max_sqr_length)
{
  Float newdist, newcenter, y, volume, rtmp1;
  bool notFound = true;

  if (svp_init_needed)
  {
    for (k = d - 1; k > kmin; k--)
    {
      cost_estimate(volume, max_sqr_length, r, k - 1);
      if (volume <= max_volume)
        break;
    }
    kmax          = k;
    svp_init_needed = false;
  }
  if (k >= d)
    return false;

  while (notFound)
  {
    // FPLLL_TRACE("Level k=" << k << " dist_k=" << dist[k] << " x_k=" << x[k]);
    y.sub(center[k], x[k]);
    newdist.mul(y, y);
    newdist.mul(newdist, r(k, k));
    newdist.add(newdist, dist[k]);

    if (newdist <= max_sqr_length)
    {
      rtmp1.sub(max_sqr_length, newdist);
      cost_estimate(volume, rtmp1, r, k - 1);
      if (k > kmin && volume >= max_volume)
      {
        k--;
        // FPLLL_TRACE("  Go down, newdist=" << newdist);

        newcenter = 0.0;
        for (int j = d - 1; j > k; j--)
          newcenter.submul(x[j], mu(j, k));

        center[k] = newcenter;
        dist[k]   = newdist;
        x[k].rnd(newcenter);
        dx[k]  = 0.0;
        ddx[k] = newcenter >= x[k] ? -1.0 : 1.0;
        continue;
      }
      sub_tree.resize(d - k);
      for (size_t j = 0; j < sub_tree.size(); j++)
        sub_tree[j]  = enumxt(x[j + k].get_d());
      // FPLLL_TRACE("  SubTree approx_size=" << volume << " coord=" << sub_tree);
      notFound = false;
    }
    else
    {
      // FPLLL_TRACE("  Go up");
      k++;
    }
    if (k < kmax)
    {
      ddx[k].neg(ddx[k]);
      dx[k].sub(ddx[k], dx[k]);
      x[k].add(x[k], dx[k]);
    }
    else
    {
      if (k >= d)
        break;
      kmax  = k;
      rtmp1 = 1.0;
      x[k].add(x[k], rtmp1);
    }
    // FPLLL_TRACE("  x[" << k << "]=" << x[k]);
  }
  return !notFound;
}

FPLLL_END_NAMESPACE
