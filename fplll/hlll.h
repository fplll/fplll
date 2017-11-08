/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.

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

#ifndef FPLLL_HLLL_H
#define FPLLL_HLLL_H

#include "householder.h"
#include <cmath>

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT> class HLLLReduction
{
public:
  /**
   * Constructor.
   * The precision of FT must be defined before creating an instance of the
   * class and must remain the same until the object is destroyed (or no longer
   * needed).
   */
  // TODO: what is c?
  HLLLReduction(MatHouseholder<ZT, FT> &arg_m, double delta, double eta, double theta, double c,
                int flags)
      : m(arg_m)
  {
    this->delta = delta;
    this->eta   = eta;
    this->theta = theta;
    this->c     = c;
    double tmp  = pow(2.0, -(double)m.get_d() * c);
    sr          = tmp;
    verbose     = flags & LLL_VERBOSE;
  }

  /**
    @brief LLL reduction.
    */

  void lll();

private:
  FT delta, eta, theta;
  MatHouseholder<ZT, FT> &m;
  // Arbitraty c > 0
  FT c;
  // Multiplicative coefficient used to check if a vector is size-reduced or not.
  FT sr;
  // Verbose mode.
  bool verbose;

  void size_reduction(int k);

  inline void print_params();
};

template <class ZT, class FT> inline void HLLLReduction<ZT, FT>::print_params()
{
  cerr << "Entering HLLL"
       << "\ndelta = " << delta << "\neta = " << eta << "\ntheta = " << theta << "\nc = " << c
       << "\nprecision = " << FT::get_prec() << endl;
}

template <class ZT, class FT>
bool is_hlll_reduced(MatHouseholder<ZT, FT> &m, double delta, double eta);

FPLLL_END_NAMESPACE

#endif /* FPLLL_HLLL_H */
