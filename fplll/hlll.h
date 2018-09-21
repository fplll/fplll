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
    dR.resize(m.get_d());
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

  /**
     @brief Size reduction.

     Perform size reduction of vector kappa.

     @param kappa index of the vector
  */
  void size_reduction(int kappa = 0);

  inline void print_params();

  // Precompute delta_ * R(k, k)^2
  vector<FT> dR;

  // Compute dR[k]
  inline void compute_dR(int k, FT delta_);
};

template <class ZT, class FT> inline void HLLLReduction<ZT, FT>::print_params()
{
  cerr << "Entering HLLL" << endl
       << "delta = " << delta << endl
       << "eta = " << eta << endl
       << "theta = " << theta << endl
       << "c = " << c << endl
       << "precision = " << FT::get_prec() << endl
       << "row_expo = " << static_cast<int>(m.is_enable_row_expo()) << endl
       << "enable_bf = " << static_cast<int>(m.is_enable_bf()) << endl
       << "long_in_size_reduction = " << static_cast<int>(m.is_row_op_force_long()) << endl;

#ifdef HOUSEHOLDER_NAIVELY
  cerr << "householder_naively = 1" << endl;
#endif  // HOUSEHOLDER_NAIVELY
}

template <class ZT, class FT> inline void HLLLReduction<ZT, FT>::compute_dR(int k, FT delta_)
{
  long expo;
  m.get_R(dR[k], k, k, expo);
  dR[k].mul(dR[k], dR[k]);
  dR[k].mul(delta_, dR[k]);  // dR[k] = delta_ * R(k, k)^2
}

template <class ZT, class FT>
bool is_hlll_reduced(MatHouseholder<ZT, FT> &m, double delta, double eta, int d, bool compute);

FPLLL_END_NAMESPACE

#endif /* FPLLL_HLLL_H */
