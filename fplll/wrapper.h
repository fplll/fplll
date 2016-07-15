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

#ifndef FPLLL_WRAPPER_H
#define FPLLL_WRAPPER_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

/* The matrix b must not be modified before calling lll().
   lll() must be called only once. */

class Wrapper {
public:
  /* u must be either empty or the identity matrix */
  Wrapper(IntMatrix& b, IntMatrix& u, IntMatrix& uInv,
          double delta, double eta, int flags);

  bool lll();

  int status;

private:
  IntMatrix& b;
  IntMatrix& u;
  IntMatrix& uInv;

#ifdef FPLLL_WITH_ZLONG
  ZZ_mat<long> bLong;
  ZZ_mat<long> uLong;    // Always empty
  ZZ_mat<long> uInvLong; // Always empty
#endif

  double delta;
  double eta;
  int goodPrec;
  bool useLong;
  int flags;

  bool little(int kappa,int precision);

  template<class Z, class F>
  int callLLL(ZZ_mat<Z>& bz, ZZ_mat<Z>& uz, ZZ_mat<Z>& uInvZ,
              LLLMethod method, int precision, double delta, double eta);

  template<class F>
  int fastLLL(double delta, double eta);

  template<class Z, class F>
  int heuristicLLL(ZZ_mat<Z>& bz, ZZ_mat<Z>& uz, ZZ_mat<Z>& uInvZ,
                   int precision, double delta, double eta);

  template<class Z, class F>
  int provedLLL(ZZ_mat<Z>& bz, ZZ_mat<Z>& uz, ZZ_mat<Z>& uInvZ,
                int precision, double delta, double eta);

  int heuristicLoop(int precision);
  int provedLoop(int precision);
  int lastLLL();

  void setUseLong(bool value);
  int increasePrec(int precision);

  int max_exponent;
  int n;
  int d;
  int lastEarlyRed;
};

FPLLL_END_NAMESPACE

#endif
