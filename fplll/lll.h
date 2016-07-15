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

#ifndef FPLLL_LLL_H
#define FPLLL_LLL_H

#include "gso.h"

FPLLL_BEGIN_NAMESPACE

/* The precision of FT must be defined before creating an instance of
   LLLReduction. The matrix b can be modified between each call to lll. */

template<class ZT, class FT>
class LLLReduction {
public:
  /**
   * Constructor.
   * The precision of FT must be defined before creating an instance of the
   * class and must remain the same until the object is destroyed (or no longer
   * needed).
   */
  LLLReduction(MatGSO<ZT, FT>& m, double delta, double eta,
               int flags);

  ~LLLReduction() {
    LDConvHelper::free();
  }
  
  bool lll(int kappaMin = 0, int kappaStart = 0, int kappaEnd = -1);
  inline bool sizeReduction(int kappaMin = 0, int kappaEnd = -1);

  int status;
  int finalKappa;
  int lastEarlyRed;
  int zeros;
  int nSwaps;

private:
  bool babai(int kappa, int nCols);
  inline bool earlyReduction(int start);
  inline void printParams();
  inline bool setStatus(int newStatus);

  MatGSO<ZT, FT>& m;
  FT delta, eta, swapThreshold;

  bool enableEarlyRed;
  bool siegel;
  bool verbose;

  vector<FT> lovaszTests;
  vector<FT> babaiMu;
  vector<long> babaiExpo;
  ZT ztmp1;
  FT muMant, ftmp1;
};

template<class ZT, class FT>
bool isLLLReduced(MatGSO<ZT, FT>& m, double delta, double eta);

template<class ZT, class FT>
inline bool LLLReduction<ZT, FT>::sizeReduction(int kappaMin, int kappaEnd) {
  if (kappaEnd == -1) kappaEnd = m.d;
  for (int k = kappaMin; k < kappaEnd; k++) {
    if ((k > 0 && !babai(k, k)) || !m.update_gso_row(k))
      return false;
  }
  return setStatus(RED_SUCCESS);
}

template<class ZT, class FT>
inline bool LLLReduction<ZT, FT>::earlyReduction(int start) {
  m.lock_cols();
  if (verbose) {
    cerr << "Early reduction start=" << start + 1 << endl;
  }
  for (int i = start; i < m.d; i++) {
    if (!babai(i, start)) return false;
  }
  m.unlock_cols();
  lastEarlyRed = start;
  return true;
}

template<class ZT, class FT>
inline void LLLReduction<ZT, FT>::printParams() {
  cerr << "Entering LLL"
       << "\ndelta = " << delta
       << "\neta = " << eta
       << "\nprecision = " << FT::getprec()
       << "\nexact_dot_product = " << static_cast<int>(m.enable_int_gram)
       << "\nrow_expo = " << static_cast<int>(m.enable_row_expo)
       << "\nearly_red = " << static_cast<int>(enableEarlyRed)
       << "\nsiegel_cond = " << static_cast<int>(siegel)
       << "\nlong_in_babai = " << static_cast<int>(m.row_op_force_long) << endl;
}

template<class ZT, class FT>
inline bool LLLReduction<ZT, FT>::setStatus(int newStatus) {
  status = newStatus;
  if (verbose) {
    if (status == RED_SUCCESS)
      cerr << "End of LLL: success" << endl;
    else
      cerr << "End of LLL: failure: " << RED_STATUS_STR[status] << endl;
  }
  return status == RED_SUCCESS;
}


FPLLL_END_NAMESPACE

#endif
