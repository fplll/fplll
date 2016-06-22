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

/* Template source file */

#include "lll.h"
#include "util.h"

FPLLL_BEGIN_NAMESPACE

static inline bool isPowerOf2(int i) {
  return (i & (i - 1)) == 0;
}

template<class ZT, class FT>
LLLReduction<ZT, FT>::LLLReduction(MatGSO<ZT, FT>& m,
        double delta, double eta, int flags) :
  status(RED_SUCCESS), lastEarlyRed(0), m(m)
{
  /* No early reduction in proved mode (i.e. enableIntGram=true).
     NOTE: To make this possible, the hypothesis "g(i, j) is valid if
     0 <= i < nKnownRows and j <= i" in gso.h should be changed and
     MatGSO<ZT, FT>::discoverRow() should be rewritten. */
  enableEarlyRed = (flags & LLL_EARLY_RED) && !m.enableIntGram;
  siegel = flags & LLL_SIEGEL;
  verbose = flags & LLL_VERBOSE;
  this->delta = delta;
  this->eta = eta;
  swapThreshold = siegel ? delta - eta * eta : delta;
  zeros = 0;
}

template<class ZT, class FT>
bool LLLReduction<ZT, FT>::lll(int kappaMin, int kappaStart, int kappaEnd) {
  if (kappaEnd == -1) kappaEnd = m.d;


  FPLLL_DEBUG_CHECK(kappaMin <= kappaStart && kappaStart < kappaEnd && kappaEnd <= m.d);
  int startTime = cputime();
  int kappa = kappaStart + 1;
  int kappaMax = 0;
  int d = kappaEnd - kappaMin;

  zeros = 0;
  nSwaps = 0;
  finalKappa = 0;
  if (verbose) printParams();
  extendVect(lovaszTests, kappaEnd);
  extendVect(babaiMu, kappaEnd);
  extendVect(babaiExpo, kappaEnd);

  for (; zeros < d && m.b[0].is_zero(); zeros++) {
    m.moveRow(kappaMin, kappaEnd - 1 - zeros);
  }

  if (zeros < d && ((kappaStart > 0 && !babai(kappaStart, kappaStart))
                    || !m.updateGSORow(kappaStart))) {
    finalKappa = kappaStart;
    return false;
  }

  long long iter, maxIter;
  maxIter = static_cast<long long>(d - 2 * d * (d + 1)
                                   * ((m.b.getMaxExp() + 3) / std::log(delta.get_d())));

  for (iter = 0; iter < maxIter && kappa < kappaEnd - zeros; iter++) {
    if (kappa > kappaMax) {
      if (verbose) {
        cerr << "Discovering vector " << kappa - kappaMin + 1 + zeros << "/"
             << d << " cputime=" << cputime() - startTime << endl;
      }
      kappaMax = kappa;
      if (enableEarlyRed && isPowerOf2(kappa) && kappa > lastEarlyRed) {
        if (!earlyReduction(kappa)) {
          finalKappa = kappa;
          return false;
        }
      }
    }

    // Lazy size reduction
    if (!babai(kappa, kappa)) {
      finalKappa = kappa;
      return false;
    }


    // Tests Lovasz's condition
    m.getGram(lovaszTests[0], kappa, kappa);
    for (int i = 1; i <= kappa; i++) {
      ftmp1.mul(m.getMuExp(kappa, i - 1), m.getRExp(kappa, i - 1));
      lovaszTests[i].sub(lovaszTests[i - 1], ftmp1);
    }
    ftmp1.mul(m.getRExp(kappa - 1, kappa - 1), swapThreshold);
    if (m.enableRowExpo) {
      ftmp1.mul_2si(ftmp1, 2 * (m.rowExpo[kappa - 1] - m.rowExpo[kappa]));
    }

    if (ftmp1 > lovaszTests[siegel ? kappa : kappa - 1]) {
      nSwaps++;
      // Failure, computes the insertion index
      int oldK = kappa;
      for (kappa--; kappa > kappaStart; kappa--) {
        ftmp1.mul(m.getRExp(kappa - 1, kappa - 1), swapThreshold);
        if (m.enableRowExpo) {
          ftmp1.mul_2si(ftmp1, 2 * (m.rowExpo[kappa - 1] - m.rowExpo[oldK]));
        }
        if (ftmp1 < lovaszTests[siegel ? kappa : kappa - 1]) break;
      }
      //FPLLL_TRACE("Lovasz's condition is not satisfied, kappa=" << kappa << " old_kappa=" << oldK);
      // Moves the vector
      if (lovaszTests[kappa] > 0) {
        m.moveRow(oldK, kappa);
      }
      else {
        zeros++;
        m.moveRow(oldK, kappaEnd - zeros);
        kappa = oldK;
        continue;
      }
    }

    m.setR(kappa, kappa, lovaszTests[kappa]);
    kappa++;
  }

  if (kappa < kappaEnd - zeros)
    return setStatus(RED_LLL_FAILURE);
  else
    return setStatus(RED_SUCCESS);
}

template<class ZT, class FT>
bool LLLReduction<ZT, FT>::babai(int kappa, int nCols) {
  //FPLLL_TRACE_IN("kappa=" << kappa);
  long maxExpo = LONG_MAX;

  for (int iter = 0;; iter++) {
    if (!m.updateGSORow(kappa, nCols - 1))
      return setStatus(RED_GSO_FAILURE);

    bool loopNeeded = false;
    for (int j = nCols - 1; j >= 0 && !loopNeeded; j--) {
      m.getMu(ftmp1, kappa, j);
      ftmp1.abs(ftmp1);
      if (ftmp1 > eta) loopNeeded = true;
    }
    if (!loopNeeded) break;

    if (iter >= 2) {
      long newMaxExpo = m.getMaxMuExp(kappa, nCols);
      if (newMaxExpo > maxExpo - SIZE_RED_FAILURE_THRESH) {
        return setStatus(RED_BABAI_FAILURE);
      }
      maxExpo = newMaxExpo;
    }

    for (int j = 0; j < nCols; j++) {
      babaiMu[j] = m.getMuExp(kappa, j, babaiExpo[j]);
    }
    m.rowOpBegin(kappa, kappa + 1);
    for (int j = nCols - 1; j >= 0; j--) {
      muMant.rnd_we(babaiMu[j], babaiExpo[j]);
      if (muMant.zero_p()) continue;
      // Approximative update of the mu_(kappa,k)'s
      for (int k = 0; k < j; k++) {
        ftmp1.mul(muMant, m.getMuExp(j, k));
        /* When enableRowExpo=true, the following line relies on the fact that
           getMuExp(a, b, expo) returns expo = rowExpo[a] - rowExpo[b]. */
        babaiMu[k].sub(babaiMu[k], ftmp1);
      }
      // Operation on the basis
      //FPLLL_TRACE("Babai : row[" << kappa << "] += " << muMant << " * 2^" << babaiExpo[j] << " * row[" << j << "]");
      muMant.neg(muMant);
      m.row_addmul_we(kappa, j, muMant, babaiExpo[j]);
    }
    m.rowOpEnd(kappa, kappa + 1);
  }
  return true;
}

template<class ZT, class FT>
bool isLLLReduced(MatGSO<ZT, FT>& m, double delta, double eta) {
  FT ftmp1;
  FT ftmp2;
  FT delta_;
  delta_.set(delta);
  m.updateGSO();
  for(int i=0; i<m.d; i++) {
    for(int j=0; j<i; j++) {
      m.getMu(ftmp1, i, j);
      ftmp1.abs(ftmp1);
      if (ftmp1 > eta)
        return false;
    }
  }
  for(int i=1; i<m.d; i++) {
    m.getMu(ftmp2, i, i-1);
    ftmp2.mul(ftmp2, ftmp2); // μ^2

    ftmp2.sub(delta_, ftmp2); // δ - μ^2
    m.getR(ftmp1, i-1, i-1);
    ftmp2.mul(ftmp1, ftmp2); // (δ - μ^2) ⋅ r_{i-1,i-1}

    m.getR(ftmp1, i, i);  // r_{i,i}

    if (ftmp1 < ftmp2)
      return false;
  }
  return true;
}

FPLLL_END_NAMESPACE
