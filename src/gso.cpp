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

#include "gso.h"

FPLLL_BEGIN_NAMESPACE

template<class ZT, class FT>
MatGSO<ZT, FT>::MatGSO(Matrix<ZT>& argB, Matrix<ZT>& argU, Matrix<ZT>& argUInvT,
                       int flags) :
  b(argB),
  enableIntGram(flags & GSO_INT_GRAM),
  enableRowExpo(flags & GSO_ROW_EXPO),
  enableTransform(argU.getRows() > 0),
  enableInvTransform(argUInvT.getRows() > 0),
  rowOpForceLong(flags & GSO_OP_FORCE_LONG),
  u(argU), uInvT(argUInvT),
  nKnownRows(0), nSourceRows(0), nKnownCols(0),
  colsLocked(false), allocDim(0)
{
  FPLLL_DEBUG_CHECK(!(enableIntGram && enableRowExpo));
  d = b.getRows();
  if (enableRowExpo) {
    tmpColExpo.resize(b.getCols());
  }
  sizeIncreased();
#ifdef DEBUG
  rowOpFirst = rowOpLast = -1;
#endif
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::invalidateGSORow(int i, int newValidCols) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows
                    && newValidCols >= 0 && newValidCols <= i + 1);
  gsoValidCols[i] = min(gsoValidCols[i], newValidCols);
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::updateBF(int i) {
  int n = max(nKnownCols, initRowSize[i]);
  if (enableRowExpo) {
    long maxExpo = LONG_MIN;
    for (int j = 0; j < n; j++) {
      b(i, j).get_f_exp(bf(i, j), tmpColExpo[j]);
      maxExpo = max(maxExpo, tmpColExpo[j]);
    }
    for (int j = 0; j < n; j++) {
      bf(i, j).mul_2si(bf(i, j), tmpColExpo[j] - maxExpo);
    }
    rowExpo[i] = maxExpo;
  }
  else {
    for (int j = 0; j < n; j++) {
      bf(i, j).set_z(b(i, j));
    }
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::invalidateGramRow(int i) {
  for (int j = 0; j <= i; j++)
    gf(i, j).set_nan();
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::rowOpEnd(int first, int last) {
#ifdef DEBUG
  FPLLL_DEBUG_CHECK(rowOpFirst == first && rowOpLast == last);
  rowOpFirst = rowOpLast = -1;
#endif
  for (int i = first; i < last; i++) {
    if (!enableIntGram) {
      updateBF(i);
      invalidateGramRow(i);
      for (int j = i + 1; j < nKnownRows; j++)
        gf(j, i).set_nan();
    }
    invalidateGSORow(i, 0);
  }
  for (int i = last; i < nKnownRows; i++) {
    invalidateGSORow(i, first);
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::discoverRow() {
  FPLLL_DEBUG_CHECK(nKnownRows < d);
  /* Early reduction (colsLocked=true) is not allowed when enableIntGram=true,
     since nKnownCols might be too small to compute all the g(i,j). */
  FPLLL_DEBUG_CHECK(!(colsLocked && enableIntGram));
  int i = nKnownRows;

  nKnownRows++;
  if (!colsLocked) {
    nSourceRows = nKnownRows;
    nKnownCols = max(nKnownCols, initRowSize[i]);
  }
  if (enableIntGram) {
    for (int j = 0; j <= i; j++)
      dotProduct(g(i, j), b[i], b[j], nKnownCols);
  }
  else {
    invalidateGramRow(i);
  }
  gsoValidCols[i] = 0;
}

template<class ZT, class FT>
long MatGSO<ZT, FT>::getMaxMuExp(int i, int nColumns) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && gsoValidCols[i] >= nColumns);
  long maxExpo = LONG_MIN, expo;
  for (int j = 0; j < nColumns; j++) {
    long expo2 = getMuExp(i, j, expo).exponent();
    maxExpo = max(maxExpo, expo + expo2);
  }
  return maxExpo;
}

template<class ZT, class FT>
bool MatGSO<ZT, FT>::updateGSORow(int i, int lastJ) {
  //FPLLL_TRACE_IN("Updating GSO up to (" << i << ", " << lastJ << ")");
  //FPLLL_TRACE("nKnownRows=" << nKnownRows << " nSourceRows=" << nSourceRows);
  if (i >= nKnownRows) {
    discoverRow();
  }
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows
                    && lastJ >= 0 && lastJ < nSourceRows);

  int j = max(0, gsoValidCols[i]);

  for (; j <= lastJ; j++) {
    getGram(ftmp1, i, j);
    FPLLL_DEBUG_CHECK(j == i || gsoValidCols[j] >= j);
    for (int k = 0; k < j; k++) {
      ftmp2.mul(mu(j, k), r(i, k));
      ftmp1.sub(ftmp1, ftmp2);
    }
    r(i, j) = ftmp1;
    if (i > j) {
      mu(i, j).div(ftmp1, r(j, j));
      if (!mu(i, j).is_finite()) return false;
    }
  }

  gsoValidCols[i] = j; // = max(0, gsoValidCols[i], lastJ + 1)
  //FPLLL_TRACE_OUT("End of GSO update");
  return true;
}


template<class ZT, class FT>
void MatGSO<ZT, FT>::row_add(int i, int j) {
  b[i].add(b[j], nKnownCols);
  if (enableTransform) {
    u[i].add(u[j]);
    if (enableInvTransform)
      uInvT[j].sub(uInvT[i]);
  }

  if (enableIntGram) {
    // g(i, i) += 2 * g(i, j) + g(j, j)
    ztmp1.mul_2si(g(i, j), 1);
    ztmp1.add(ztmp1, g(j, j));
    g(i, i).add(g(i, i), ztmp1);

    for (int k = 0; k < nKnownRows; k++)
      if (k != i)
        symG(i, k).add(symG(i, k), symG(j, k));
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::row_sub(int i, int j) {
  b[i].sub(b[j], nKnownCols);
  if (enableTransform) {
    u[i].sub(u[j]);
    if (enableInvTransform)
      uInvT[j].add(uInvT[i]);
  }

  if (enableIntGram) {
    // g(i, i) += g(j, j) - 2 * g(i, j)
    ztmp1.mul_2si(g(i, j), 1);
    ztmp1.sub(g(j, j), ztmp1);
    g(i, i).add(g(i, i), ztmp1);

    for (int k = 0; k < nKnownRows; k++)
      if (k != i)
        symG(i, k).sub(symG(i, k), symG(j, k));
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_si(int i, int j, long x) {
  b[i].addmul_si(b[j], x, nKnownCols);
  if (enableTransform) {
    u[i].addmul_si(u[j], x);
    if (enableInvTransform)
      uInvT[j].addmul_si(uInvT[i], -x);
  }

  if (enableIntGram) {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul_si(g(i, j), x);
    ztmp1.mul_2si(ztmp1, 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul_si(g(j, j), x);
    ztmp1.mul_si(ztmp1, x);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < nKnownRows; k++) {
      if (k == i) continue;
      ztmp1.mul_si(symG(j, k), x);
      symG(i, k).add(symG(i, k), ztmp1);
    }
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_si_2exp(int i, int j,
        long x, long expo) {
  b[i].addmul_si_2exp(b[j], x, expo, nKnownCols, ztmp1);
  if (enableTransform) {
    u[i].addmul_si_2exp(u[j], x, expo, ztmp1);
    if (enableInvTransform)
      uInvT[j].addmul_si_2exp(uInvT[i], -x, expo, ztmp1);
  }

  if (enableIntGram) {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul_si(g(i, j), x);
    ztmp1.mul_2si(ztmp1, expo + 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul_si(g(j, j), x);
    ztmp1.mul_si(ztmp1, x);
    ztmp1.mul_2si(ztmp1, 2 * expo);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < nKnownRows; k++) {
      if (k == i) continue;
      ztmp1.mul_si(symG(j, k), x);
      ztmp1.mul_2si(ztmp1, expo);
      symG(i, k).add(symG(i, k), ztmp1);
    }
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_2exp(int i, int j, const ZT& x, long expo) {
  b[i].addmul_2exp(b[j], x, expo, nKnownCols, ztmp1);
  if (enableTransform) {
    u[i].addmul_2exp(u[j], x, expo, ztmp1);
    if (enableInvTransform) {
      ZT minusX;
      minusX.neg(x);
      uInvT[j].addmul_2exp(uInvT[i], minusX, expo, ztmp1);
    }
  }

  if (enableIntGram) {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul(g(i, j), x);
    ztmp1.mul_2si(ztmp1, expo + 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul(g(j, j), x);
    ztmp1.mul(ztmp1, x);
    ztmp1.mul_2si(ztmp1, 2 * expo);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < nKnownRows; k++) {
      if (k == i) continue;
      ztmp1.mul(symG(j, k), x);
      ztmp1.mul_2si(ztmp1, expo);
      symG(i, k).add(symG(i, k), ztmp1);
    }
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_we(int i, int j, const FT& x, long expoAdd) {
  FPLLL_DEBUG_CHECK(j >= 0 && /*i > j &&*/ i < nKnownRows && j < nSourceRows);
  long expo;
  long lx = x.get_si_exp_we(expo, expoAdd);

  if (expo == 0) {
    if (lx == 1)
      row_add(i, j);
    else if (lx == -1)
      row_sub(i, j);
    else if (lx != 0)
      row_addmul_si(i, j, lx);
  }
  else if (rowOpForceLong) {
    row_addmul_si_2exp(i, j, lx, expo);
  }
  else {
    x.get_z_exp_we(ztmp2, expo, expoAdd);
    row_addmul_2exp(i, j, ztmp2, expo);
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::rowSwap(int i, int j) {
  FPLLL_DEBUG_CHECK(!enableInvTransform);
  b.swapRows(i, j);
  if (enableTransform) {
    u.swapRows(i, j);
  }

  if (enableIntGram) {
    for (int k = 0; k < i; k++)
      g(i, k).swap(g(j, k));
    for (int k = i + 1; k < j; k++)
      g(k, i).swap(g(j, k));
    for (int k = j + 1; k < nKnownRows; k++)
      g(k, i).swap(g(k, j));
    g(i, i).swap(g(j, j));
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::moveRow(int oldR, int newR) {
  FPLLL_DEBUG_CHECK(!colsLocked);
  if (newR < oldR) {
    FPLLL_DEBUG_CHECK(oldR < nKnownRows && !colsLocked);
    for (int i = newR; i < nKnownRows; i++) {
      invalidateGSORow(i, newR);
    }
    rotate(gsoValidCols.begin() + newR, gsoValidCols.begin() + oldR,
           gsoValidCols.begin() + oldR + 1);
    mu.rotateRight(newR, oldR);
    r.rotateRight(newR, oldR);
    b.rotateRight(newR, oldR);
    if (enableTransform) {
      u.rotateRight(newR, oldR);
      if (enableInvTransform)
        uInvT.rotateRight(newR, oldR);
    }
    if (enableIntGram)
      g.rotateGramRight(newR, oldR, nKnownRows);
    else {
      gf.rotateGramRight(newR, oldR, nKnownRows);
      bf.rotateRight(newR, oldR);
    }
    if (enableRowExpo)
      rotate(rowExpo.begin() + newR, rowExpo.begin() + oldR,
             rowExpo.begin() + oldR + 1);
  }
  else if (newR > oldR) {
    for (int i = oldR; i < nKnownRows; i++) {
      invalidateGSORow(i, oldR);
    }
    rotate(gsoValidCols.begin() + oldR, gsoValidCols.begin() + oldR + 1,
           gsoValidCols.begin() + newR + 1);
    mu.rotateLeft(oldR, newR);
    r.rotateLeft(oldR, newR);
    b.rotateLeft(oldR, newR);
    if (enableTransform) {
      u.rotateLeft(oldR, newR);
      if (enableInvTransform)
        uInvT.rotateLeft(oldR, newR);
    }
    if (enableIntGram) {
      if (oldR < nKnownRows - 1)
        g.rotateGramLeft(oldR, min(newR, nKnownRows - 1), nKnownRows);
    }
    else {
      if (oldR < nKnownRows - 1)
        gf.rotateGramLeft(oldR, min(newR, nKnownRows - 1), nKnownRows);
      bf.rotateLeft(oldR, newR);
    }
    if (enableRowExpo)
      rotate(rowExpo.begin() + oldR, rowExpo.begin() + oldR + 1,
             rowExpo.begin() + newR + 1);
    if (newR >= nKnownRows) {
      rotate(initRowSize.begin() + oldR, initRowSize.begin() + oldR + 1,
             initRowSize.begin() + newR + 1);
      if (oldR < nKnownRows) {
        nKnownRows--;
        nSourceRows = nKnownRows;
        initRowSize[newR] = max(b[newR].sizeNZ(), 1);
      }
    }
  }
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::lockCols() {
  colsLocked = true;
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::unlockCols() {
  nKnownRows = nSourceRows;
  colsLocked = false;
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::applyTransform(const Matrix<FT>& transform, int srcBase, int targetBase) {
  int targetSize = transform.getRows(), srcSize = transform.getCols();
  int oldD = d;
  createRows(targetSize);
  for (int i = 0; i < targetSize; i++) {
    for (int j = 0; j < srcSize; j++) {
      row_addmul(oldD + i, srcBase + j, transform(i, j));
    }
  }
  rowOpBegin(targetBase, targetBase + targetSize);
  for (int i = 0; i < targetSize; i++) {
    rowSwap(targetBase + i, oldD + i);
  }
  rowOpEnd(targetBase, targetBase + targetSize);
  removeLastRows(targetSize);
}

template<class ZT, class FT>
void MatGSO<ZT, FT>::sizeIncreased() {
  int oldD = mu.getRows();

  if (d > allocDim) {
    if (enableIntGram)
      g.resize(d, d);
    else {
      bf.resize(d, b.getCols());
      gf.resize(d, d);
    }
    mu.resize(d, d);
    r.resize(d, d);
    gsoValidCols.resize(d);
    initRowSize.resize(d);
    if (enableRowExpo) {
      rowExpo.resize(d);
    }
    allocDim = d;
  }

  for (int i = oldD; i < d; i++) {
    initRowSize[i] = max(b[i].sizeNZ(), 1);
    if (!enableIntGram) {
      bf[i].fill(0); // updateBF might not copy all the zeros of b[i]
      updateBF(i);
    }
  }
}

FPLLL_END_NAMESPACE
