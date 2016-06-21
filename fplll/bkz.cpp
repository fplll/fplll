/* Copyright (C) 2011 Xavier Pujol
   (C) 2014 Martin Albrecht.

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

#include <iomanip>

/* Template source file */
#include "bkz.h"
#include "bkz_params.h"
#include "enum/enumerate.h"
#include <iomanip>

FPLLL_BEGIN_NAMESPACE

template <class FT>
BKZReduction<FT>::BKZReduction(MatGSO<Integer, FT> &m, LLLReduction<Integer, FT> &lllObj,
                               const BKZParam &param)
    : status(RED_SUCCESS), param(param), m(m), lllObj(lllObj) {
  for (numRows = m.d; numRows > 0 && m.b[numRows - 1].is_zero(); numRows--) {
  }
  this->delta = param.delta;
}

template <class FT> BKZReduction<FT>::~BKZReduction() {}

template <class FT> double getCurrentSlope(MatGSO<Integer, FT> &m, int startRow, int stopRow) {
  FT f, logF;
  long expo;
  vector<double> x;
  x.resize(stopRow);
  for (int i = startRow; i < stopRow; i++) {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    x[i] = logF.get_d() + expo * std::log(2.0);
  }
  int n        = stopRow - startRow;
  double iMean = (n - 1) * 0.5 + startRow, xMean = 0, v1 = 0, v2 = 0;
  for (int i = startRow; i < stopRow; i++) {
    xMean += x[i];
  }
  xMean /= n;
  for (int i = startRow; i < stopRow; i++) {
    v1 += (i - iMean) * (x[i] - xMean);
    v2 += (i - iMean) * (i - iMean);
  }
  return v1 / v2;
}

template <class FT>
void computeGaussHeurDist(MatGSO<Integer, FT> &m, FT &maxDist, long maxDistExpo, int kappa,
                          int blockSize, double ghFactor) {
  double t = (double)blockSize / 2.0 + 1;
  t        = tgamma(t);
  t        = pow(t, 2.0 / (double)blockSize);
  t        = t / M_PI;
  FT f, g, h;
  f = t;
  m.getR(h, kappa, kappa);
  g.log(h);
  for (int i = kappa + 1; i < kappa + blockSize; i++) {
    m.getR(h, i, i);
    h.log(h);
    g.add(g, h);
  }
  h = (double)blockSize;
  g.div(g, h);
  g.exponential(g);
  f.mul(f, g);
  f.mul_2si(f, -maxDistExpo);
  h = ghFactor;
  f.mul(f, h);
  if (f < maxDist) {
    maxDist = f;
  }
}

/** Randomize basis between from ``min_row`` and ``max_row`` (exclusive)

    1. permute rows
    2. apply lower triangular matrix with coefficients in -1,0,1
    3. LLL reduce result

    :param min_row: start in this row
    :param max_row: stop at this row (exclusive)
    :param stats: object for maintaining statistics
    :param triangular_density: number of non-zero coefficients in lower triangular
    transformation matrix
**/

template <class FT>
void BKZReduction<FT>::rerandomizeBlock(int minRow, int maxRow, int density) {
  if (maxRow - minRow < 2)
    return;

  // 1. permute rows
  size_t niter = 4 * (maxRow-minRow); // some guestimate

  for(size_t i=0; i<niter; ++i) {
    size_t a = gmp_urandomm_ui(RandGen::getGMPState(), maxRow-minRow-1) + minRow;
    size_t b = a;
    while (b == a) {
      b = gmp_urandomm_ui(RandGen::getGMPState(), maxRow-minRow-1) + minRow;
    }
    m.moveRow(b, a);
  }

  // 2. triangular transformation matrix with coefficients in -1,0,1
  m.rowOpBegin(minRow, maxRow);
  FT x;
  for(long a=minRow; a<maxRow-2; ++a) {
    for(long i=0; i<density; i++) {
      size_t b = gmp_urandomm_ui(RandGen::getGMPState(), maxRow-(a+1)-1) + a+1;
      if (gmp_urandomm_ui(RandGen::getGMPState(), 2))
        m.row_add(a, b);
      else
        m.row_sub(a, b);
    }
  }
  m.rowOpEnd(minRow, maxRow);

  // 3. LLL reduce
  if (!lllObj.lll(0, minRow, maxRow))
    throw lllObj.status;
  return;
}


template <class FT>
bool BKZReduction<FT>::svpPreprocessing(int kappa, int blockSize, const BKZParam &param) {
  bool clean = true;

  FPLLL_DEBUG_CHECK(param.strategies.size() > blockSize);

  auto &preproc = param.strategies[blockSize].preprocessing_blocksizes;
  for(auto it = preproc.begin(); it != preproc.end(); ++it) {
    int dummyKappaMax = numRows;
    BKZParam prepar = BKZParam(*it, param.strategies);
    clean &= tour(0, dummyKappaMax, prepar, kappa, kappa + blockSize);
  }

  return clean;
}

template <class FT>
const Pruning &BKZReduction<FT>::getPruning(int kappa, int blockSize, const BKZParam &par) const {

  FPLLL_DEBUG_CHECK(param.strategies.size() > blockSize);

  Strategy &strat = par.strategies[blockSize];

  long maxDistExpo;
  FT maxDist = m.getRExp(kappa, kappa, maxDistExpo);

  FT ghMaxDist;
  computeGaussHeurDist(m, ghMaxDist, maxDistExpo, kappa, blockSize, 1.0);
  return strat.getPruning(maxDist.get_d() * pow(2, maxDistExpo),
                          ghMaxDist.get_d() * pow(2, maxDistExpo));
}

template <class FT>
bool BKZReduction<FT>::svpPostprocessing(int kappa, int blockSize, const vector<FT> &solution) {
  // Is it already in the basis ?
  int nzVectors = 0, iVector = -1;
  for (int i = 0; i < blockSize; i++) {
    if (!solution[i].is_zero()) {
      nzVectors++;
      if (iVector == -1 && fabs(solution[i].get_d()) == 1)
        iVector = i;
    }
  }

  FPLLL_DEBUG_CHECK(nzVectors > 0);

  if (nzVectors == 1) {
    // Yes, it is another vector
    FPLLL_DEBUG_CHECK(iVector != -1 && iVector != 0);
    m.moveRow(kappa + iVector, kappa);
    if (!lllObj.sizeReduction(kappa, kappa + 1))
      throw lllObj.status;
  } else {
    // No, general case
    int d = m.d;
    m.createRow();
    m.rowOpBegin(d, d + 1);
    for (int i = 0; i < blockSize; i++) {
      m.row_addmul(d, kappa + i, solution[i]);
    }
    m.rowOpEnd(d, d + 1);
    m.moveRow(d, kappa);
    if (!lllObj.lll(kappa, kappa, kappa + blockSize + 1))
      throw lllObj.status;
    FPLLL_DEBUG_CHECK(m.b[kappa + blockSize].is_zero());
    m.moveRow(kappa + blockSize, d);
    m.removeLastRow();
  }
  return false;
}

template <class FT>
bool BKZReduction<FT>::svpReduction(int kappa, int blockSize, const BKZParam &par) {
  bool clean = true;

  int lllStart = (par.flags & BKZ_BOUNDED_LLL) ? kappa : 0;

  if (!lllObj.lll(lllStart, kappa, kappa + blockSize)) {
    throw lllObj.status;
  }

  clean &= (lllObj.nSwaps == 0);

  size_t trial = 0;
  double remaining_probability = 1.0;

  while (remaining_probability > 1. - par.minSuccessProbability) {
    if (trial > 0) {
      rerandomizeBlock(kappa+1, kappa+blockSize, par.rerandomizationDensity);
    }

    clean &= svpPreprocessing(kappa, blockSize, par);

    long maxDistExpo;
    FT maxDist = m.getRExp(kappa, kappa, maxDistExpo);
    FT deltaMaxDist;
    deltaMaxDist.mul(delta, maxDist);

    if ((par.flags & BKZ_GH_BND) && blockSize > 30) {
      computeGaussHeurDist(m, maxDist, maxDistExpo, kappa, blockSize, par.ghFactor);
    }

    const Pruning &pruning = getPruning(kappa, blockSize, par);

    vector<FT> &solCoord = evaluator.solCoord;
    solCoord.clear();
    Enumeration::enumerate(m, maxDist, maxDistExpo, evaluator, emptySubTree, emptySubTree, kappa,
                           kappa + blockSize, pruning.coefficients);

    nodes += Enumeration::getNodes();

    if (solCoord.empty()) {
      if (par.flags & BKZ_GH_BND) {
        return clean;
      } else {
        throw RED_ENUM_FAILURE;
      }
    }

    if (maxDist < deltaMaxDist) {
      clean &= svpPostprocessing(kappa, blockSize, solCoord);
    }
    remaining_probability *= (1- pruning.probability);
    trial += 1;
  }
  return clean;
}

template <class FT>
bool BKZReduction<FT>::tour(const int loop, int &kappaMax, const BKZParam &par, int minRow,
                            int maxRow) {
  bool clean = true;

  for (int kappa = minRow; kappa < maxRow - 1; kappa++) {
    int blockSize = min(par.blockSize, maxRow - kappa);
    clean &= svpReduction(kappa, blockSize, par);
    if ((par.flags & BKZ_VERBOSE) && kappaMax < kappa && clean) {
      cerr << "Block [1-" << setw(4) << kappa + 1 << "] BKZ-" << setw(0) << par.blockSize
           << " reduced for the first time" << endl;
      kappaMax = kappa;
    }
  }

  if (par.flags & BKZ_VERBOSE) {
    printTour(loop, minRow, maxRow);
  }

  if (par.flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "End of BKZ loop " << std::setw(4) << loop;
    prefix << " (" << std::fixed << std::setw(9) << std::setprecision(3)
           << (cputime() - cputimeStart) * 0.001 << "s)";
    dumpGSO(par.dumpGSOFilename, prefix.str());
  }

  return clean;
}

template <class FT> bool BKZReduction<FT>::bkz() {
  int flags       = param.flags;
  int finalStatus = RED_SUCCESS;
  nodes           = 0;

  if (flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "Input";
    dumpGSO(param.dumpGSOFilename, prefix.str(), false);
  }

  if (param.blockSize < 2)
    return setStatus(RED_SUCCESS);

  int i = 0;
  BKZAutoAbort<FT> autoAbort(m, numRows);

  if (flags & BKZ_VERBOSE) {
    cerr << "Entering BKZ:" << endl;
    printParams(param, cerr);
    cerr << endl;
  }
  cputimeStart = cputime();

  m.discoverAllRows();

  int kappaMax;
  bool clean = true;
  for (i = 0;; i++) {
    if ((flags & BKZ_MAX_LOOPS) && i >= param.maxLoops) {
      finalStatus = RED_BKZ_LOOPS_LIMIT;
      break;
    }
    if ((flags & BKZ_MAX_TIME) && (cputime() - cputimeStart) * 0.001 >= param.maxTime) {
      finalStatus = RED_BKZ_TIME_LIMIT;
      break;
    }
    if ((flags & BKZ_AUTO_ABORT) &&
        autoAbort.testAbort(param.autoAbort_scale, param.autoAbort_maxNoDec)) {
      break;
    }

    try {
      clean = tour(i, kappaMax, param, 0, numRows);
    } catch (RedStatus &e) {
      return setStatus(e);
    }

    if (clean || param.blockSize >= numRows)
      break;
  }

  if (flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "Output ";
    prefix << " (" << std::fixed << std::setw(9) << std::setprecision(3)
           << (cputime() - cputimeStart) * 0.001 << "s)";
    dumpGSO(param.dumpGSOFilename, prefix.str());
  }
  return setStatus(finalStatus);
}

template <class FT> void BKZReduction<FT>::printTour(const int loop, int minRow, int maxRow) {
  FT r0;
  Float fr0;
  long expo;
  r0  = m.getRExp(minRow, minRow, expo);
  fr0 = r0.get_d();
  fr0.mul_2si(fr0, expo);
  cerr << "End of BKZ loop " << std::setw(4) << loop << ", time = " << std::fixed << std::setw(9)
       << std::setprecision(3) << (cputime() - cputimeStart) * 0.001 << "s";
  cerr << ", r_" << minRow << " = " << fr0;
  cerr << ", slope = " << std::setw(9) << std::setprecision(6)
       << getCurrentSlope(m, minRow, maxRow);
  cerr << ", log2(nodes) = " << std::setw(9) << std::setprecision(6) << log2(nodes) << endl;
}

template <class FT> void BKZReduction<FT>::printParams(const BKZParam &param, ostream &out) {
  out << "blocksize: " << std::setw(3) << param.blockSize << ", ";
  out << "flags: 0x" << std::setw(4) << setfill('0') << std::hex << param.flags << ", " << std::dec
      << std::setfill(' ');
  out << "maxLoops: " << std::setw(3) << param.maxLoops << ", ";
  out << "maxTime: " << std::setw(0) << std::fixed << std::setprecision(1) << param.maxTime << ", ";
  if (param.flags & BKZ_AUTO_ABORT) {
    out << "autoAbort: (" << std::setw(0) << std::fixed << std::setprecision(4)
        << param.autoAbort_scale;
    out << ", " << std::setw(2) << param.autoAbort_maxNoDec << "), ";
  } else {
    out << "autoAbort: (     -,  -), ";
  }
  out << endl;
}

template <class FT> bool BKZReduction<FT>::setStatus(int newStatus) {
  status = newStatus;
  if (param.flags & BKZ_VERBOSE) {
    if (status == RED_SUCCESS)
      cerr << "End of BKZ: success" << endl;
    else
      cerr << "End of BKZ: failure: " << RED_STATUS_STR[status] << endl;
  }
  return status == RED_SUCCESS;
}

template <class FT>
void BKZReduction<FT>::dumpGSO(const std::string filename, const std::string prefix, bool append) {
  ofstream dump;
  if (append)
    dump.open(filename.c_str(), std::ios_base::app);
  else
    dump.open(filename.c_str());
  dump << std::setw(4) << prefix << ": ";
  FT f, logF;
  long expo;
  for (int i = 0; i < numRows; i++) {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    dump << std::setprecision(8) << logF.get_d() + expo * std::log(2.0) << " ";
  }
  dump << std::endl;
  dump.close();
}

/**
   Force instantiation of templates
*/

template <class FT> bool BKZAutoAbort<FT>::testAbort(double scale, int maxNoDec) {
  double newSlope = -getCurrentSlope(m, startRow, numRows);
  if (noDec == -1 || newSlope < scale * oldSlope)
    noDec = 0;
  else
    noDec++;
  oldSlope = min(oldSlope, newSlope);
  return noDec >= maxNoDec;
}

template class BKZReduction<FP_NR<double> >;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class BKZReduction<FP_NR<long double> >;
#endif

#ifdef FPLLL_WITH_QD
template class BKZReduction<FP_NR<dd_real> >;

template class BKZReduction<FP_NR<qd_real> >;
#endif

#ifdef FPLLL_WITH_DPE
template class BKZReduction<FP_NR<dpe_t> >;
#endif

template class BKZReduction<FP_NR<mpfr_t> >;

FPLLL_END_NAMESPACE
