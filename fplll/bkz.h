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

#ifndef FPLLL_BKZ_H
#define FPLLL_BKZ_H

#include "lll.h"
#include "enum/evaluator.h"
#include "bkz_params.h"

FPLLL_BEGIN_NAMESPACE

template <class FT> class BKZAutoAbort {
public:
  BKZAutoAbort(MatGSO<Integer, FT> &m, int numRows, int startRow = 0)
      : m(m), oldSlope(numeric_limits<double>::max()), noDec(-1), numRows(numRows),
        startRow(startRow) {}
  bool testAbort(double scale = 1.0, int maxNoDec = 5);

private:
  MatGSO<Integer, FT> &m;
  double oldSlope;
  int noDec;
  int numRows;
  int startRow;
};

/** Finds the slope of the curve fitted to the lengths of the vectors from
    startRow to stopRow. The slope gives an indication of the quality of the
    LLL-reduced basis.
*/
template <class FT> double getCurrentSlope(MatGSO<Integer, FT> &m, int startRow, int stopRow);

/** Uses the Gaussian Heuristic Distance to compute a bound on the length of the
    shortest vector.
*/
template <class FT>
void computeGaussHeurDist(MatGSO<Integer, FT> &m, FT &maxDist, long maxDistExpo, int kappa,
                          int blockSize, double ghFactor);

/* The matrix must be LLL-reduced */
template <class FT> class BKZReduction {
public:
  BKZReduction(MatGSO<Integer, FT> &m, LLLReduction<Integer, FT> &lllObj, const BKZParam &param);
  ~BKZReduction();

  bool  svpPreprocessing(int kappa, int blockSize, const BKZParam &param);

  bool svpPostprocessing(int kappa, int blockSize, const vector<FT> &solution);

  /**
     Run enumeration to find a new shortest vector in the sublattice B[kappa,kappa+blockSize]

     @param kappa Start row
     @param blockSize Block size to use, this may be < param.blockSize
     @param param Parameters to use for this enumeration (blockSize is ignored)

     ``_ex`` variant is exception handling.
  */

  bool svpReduction(int kappa, int blockSize, const BKZParam &param);

  bool svpReduction_ex(int kappa, int blockSize, const BKZParam &param, bool &clean) {
    try {
      clean = svpReduction(kappa, blockSize, param);
      return true;
    } catch (RedStatus & e) {
      return setStatus(e);
    }
  }

  bool tour(const int loop, int &kappaMax, const BKZParam &param, int minRow, int maxRow);

  bool tour_ex(const int loop, int &kappaMax, const BKZParam &param, int minRow, int maxRow, bool &clean) {
    try {
      clean = tour(loop, kappaMax, param, minRow, maxRow);
      return true;
    } catch (RedStatus & e) {
      return setStatus(e);
    }
  }

  bool bkz();

  void rerandomizeBlock(int minRow, int maxRow, int density);

  /** I/O **/

  void dumpGSO(const std::string filename, const std::string prefix, bool append = true);

  int status;

  /**
      Number of nodes visited during enumeration.
  */

  long nodes;

private:
  void printTour(const int loop, int minRow, int maxRow);
  void printParams(const BKZParam &param, ostream &out);

  bool setStatus(int newStatus);

  const Pruning &getPruning(int kappa, int blockSize, const BKZParam &par) const;

  const BKZParam &param;
  int numRows;
  MatGSO<Integer, FT> &m;
  LLLReduction<Integer, FT> &lllObj;
  FastEvaluator<FT> evaluator;
  FT delta;

  // Temporary data
  const vector<FT> emptyTarget, emptySubTree;
  FT maxDist, deltaMaxDist;
  double cputimeStart;



};


FPLLL_END_NAMESPACE

#endif
