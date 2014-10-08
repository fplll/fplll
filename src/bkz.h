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
#include "evaluator.h"

FPLLL_BEGIN_NAMESPACE

class BKZParam {
public:

  BKZParam(int blockSize=0, double delta=LLL_DEF_DELTA, int flags=BKZ_DEFAULT,
           int maxLoops=0, double maxTime=0, int linearPruningLevel=0,
           double autoAbort_scale=1.0, int autoAbort_maxNoDec=5) :
  blockSize(blockSize), delta(delta), flags(flags),
  maxLoops(maxLoops), maxTime(maxTime),
  autoAbort_scale(autoAbort_scale), autoAbort_maxNoDec(autoAbort_maxNoDec),
  dumpGSOFilename("gso.log"), preprocessing(NULL) {
    if (linearPruningLevel > 0) {
      enableLinearPruning(linearPruningLevel);
    }
  }

  /** Block size used for enumeration **/
  int    blockSize;

  /** LLL parameter delta **/
  double delta;

  /** See BKZFlags **/
  int    flags;

  /** Maximum number of loops to execute **/
  int    maxLoops;

  /** Maximum time to spend **/
  double maxTime;

  /** If BKZ_AUTOABORT is set, We abort if newSlope < autoAbort_scale * oldSlope
      is true for autoAbort_maxNoDec loops.
   */
  double autoAbort_scale;
  int autoAbort_maxNoDec;

  /** If not empty these are the prunning coefficients used for prunned
      enumeration. If linearPruningLevel > 0, the first blockSize -
      linearPruningLevel coefficients will be one. Afterwards the coefficients
      drop linearly with slope -1/blockSize.
  */
  
  vector<double> pruning;

  /** If BKZ_DUMP_GSO us set, the norms of the GSO matrix are written to this
      file after each complete round.
  */
  
  string dumpGSOFilename;

  /** If not NULL, these parameters are used for BKZ preprocessing. It is
      allowed to nest these preprocessing parameters
  */
  
  BKZParam *preprocessing;

  inline void enableLinearPruning(int level) {
    int startDescent = blockSize - level;
    
    if (startDescent > blockSize)
      startDescent = blockSize;
    
    pruning.resize(blockSize);
    for(int k=0; k< startDescent; k++)
      pruning[k] = 1.0;
    for(int k=startDescent; k<blockSize; k++)
      pruning[k] = ((double)(blockSize-k-startDescent))/blockSize;
  }
                     
};

template<class FT>
static double getCurrentSlope(MatGSO<Integer, FT>& m, int startRow, int stopRow);


/* The matrix must be LLL-reduced */
template<class FT>
class BKZReduction {
public:
  BKZReduction(MatGSO<Integer, FT>& m, LLLReduction<Integer, FT>& lllObj,
               const BKZParam& param);
  ~BKZReduction();

  /**
     Run enumeration to find a new shortest vecto in the sublattice B[kappa,kappa+blockSize]

     @param kappa Start row
     @param blockSize Block size to use, this may be < param.blockSize
     @param param Parameters to use for this enumeration (blockSize is ignored)
     @param clean Did we change anything?
  */
  
  bool svpReduction(int kappa, int blockSize, const BKZParam &param, bool& clean);
  bool bkzLoop(const int loop, int& kappaMax, const BKZParam &param, int minRow, int maxRow, bool& clean);
  bool bkz();
  void dumpGSO(const std::string filename, const std::string prefix, bool append = true);
 
  int status;

private:
  void printParams(const BKZParam &param, ostream &out);
  bool setStatus(int newStatus);

  const BKZParam& param;
  int numRows;
  MatGSO<Integer, FT>& m;
  LLLReduction<Integer, FT>& lllObj;
  FastEvaluator<FT> evaluator;
  FT delta;

  // Temporary data
  const vector<FT> emptyTarget, emptySubTree;
  FT maxDist, deltaMaxDist;
  double cputimeStart;
};

template<class FT>
class BKZAutoAbort {
public:
  BKZAutoAbort(MatGSO<Integer, FT>& m, int numRows, int startRow = 0): m(m),
    oldSlope(numeric_limits<double>::max()), noDec(-1), numRows(numRows), startRow(startRow) {}
  bool testAbort(double scale=1.0, int maxNoDec=5);

private:
  MatGSO<Integer, FT>& m;
  double oldSlope;
  int noDec;
  int numRows;
  int startRow;
};

FPLLL_END_NAMESPACE

#endif
