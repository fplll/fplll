#ifndef BKZ_PARAMS_H
#define BKZ_PARAMS_H

/* (C) 2014-2016 Martin Albrecht.

   This file is part of fplll. fplll is free software: you can
   redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software
   Foundation, either version 2.1 of the License, or (at your option)
   any later version.

   fplll is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with fplll. If not, see
   <http://www.gnu.org/licenses/>.

*/

#include <string>
#include <vector>
#include "defs.h"

FPLLL_BEGIN_NAMESPACE

class Pruning {

public:

  double radiusFactor;              //< radius/Gaussian heuristic
  double probability;               //< success probability
  std::vector<double> coefficients;  //< pruning coefficients

  Pruning() : radiusFactor(1.), probability(1.) {};

  /** Sets all pruning coefficients to 1, except the last <level>
      coefficients, these will be linearly with slope -1 / blockSize.

      @param level number of levels in linear descent
  */

  static Pruning LinearPruning(int blockSize, int level) {

    Pruning pruning = Pruning();
    int startDescent = blockSize - level;

    if (startDescent > blockSize)
      startDescent = blockSize;

    if (startDescent < 1)
      startDescent = 1;

    pruning.coefficients.resize(blockSize);
    for (int k = 0; k < startDescent; k++) {
      pruning.coefficients[k] = 1.0;
    }
    for (int k = 0; k < blockSize - startDescent; k++) {
      pruning.coefficients[startDescent + k] = ((double)(blockSize - k - 1)) / blockSize;
    }
    // TODO: need to adapt probability
    pruning.radiusFactor = 1.0;
    pruning.probability = 1.0;

    return pruning;
  }
};


class Strategy {
public:
  vector<Pruning> pruning_parameters;
  vector<int> preprocessing_blocksizes;

  static Strategy EmptyStrategy() {
    Strategy strat;
    strat.pruning_parameters.emplace_back(Pruning());
    return strat;
  };

  const Pruning &getPruning(double radius, double gh) const;
};


class BKZParam {
public:
  BKZParam(int blockSize, vector<Strategy> &strategies,
           double delta = LLL_DEF_DELTA, int flags = BKZ_DEFAULT,
           int maxLoops = 0, double maxTime = 0, double autoAbort_scale = 1.0,
           int autoAbort_maxNoDec = 5, double ghFactor = 1.1,
           double minSuccessProbability = 0.5,
           int rerandomizationDensity = 3)
    : blockSize(blockSize), strategies(strategies), delta(delta), flags(flags), maxLoops(maxLoops), maxTime(maxTime),
      autoAbort_scale(autoAbort_scale), autoAbort_maxNoDec(autoAbort_maxNoDec),
      ghFactor(ghFactor), dumpGSOFilename("gso.log"),
      minSuccessProbability(minSuccessProbability),
      rerandomizationDensity(rerandomizationDensity) {
    if (strategies.empty()) {
      strategies = vector<Strategy>();
      for (long b = 0; b <= blockSize; ++b) {
        strategies.emplace_back(std::move(Strategy::EmptyStrategy()));
      }
    }
  };

  /** Block size used for enumeration **/
  int blockSize;

  /** Strategies (pruning coefficients, preprocessing)  */

  vector<Strategy> &strategies;

  /** LLL parameter delta **/
  double delta;

  /** See BKZFlags **/
  int flags;

  /** Maximum number of loops to execute **/
  int maxLoops;

  /** Maximum time to spend **/
  double maxTime;

  /** If BKZ_AUTOABORT is set, We abort if newSlope < autoAbort_scale * oldSlope
      is true for autoAbort_maxNoDec loops.
   */
  double autoAbort_scale;
  int    autoAbort_maxNoDec;

  /** If BKZ_GH_BND is set, the enumeration bound will be set to ghFactor times
      the Gaussian Heuristic
  */

  double ghFactor;

  /** If BKZ_DUMP_GSO is set, the norms of the GSO matrix are written to this
      file after each complete round.
  */

  string dumpGSOFilename;

  /** minimum success probability when using extreme pruning */

  double minSuccessProbability;

  /** density of rerandomization operation when using extreme pruning **/

  int rerandomizationDensity;

};

FPLLL_END_NAMESPACE
#endif /* BKZ_PARAMS_H */
