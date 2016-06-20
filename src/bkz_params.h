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
  double radiusFactor;               //< radius/Gaussian heuristic
  double probability;                //< success probability
  std::vector<double> coefficients;  //< pruning coefficients

  /** Sets all pruning coefficients to 1, except the last <level> coefficients,
      these will be linearly with slope -1 / blockSize.

      @param level number of levels in linear descent
  */

  inline Pruning() : radiusFactor(1.),probability(1.) {};
  Pruning(int blockSize, int level);
};

class BKZParam {
public:

  BKZParam(int blockSize = 0, double delta = LLL_DEF_DELTA, int flags = BKZ_DEFAULT,
           int maxLoops = 0, double maxTime = 0, int linearPruningLevel = 0,
           double autoAbort_scale = 1.0, int autoAbort_maxNoDec = 5, double ghFactor = 1.1)
      : blockSize(blockSize), delta(delta), flags(flags), maxLoops(maxLoops), maxTime(maxTime),
        autoAbort_scale(autoAbort_scale), autoAbort_maxNoDec(autoAbort_maxNoDec),
        ghFactor(ghFactor), dumpGSOFilename("gso.log"), preprocessing(NULL) {

    if (linearPruningLevel > 0) {
      pruning.emplace_back(blockSize, linearPruningLevel);
    }
  }

  /** Block size used for enumeration **/
  int blockSize;

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
  int autoAbort_maxNoDec;

  /** If not empty these are the prunning coefficients used for prunned
      enumeration.
  */

  std::vector<Pruning> pruning;

  /** If BKZ_GH_BND is set, the enumeration bound will be set to ghFactor times
      the Gaussian Heuristic
  */

  double ghFactor;

  /** If BKZ_DUMP_GSO is set, the norms of the GSO matrix are written to this
      file after each complete round.
  */

  string dumpGSOFilename;

  /** If not NULL, these parameters are used for BKZ preprocessing. It is
      allowed to nest these preprocessing parameters
  */

  BKZParam *preprocessing;

  const Pruning &getPruning(double radius, double gh) const;
};

FPLLL_END_NAMESPACE
#endif /* BKZ_PARAMS_H */
