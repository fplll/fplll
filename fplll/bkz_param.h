#ifndef BKZ_PARAM_H
#define BKZ_PARAM_H

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

#include "defs.h"
#include <string>
#include <vector>

FPLLL_BEGIN_NAMESPACE

/**
   Pruning parameters for one radius (expressed as a ratio to the Gaussian heuristic)
 */

class Pruning
{

public:
  double radius_factor;              //< radius/Gaussian heuristic
  std::vector<double> coefficients;  //< pruning coefficients
  double expectation;                //< either success probability or number of solutions
  PrunerMetric
      metric;  //< metric used for optimisation (success probability or number of solutions)
  std::vector<double> detailed_cost;  //< Expected nodes per level

  /**
     The default constructor means no pruning.
  */

  Pruning() : radius_factor(1.), expectation(1.), metric(PRUNER_METRIC_PROBABILITY_OF_SHORTEST){};

  /** Set all pruning coefficients to 1, except the last <level>
      coefficients, these will be linearly with slope `-1 /
      block_size`.

      @param level number of levels in linear descent
  */

  static Pruning LinearPruning(int block_size, int level);
};

/**
   A strategy covers pruning parameters and preprocessing block_sizes
*/

class Strategy
{
public:
  size_t block_size;                         //< block size
  vector<Pruning> pruning_parameters;        //< Pruning parameters
  vector<size_t> preprocessing_block_sizes;  //< For each block size we run one tour

  /** Construct an empty strategy

      @note Use this instead of the default constructor. The default
      constructor does not add default pruning parameters.

   */

  static Strategy EmptyStrategy(size_t block_size)
  {
    Strategy strat;
    strat.block_size = block_size;
    strat.pruning_parameters.emplace_back(Pruning());
    return strat;
  };

  /** Select the best pruning parameters for the input `radius`. The
      parameter `gh` is used to establish the ratio between `radius`
      and the Gaussian heuristic, which is used for sizes.

     @param radius radius of the currently shortest vector
     @param gh Gaussian heuristic prediction for radius

   */

  const Pruning &get_pruning(double radius, double gh) const;
};

class BKZParam
{
public:
  /**
     @brief Create BKZ parameters

     @param block_size
        block size for the reduction
     @param strategies
        vector of strategies used for pruning and preprocessing
     @param delta
        LLL parameter delta
     @param flags
        various flags that can be arbitrarily combined (using |):
          - BKZ_VERBOSE       print additional information during reduction
          - BKZ_NO_LLL        do not run LLL before block reduction (use at your own risk)
          - BKZ_MAX_LOOPS     terminate after max_loops iterations
          - BKZ_MAX_TIME      terminate after max_time time
          - BKZ_BOUNDED_LLL   only run LLL in current block during SVP preprocessing (use at your
     own
     risk)
          - BKZ_AUTO_ABORT    heuristically terminate the reduction if progress stalls
          - BKZ_DUMP_GSO      after every iteration write the shape of the current basis to a file
          - BKZ_GH_BND        use the Gaussian heuristic to reduce the enumeration bound of possible
          - BKZ_SD_VARIANT    run SD-BKZ
          - BKZ_SLD_RED       run slide reduction
     @param max_loops
        maximum number of loops (or zero to disable this)
     @param max_time
        maximum number of time  (or zero to disable this)
     @param auto_abort_scale
        auto abort when next tour does not improve slope over `scale`* previous tour
     @param auto_abort_max_no_dec
        auto abort when next tour does not improve slope `no_dec` times
     @param gh_factor
        set enumeration bound to Gaussian heuristic times `gh_factor`
     @param min_success_probability
        minimum success probability in an SVP reduction (when using pruning)
     @param rerandomization_density
        the heavier rerandomization, the better our guarantees and costs
  */

  BKZParam(int block_size, vector<Strategy> &strategies, double delta = LLL_DEF_DELTA,
           int flags = BKZ_DEFAULT, int max_loops = 0, double max_time = 0,
           double auto_abort_scale        = BKZ_DEF_AUTO_ABORT_SCALE,
           int auto_abort_max_no_dec      = BKZ_DEF_AUTO_ABORT_MAX_NO_DEC,
           double gh_factor               = BKZ_DEF_AUTO_ABORT_MAX_NO_DEC,
           double min_success_probability = BKZ_DEF_MIN_SUCCESS_PROBABILITY,
           int rerandomization_density    = BKZ_DEF_RERANDOMIZATION_DENSITY)
      : block_size(block_size), strategies(strategies), delta(delta), flags(flags),
        max_loops(max_loops), max_time(max_time), auto_abort_scale(auto_abort_scale),
        auto_abort_max_no_dec(auto_abort_max_no_dec), gh_factor(gh_factor),
        dump_gso_filename("gso.json"), min_success_probability(min_success_probability),
        rerandomization_density(rerandomization_density)
  {

    // we create dummy strategies
    if (strategies.empty())
    {
      strategies = vector<Strategy>();
      for (long b = 0; b <= block_size; ++b)
      {
        strategies.emplace_back(Strategy::EmptyStrategy(b));
      }
    }
  };

  /** Block size used for enumeration **/
  int block_size;

  /** Strategies (pruning coefficients, preprocessing)  */
  vector<Strategy> &strategies;

  /** LLL parameter delta **/
  double delta;

  /** See BKZFlags **/
  int flags;

  /** Maximum number of loops to execute **/
  int max_loops;

  /** Maximum time to spend **/
  double max_time;

  /** If BKZ_AUTOABORT is set, We abort if `new_slope < auto_abort_scale * old_slope`
      is true for `auto_abort_max_no_dec` loops.
   */
  double auto_abort_scale;
  int auto_abort_max_no_dec;

  /** If BKZ_GH_BND is set, the enumeration bound will be set to gh_factor times
      the Gaussian Heuristic
  */
  double gh_factor;

  /** If BKZ_DUMP_GSO is set, the norms of the GSO matrix are written to this
      file after each complete round.
  */
  string dump_gso_filename;

  /** minimum success probability when using extreme pruning */

  double min_success_probability;

  /** density of rerandomization operation when using extreme pruning **/

  int rerandomization_density;
};

/**
   Load BKZ pruning and preprocessing strategies from a json file.

   @note All parameters except pruning and preprocessing are silently ignored.
*/

vector<Strategy> load_strategies_json(const std::string &filename);

/**
   Return default directory to search for strategies (if any)
*/

const std::string &default_strategy_path();

/**
   Return default strategy (if any)
*/

const std::string &default_strategy();

/**
   Attempt to expand strategy path to full path
*/

const std::string strategy_full_path(const string &strategy_path);

FPLLL_END_NAMESPACE
#endif /* BKZ_PARAM_H */
