#ifndef FPLLL_PRUNER_H
#define FPLLL_PRUNER_H

/* Copyright (C) 2015-2017 Martin Albrecht, Leo Ducas.

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTAbILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#include <vector>
#include "defs.h"

FPLLL_BEGIN_NAMESPACE

#define PRUNER_MAX_N 2047

/**
   Pruning parameters for one radius (expressed as a ratio to the Gaussian heuristic)
 */

class PruningParams
{

public:
  double gh_factor;                  //< radius^2/Gaussian heuristic^2
  std::vector<double> coefficients;  //< pruning coefficients
  double expectation;                //< either expected success probability or number of solutions
  /**
      metric used for optimisation (success probability or number of solutions)
   */
  PrunerMetric metric;
  std::vector<double> detailed_cost;  //< Expected nodes per level

  /**
     The default constructor means no pruning.
  */

  PruningParams() : gh_factor(1.), expectation(1.), metric(PRUNER_METRIC_PROBABILITY_OF_SHORTEST){};

  /** Set all pruning coefficients to 1, except the last <level>
      coefficients, these will be linearly with slope `-1 /
      block_size`.

      @param level number of levels in linear descent
  */

  static PruningParams LinearPruningParams(int block_size, int level);
};



/**
   @brief Search for optimal pruning parameters.

   Wrapping function for the Pruner class. The function is templated by a Floating Point Type, which
   are used for internal computations.

   Calling M the metric function (Success proba or expectation of number of results), C the cost
   function of enumeration, P the preproc_cost, this function tries to provide a pruning vector p
   essentially minimizing (C(p) + P) / M(p). The pruning vector is expressed as relative (squared)
   length with respect to the enumeration radius, i.e. the entries of p are in the interval (0,1),
   and are always decreasing.

   Small amendents are made to also take account for the desired target value of the metric after
   repeating several enumeration. In particular the chosen pruning value should not excess the
   target, and a bit more enumeration may be chosen if it get is to the the desired result without
   any retrial.

   The cost C of enumeration under the gaussian heuristic, and depends on gso_r. The function S is
   define by the metric parameter (and depend on gso_r for EXPECTED_SOLUTIONS)

   This function provides access to three algorithm to optimize those parameter

   - greedy                (deactivated by setting flag PRUNER_START_FROM_INPUT)
   - gradient decent       (activated by setting flag PRUNER_GRADIENT)
   - Nelder-Mead Algorithm (activated by setting flag PRUNER_NELDER_MEAD)

   If several algortihm are activated, they will all be ran, using the output of the previous as
   there starting point.

   Greedy shoud be very fast, and essentially tries to produce pruning parameters generating a flat
   enumeration tree. It is meant to quickly provide pruning parameters, say up to block 40. It is
   also used as a decent starting point for the gradient descent to be quicker. It can be argued
   that this strategy is a factor O(n) away from optimal, probably less in practice.

   The gradient descent is rather standard. The gradient is computed numerically, and then several
   steps are made keeping the same direction until a local minima in this direction is found (the
   rational is that computing the gradient is 2n times more expensive than test a point). It is to
   be expected that this method suffers from convergence and numerical instability especially due
   to numerical differentiation. It remains rather fast and seems to properly handle basis up to
   dimension LLL reduced-basis up to dim ~70, and more if the basis is strongly reduced.

   The Nelder-Mead algorithm is a descent algorithm praised for its robustness. It's pretty cool,
   you should really read about it: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method . In
   brief, it maintains a simplex of points and each step consist in improving the worst of them by
   reflexion expansion or contraction, or shrinking the all simplex if the above failed. This
   method is somehow aware not only of the gradient, but also the curvature of the optimized
   function. Thanks to contraction and expansion, it can slowly "pass through the eye of a needle,"
   and re-accelerate the descent after that.

   Maybe you don't care so much of all this... I'm just procratinating low-level documentation.

   @param pruning Output of the function. Also used as an input if PRUNER_START_FROM_INPUT is on.
   @param enumeration_radius radius of the desired enumeration
   @param preproc_cost cost of preprocessing (i.e. additive cost for a retrying an enumeration)
   @param gso_r vector of Gram-Schmidt length (squared) of the enumerated block
   @param target desired target success probability/expected solutions after all retrial.
   @param metric metric is to be optimized : PRUNER_METRIC_PROBABILITY_OF_SHORTEST or
   PRUNER_METRIC_EXPECTED_SOLUTIONS
   @param flags complementary parameters : PRUNER_CVP  PRUNER_START_FROM_INPUT  PRUNER_GRADIENT
   PRUNER_NELDER_MEAD  PRUNER_VERBOSE

   @return Return value is done by reference on parameter pruning.
*/

template <class FT>
void prune(/*(input)output*/ PruningParams &pruning,
           /*inputs*/
           const double enumeration_radius, const double preproc_cost, const vector<double> &gso_r,
           const double target       = .9,
           const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
           const int flags           = PRUNER_GRADIENT);

/**
   @brief Search for optimal Pruning parameters, averaging over several basis

   Wrapping function for the Pruner class.

   Same as prune(), but averaging cost over several bases (without slowing down the search)

   @param pruning Output of the function. Also used as an input if PRUNER_START_FROM_INPUT is on.
   @param enumeration_radius radius of the desired enumeration
   @param preproc_cost cost of preprocessing (i.e. additive cost for a retrying an enumeration)
   @param gso_rs vector of vector Gram-Schmidt length (squared) of the enumerated block over
   several bases
   @param target desired target success probability/expected solutions after all retrial.
   @param metric metric is to be optimized : PRUNER_METRIC_PROBABILITY_OF_SHORTEST or
   PRUNER_METRIC_EXPECTED_SOLUTIONS
   @param flags complementary parameters : PRUNER_CVP  PRUNER_START_FROM_INPUT  PRUNER_GRADIENT
   PRUNER_NELDER_MEAD  PRUNER_VERBOSE

*/
template <class FT>
void prune(/*output*/ PruningParams &pruning,
           /*inputs*/
           double enumeration_radius, const double preproc_cost,
           const vector<vector<double>> &gso_rs, const double target = .9,
           const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
           const int flags           = PRUNER_GRADIENT);

/**
   @brief Probability that a pruned enumeration indeed returns the shortest vector.

   Probability that a pruned enumeration indeed returns the shortest vector. Wrapping function for
   the Pruner class.

   @param pruning A pruning object. Only uses the coefficient field. (coefficients must start with
   two 1s, be decreasing and in the interval (0,1))

   @return probability
*/
template <class FT> FT svp_probability(const PruningParams &pruning);

/**
   @brief Probability that a pruned enumeration indeed returns the shortest vector.

   @param pr A vector of pruning bounds (coefficients must start with two 1s, be decreasing and in
   the interval (0,1))

   @return probability
*/
template <class FT> FT svp_probability(const vector<double> &pr);

/**
   @brief Pruner class, to compute and optimize cost and success probability of pruned enumeration.

   The class is templated by a Floating Point Type, which are used for internal computations.

   This class provides an implementation of the numerically optimized pruning strategy of [Gama
   Nguyen Regev 2010].

   Many details for implementation follows from the thesis [Chen 2013] Some simplifications have
   been made, for example we restrict ourselves to ``even vector bounds'' i.e., the bound at indices
   2i and 2i+1 are kept equal, as to allow volume estimation by symbolic integration (see
   [GNR10,Chen13])

   naming conventions:
   - b is for bound (squared)
   - ipv is for inverse partial volumes (NOT squared)
   - r is for gram schmidt length (squared). Automatically renormalized to avoid
   overflowing partial volumes
   - p is for polynomial

   inside this code, b,pv,and r are in reversed order as to conform with the algorithm description
   of [Chen13] reversing the order is handled by the load and save methods.

   - n is for the dimension of the basis to prune
   - d is for degrees of polynomials. md is for max_degree
   - d is floor(n/2 at most). Odd n are dealt with by ignoring the first component
*/
template <class FT> class Pruner
{
public:
  class TestPruner;
  friend class TestPruner;

  /**
   @brief (limited) Constructor for Pruner. Only for svp_probability.
   @param n dimension of the pruned block
  */
  explicit Pruner(const int n) : metric(PRUNER_METRIC_PROBABILITY_OF_SHORTEST), flags(0), n(n)
  {
    verbosity = flags & PRUNER_VERBOSE;
    set_tabulated_consts();
    d = n / 2;
    min_pruning_coefficients.resize(d);
    fill(min_pruning_coefficients.begin(), min_pruning_coefficients.end(), 0.);
  }

  /**
    @brief Constructor for Pruner.

    Algorithmic details are given for the wrapping function prune().

    @param enumeration_radius radius of the desired enumeration
    @param preproc_cost cost of preprocessing (i.e. additive cost for a retrying an enumeration)
    @param gso_r  vector Gram-Schmidt length (squared) of the enumerated block over several bases
    @param target desired target success probability/expected solutions after all retrial.
    @param metric metric is to be optimized : PRUNER_METRIC_PROBABILITY_OF_SHORTEST or
    PRUNER_METRIC_EXPECTED_SOLUTIONS
    @param flags complementary parameters : PRUNER_CVP  PRUNER_START_FROM_INPUT  PRUNER_GRADIENT
    PRUNER_NELDER_MEAD  PRUNER_VERBOSE
  */

  Pruner(const FT enumeration_radius, const FT preproc_cost, const vector<double> &gso_r,
         const FT target = 0.9, const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
         int flags = PRUNER_GRADIENT)
      : enumeration_radius(enumeration_radius), preproc_cost(preproc_cost), target(target),
        metric(metric), flags(flags)
  {
    verbosity = flags & PRUNER_VERBOSE;

    n = gso_r.size();
    d = n / 2;
    if (flags & PRUNER_CVP)
    {
      symmetry_factor = 1;
    }
    min_pruning_coefficients.resize(d);
    fill(min_pruning_coefficients.begin(), min_pruning_coefficients.end(), 0.);
    set_tabulated_consts();
    load_basis_shape(gso_r);
  }

  /**
     @brief Constructor for Pruner with multiple basis (averaging the cost function over all of
     them)

     Algorithmic details are given for the wrapping function prune().

     @param enumeration_radius radius of the desired enumeration
     @param preproc_cost cost of preprocessing (i.e. additive cost for a retrying an enumeration)
     @param gso_rs vector of vector Gram-Schmidt length (squared) of the enumerated block over
     several bases
     @param target desired target success probability/expected solutions after all retrial.
     @param metric metric is to be optimized : PRUNER_METRIC_PROBABILITY_OF_SHORTEST or
     PRUNER_METRIC_EXPECTED_SOLUTIONS
     @param flags complementary parameters : PRUNER_CVP  PRUNER_START_FROM_INPUT  PRUNER_GRADIENT
     PRUNER_NELDER_MEAD  PRUNER_VERBOSE
  */

  Pruner(const FT enumeration_radius, const FT preproc_cost, const vector<vector<double>> &gso_rs,
         const FT target = 0.9, const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
         int flags = PRUNER_GRADIENT)
      : enumeration_radius(enumeration_radius), preproc_cost(preproc_cost), target(target),
        metric(metric), flags(flags)
  {
    verbosity = flags & PRUNER_VERBOSE;
    n         = gso_rs[0].size();
    d         = n / 2;
    if (flags & PRUNER_CVP)
    {
      symmetry_factor = 1;
    }
    min_pruning_coefficients.resize(d);
    fill(min_pruning_coefficients.begin(), min_pruning_coefficients.end(), 0.);

    set_tabulated_consts();
    load_basis_shapes(gso_rs);
  }

  /**
     Run the optimization process, successively using the algorithm activated

     @brief run the optimization process
  */
  void optimize_coefficients(/*io*/ vector<double> &pr);

  /** @brief Compute the cost of a single enumeration */

  double single_enum_cost(/*i*/ const vector<double> &pr, vector<double> *detailed_cost = nullptr)
  {
    evec b(n);
    load_coefficients(b, pr);
    return single_enum_cost(b, detailed_cost).get_d();
  }

  /**
     @brief Cost of repeating enumeration and preprocessing until reaching target

     Compute the cost of r enumeration and (r-1) preprocessing, where r is the required number of
     retrials to reach target.
  */
  double repeated_enum_cost(/*i*/ const vector<double> &pr)
  {
    evec b(n);
    load_coefficients(b, pr);
    return repeated_enum_cost(b).get_d();
  }

  /**
     @brief Compute the success proba of a single enumeration
  */
  double measure_metric(/*i*/ const vector<double> &pr)
  {
    evec b(n);
    load_coefficients(b, pr);
    return measure_metric(b).get_d();
  }

  /**
     @brief gaussian heuristic for the shortest vector length of the loaded basis.

     @note If multiple loaded then it is the gaussian heuristic of the first one.
  */
  FT gaussian_heuristic();

private:
  FT enumeration_radius;
  FT preproc_cost;
  FT target;
  PrunerMetric metric;
  bool shape_loaded = false;  ///< Keep track of whether a basis was loaded. Only allows
                              /// single_enum_cost() with metric PROBABILITY_OF_SHORTEST if not
  int flags;
  int n;  ///< Dimension of the block to be enumerated
  int d;  ///< Half dimension: most function only use half dimension, as our Pruner implcitly
          /// enforce parity
  vector<FT> min_pruning_coefficients;

  double descent_starting_clock;
  static bool tabulated_values_imported;  ///< (static) has tabulated constant been imported
  static FT tabulated_factorial[PRUNER_MAX_N];
  static FT tabulated_ball_vol[PRUNER_MAX_N];

  FT epsilon  = std::pow(2., -7);  ///< Epsilon to use for numerical differentiation
  FT min_step = std::pow(2., -6);  ///< Minimal step in a given direction
  FT min_cf_decrease =
      .995;  //< Maximal ratio of two consectuive cost_factor in the descent before stopping
  FT step_factor     = std::pow(2, .5);  ///< Increment factor for steps in a given direction
  FT shell_ratio     = .995;             ///< Shell thickness Ratio when evaluating svp proba
  FT symmetry_factor = .5;               ///< Account for symmetry in costs: 0.5 for SVP, 1 for CVP.

  // Following typedefs are given for readability
  using vec  = vector<FT>;  ///< Vectors of dimension n
  using evec = vector<FT>;  ///< Vectors of dimension d=n/2
  using poly = vector<FT>;  ///< Polynomials of degree at most d, as vectors of dimension d+1

  vec r;                    ///< Gram-Schmidt length (squared, and reversed ordered)
  vec ipv;                  ///< Inverse Partial Volumes of projected sublattices (not squared !)
  FT normalization_factor;  ///< internal renormalization factor to avoid over/underflows
  FT normalized_radius;     ///< renormalized enumeration radius
  int verbosity = 0;

  /** @brief Static of load tabulated constants (done once over all instances). */
  void set_tabulated_consts();

  /** @brief load the shape of a basis from vector<double>.
      @param gso_r vector of gso coefficients
      @param reset_normalization boolean triggering computation of the renormalization factor
    */
  void load_basis_shape(const vector<double> &gso_r, bool reset_normalization = true);

  /**
      @brief Load the shapes of many bases from vector<vector<double>>.

      Costs will be averaged over all bases. Normalization is done according to the first basis.

      @param gso_r vector of vector of gso coefficients
    */
  void load_basis_shapes(const vector<vector<double>> &gso_rs);

  /**
     @brief convert pruning coefficient from external to internal format.

     Convert type, reverse the ordering, and only select coefficent at even position

     @param b Output in internal format
     @param pr Input in external format
  */
  void load_coefficients(/*o*/ evec &b, /*i*/ const vector<double> &pr);

  /**
     @brief Convert pruning coefficient from internal to external.

     Convert type, reverse the ordering, and repeat each coefficent twice

     @param pr Output in external format
     @param b Input in internal format
  */
  void save_coefficients(/*o*/ vector<double> &pr, /*i*/ const evec &b);

  /**
     @brief Enforce contraints on pruning coefficients

     Enforce that pruning coeffient terminate with a 1, are decreasing, and are not too close to 0
     (for numerical stability).

     @param b input/output
     @param j Keeps index j unchanged when possible
     @return was a constraint violated ?
  */

  inline bool enforce(/*io*/ evec &b, /*opt i*/ const int j = 0);

  /**
     @brief Evaluate a polynomial
     @param ld degree of the polynomial p
     @param p the polynomial
     @param x value at which to evaluate
     @return p evaluated at x
  */
  inline FT eval_poly(const int ld, /*i*/ const poly &p, const FT x);

  /**
     @brief Integrate a polynomial (in place)
     @param ld degree of the polynomial p
     @param p the polynomial
     @return p evaluated at x
  */
  inline void integrate_poly(const int ld, /*io*/ poly &p);

  /**
     @brief Compute the relative volume of a cylinder intersection to the unit sphere (of a given
     subdimension)
     @param rd sub-dimension
     @param b bounds of the cylinder intersection
     @return relative volume
     @note Implicitly renormalized as if b[rd-1] = 1
  */
  inline FT relative_volume(/*i*/ const int rd, const evec &b);

  /**
     @brief Compute the cost of a pruned enumeration
     @param b pruning bounds
     @param detailed_cost (optional) store the details of node per level of the enum tree there
     @return cost, in enumeration node
  */
  FT single_enum_cost(/*i*/ const evec &b, vector<double> *detailed_cost = nullptr);

  /**
     @brief Compute the success probability for SVP/CVP of a single enumeration
     @param b pruning bounds
     @return success probability
  */

  FT svp_probability(/*i*/ const evec &b);
  /**
     @brief Compute the expected number of solution of a single of a single enumeration
     @param b pruning bounds
     @return expected number of solutions
  */
  FT expected_solutions(/*i*/ const evec &b);

  /**
     @brief links to either svp_probability() or expected_solution() depending on metric
     @param b pruning bounds
     @return expected number of solutions or success probability
  */
  FT measure_metric(/*i*/ const evec &b);

  /**
     Compute the cost of r enumeration and (r-1) preprocessing,
     where r is the required number of retrials to reach target/target_solution
     @brief cost of repeating enumeration and preprocessing
     @param b pruning bounds
     @return cost
  */

  FT repeated_enum_cost(/*i*/ const evec &b);

  /**
     Compute the gradient cost of r enumeration and (r-1) preprocessing,
     where r is the required number of retrials to reach target/target_solution
     @brief gradient of the cost of repeating enumeration and preprocessing
     @param b pruning bounds
     @param res reference for output
     @return cost
  */
  void repeated_enum_cost_gradient(/*i*/ const evec &b, /*o*/ evec &res);

  /**
     @brief gradient of the cost of repeating enumeration and preprocessing
     @param b reference for output
     @note Greedy (with less cost) is also used as a heuristic to lower bound
     pruning coeffients during the descent and avoid numerical stability issues.
  */
  void greedy(evec &b);

  /**
     @brief Perform one step of the gradient descent, in place
     @param b input/output
  */
  int gradient_descent_step(/*io*/ evec &b);

  /**
     @brief Runs nelder-mead algorithm
     @param b input/output
     @note As recommended on Wikipedia, the algorithm is re-ran several time
  */
  int nelder_mead_step(/*io*/ evec &b);
};

template <class FT>
bool Pruner<FT>::tabulated_values_imported = false;  ///< (static) tabulated value not loaded yet
template <class FT> FT Pruner<FT>::tabulated_factorial[PRUNER_MAX_N];
template <class FT> FT Pruner<FT>::tabulated_ball_vol[PRUNER_MAX_N];

FPLLL_END_NAMESPACE

#endif /* FPLLL_PRUNER_H */
