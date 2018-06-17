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

#ifndef FPLLL_PRUNER_H
#define FPLLL_PRUNER_H

#include "fplll/defs.h"
#include "fplll/lll.h"
#include <vector>

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
 * @brief Performs pruning using PruningParams object with specific float type FT.
 *
 * @param B
 *    basis of the lattice to be reduced
 * @param param
 *    parameter object
 * @param sel_ft
 *    specifies the float type used for GSO computations
 * @param precision
 *    specifies the precision if sel_ft=FT_MPFR (and needs to be > 0 in that case)
 * @param prune_start
 *    start index of pruning (first index being 0)
 * @param prune_end
 *    (prune_end-1) is end index of pruning
 * @param prune_pre_nodes
 *     preprocessing cost in the number of nodes
 * @param prune_min_prob
 *     target probability. If it is -1, it will optimize
 *       single_enum_cost/succ. prob while not fixing the succ. prob.
 *       If it is > 0, it will fix the succ. prob and optimize the
 *       single_enum_cost.
 * @param gh_factor
 *      input GH factor to compute the enumeration radius. The computed
 *        enumeration radius will be min(GH*gh_factor, |b_i*|).
 * @return
 *    the status of the pruning
 */
template <class FT>
int run_pruner_f(ZZ_mat<mpz_t> &B, const PruningParams &param, int sel_ft, int precision = 0,
                 int prune_start = 0, int prune_end = 1, double prune_pre_nodes = 1e6,
                 double prune_min_prob = -1, double gh_factor = 1.0);

/**
 * @brief Performs pruning using PruningParams object.
 *
 * @param B
 *    basis of the lattice to be reduced
 * @param param
 *    parameter object
 * @param float_type
 *    specifies the data type used for GSO computations (see defs.h for options)
 * @param precision
 *    specifies the precision if float_type=FT_MPFR (and needs to be > 0 in that case)
 *    ignored otherwise
 * @param prune_start
 *    start index of pruning (first index being 0)
 * @param prune_end
 *    (prune_end-1) is end index of pruning
 * @param prune_pre_nodes
 *     preprocessing cost in the number of nodes
 * @param prune_min_prob
 *     target probability. If it is -1, it will optimize
 *       single_enum_cost/succ. prob while not fixing the succ. prob.
 *       If it is > 0, it will fix the succ. prob and optimize the
 *       single_enum_cost.
 * @param gh_factor
 *      input GH factor to compute the enumeration radius. The computed
 *        enumeration radius will be min(GH*gh_factor, |b_i*|).
 * @return
 *    the status of the prunign
 */
int run_pruner(ZZ_mat<mpz_t> &B, const PruningParams &param, FloatType float_type = FT_DEFAULT,
               int precision = 0, int prune_start = 0, int prune_end = 1,
               double prune_pre_nodes = 1e6, double prune_min_prob = -1, double gh_factor = 1.0);

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
    btmp.resize(d);
    bftmp.resize(n);
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
    btmp.resize(d);
    bftmp.resize(n);
    fill(min_pruning_coefficients.begin(), min_pruning_coefficients.end(), 0.);
    set_tabulated_consts();

    // need to fix target if possible
    if (metric == PRUNER_METRIC_PROBABILITY_OF_SHORTEST)
    {
      // if target > 1 or < 0, optimize overall cost, use
      //   0.9 as target in repeated_enum_cost().
      if (this->target > 1.0 || this->target < 0.0)
      {
        this->target = 0.99;
        opt_overall  = true;
      }
    }
    else if (metric == PRUNER_METRIC_EXPECTED_SOLUTIONS)
    {
      if (this->target < 0.0)
      {
        this->target = 0.99;
        opt_overall  = true;
      }
    }
    else
    {
      throw std::invalid_argument("Pruner was set to an unknown metric");
    }
    // normalization on |b_i^*|
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
    btmp.resize(d);
    bftmp.resize(n);
    fill(min_pruning_coefficients.begin(), min_pruning_coefficients.end(), 0.);
    set_tabulated_consts();
    // need to fix target if possible
    if (metric == PRUNER_METRIC_PROBABILITY_OF_SHORTEST)
    {
      // if target > 1 or < 0, optimize overall cost, use
      //   0.99 as target for the computation in repeated_enum_cost().
      if (this->target > 1.0 || this->target < 0.0)
      {
        this->target = 0.99;
        opt_overall  = true;
      }
    }
    else if (metric == PRUNER_METRIC_EXPECTED_SOLUTIONS)
    {
      // if target < 0, optimize overall cost, use
      //   0.99 as target for the computation in repeated_enum_cost().
      if (this->target < 0.0)
      {
        this->target = 0.99;
        opt_overall  = true;
      }
    }
    else
    {
      throw std::invalid_argument("Pruner was set to an unknown metric");
    }
    load_basis_shapes(gso_rs);
  }

  /**
     @brief Main interface to optimize pruning coefficients

     Main interface to optimize pruning coefficients.
     It will in turn invoke either of the two functions:
      (1)   optimize_coefficients_cost() or
      (2)   optimize_coefficients_prob()
     depending on the input "target".

     If the target is negative (e.g. -1), it calls function (1)
     such that goal is to optimize the
     single_enum_cost(pr) * trials + preproc_cost * (trials - 1.0).

     If the target is > 0, it calls function (2) such that the
     goal is to optimize the single_enum_cost(pr) while fixing
     the succ. probability == target.
  */
  void optimize_coefficients(/*io*/ vector<double> &pr);

  /**
     @brief main interface to optimize the pruning coefficients with
     respect to the overall enumeraiton time.

     Main interface to optimize the overall enumeraiton time where the
     target function is:

     single_enum_cost(pr) * trials + preproc_cost * (trials - 1.0);
  */
  void optimize_coefficients_cost(/*io*/ vector<double> &pr);

  /**
     @brief main interface to optimize the pruning coefficients with
     repeect to the single enumeraiton time while fixing the succ. prob
     or expected number of solutions.

     Main interface to optimize the single enumeration time with the
     constraint such that the succ. prob (or expected solutions) is
     fixed (and given) from input to the Pruner constructor.
  */
  void optimize_coefficients_prob(/*io*/ vector<double> &pr);

  /**
     @brief run the optimization process using 'even' coefficients. Note
     the optimization only applies to the pruning coefficients indexed
     by (0, 2, 4, ... n). It thus uses half of the coefficients.

     Run the optimization process, successively using the algorithm activated
     using using half coefficients: the input pr has length n; but only the
     even indices in the vector will be used in the optimization.
     In the end, we have pr_i = pr_{i+1}.
     Note it only optimize the overall enumeraiton time where the
     target function is:
     single_enum_cost(pr) * trials + preproc_cost * (trials - 1.0);
  */
  void optimize_coefficients_evec(/*io*/ vector<double> &pr);

  /**
     @brief run the optimization process using all the coefficients.

     Run the optimization process, successively using the algorithm activated
     using using full coefficients. That is, we don't have the
     constraint pr_i = pr_{i+1} anymore in this function.
     Note it only optimize the overall enumeraiton time where the
     target function is:
     single_enum_cost(pr) * trials + preproc_cost * (trials - 1.0);

  */
  void optimize_coefficients_full(/*io*/ vector<double> &pr);

  /**
      @brief Compute the cost of a single enumeration
  */
  double single_enum_cost(/*i*/ const vector<double> &pr, vector<double> *detailed_cost = nullptr)
  {
    // cout << "# pr is --> " << pr[0] << endl;
    evec b(d);
    load_coefficients(b, pr);

    return single_enum_cost(b, detailed_cost).get_d();
  }

  /**
     Compute the cost of r enumeration and (r-1) preprocessing,
     where r is the required number of retrials to reach target/target_solution
     @brief cost of repeating enumeration and preprocessing
     @param b pruning bounds
     @return cost
  */
  double repeated_enum_cost(/*i*/ const vector<double> &pr)
  {
    vec b(n);
    load_coefficients(b, pr);
    return repeated_enum_cost(b).get_d();
  }

  /**
     @brief Compute the success proba of a single enumeration
  */
  double measure_metric(/*i*/ const vector<double> &pr)
  {
    vec b(n);
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
  bool opt_overall = false;

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

  /* normalised by 1/det^(1/n) */
  vec r;                    ///< Gram-Schmidt length (squared, and reversed ordered)
  vec ipv;                  ///< Inverse Partial Volumes of projected sublattices (not squared !)
  FT normalization_factor;  ///< internal renormalization factor to avoid over/underflows
  FT normalized_radius;     ///< renormalized enumeration radius
  int verbosity = 0;
  vec r_old;  ///< a copy of gso_r)
  FT logvol;

  /* temp. variables used internally */
  vector<FT> btmp;
  vector<FT> bftmp;

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
     @brief auxiliary function to print the coefficients.
  */
  void print_coefficients(/*i*/ const vector<double> &pr);

  /**
     @brief auxiliary function to print the coefficients.
  */
  void print_coefficients(/*i*/ const evec &pr);

  /**
     @brief Convert pruning coefficient from internal to external.

     Convert type, reverse the ordering, and repeat each coefficent twice

     @param pr Output in external format
     @param b Input in internal format
  */
  void save_coefficients(/*o*/ vector<double> &pr, /*i*/ const evec &b);

  /**
     @brief Enforce contraints on pruning coefficients

     Enforce that pruning coeffient terminate with a 1, are decreasing, and
     are not too close to 0 (for numerical stability).

     @param b input/output
     @param j Keeps index j unchanged when possible
     @return was a constraint violated ?
  */
  inline bool enforce(/*io*/ vec &b, /*opt i*/ const int j = 0);

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

     Function single_enum_cost() is the main interface to compute single enumeration cost.
     Function single_enum_cost_lower() is an auxiliary function to compute the
               lower bound of single enumation cost.
     Function single_enum_cost_upper() is an auxiliary function to compute the
               upper bound of single enumation cost.
     Function single_enum_cost_evec() is the underlying backend for the above two.
  */
  FT single_enum_cost(/*i*/ const vec &b, vector<double> *detailed_cost = nullptr);
  FT single_enum_cost_evec(/*i*/ const evec &b, vector<double> *detailed_cost = nullptr);
  FT single_enum_cost_lower(/*i*/ const vec &b, vector<double> *detailed_cost = nullptr);
  FT single_enum_cost_upper(/*i*/ const vec &b, vector<double> *detailed_cost = nullptr);

  /**
     @brief Compute the success probability for SVP/CVP of a single enumeration
     @param b pruning bounds
     @return success probability

     Function svp_probability() is the main interface.
     Function svp_probability_lower() is an auxiliary function to compute the
                                      approx. lower bound of succ. probability
     Function svp_probability_upper() is an auxiliary function to compute the
                                      upper bound of succ. probability
     The svp_probability_evec() is the underlying backend for above.
  */
  FT svp_probability(/*i*/ const vec &b);
  FT svp_probability_evec(/*i*/ const evec &b);
  FT svp_probability_lower(/*i*/ const vec &b);
  FT svp_probability_upper(/*i*/ const vec &b);

  /**
     @brief Compute the expected number of solution of a single of a single enumeration
     @param b pruning bounds
     @return expected number of solutions

     The naming conventions is the same as those in svp_probability_*()
  */
  FT expected_solutions(/*i*/ const vec &b);
  FT expected_solutions_evec(/*i*/ const evec &b);
  FT expected_solutions_lower(/*i*/ const vec &b);
  FT expected_solutions_upper(/*i*/ const vec &b);

  /**
     @brief links to either svp_probability() or expected_solution() depending on metric
     @param b pruning bounds
     @return expected number of solutions or success probability
  */
  FT measure_metric(/*i*/ const vec &b);

  /**
     @brief Compute the target function to be optimized. This could be the
     cost of repeating enumeration and preprocessing until reaching target
     or the cost of a single enumeration.

     Note this function does not do optimization. Depending on the
     optimization goal, it will return the cost based on different
     models/modes. There are two modes depending on the target.

     (1) if target == -1, it will compute the cost of r * single_enum_cost()
         plus (r-1) * preprocessing where r is the expected number of
         trials.

     (2) if target != -1 (target must also be given), it will compute
         the cost of a single_enum_cost(). Note in such case the input
         target prob is fixed. Whether we could achieve this input
         target prob will be determined in the optimization procedure.
         See also below.

     Note one should set the radius with care, since if the radius is too
     small, perhaps with all 1's the pruning coefficients are not able
     to support a target expectation/probability. Consider an example
     if using PRUNER_METRIC_EXPECTED_SOLUTIONS:

     (a) If the radius is large enough to support target many solutions,
     the optimizer will optimize the single_enum_cost while fixing
     the target.

     (b) if the radius is too small to support target many solutions,
     it will likely return all 1's in the coefficients (attempt
     to achieve the target but which is not possible)

     Note the radius is either from the first b_i^* of the block lattice;
     or set as min(GH * gh_ratio, |b_i^*|). Hence it will not be arbitrarily
     large.
  */
  FT target_function(/*i*/ const evec &b);

  /**
     @brief Compute the gradient cost of the target function to be optimized
     @param b pruning bounds
     @param res reference for output
     @return cost
  */
  void target_function_gradient(/*i*/ const vec &b, /*o*/ vec &res);

  /**
     Compute the cost of r enumeration and (r-1) preprocessing,
     where r is the required number of retrials to reach target/target_solution
     @brief cost of repeating enumeration and preprocessing
     @param b pruning bounds
     @return cost
  */
  FT repeated_enum_cost(/*i*/ const evec &b);

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
  int gradient_descent_step(/*io*/ vec &b);

  /**
     @brief Perform several steps of the gradient descent, in place
            It calls the gradient_descent_step() for several times.
     @param b input/output
  */
  int gradient_descent(/*io*/ vec &b);

  /**
     @brief Runs nelder-mead algorithm
     @param b input/output
     @note As recommended on Wikipedia, the algorithm is re-ran several time
  */
  int nelder_mead_step(/*io*/ evec &b);

  /**
     @brief tune the pruning parameter to reduce single enumeration time.
     The purpose is to reduce the single enumeration time while not
     decreasing the succ. probability too much.

     Optimization process to the pruning parameters to reduce single enumeration
     time with the hope that it also reduces the overall enumeration time. Local
     implies the code will hopefully only do small changes to the pruning
     coefficients.
  */
  void optimize_coefficients_local_tune_single_enum(/*io*/ vector<double> &pr);

  /**
     @brief tune the pruning parameter to increase succ. probability.
     The purpose is to increase the succ. probability while not incresing
     the single enumeration time too much.

     Optimization process to the pruning parameters to increase succ.
     probability with the restriction that the single enumeration time does
     not increase significantly.
  */
  void optimize_coefficients_local_tune_succ_prob(/*io*/ vector<double> &pr);

  /**
     @brief tune the pruning parameter to make the curve of the pruning
     parameter more smooth.

     Optimization process to the pruning parameters to make the curve
     more smooth if the input has obvious discountinuties on consecutive
     pruning parameters.
  */
  void optimize_coefficients_local_tune_smooth(/*io*/ vector<double> &pr);

  /**
     @brief auxiliary function in optimizing the single enumeraiton time
     fixing the succ. prob. The purpose is to increase the probability
     to be close to the target (as much as it can). This is used
     to make sure the probability is ''somewhat'' close to the target.

     Auxiliary function to optimize the single enumeration time with the
     constraint that the succ. prob (or expected solutions) is fixed. It
     is used if the given targeted probaility is larger than the one
     computed using pruning parameters. Then one increases the succ.
     probability by increasing the pruning parameters. Note since the
     radius is fixed as input, it may not be possible to achieve the
     given targeted probability even if all coefficients are (1, ..., 1).
     In such case, it will try to increase the succ. prob as much as
     it can and possibly return (1, ..., 1).
  */
  void optimize_coefficients_incr_prob(/*io*/ vector<double> &pr);

  /**
     @brief auxiliary function in optimizing the single enumeraiton time
     fixing the succ. prob. The purpose is to decrease the probability
     to be close to the target. This is used to make sure the
     probability is ''somewhat'' close to the target.

     Auxiliary function to optimize the single enumeration time with the
     constraint that the succ. prob (or expected solutions) is fixed. It
     is used if the given targeted probaility is smaller than the one
     computed using pruning parameters. Then one decreases the succ.
     probability by descreasing the pruning parameters.
  */
  void optimize_coefficients_decr_prob(/*io*/ vector<double> &pr);

  /**
     @brief auxiliary function in optimizing the single enumeraiton
     time fixing succ. prob. This is used to make sure the probability
     is ''sufficiently'' close to the target (if possible).

     Heuristic tuning procedure which seems to be useful. This is used to
     make the ratio between the succ. prob and the target succ. prob are
     sufficiently close. Depending on whether the succ. prob is larger
     (or smaller), it will try to reduce the pruning coefficients
     (or increase) to make succ. prob \approx the target succ. prob.
     Sometimes the succ. prob can not be modified closer to target
     succ. prob due to the contraints (in such case it just returns).
  */
  void optimize_coefficients_tune_prob(/*io*/ vector<double> &pr);
};

template <class FT>
bool Pruner<FT>::tabulated_values_imported = false;  ///< (static) tabulated value not loaded yet
template <class FT> FT Pruner<FT>::tabulated_factorial[PRUNER_MAX_N];
template <class FT> FT Pruner<FT>::tabulated_ball_vol[PRUNER_MAX_N];

/**
    enforce (and possible fix) for half/full coefficients
 */
template <class FT> inline bool Pruner<FT>::enforce(/*io*/ vec &b, /*opt i*/ const int j)
{
  int dn      = b.size();
  int c       = (dn == d) ? 1 : 2;
  bool status = false;

  // the last coefficient being 1
  if ((b[dn - 1] < .999) & (j != dn - 1))
  {
    status    = 1;
    b[dn - 1] = 1.;
  }

  for (int i = 0; i < dn; ++i)
  {
    status |= (b[i] > 1.0001);
    b[i] = b[i] > 1 ? 1. : b[i];

    // note min_pruning_coefficients always has length n
    if (b[i] <= min_pruning_coefficients[i / c])
      b[i] = min_pruning_coefficients[i / c];
  }

  for (int i = j; i < dn - 1; ++i)
  {
    if (b[i + 1] < b[i])
    {
      status |= (b[i + 1] + .000001 < b[i]);
      b[i + 1] = b[i];
    }
  }

  for (int i = j - 1; i >= 0; --i)
  {
    if (b[i + 1] < b[i])
    {
      status |= (b[i + 1] + .000001 < b[i]);
      b[i] = b[i + 1];
    }
  }
  return status;
}

// put here to avoid warnings on inlines function should be in the header.
#include "pruner_simplex.h"

FPLLL_END_NAMESPACE

#endif /* FPLLL_PRUNER_H */
