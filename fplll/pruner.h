#ifndef FPLLL_PRUNER_H
#define FPLLL_PRUNER_H

/* Copyright (C) 2015-2016 Martin Albrecht, Leo Ducas.

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

#include "bkz_param.h"
#include "defs.h"
#include "gso.h"
#include <array>

FPLLL_BEGIN_NAMESPACE

/**
   This file provides an implementation of the numerically optimized pruning
   strategy of [GNR10].

   Many details for implementation follows from the thesis [Chen13] Some
   simplifications have been made, for example we restrict ourselves to ``even
   vector bounds'' i.e., the bound at indices 2i and 2i+1 are kept equal, as to
   allow volume estimation by pure integration (see [GNR10,Chen13])

   The current descent method use gradients, but alternative algorithms
   (Nelder-Mead) May be added soon.

   naming conventions:

   - b is for bound (squared)

   - ipv is for inverse partial volumes (NOT squared)

   - r is for gram schmidt length (squared). Automatically renormalized to avoid
   overflowing partial volumes

   - p is for polynomial

   inside this code, b,pv,and R are in reversed order as to conform with the
   algorithm desciption of [Chen13] reversing the order is handled by the load
   and save methods.

   - n is for the dimension of the basis to prune

   - d is for degrees of polynomials. md is for max_degree

   - d is floor(n/2 at most). Odd n are dealt with by ignoring the first component
*/

#define PRUNER_MAX_PREC 1000
#define PRUNER_MAX_D 1023
#define PRUNER_MAX_N 2047

#define PRUNER_MIN_BOUND .1

/**
   @brief prune function, hiding the Pruner class

   @param pruning store output here
   @param preproc_cost cost of preprocessing
   @param target overall target success probability/expected solutionss
   @param m GSO matrix
   @param method for the descent (gradient, NM, both)
   @param metric either success probability or expected number of solutions
   @param reset reset pruning coefficients

*/

template <class FT>
void prune(/*output*/ Pruning &pruning,
           /*inputs*/ double &enumeration_radius, const double preproc_cost, const double target,
           vector<double> &r, const PrunerMethod method = PRUNER_METHOD_HYBRID,
           const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST, bool reset = true);

/**
   @brief prune function averaging over several bases

   @param pruning store output here
   @param preproc_cost cost of preprocessing
   @param target overall target success probability/expected solutionss
   @param m GSO matrices
   @param method for the descent (gradient, NM, both)
   @param metric either success probability or expected number of solutions
   @param reset reset pruning coefficients

*/

template <class FT>
void prune(/*output*/ Pruning &pruning,
           /*inputs*/ double &enumeration_radius, const double preproc_cost, const double target,
           vector<vector<double>> &rs, const PrunerMethod method = PRUNER_METHOD_HYBRID,
           const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST, bool reset = true);

/**
   @brief svp_probability function, hiding the Pruner class

   @param pruning
   @return pruning probability
*/

template <class FT> FT svp_probability(const Pruning &pruning);
template <class FT> FT svp_probability(const vector<double> &pr);

template <class FT> class Pruner
{
public:
  class TestPruner;
  friend class TestPruner;

  /** @brief enumeration radius (squared) */
  FT enumeration_radius;

  /** @brief cost of pre-processing a basis for a retrial

      This cost should be expressed in terms of ``nodes'' in an enumeration.
      Roughly, a node is equivalent to 100 CPU cycles.
  */
  FT preproc_cost;

  /** @brief desired success probability after several retrial

      @note one can try to force probability >= target by setting
      a prohibitive preproc_cost. But beware: this may induces numerical
      stability issue, especially with the gradient method.
  */

  FT target;

  int verbosity = 0;

  PrunerMethod method;
  PrunerMetric metric;

  Pruner(FT enumeration_radius = 0.0, FT preproc_cost = 0.0, FT target = 0.9,
         PrunerMethod method = PRUNER_METHOD_HYBRID,
         PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST, size_t n = 0, size_t d = 0)
      : enumeration_radius(enumeration_radius), preproc_cost(preproc_cost), target(target),
        method(method), metric(metric), n(n), d(d)
  {

    if (!tabulated_value_imported)
    {
      set_tabulated_consts();
      tabulated_value_imported = true;
    }
    epsilon     = std::pow(2., -7);  // Guesswork. Will become obsolete with Nelder-Mead
    min_step    = std::pow(2., -6);  // Guesswork. Will become obsolete with Nelder-Mead
    step_factor = std::pow(2, .5);   // Guesswork. Will become obsolete with Nelder-Mead
    shell_ratio = .995;  // This approximation means that SVP will in fact be approx-SVP with factor
                         // 1/.995. Sounds fair.
    min_cf_decrease = .995;  // We really want the gradient descent to reach the minima
    symmetry_factor = 2;     // For now, we are just considering SVP
  }

  /** @brief load the shape of a basis from vector<double>. Mostly for testing purposes */
  void load_basis_shape(const vector<double> &gso_sq_norms, bool reset_renorm = true);

  /** @brief load the shapes of may bases from vector<vector<double>> . Cost are average over all
   * bases. Mostly for testing purposes */

  void load_basis_shapes(const vector<vector<double>> &gso_sq_norms_vec);

  /** @brief optimize pruning coefficients

      @note Basis Shape and other parameters must have been set beforehand. See
      auto_prune for an example of proper usage.
  */
  void optimize_coefficients(/*io*/ vector<double> &pr, /*i*/ const bool reset = true);

  /** @brief Compute the cost of a single enumeration */

  double single_enum_cost(/*i*/ const vector<double> &pr, vector<double> *detailed_cost = nullptr)
  {
    evec b;
    load_coefficients(b, pr);
    return single_enum_cost(b, detailed_cost).get_d();
  }

  /** @brief Compute the cost of r enumeration and (r-1) preprocessing, where r
      is the required number of retrials to reach target
  */
  double repeated_enum_cost(/*i*/ const vector<double> &pr)
  {
    evec b;
    load_coefficients(b, pr);
    return repeated_enum_cost(b).get_d();
  }

  /**
     @brief Compute the success proba of a single enumeration
  */
  double measure_metric(/*i*/ const vector<double> &pr)
  {
    if (!n)
    {  // Can be called even if no basis has been loaded. In that case, set the dims
      n = pr.size();
      d = n / 2;
    }
    evec b;
    load_coefficients(b, pr);
    return measure_metric(b).get_d();
  }

private:
  using vec  = array<FT, PRUNER_MAX_N>;
  using evec = array<FT, PRUNER_MAX_D>;
  // Even vectors, i.e. only one every two entry is stored: V[2i] = V[2i+1] =E[i]
  using poly = array<FT, PRUNER_MAX_D + 1>;

  // Load the constants for factorial and ball-volumes
  void set_tabulated_consts();
  size_t n;  // Dimension of the (sub)-basis
  size_t d;  // Degree d = floor(n/2)

  vec r;                      // Gram-Schmidt length (squared, inverted ordering)
  vec ipv;                    // Partial volumes (inverted ordering)
  FT renormalization_factor;  // internal renormalization factor to avoid over/underflows

  // Sanity check: has a basis indeed been loaded ?
  bool check_basis_loaded();
  // Initialize pruning coefficients (linear pruning)
  void init_coefficients(evec &b);
  // Load pruning coefficient from double*
  void load_coefficients(/*o*/ evec &b, /*i*/ const vector<double> &pr);
  // Save pruning coefficients to double*
  void save_coefficients(/*o*/ vector<double> &pr, /*i*/ const evec &b);
  // Enforce reasonable contraints on pruning bounds (inside [0,1], increasing).
  // Keeps index j unchanged when possible
  bool enforce_bounds(/*io*/ evec &b, /*opt i*/ const int j = 0);
  // Evaluate a polynomial
  FT eval_poly(const int ld, /*i*/ const poly &p, const FT x);
  // Integrate a polynomial
  void integrate_poly(const int ld, /*io*/ poly &p);
  // Compute the relative volume of a cylinder interesection of dim rd, and bounds b[0:rd]
  FT relative_volume(/*i*/ const int rd, const evec &b);
  // Compute the cost of a single enumeration
  FT single_enum_cost(/*i*/ const evec &b, vector<double> *detailed_cost = nullptr);
  // Compute the success probability for SVP/CVP of a single enumeration
  FT svp_probability(/*i*/ const evec &b);
  // Compute the expected nmber of solution of a single of a single enumeration
  FT expected_solutions(/*i*/ const evec &b);
  // One of the above depending on metric
  FT measure_metric(/*i*/ const evec &b);
  // Compute the cost of r enumeration and (r-1) preprocessing,
  // where r is the required number of retrials to reach target/target_solution
  FT repeated_enum_cost(/*i*/ const evec &b);
  // Compute the gradient of the above function
  void repeated_enum_cost_gradient(/*i*/ const evec &b, /*o*/ evec &res);
  // Choose pruning parameter to have a near-constant width enum tree
  void greedy(evec &b);
  // Improve the pruning bounds a bit,  using one gradient step
  int improve(/*io*/ evec &b);
  // Improve the pruning bounds substantially, using Nelder-Mead method
  int nelder_mead(/*io*/ evec &b);
  // Run the whole escent to optimize pruning bounds
  void descent(/*io*/ evec &b);

  static FT tabulated_factorial[PRUNER_MAX_N];
  static FT tabulated_ball_vol[PRUNER_MAX_N];
  static bool tabulated_value_imported;

  FT epsilon;          //< Epsilon to use for numerical differentiation
  FT min_step;         //< Minimal step in a given direction
  FT min_cf_decrease;  //< Maximal ratio of two consectuive cf in the descent before stopping

  FT step_factor;      //< Increment factor for steps in a given direction
  FT shell_ratio;      //< Shell thickness Ratio when evaluating svp proba
  FT symmetry_factor;  //< Set at 2 for SVP enumeration assuming the implem only explore half the
                       //< space
};

template <class FT> bool Pruner<FT>::tabulated_value_imported = false;
template <class FT> FT Pruner<FT>::tabulated_factorial[PRUNER_MAX_N];
template <class FT> FT Pruner<FT>::tabulated_ball_vol[PRUNER_MAX_N];

FPLLL_END_NAMESPACE

#endif /* FPLLL_PRUNER_H */
