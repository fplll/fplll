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

#define PRUNER_MAX_N 2047

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
           /*inputs*/ 
const double enumeration_radius, const double preproc_cost, 
vector<double> &gso_r, const double target = .9,
const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST, 
const int flags=PRUNER_DEFAULT, const double timeout = -1.);

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
           /*inputs*/ double enumeration_radius, const double preproc_cost,
vector<vector<double>> &gso_rs, const double target = .9, 
const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST, 
const int flags=PRUNER_DEFAULT, const double timeout = -1.);

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

  // TODO : shall those be public ?
  FT enumeration_radius;
  FT preproc_cost;
  FT target;
  PrunerMetric metric;
  int shapes_loaded = 0;
  int flags;
  size_t n;  // Dimension of the (sub)-basis
  size_t d;  // Degree d = floor(n/2)
  double timeout;
  double min_pruning_bound;

  void import_tabulated_values()
  {
    if (!tabulated_values_imported)
    {
      set_tabulated_consts();
      tabulated_values_imported = true;
    }
  }
  
  // TODO use channels to avoid explicit vebosity conditions
  // See https://stackoverflow.com/questions/11826554/standard-no-op-output-stream
  // ostream &channel1; // Channel for important messages
  // ostream &channel2; // Channel for less important message

  Pruner(const int n): n(n)
  {
    import_tabulated_values();
  }

  Pruner(const FT enumeration_radius, const FT preproc_cost,
   const vector<double> &gso_r, const FT target = 0.9,
   const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST, 
   int flags = PRUNER_DEFAULT, double timeout = -1): 
  enumeration_radius(enumeration_radius), 
  preproc_cost(preproc_cost), target(target),
  metric(metric),flags(flags), timeout(timeout)
  {
    import_tabulated_values();
    load_basis_shape(gso_r);
    set_min_pruning_bound();
    if (timeout<0) 
    {
      timeout = PRUNER_DEFAULT_TIMEOUT_CONST * n * n;
      flags |= PRUNER_TIMOUT_WARNING;
    }
    if (timeout==0)
    {
      // TODO : actually implement timeout
      timeout = 4.2e17; // The age of the universe in seconds
    }
  // TODO connect channels to cerr or   
  }


  Pruner(const FT enumeration_radius, const FT preproc_cost,
   const vector<vector<double>> &gso_rs, const FT target = 0.9,
   const PrunerMetric metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST, 
   int flags = PRUNER_DEFAULT, double timeout = -1): 
  enumeration_radius(enumeration_radius), 
  preproc_cost(preproc_cost), target(target),
  metric(metric),flags(flags), timeout(timeout)
  {
    import_tabulated_values();
    load_basis_shapes(gso_rs);
    set_min_pruning_bound();
    if (timeout<0) 
    {
      timeout = PRUNER_DEFAULT_TIMEOUT_CONST * n * n;
      flags |= PRUNER_TIMOUT_WARNING;
    }
    if (timeout==0)
    {
      // TODO : actually implement timeout
      timeout = 4.2e17; // The age of the universe in seconds
    }
  // TODO connect channels to cerr or   
  }


  /** @brief optimize pruning coefficients
      @note Basis Shape and other parameters must have been set beforehand. See
      auto_prune for an example of proper usage.
  */
  void optimize_coefficients(/*io*/ vector<double> &pr);

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
  double descent_starting_clock;
  static FT tabulated_factorial[PRUNER_MAX_N];
  static FT tabulated_ball_vol[PRUNER_MAX_N];
  static bool tabulated_values_imported;

  FT epsilon = std::pow(2., -7);    //< Epsilon to use for numerical differentiation
  FT min_step = std::pow(2., -6);   //< Minimal step in a given direction
  FT min_cf_decrease  = .995;       //< Maximal ratio of two consectuive cost_factor in the descent before stopping
  FT step_factor = std::pow(2, .5); //< Increment factor for steps in a given direction
  FT shell_ratio  = .995;           //< Shell thickness Ratio when evaluating svp proba
  FT symmetry_factor = 2;           //< 2 for SVP, 1 for CVP.

  // Following typedefs are given for readability
  using vec  = vector<FT>; // Those have dimension n
  using evec = vector<FT>; // Those have dimension d
  using poly = vector<FT>; // Those have dimension d+1 (or less)


  vec r;                      // Gram-Schmidt length (squared, inverted ordering)
  vec ipv;                    // Partial volumes (inverted ordering)
  FT normalization_factor;  // internal renormalization factor to avoid over/underflows
  FT normalized_radius;  // internal renormalization factor to avoid over/underflows
  int verbosity = 0;


  // Load the constants for factorial and ball-volumes
  void set_tabulated_consts();
  // Has the descent exceeeded the timeout
  bool timeouted();
  /** @brief load the shape of a basis from vector<double>.  */
  void load_basis_shape(const vector<double> &gso_r, bool reset_normalization = true);
  /** @brief load the shapes of many bases from vector<vector<double>>. 
      Costs are average over all bases.  */
  void load_basis_shapes(const vector<vector<double>> &gso_rs);
  // Set the min_puning_bound
  void set_min_pruning_bound();
  // Removed : now just use greedy as the starting point
  // Initialize pruning coefficients (linear pruning)
  // void init_coefficients(evec &b);

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
  void repeated_enum_cost_gradient(/*i*/ const evec &b, /*o*/ evec &gso_res);
  // Choose pruning parameter to have a near-constant width enum tree
  void greedy(evec &b);
  // Improve the pruning bounds a bit,  using one gradient step
  int gradient_step(/*io*/ evec &b);
  // Improve the pruning bounds substantially, using Nelder-Mead method
  int nelder_mead_step(/*io*/ evec &b);
  // Run the whole escent to optimize pruning bounds
  void descent(/*io*/ evec &b);
};

template <class FT> bool Pruner<FT>::tabulated_values_imported = false;
template <class FT> FT Pruner<FT>::tabulated_factorial[PRUNER_MAX_N];
template <class FT> FT Pruner<FT>::tabulated_ball_vol[PRUNER_MAX_N];

FPLLL_END_NAMESPACE

#endif /* FPLLL_PRUNER_H */
