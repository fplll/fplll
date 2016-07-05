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


#include <array>
#include "factorial.const"
#include "ballvol.const"
#include "defs.h"
#include "gso.h"

FPLLL_BEGIN_NAMESPACE

// This file provides an implementation of the numerically optimized
// pruning strategy of [GNR10]. 

// Many details for implementation follows from the thesis [Chen13]
// Some simplifications have been made, for example we restrict ourselves
// to ``even vector bounds'' i.e., the bound at indices 2i and 2i+1 are
// kept equal, as to allow volume estimation by pure integration 
// (see [GNR10,Chen13])

// The current descent method use gradients, but alternative algorithms (Nelder-Mead)
// May be added soon.

// naming conventions:
// b is for bound (squared)
// pv is for partial volumes (NOT squared)
// r is for gram schmidt length (squared). Automatically renormalized
// to avoid overflowing partial volumes
// p is for polynomial

// inside this code, b,pv,and R are in reversed order
// as to conform with the algorithm desciption of [Chen13]
// reversing the order is handled by the load and save methods.


// n is for the dimension of the basis to prune
// d is for degrees of polynomials. md is for max_degree
// d is floor(n/2 at most). Odd n are dealt with by ignoring the first component



#define PRUNER_MAX_PREC 1000
#define PRUNER_MAX_D 1023
#define PRUNER_MAX_N 2047



template<class FT> 
class Pruner{
  public:
    class TestPruner;
    friend class TestPruner;

    // Defines the cost of re-processing a basis for a retrial
    // This cost should be expressed in terms of ``Nodes'' in an enumeration
    // Roughly, a Node is equivalent to 100 CPU cycles
    FT preproc_cost;
    // Defines the desired success probability after several retrial
    FT target_success_proba;
    // Defines the enumeration radius (squared)
    FT enumeration_radius;

    // Note: one can try to force success_proba = target_success_proba by
    // setting a prohibitive preproc_cost. But beware: this may induces 
    // numerical stability issue, especially with the gradient method.
    // Melder-Mead should be more robust.

    
    Pruner();

    // Load the shape of a basis from a MatGSO object. Can select a projected sub-lattice [beginning,end-1]
    template<class GSO_ZT,class GSO_FT>
    void load_basis_shape(const MatGSO<GSO_ZT, GSO_FT>& gso, const int beginning = 0, const int end = 0);
    // Load the shape of a basis from a double*. Mostly for testing purposes ?
    void load_basis_shape(const int dim, const double* gso_sq_norms);
    


    
    // Optimize pruning coefficients. 
    // Basis Shape and other parameters must have been set beforehands. 
    // See auto_prune for an example of proper usage.
    void optimize_pruning_coeffs(/*io*/double* pr, /*i*/const int reset = 1);
    // Compute the cost of a single enumeration
    double get_single_enum_cost(/*i*/const double* pr);
    // Compute the cost of r enumeration and (r-1) preprocessing, 
    // where r is the required number of retrials to reach target_success_proba
    double get_repeated_enum_cost(/*i*/const double* pr);
    // Compute the success proba of a single enumeration
    double get_svp_success_proba(/*i*/const double* pr);



  private:
    using vec = array<FT, PRUNER_MAX_N>;
    using evec = array<FT, PRUNER_MAX_D>; 
          // Even vectors, i.e. only one every two entry is stored: V[2i] = V[2i+1] =E[i]
    using poly = array<FT, PRUNER_MAX_D + 1>;

    // Load the constants for factorial and ball-volumes
    void set_tabulated_consts();
    int n;  // Dimension of the (sub)-basis
    int d;  // Degree d = floor(n/2)

    vec r;  // Gram-Schmidt length (squared, inverted ordering)
    vec pv; // Partial volumes (inverted ordering)
    FT renormalization_factor; // internal renormalization factor to avoid over/underflows


    // Sanity check: has a basis indeed been loaded ?
    int check_loaded_basis();
    // Initialize pruning coefficients (linear pruning)
    void init_prunning_coeffs(evec &b);
    // Load pruning coefficient from double*
    void load_prunning_coeffs(/*i*/const double* pr, /*o*/ evec& b);
    // Save pruning coefficients to double*
    void save_prunning_coeffs(/*o*/double* pr, /*i*/const evec& b);
    // Enforce reasonable contraints on pruning bounds (inside [0,1], increasing). 
    // Keeps index j unchanged when possible
    inline int enforce_bounds(/*io*/ evec &b, /*opt i*/const int j = 0);
    // Evaluate a polynomial
    inline FT eval_poly(const int ld,/*i*/const poly *p,const FT x);
    // Integrate a polynomial
    inline void integrate_poly(const int ld,/*io*/ poly *p);
    // Compute the relative volume of a cylinder interesection of dim rd, and bounds b[0:rd]
    inline FT relative_volume(/*i*/const int rd, const evec &b);
    // Compute the cost of a single enumeration
    inline FT single_enum_cost(/*i*/const evec &b);
    // Compute the success probability of a single enumeration
    inline FT svp_success_proba(/*i*/const evec &b);
    // Compute the cost of r enumeration and (r-1) preprocessing, 
    // where r is the required number of retrials to reach target_success_proba
    inline FT repeated_enum_cost(/*i*/const evec &b);
    // Compute the gradient of the above function
    void repeated_enum_cost_gradient(/*i*/const evec &b, /*o*/ evec &res);
    // Improve the pruning bounds a bit,  using one gradient step
    int improve(/*io*/ evec &b);
    // Run the whole escent to optimize pruning bounds
    void descent(/*io*/ evec &b);

    FT tabulated_factorial[PRUNER_MAX_N];
    FT tabulated_ball_vol[PRUNER_MAX_N];
    

    FT epsilon;    // Epsilon to use for numerical differentiation
    FT min_step;    // Minimal step in a given direction 
    FT min_cf_decrease;    // Maximal ratio of two consectuive cf in the descent before stopping

    FT step_factor; // Increment factor for steps in a given direction
    FT shell_ratio; // Shell thickness Ratio when evaluating svp proba
    FT symmetry_factor; // Set at 2 for SVP enumeration assuming the implem only explore half the space
};

template<class FT>
Pruner<FT>::Pruner(){
  n = 0;
  d = 0;
  set_tabulated_consts();

  epsilon = std::pow(2., -13);    // Guesswork. Will become obsolete with Nelder-Mead
  min_step = std::pow(2., -12);   // Guesswork. Will become obsolete with Nelder-Mead  
  step_factor = std::pow(2, .5);  // Guesswork. Will become obsolete with Nelder-Mead 
  shell_ratio = .995;             // This approximation means that SVP will in fact be approx-SVP with factor 1/.995. Sounds fair.
  min_cf_decrease = .9999;        // We really want the gradient descent to reach the minima
  symmetry_factor = 2;            // For now, we are just considering SVP

  preproc_cost = 0.0;              // Please set your own value before running.
  enumeration_radius = 0.0;        // Please set your own value before running.
  target_success_proba = .90;     // Please set your own value before running.
  preproc_cost = 0.0;               // Please set your own value before running.
}


template<class FT>
void Pruner<FT>::set_tabulated_consts(){
  mpfr_t tmp;
  mpfr_init2(tmp, PRUNER_MAX_PREC);
  for (int i = 0; i < PRUNER_MAX_N; ++i)
  {
    mpfr_set_str(tmp, pre_factorial[i], 10, MPFR_RNDN);
    tabulated_factorial[i] = tmp;
    mpfr_set_str(tmp, pre_ball_vol[i], 10, MPFR_RNDN);
    tabulated_ball_vol[i] = tmp;
  }
  return;
}



/// PUBLIC METHODS

template<class FT>
template<class GSO_ZT,class GSO_FT>
void Pruner<FT>::load_basis_shape(const MatGSO<GSO_ZT, GSO_FT>& gso,const int beginning,const int end){
  if (!end){
    end = gso.d;
  }
  n = end - beginning;
  d = n/2;
  if (!d){
    throw std::runtime_error("Inside Pruner : Needs a dimension n>1");
  }
  GSO_FT f;
  FT logvol,tmp;
  logvol = 0.0;
  for (int i = 0; i < n; ++i)
  {
    gso.getR(f, end - 1 - i, end - 1 - i);
    r[i] = f;
    logvol += log(f);
  }
  tmp = - n;
  renormalization_factor = exp(logvol / tmp);
  for (int i = 0; i < n; ++i)
  {
    r[i] *= renormalization_factor;
  }

  tmp = 1.;
  for (int i = 0; i < 2 * d; ++i) {
    tmp *= sqrt(r[i]);
    pv[i] = tmp;
  }
}

template<class FT>
void Pruner<FT>::load_basis_shape(const int dim,const double* gso_sq_norms){
  n = dim;
  d = n/2;
  if (!d){
    throw std::runtime_error("Inside Pruner : Needs a dimension n>1");
  }
  FT logvol,tmp;
  logvol = 0.0;
  for (int i = 0; i < n; ++i)
  {
    r[i] = gso_sq_norms[n - 1 - i];

    logvol += log(r[i]);
  }
  tmp = - n;
  renormalization_factor = exp(logvol / tmp);

  for (int i = 0; i < n; ++i)
  {
    r[i] *= renormalization_factor;

  }
  tmp = 1.;
  for (int i = 0; i < 2 * d; ++i) {
    tmp *= sqrt(r[i]);
    pv[i] = tmp;
  }
}


template<class FT>
double Pruner<FT>::get_svp_success_proba(/*i*/const double* pr){
  evec b;
  load_prunning_coeffs(pr, b);
  return svp_success_proba(b).get_d();
}

template<class FT>
double Pruner<FT>::get_single_enum_cost(/*i*/const double* pr){
  evec b;
  load_prunning_coeffs(pr, b);
  return single_enum_cost(b).get_d();
}

template<class FT>
double Pruner<FT>::get_repeated_enum_cost(/*i*/const double* pr){
  evec b;
  load_prunning_coeffs(pr, b);
  return repeated_enum_cost(b).get_d();
}


template<class FT>
void Pruner<FT>::optimize_pruning_coeffs(/*io*/double* pr, /*i*/const int reset){
  evec b;
  if (reset){
    init_prunning_coeffs(b);
  }
  else{
    load_prunning_coeffs(pr, b);
  }
  descent(b);
  save_prunning_coeffs(pr, b);
}

// PRIVATE METHODS

template<class FT>
void Pruner<FT>::load_prunning_coeffs(/*i*/const double* pr, /*o*/ evec& b){
  for (int i = 0; i < d; ++i) {
    b[i] = pr[n - 1 - 2 * i];
  }
  if (enforce_bounds(b)){
    throw std::runtime_error(
      "Inside Pruner : Ill formed pruning coefficients (must be decreasing, starting with two 1.0)");
  }
}

template<class FT>
int Pruner<FT>::check_loaded_basis(){
  if (d){
      return 0;
    }
  throw std::runtime_error("Inside Pruner : No basis loaded");
  return 1;
}

template<class FT>
void Pruner<FT>::save_prunning_coeffs(/*o*/double* pr, /*i*/const evec& b){
  for (int i = 0; i < d; ++i) {
    pr[n - 1 - 2 * i] = b[i].get_d();
    pr[n - 2 - 2 * i] = b[i].get_d();
  }
  pr[0] = 1.;
}

template<class FT>
inline int Pruner<FT>::enforce_bounds(/*io*/ evec &b, /*opt i*/const int j){
  int status = 0;
  if (b[d - 1] < 1){
    status = 1;
  }
  b[d - 1] = 1;
  for (int i = 0; i < d; ++i){
    if (b[i] > 1) {b[i] = 1.0; status = 1;}
    if (b[i] <= .1) b[i] = .1;
  }
  for (int i = j; i < d - 1; ++i){
    if (b[i + 1] < b[i]) {b[i + 1] = b[i]; status = 1;}
  }
  for (int i = j - 1; i >= 0; --i){
    if (b[i + 1] < b[i]) {b[i] = b[i + 1]; status = 1;}
  }  
  return status;
}

template<class FT>
inline FT Pruner<FT>::eval_poly(const int ld,/*i*/const poly *p,const FT x){
  FT acc;
  acc = 0.0;
  for (int i = ld; i >= 0; --i) {
    acc = acc * x;
    acc = acc + (*p)[i];
  }
  return acc;
}

template<class FT>
inline void Pruner<FT>::integrate_poly(const int ld,/*io*/ poly *p){
  for (int i = ld; i >= 0; --i) {
    FT tmp;
    tmp = i + 1.;
    (*p)[i + 1] = (*p)[i] / tmp;
  }
  (*p)[0] = 0.0;
}



template<class FT>
inline FT Pruner<FT>::relative_volume(const int rd, /*i*/const evec &b){
  poly P;
  P[0] = 1;
  int ld = 0;
  for (int i = rd - 1; i >= 0; --i) {
    integrate_poly(ld, &P);
    ld++;
    P[0] = -1.0 * eval_poly(ld, &P, b[i] / b[rd - 1]);
  }
  if (rd % 2) {
    return -1.0 * P[0] * tabulated_factorial[rd];
  } else {
    return P[0] * tabulated_factorial[rd];
  }
}

template<class FT>
inline FT Pruner<FT>::single_enum_cost(/*i*/const evec &b){
  vec rv; // Relative volumes at each level

  for (int i = 0; i < d; ++i) {

    rv[2 * i + 1] = relative_volume(i + 1, b);
  }

  rv[0] = 1;
  for (int i = 1; i < d; ++i) {
    rv[2 * i] =
        sqrt(rv[2 * i - 1] * rv[2 * i + 1]); // Interpolate even values
  }

  FT total;
  total = 0.0;
  FT normalized_radius;
  normalized_radius = sqrt(enumeration_radius * renormalization_factor);
  
  for (int i = 0; i < 2 * d; ++i) {
    FT tmp;
    tmp = pow_si(normalized_radius, 1 + i) *
          rv[i] * tabulated_ball_vol[i + 1] *
          sqrt(pow_si(b[i / 2], 1 + i)) / pv[i];
    total += tmp;
  }
  total /= symmetry_factor;  
  //exit(1);
  return total;
}


template<class FT>
inline FT Pruner<FT>::svp_success_proba(/*i*/const evec &b){

  evec b_minus_db;
  FT dx = shell_ratio;

  for (int i = 0; i < d; ++i) {
    b_minus_db[i] = b[i] / (dx * dx);
    if (b_minus_db[i] > 1)
      b_minus_db[i] = 1;
  }

  FT vol = relative_volume(d, b);
  FT dxn = pow_si(dx, 2 * d);
  FT dvol =  dxn * relative_volume(d, b_minus_db) - vol;
  return dvol / (dxn - 1.);

}



template<class FT>
inline FT Pruner<FT>::repeated_enum_cost(/*i*/const evec &b){

  FT success_proba = svp_success_proba(b);

  if (success_proba >= target_success_proba)
    return single_enum_cost(b);

  FT trials =  log(1.0 - target_success_proba) / log(1.0 - success_proba);
  return single_enum_cost(b) * trials + preproc_cost * (trials-1);
}



template<class FT>
void Pruner<FT>::repeated_enum_cost_gradient(/*i*/const evec &b, /*o*/ evec &res){
  evec bpDb;
  res[d - 1] = 0.0;
  for (int i = 0; i < d-1; ++i) {
    bpDb = b;
    bpDb[i] *= (1.0 - epsilon);
    enforce_bounds(bpDb, i);
    FT X = repeated_enum_cost(bpDb);

    bpDb = b;
    bpDb[i] *= (1.0 + epsilon);
    enforce_bounds(bpDb, i);
    FT Y = repeated_enum_cost(bpDb);
    res[i] = (log(X) - log(Y)) / epsilon;
  }
}



template<class FT>
int Pruner<FT>::improve(/*io*/ evec &b){

  FT cf = repeated_enum_cost(b);
  FT old_cf = cf;
  evec newb;
  evec gradient;
  repeated_enum_cost_gradient(b, gradient);
  FT norm = 0;

  // normalize the gradient
  for (int i = 0; i < d; ++i) {
    norm += gradient[i] * gradient[i];
    newb[i] = b[i];
  }

  norm = sqrt(norm /  (1.0 * (1. * d)) );
  if (norm <= 0.)
    return 0;

  for (int i = 0; i < d; ++i) {
    gradient[i] /= norm;
  }
  FT new_cf;

  FT step = min_step;
  int i;

  for (i = 0;; ++i) {
    for (int i = 0; i < d; ++i) {
      newb[i] = newb[i] + step * gradient[i];
    }

    enforce_bounds(newb);
    new_cf = repeated_enum_cost_factor(newb);



    if (new_cf >= cf){
      break;
    }
    b = newb;
    cf = new_cf;
    step *= step_factor;
  }

  if (cf > old_cf * min_cf_decrease){
    return 0;
  }
  return i;
}

template<class FT>
void Pruner<FT>::descent(/*io*/ evec &b){
  while (improve(b)) { };
}


template<class FT>
void Pruner<FT>::init_prunning_coeffs(evec &b) {
  for (int i = 0; i < d; ++i) {
    b[i] = .1 + ((1.*i) / d);
  }
  enforce_bounds(b);
}





/// Autoprune function, hiding the Pruner class


template<class FT, class GSO_ZT, class GSO_FT> 
void auto_prune(/*output*/ double* pr, double& success_proba,
                /*inputs*/const double enumeration_radius,const double preproc_cost,
                const double target_success_proba,
                const MatGSO<GSO_ZT, GSO_FT>& gso, int beginning = 0, int end = 0){

  Pruner<FP_NR<double>> pru;

  pru.enumeration_radius = enumeration_radius;
  pru.target_success_proba = target_success_proba;
  pru.preproc_cost = preproc_cost;
  load_basis_shape(gso,beginning, end);
  pru.optimize_pruning_coeffs(pr);
  success_proba = pru.get_svp_success_proba(pr);
}


FPLLL_END_NAMESPACE

#endif /* FPLLL_PRUNER_H */
