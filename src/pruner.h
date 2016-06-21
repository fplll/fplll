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
#include "fplll.h"

FPLLL_BEGIN_NAMESPACE


// naming conventions:
// b is for bound (squared)
// pv is for partial volumes (NOT squared)
// r is for gram schmidt length (squared). Must be renormaliuzed by the
// enumeration radius
// p is for polynomial

// inside this code, b,pv,and R are in reversed order
// as to conform with the algorithm desciption of [ChenThesis]
// reversing output and ouput is done by the C extern function

// n is for the dimension of the basis to prune
// d is for degrees of polynomials. md is for max_degree
// d is floor(n/2 at most). Odd n are dealt with by ignoring the first component



#define PRUNER_MAX_PREC 1000
#define PRUNER_MAX_D 1023
#define PRUNER_MAX_N 2047


template<class FT> 
class Pruner{
  public:
    
    Pruner();

    template<class ZT2,class FT2>
    void load_basis_shape(MatGSO<ZT2, FT2>& gso, int beginning = 0, int end = 0);
    void set_parameters(FT preproc_cost, FT target_sucess_proba);
    
    void optimize_pruning_coeffs(/*io*/double* pr, /*i*/ int reset = 1);
    void get_cost(/*i*/double* pr);
    FT get_svp_success_proba(/*i*/double* pr);

  private:
    using vec = array<FT, PRUNER_MAX_N>;
    using evec = array<FT, PRUNER_MAX_D>; 
          // Even vectors, i.e. only one every two entry is stored: V[2i] = V[2i+1] =E[i]
    using poly = array<FT, PRUNER_MAX_D + 1>;


    void set_tabulated_consts();
    int n;  // Dimension of the (sub)-basis
    int d;  // Degree d = floor(n/2)

    vec r;
    vec pv;
    FT preproc_cost;
    FT target_sucess_proba;

    inline void enforce_constraints(/*io*/ evec &b, /*opt i*/ int j = 0);
    inline FT eval_poly(int ld,/*i*/ poly *p, FT x);
    inline void integrate_poly(int ld,/*io*/ poly *p);
    inline FT relative_volume(/*i*/int rd, evec &b);
    inline FT node_count_predict(/*i*/ evec &b);
    inline FT svp_success_proba(/*i*/ evec &b);
    inline FT cost_factor(/*i*/ evec &b);
    FT cost_factor_derivative(/*i*/ evec &b, /*o*/ evec &res);
    int improve(/*io*/ evec &b);
    int descent(/*io*/ evec &b);

    FT tabulated_factorial[PRUNER_MAX_N];
    FT tabulated_ball_vol[PRUNER_MAX_N];
    
    FT epsilon;    // Epsilon to use for numerical differentiation
    FT min_step;    // Minimal step in a given direction 
    FT step_factor; // Increment factor for steps in a given direction
    FT shell_ratio; // Shell thickness Ratio when evaluating svp probability
};

template<class FT>
Pruner<FT>::Pruner(){
  set_tabulated_consts();
  epsilon = pow(2., -13);
  min_step = pow(2., -9);
  step_factor = pow(2, .5);
  shell_ratio = .995;
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
template<class ZT2,class FT2>
void Pruner<FT>::load_basis_shape(MatGSO<ZT2, FT2>& gso, int beginning, int end){
  if (!end){
    end = gso.d;
  }
  n = end - beginning;
  d = n/2;
  FT2 f;
  for (int i = 0; i < n; ++i)
  {
    gso.getR(f, beginning + i, beginning + i);
    r[i] = f;
    std::cerr << "HAHAAHA" << i << " BEBEBBE " << r[i] << endl;
  }
}

template<class FT>
FT Pruner<FT>::get_svp_success_proba(/*i*/double* pr){
  vec b;
  for (int i = 0; i < d; ++i) {
    b[i] = pr[n - 1 - 2 * i];
  }
  return svp_success_proba(b);

}




// PRIVATE METHODS











template<class FT>
inline void Pruner<FT>::enforce_constraints(/*io*/ evec &b, /*opt i*/ int j){
  b[d - 1] = 1;
  for (int i = 0; i < d; ++i){
    if (b[i] > 1) b[i] = 1;
    if (b[i] <= .1) b[i] = .1;
  }
  for (int i = j; i < d - 1; ++i){
    if (b[i + 1] < b[i]) b[i + 1] = b[i];
  }
  for (int i = j - 1; i >= 0; --i){
    if (b[i + 1] < b[i]) b[i] = b[i + 1];
  }  
}

template<class FT>
inline FT Pruner<FT>::eval_poly(int ld,/*i*/ poly *p, FT x){
  FT acc = 0.0;
  for (int i = ld; i >= 0; --i) {
    acc *= x;
    acc += (*p)[i];
  }
  return acc;
}

template<class FT>
inline void Pruner<FT>::integrate_poly(int ld,/*io*/ poly *p){
  for (int i = ld; i >= 0; --i) {
    (*p)[i + 1] = (*p)[i] / (i + 1.);
  }
  (*p)[0] = 0;
}



template<class FT>
inline FT Pruner<FT>::relative_volume(int rd, /*i*/ evec &b){
  poly P;
  P[0] = 1;
  int ld = 0;
  for (int i = rd - 1; i >= 0; --i) {
    integrate_poly(ld, &P);
    ld++;
    P[0] = -eval_poly(ld, &P, b[i] / b[rd - 1]);
  }
  if (rd % 2) {
    return -P[0] * tabulated_factorial[rd];
  } else {
    return P[0] * tabulated_factorial[rd];
  }
}

template<class FT>
inline FT Pruner<FT>::node_count_predict(/*i*/ evec &b){
  vec rv; // Relative volumes at each level

  for (int i = 0; i < d; ++i) {
    rv[2 * i + 1] = relative_volume(i + 1, b);
  }

  rv[0] = 1;
  for (int i = 1; i < d; ++i) {
    rv[2 * i] =
        sqrt(rv[2 * i - 1] * rv[2 * i + 1]); // Interpolate even values
  }

  FT total = 0;
  for (int i = 0; i < 2 * d; ++i) {
    FT tmp;
    tmp = rv[i] * tabulated_ball_vol[i + 1] *
             pow(b[i / 2], (1. + i) / 2.) / pv[i];
    total += tmp;
  }
  return total;
}


template<class FT>
inline FT Pruner<FT>::svp_success_proba(/*i*/ evec &b){

  evec b_minus_db;
  FT dx = shell_ratio;

  for (int i = 0; i < d; ++i) {
    b_minus_db[i] = b[i] / (dx * dx);
    if (b_minus_db[i] > 1)
      b_minus_db[i] = 1;
  }

  FT vol = relative_volume(d, b);
  FT dxn = pow(dx, 2. * d);
  FT dvol = vol - dxn * relative_volume(d, b_minus_db);
  return dvol / (1 - dxn);

}



// template<class FT>
// inline FT Pruner<FT>::cost_factor(/*i*/ evec &b);
// template<class FT>
// FT Pruner<FT>::cost_factor_derivative(/*i*/ evec &b, /*o*/ evec &res);
// template<class FT>
// int Pruner<FT>::improve(/*io*/ evec &b);
// template<class FT>
// int Pruner<FT>::descent(/*io*/ evec &b);





FPLLL_END_NAMESPACE

