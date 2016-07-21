/* Copyright (C) 2015-2016 Martin Albrecht, Leo Ducas.

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


#include "ballvol.const"
#include "factorial.const"
#include "pruner.h"
#include "gso.h"

FPLLL_BEGIN_NAMESPACE

template <class FT, class GSO_ZT, class GSO_FT>
void prune(/*output*/ vector<double> &pr, double &probability,
           /*inputs*/ const double enumeration_radius, const double preproc_cost,
           const double target_probability, const MatGSO<GSO_ZT, GSO_FT> &m,
           int start_row, int end_row)
{

  Pruner<FT> pruner(enumeration_radius, preproc_cost, target_probability);
  pruner.load_basis_shape(m, start_row, end_row);
  pruner.optimize_coefficients(pr);
  probability = pruner.svp_probability(pr);
}

template <class FT, class GSO_ZT, class GSO_FT>
void prune(Pruning &pruning,
           /*inputs*/ const double enumeration_radius, const double preproc_cost,
           const double target_probability, MatGSO<GSO_ZT, GSO_FT> &m,
           int start_row, int end_row)
{
  if (!end_row)
    end_row = m.d;
  Pruner<FT> pruner(enumeration_radius, preproc_cost, target_probability);
  pruner.load_basis_shape(m, start_row, end_row);

  long expo = 0;
  FT gh_radius = m.get_r_exp(start_row, start_row, expo);
  FT root_det = m.get_root_det(start_row, end_row);
  gaussian_heuristic(gh_radius, expo, end_row - start_row, root_det, 1.0);

  pruner.optimize_coefficients(pruning.coefficients);
  pruning.probability = pruner.svp_probability(pruning.coefficients);
  pruning.radius_factor = enumeration_radius/(gh_radius.get_d() * pow(2,expo) );
}


template <class FT>
double svp_probability(const Pruning &pruning){
  Pruner<FT> pru;
  return pru.svp_probability(pruning.coefficients);
}

template <class FT>
double svp_probability(const vector<double> &pr){
  Pruner<FT> pru;
  return pru.svp_probability(pr);
}

template <class FT> void Pruner<FT>::set_tabulated_consts()
{
  for (int i = 0; i < PRUNER_MAX_N; ++i)
  {
    tabulated_factorial[i] = pre_factorial[i];
    tabulated_ball_vol[i]  = pre_ball_vol[i];
  }
  return;
}

/// PUBLIC METHODS

template <class FT>
template <class GSO_ZT, class GSO_FT>
void Pruner<FT>::load_basis_shape(MatGSO<GSO_ZT, GSO_FT> &m, int start_row, int end_row, int /*reset_renorm*/)
{
  if (!end_row)
  {
    end_row = m.d;
  }
  n = end_row - start_row;
  d = n / 2;
  if (!d)
  {
    throw std::runtime_error("Inside Pruner : Needs a dimension n>1");
  }
  vector<double> gso_sq_norms;
  GSO_FT f;
  for (size_t i = 0; i < n; ++i)
  {
    m.get_r(f, start_row+i, start_row+i);
    gso_sq_norms.emplace_back(f.get_d());
  }
  load_basis_shape(gso_sq_norms);
}

template <class FT> void Pruner<FT>::load_basis_shape(const vector<double> &gso_sq_norms, int reset_renorm)
{
  n = gso_sq_norms.size();
  d = n / 2;
  if (!d)
  {
    throw std::runtime_error("Inside Pruner : Needs a dimension n>1");
  }
  FT logvol, tmp;
  logvol = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    r[i] = gso_sq_norms[n - 1 - i];

    logvol += log(r[i]);
  }
  if (reset_renorm){
    renormalization_factor = exp(logvol / (-1.0 * n));
  }

  for (size_t i = 0; i < n; ++i)
  {
    r[i] *= renormalization_factor;
  }
  tmp = 1.;
  for (size_t i = 0; i < 2 * d; ++i)
  {
    tmp *= sqrt(r[i]);
    ipv[i] = 1.0 / tmp;
  }
}

template <class FT> 
void Pruner<FT>::load_basis_shapes(const vector<vector<double> *> &gso_sq_norms_vec)
{
  vec sum_ipv;
  n = (*gso_sq_norms_vec[0]).size();
  for (size_t i = 0; i < n; ++i)
  {
    sum_ipv[i] = 0.;
  }
  int count = gso_sq_norms_vec.size();
  for (int k = 0; k < count; ++k)
  {
    if ((*gso_sq_norms_vec[k]).size() != n)
    {
      throw std::runtime_error("Inside Pruner : loading several bases with different dimensions");  
    }
    int reset_renorm = (k==0);
    load_basis_shape(*gso_sq_norms_vec[k], reset_renorm);
    for (size_t i = 0; i < n; ++i)
    {
      sum_ipv[i] += ipv[i];
    }
  }
  for (size_t i = 0; i < n; ++i){
    ipv[i] = sum_ipv[i] / (1.0 * count);
  }
}

template <class FT> 
template <class GSO_ZT, class GSO_FT>
void Pruner<FT>::load_basis_shapes(vector<MatGSO<GSO_ZT, GSO_FT> *>&m_vec, int start_row, int end_row)
{
  vec sum_ipv;
  if (!end_row)
  {
    end_row = m_vec[0]->d;
  }
  n = end_row - start_row;
  d = n / 2;

  for (size_t i = 0; i < n; ++i)
  {
    sum_ipv[i] = 0.;
  }
  int count = m_vec.size();
  for (int k = 0; k < count; ++k)
  {
    int reset_renorm = (k==0);    
    load_basis_shape(*m_vec[k], start_row, end_row, reset_renorm);
    for (size_t i = 0; i < n; ++i)
    {
      sum_ipv[i] += ipv[i];
    }
  }
  for (size_t i = 0; i < n; ++i){
    ipv[i] = sum_ipv[i] / (1.0 * count);
  }
}


template <class FT>
void Pruner<FT>::optimize_coefficients(/*io*/ vector<double> &pr, /*i*/ const int reset)
{
  evec b;
  if (reset)
  {
    init_coefficients(b);
  }
  else
  {
    load_coefficients(b, pr);
  }
  descent(b);
  save_coefficients(pr, b);
}

// PRIVATE METHODS

template <class FT> void Pruner<FT>::load_coefficients(/*o*/ evec &b, /*i*/ const vector<double> &pr)
{
  for (size_t i = 0; i < d; ++i)
  {
    b[i] = pr[n - 1 - 2 * i];
  }
  if (enforce_bounds(b))
  {
    throw std::runtime_error("Inside Pruner : Ill formed pruning coefficients (must be decreasing, "
                             "starting with two 1.0)");
  }
}

template <class FT> int Pruner<FT>::check_basis_loaded()
{
  if (d)
  {
    return 0;
  }
  throw std::runtime_error("Inside Pruner : No basis loaded");
  return 1;
}

template <class FT> void Pruner<FT>::save_coefficients(/*o*/ vector<double> &pr, /*i*/ const evec &b)
{
  pr.resize(n);
  for (size_t i = 0; i < d; ++i)
  {
    pr[n - 1 - 2 * i] = b[i].get_d();
    pr[n - 2 - 2 * i] = b[i].get_d();
  }
  pr[0] = 1.;
}

template <class FT> inline int Pruner<FT>::enforce_bounds(/*io*/ evec &b, /*opt i*/ const int j)
{
  int status = 0;
  if (b[d - 1] < 1)
  {
    status = 1;
  }
  b[d - 1] = 1;
  for (size_t i = 0; i < d; ++i)
  {
    if (b[i] > 1)
    {
      b[i]   = 1.0;
      status = 1;
    }
    if (b[i] <= .1)
      b[i] = .1;
  }
  for (size_t i = j; i < d - 1; ++i)
  {
    if (b[i + 1] < b[i])
    {
      b[i + 1] = b[i];
      status   = 1;
    }
  }
  for (int i = j - 1; i >= 0; --i)
  {
    if (b[i + 1] < b[i])
    {
      b[i]   = b[i + 1];
      status = 1;
    }
  }
  return status;
}

template <class FT> inline FT Pruner<FT>::eval_poly(const int ld, /*i*/ const poly &p, const FT x)
{
  FT acc;
  acc = 0.0;
  for (int i = ld; i >= 0; --i)
  {
    acc = acc * x;
    acc = acc + p[i];
  }
  return acc;
}

template <class FT> inline void Pruner<FT>::integrate_poly(const int ld, /*io*/ poly &p)
{
  for (int i = ld; i >= 0; --i)
  {
    FT tmp;
    tmp         = i + 1.;
    p[i + 1] = p[i] / tmp;
  }
  p[0] = 0.0;
}

template <class FT> inline FT Pruner<FT>::relative_volume(const int rd, /*i*/ const evec &b)
{
  poly P;
  P[0]   = 1;
  int ld = 0;
  for (int i = rd - 1; i >= 0; --i)
  {
    integrate_poly(ld, P);
    ld++;
    P[0] = -1.0 * eval_poly(ld, P, b[i] / b[rd - 1]);
  }
  if (rd % 2)
  {
    return -1.0 * P[0] * tabulated_factorial[rd];
  }
  else
  {
    return P[0] * tabulated_factorial[rd];
  }
}

template <class FT> inline FT Pruner<FT>::single_enum_cost(/*i*/ const evec &b)
{
  vec rv;  // Relative volumes at each level

  for (size_t i = 0; i < d; ++i)
  {

    rv[2 * i + 1] = relative_volume(i + 1, b);
  }

  rv[0] = 1;
  for (size_t i = 1; i < d; ++i)
  {
    rv[2 * i] = sqrt(rv[2 * i - 1] * rv[2 * i + 1]);  // Interpolate even values
  }

  FT total;
  total = 0.0;
  FT normalized_radius;
  normalized_radius = sqrt(enumeration_radius * renormalization_factor);

  for (size_t i = 0; i < 2 * d; ++i)
  {
    FT tmp;
    tmp = pow_si(normalized_radius, 1 + i) * rv[i] * tabulated_ball_vol[i + 1] *
          sqrt(pow_si(b[i / 2], 1 + i)) * ipv[i];
    total += tmp;
  }
  total /= symmetry_factor;
  // exit(1);
  return total;
}

template <class FT> inline FT Pruner<FT>::svp_probability(/*i*/ const evec &b)
{

  evec b_minus_db;
  FT dx = shell_ratio;

  for (size_t i = 0; i < d; ++i)
  {
    b_minus_db[i] = b[i] / (dx * dx);
    if (b_minus_db[i] > 1)
      b_minus_db[i] = 1;
  }

  FT vol  = relative_volume(d, b);
  FT dxn  = pow_si(dx, 2 * d);
  FT dvol = dxn * relative_volume(d, b_minus_db) - vol;
  return dvol / (dxn - 1.);
}

template <class FT> inline FT Pruner<FT>::repeated_enum_cost(/*i*/ const evec &b)
{

  FT probability = svp_probability(b);

  if (probability >= target_probability)
    return single_enum_cost(b);

  FT trials = log(1.0 - target_probability) / log(1.0 - probability);
  return single_enum_cost(b) * trials + preproc_cost * (trials - 1.0);
}

template <class FT>
void Pruner<FT>::repeated_enum_cost_gradient(/*i*/ const evec &b, /*o*/ evec &res)
{
  evec bpDb;
  res[d - 1] = 0.0;
  for (size_t i = 0; i < d - 1; ++i)
  {
    bpDb = b;
    bpDb[i] *= (1.0 - epsilon);
    enforce_bounds(bpDb, i);
    FT X = repeated_enum_cost(bpDb);

    bpDb = b;
    bpDb[i] *= (1.0 + epsilon);
    enforce_bounds(bpDb, i);
    FT Y   = repeated_enum_cost(bpDb);
    res[i] = (log(X) - log(Y)) / epsilon;
  }
}

template <class FT> int Pruner<FT>::improve(/*io*/ evec &b)
{

  FT cf     = repeated_enum_cost(b);
  FT old_cf = cf;
  evec newb;
  evec gradient;
  repeated_enum_cost_gradient(b, gradient);
  FT norm = 0.0;

  // normalize the gradient
  for (size_t i = 0; i < d; ++i)
  {
    norm += gradient[i] * gradient[i];
    newb[i] = b[i];
  }

  norm /= (1.0 * (1. * d));
  norm = sqrt(norm);
  if (norm <= 0.)
    return 0;

  for (size_t i = 0; i < d; ++i)
  {
    gradient[i] /= norm;
  }
  FT new_cf;

  FT step = min_step;
  size_t i;

  for (i = 0;; ++i)
  {
    for (size_t i = 0; i < d; ++i)
    {
      newb[i] = newb[i] + step * gradient[i];
    }

    enforce_bounds(newb);
    new_cf = repeated_enum_cost(newb);

    if (new_cf >= cf)
    {
      break;
    }
    b  = newb;
    cf = new_cf;
    step *= step_factor;
  }
  if (cf > old_cf * min_cf_decrease)
  {
    return 0;
  }
  return i;
}

template <class FT> void Pruner<FT>::descent(/*io*/ evec &b)
{
  while (improve(b))
  {
  };
}

template <class FT> void Pruner<FT>::init_coefficients(evec &b)
{
  for (size_t i = 0; i < d; ++i)
  {
    b[i] = .1 + ((1. * i) / d);
  }
  enforce_bounds(b);
}

template <class FT, class GSO_ZT, class GSO_FT>
void prune(/*output*/ vector<double> &pr, double &probability,
           /*inputs*/ const double enumeration_radius, const double preproc_cost,
           const double target_probability, MatGSO<GSO_ZT, GSO_FT> &m,
           int start_row, int end_row)
{

  Pruner<FP_NR<double> > pruner(enumeration_radius, preproc_cost, target_probability);
  pruner.load_basis_shape(m, start_row, end_row);
  pruner.optimize_coefficients(pr);
  probability = pruner.svp_probability(pr);
}

template <class FT, class GSO_ZT, class GSO_FT>
Pruning prune(/*inputs*/ const double enumeration_radius, const double preproc_cost,
              const double target_probability, MatGSO<GSO_ZT, GSO_FT> &m,
              int start_row, int end_row)
{
  Pruning pruning;
  if (!end_row)
    end_row = m.d;
  Pruner<FP_NR<double> > pruner(enumeration_radius, preproc_cost, target_probability);
  pruner.load_basis_shape(m, start_row, end_row);

  long expo = 0;
  FT gh_radius = m.get_r_exp(start_row, start_row, expo);
  FT root_det = m.get_root_det(start_row, end_row);
  gaussian_heuristic(gh_radius, expo, end_row - start_row, root_det, 1.0);

  pruner.optimize_coefficients(pruning.coefficients);
  pruning.probability = pruner.svp_probability(pruning.coefficients);
  pruning.radius_factor = enumeration_radius/(gh_radius.get_d() * pow(2,expo) );
  return pruning;
}


template <class FT, class GSO_ZT, class GSO_FT>
Pruning prune(/*inputs*/ const double enumeration_radius, const double preproc_cost,
              const double target_probability, vector<MatGSO<GSO_ZT, GSO_FT> *> &m,
              int start_row, int end_row)
{
  Pruning pruning;
  if (!end_row)
    end_row = m[0]->d;
  Pruner<FP_NR<double> > pruner(enumeration_radius, preproc_cost, target_probability);
  pruner.load_basis_shapes(m, start_row, end_row);


  FT gh_radius = 0.0;
  FT root_det = 0.0;
  for(auto it = m.begin(); it != m.end(); ++it)
  {
    FT tmp;
    (*it)->get_r(tmp, start_row, start_row);
    gh_radius += tmp;
    root_det += (*it)->get_root_det(start_row, end_row);
  }
  gh_radius /= m.size();
  root_det /= m.size();

  int expo = 0;
  gaussian_heuristic(gh_radius, expo, end_row - start_row, root_det, 1.0);

  pruner.optimize_coefficients(pruning.coefficients);
  pruning.probability = pruner.svp_probability(pruning.coefficients);
  pruning.radius_factor = enumeration_radius/(gh_radius.get_d() * pow(2,expo) );
  return pruning;
}


/** instantiate functions **/

template class Pruner<FP_NR<double> >;
template void prune<FP_NR<double>, Z_NR<mpz_t>, FP_NR<double> > (vector<double> &, double &, const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<double> >&, int, int);
template Pruning prune<FP_NR<double>, Z_NR<mpz_t>, FP_NR<double> > (const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<double> >&, int, int);
template Pruning prune<FP_NR<double>, Z_NR<mpz_t>, FP_NR<double> > (const double, const double, const double, vector<MatGSO<Z_NR<mpz_t>, FP_NR<double> >* >&, int, int);
template double svp_probability<FP_NR<double> >(const Pruning&);
template double svp_probability<FP_NR<double> >(const vector<double> &);

#ifdef FPLLL_WITH_LONG_DOUBLE
template class Pruner<FP_NR<long double> >;
template void prune<FP_NR<long double>, Z_NR<mpz_t>, FP_NR<long double> > (vector<double> &, double &, const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<long double> >&, int, int);
template Pruning prune<FP_NR<long double>, Z_NR<mpz_t>, FP_NR<long double> > (const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<long double> >&, int, int);
template Pruning prune<FP_NR<long double>, Z_NR<mpz_t>, FP_NR<long double> > (const double, const double, const double, vector<MatGSO<Z_NR<mpz_t>, FP_NR<long double> >* >&, int, int);
template double svp_probability<FP_NR<long double> >(const Pruning&);
template double svp_probability<FP_NR<long double> >(const vector<double> &);
#endif

#ifdef FPLLL_WITH_QD
template class Pruner<FP_NR<dd_real> >;
template void prune<FP_NR<dd_real>, Z_NR<mpz_t>, FP_NR<dd_real> > (vector<double> &, double &, const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<dd_real> >&, int, int);
template Pruning prune<FP_NR<dd_real>, Z_NR<mpz_t>, FP_NR<dd_real> > (const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<dd_real> >&, int, int);
template Pruning prune<FP_NR<dd_real>, Z_NR<mpz_t>, FP_NR<dd_real> > (const double, const double, const double, vector<MatGSO<Z_NR<mpz_t>, FP_NR<dd_real> >* >&, int, int);
template double svp_probability<FP_NR<dd_real> >(const Pruning&);
template double svp_probability<FP_NR<dd_real> >(const vector<double> &);

template class Pruner<FP_NR<qd_real> >;
template void prune<FP_NR<qd_real>, Z_NR<mpz_t>, FP_NR<qd_real> > (vector<double> &, double &, const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<qd_real> >&, int, int);
template Pruning prune<FP_NR<qd_real>, Z_NR<mpz_t>, FP_NR<qd_real> > (const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<qd_real> >&, int, int);
template Pruning prune<FP_NR<qd_real>, Z_NR<mpz_t>, FP_NR<qd_real> > (const double, const double, const double, vector<MatGSO<Z_NR<mpz_t>, FP_NR<qd_real> >* >&, int, int);
template double svp_probability<FP_NR<qd_real> >(const Pruning&);
template double svp_probability<FP_NR<qd_real> >(const vector<double> &);
#endif

#ifdef FPLLL_WITH_DPE
template class Pruner<FP_NR<dpe_t> >;
template void prune<FP_NR<dpe_t>, Z_NR<mpz_t>, FP_NR<dpe_t> > (vector<double> &, double &, const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<dpe_t> >&, int, int);
template Pruning prune<FP_NR<dpe_t>, Z_NR<mpz_t>, FP_NR<dpe_t> > (const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<dpe_t> >&, int, int);
template Pruning prune<FP_NR<dpe_t>, Z_NR<mpz_t>, FP_NR<dpe_t> > (const double, const double, const double, vector<MatGSO<Z_NR<mpz_t>, FP_NR<dpe_t> >* >&, int, int);
template double svp_probability<FP_NR<dpe_t> >(const Pruning&);
template double svp_probability<FP_NR<dpe_t> >(const vector<double> &);
#endif

template class Pruner<FP_NR<mpfr_t> >;
template void prune<FP_NR<mpfr_t>, Z_NR<mpz_t>, FP_NR<mpfr_t> > (vector<double> &, double &, const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t> >&, int, int);
template Pruning prune<FP_NR<mpfr_t>, Z_NR<mpz_t>, FP_NR<mpfr_t> > (const double, const double, const double, MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t> >&, int, int);
template Pruning prune<FP_NR<mpfr_t>, Z_NR<mpz_t>, FP_NR<mpfr_t> > (const double, const double, const double, vector<MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t> >* >&, int, int);
template double svp_probability<FP_NR<mpfr_t> >(const Pruning&);
template double svp_probability<FP_NR<mpfr_t> >(const vector<double> &);

FPLLL_END_NAMESPACE
