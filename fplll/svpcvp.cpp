/* Copyright (C) 2008-2011 Xavier Pujol.
    (C) 2015 Michael Walter.

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

#include "svpcvp.h"
#include "enum/enumerate.h"
#include "enum/topenum.h"

FPLLL_BEGIN_NAMESPACE

/* Shortest vector problem
   ======================= */

/* Returns i such that the shortest vector of L(b) belongs to
   L(b_0,...,b_(i-1)), assuming that the error on rdiag's is less than 100%.
   If b is LLL-reduced, then for any reasonnable dimension,
   max(rdiag[0],...,rdiag[i-1]) / min(rdiag[0],...,rdiag[i-1])
   is much smaller than numeric_limits<double>::max */
static int last_useful_index(const Matrix<FP_NR<mpfr_t>> &r)
{
  int i;
  FP_NR<mpfr_t> rdiag_min_value;
  rdiag_min_value.mul_2si(r(0, 0), 1);
  for (i = r.get_rows() - 1; i > 0; i--)
  {
    if (r(i, i) <= rdiag_min_value)
      break;
  }
  return i + 1;
}

/* Finds the shortest vector of the basis b and returns its squared norm in
   basisMin */
static void get_basis_min(Z_NR<mpz_t> &basis_min, const ZZ_mat<mpz_t> &b, int first, int last)
{
  Z_NR<mpz_t> sq_norm;
  int n = b.get_cols();
  b[first].dot_product(basis_min, b[first], n);

  for (int i = first + 1; i < last; i++)
  {
    b[i].dot_product(sq_norm, b[i], n);
    if (sq_norm < basis_min)
      basis_min = sq_norm;
  }
}

static void get_basis_min(Z_NR<mpz_t> &basis_min, MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso,
                          int first, int last)
{
  Z_NR<mpz_t> sq_norm;
  gso.get_int_gram(basis_min, first, first);
  for (int i = first + 1; i < last; i++)
  {
    gso.get_int_gram(sq_norm, i, i);
    if (sq_norm < basis_min)
      basis_min = sq_norm;
  }
}

static bool enumerate_svp(int d, MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso,
                          FP_NR<mpfr_t> &max_dist, ErrorBoundedEvaluator &evaluator,
                          const vector<enumf> &pruning, int flags)
{
  Enumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>> enumobj(gso, evaluator);
  bool dual = (flags & SVP_DUAL);
  if (d == 1 || !pruning.empty() || dual)
  {
    enumobj.enumerate(0, d, max_dist, 0, vector<FP_NR<mpfr_t>>(), vector<enumxt>(), pruning, dual);
  }
  else
  {
    Enumerator enumerator(d, gso.get_mu_matrix(), gso.get_r_matrix());
    FP_NR<mpfr_t> bestdist = -1;
    while (enumerator.enum_next(max_dist))
    {
      if (flags & SVP_VERBOSE)
      {
        cerr << enumerator.get_sub_tree();
        if (evaluator.eval_mode != EVALMODE_SV)
          cerr << " (count=2*" << evaluator.size() << ")";
      }

      /* Enumerates short vectors only in enumerator.get_sub_tree()
        (about maxVolume iterations or less) */
      enumobj.enumerate(0, d, max_dist, 0, vector<FP_NR<mpfr_t>>(), enumerator.get_sub_tree(),
                        pruning);

      if (flags & SVP_VERBOSE)
      {
        cerr << "\r" << (char)27 << "[K";
        if (evaluator.eval_mode == EVALMODE_SV && !evaluator.empty() &&
            evaluator.begin()->first != bestdist)
        {
          bestdist = evaluator.begin()->first;
          cerr << "Solution norm^2=" << bestdist << " value=" << evaluator.begin()->second << endl;
        }
      }
    }
  }
  return !evaluator.empty();
}

static int shortest_vector_ex(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord, SVPMethod method,
                              const vector<double> &pruning, int flags, EvaluatorMode eval_mode,
                              long long &sol_count,
                              vector<vector<Z_NR<mpz_t>>> *subsol_coord = nullptr,
                              vector<enumf> *subsol_dist                = nullptr,
                              vector<vector<Z_NR<mpz_t>>> *auxsol_coord = nullptr,
                              vector<enumf> *auxsol_dist = nullptr, int max_aux_sols = 0)
{
  bool findsubsols = (subsol_coord != nullptr) && (subsol_dist != nullptr);
  bool findauxsols = (auxsol_coord != nullptr) && (auxsol_dist != nullptr) && (max_aux_sols != 0);

  // d = lattice dimension (note that it might decrease during preprocessing)
  int d = b.get_rows();
  // n = dimension of the space
  int n = b.get_cols();

  FPLLL_CHECK(d > 0 && n > 0, "shortestVector: empty matrix");
  FPLLL_CHECK(d <= n, "shortestVector: number of vectors > size of the vectors");

  // Sets the floating-point precision
  // Error bounds on GSO are valid if prec >= minprec
  double rho;
  int min_prec = gso_min_prec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec     = max(53, min_prec + 10);
  int old_prec = FP_NR<mpfr_t>::set_prec(prec);

  // Allocates space for vectors and matrices in constructors
  ZZ_mat<mpz_t> empty_mat;
  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> gso(b, empty_mat, empty_mat, GSO_INT_GRAM);
  FP_NR<mpfr_t> max_dist;
  Z_NR<mpz_t> int_max_dist;
  Z_NR<mpz_t> itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.update_gso();
  gen_zero_vect(sol_coord, d);

  // If the last b_i* are too large, removes them to avoid an underflow
  int new_d = last_useful_index(gso.get_r_matrix());
  if (new_d < d)
  {
    // FPLLL_TRACE("Ignoring the last " << d - new_d << " vector(s)");
    d = new_d;
  }

  if (flags & SVP_DUAL)
  {
    max_dist = 1.0 / gso.get_r_exp(d - 1, d - 1);
    if (flags & SVP_VERBOSE)
    {
      cout << "max_dist = " << max_dist << endl;
    }
  }
  else
  {
    /* Computes a bound for the enumeration. This bound would work for an
       exact algorithm, but we will increase it later to ensure that the fp
       algorithm finds a solution */
    get_basis_min(int_max_dist, b, 0, d);
    max_dist.set_z(int_max_dist, GMP_RNDU);
  }

  // Initializes the evaluator of solutions
  ErrorBoundedEvaluator *evaluator;
  if (method == SVPM_FAST)
  {
    evaluator =
        new FastErrorBoundedEvaluator(d, gso.get_mu_matrix(), gso.get_r_matrix(), eval_mode,
                                      max_aux_sols + 1, EVALSTRATEGY_BEST_N_SOLUTIONS, findsubsols);
  }
  else if (method == SVPM_PROVED)
  {
    ExactErrorBoundedEvaluator *p = new ExactErrorBoundedEvaluator(
        // d, b, gso.get_mu_matrix(), gso.get_r_matrix()
        gso, eval_mode, max_aux_sols + 1, EVALSTRATEGY_BEST_N_SOLUTIONS, findsubsols);
    p->int_max_dist = int_max_dist;
    evaluator       = p;
  }
  else
  {
    FPLLL_ABORT("shortestVector: invalid evaluator type");
  }
  evaluator->init_delta_def(prec, rho, true);

  if (!(flags & SVP_OVERRIDE_BND) && (eval_mode == EVALMODE_SV || method == SVPM_PROVED))
  {
    FP_NR<mpfr_t> ftmp1;
    bool result = evaluator->get_max_error_aux(max_dist, true, ftmp1);
    FPLLL_CHECK(result, "shortestVector: cannot compute an initial bound");
    max_dist.add(max_dist, ftmp1, GMP_RNDU);
  }

  // Main loop of the enumeration
  enumerate_svp(d, gso, max_dist, *evaluator, pruning, flags);

  int result = RED_ENUM_FAILURE;
  if (eval_mode != EVALMODE_SV)
  {
    result    = RED_SUCCESS;
    sol_count = evaluator->sol_count * 2;
  }
  else if (!evaluator->empty())
  {
    /*FP_NR<mpfr_t> fMaxError;
    validMaxError = evaluator->get_max_error(fMaxError);
    max_error = fMaxError.get_d(GMP_RNDU);*/
    for (int i = 0; i < d; i++)
    {
      itmp1.set_f(evaluator->begin()->second[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
    }
    result = RED_SUCCESS;
  }

  if (findsubsols)
  {
    subsol_coord->clear();
    subsol_dist->clear();
    subsol_dist->resize(evaluator->sub_solutions.size());
    for (size_t i = 0; i < evaluator->sub_solutions.size(); ++i)
    {
      (*subsol_dist)[i] = evaluator->sub_solutions[i].first.get_d();

      vector<Z_NR<mpz_t>> ss_c;
      for (size_t j = 0; j < evaluator->sub_solutions[i].second.size(); ++j)
      {
        itmp1.set_f(evaluator->sub_solutions[i].second[j]);
        ss_c.emplace_back(itmp1);
      }
      subsol_coord->emplace_back(std::move(ss_c));
    }
  }
  if (findauxsols)
  {
    auxsol_coord->clear();
    auxsol_dist->clear();
    // iterators over all solutions
    auto it = evaluator->begin(), itend = evaluator->end();
    // skip shortest solution
    ++it;
    for (; it != itend; ++it)
    {
      auxsol_dist->push_back(it->first.get_d());

      vector<Z_NR<mpz_t>> as_c;
      for (size_t j = 0; j < it->second.size(); ++j)
      {
        itmp1.set_f(it->second[j]);
        as_c.emplace_back(itmp1);
      }
      auxsol_coord->emplace_back(std::move(as_c));
    }
  }

  delete evaluator;
  FP_NR<mpfr_t>::set_prec(old_prec);
  return result;
}

int shortest_vector(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord, SVPMethod method, int flags)
{
  long long tmp;
  return shortest_vector_ex(b, sol_coord, method, vector<double>(), flags, EVALMODE_SV, tmp);
}

int shortest_vector_pruning(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord,
                            const vector<double> &pruning, int flags)
{
  long long tmp;
  return shortest_vector_ex(b, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp);
}

int shortest_vector_pruning(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &subsol_coord, vector<enumf> &subsol_dist,
                            const vector<double> &pruning, int flags)
{
  long long tmp;
  return shortest_vector_ex(b, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp,
                            &subsol_coord, &subsol_dist);
}

int shortest_vector_pruning(ZZ_mat<mpz_t> &b, vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &auxsol_coord, vector<enumf> &auxsol_dist,
                            const int max_aux_sols, const vector<double> &pruning, int flags)
{
  long long tmp;
  return shortest_vector_ex(b, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp, nullptr,
                            nullptr, &auxsol_coord, &auxsol_dist, max_aux_sols);
}

//////////////////////////////////////////
////// SVP FOR GSO OBJECT  ///////////////
//////////////////////////////////////////

static int shortest_vector_ex(
    MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso, vector<Z_NR<mpz_t>> &sol_coord,
    SVPMethod method, const vector<double> &pruning, int flags, EvaluatorMode eval_mode,
    long long &sol_count, vector<vector<Z_NR<mpz_t>>> *subsol_coord = nullptr,
    vector<enumf> *subsol_dist = nullptr, vector<vector<Z_NR<mpz_t>>> *auxsol_coord = nullptr,
    vector<enumf> *auxsol_dist = nullptr, int max_aux_sols = 0, bool merge_sol_in_aux = false)
{
  bool findsubsols = (subsol_coord != nullptr) && (subsol_dist != nullptr);
  bool findauxsols = (auxsol_coord != nullptr) && (auxsol_dist != nullptr) && (max_aux_sols != 0);

  // d = lattice dimension (note that it might decrease during preprocessing)
  // int d = b.get_rows();
  int d = gso.d;  // Number of rows of b in the GSO

  // n = dimension of the space
  // int n = b.get_cols();
  int n = gso.get_cols_of_b();  // number of columns of b in the GSO

  FPLLL_CHECK(d > 0 && n > 0, "shortestVector: empty matrix");
  FPLLL_CHECK(d <= n, "shortestVector: number of vectors > size of the vectors");

  // Sets the floating-point precision
  // Error bounds on GSO are valid if prec >= minprec
  double rho;
  int min_prec = gso_min_prec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec     = max(53, min_prec + 10);
  int old_prec = FP_NR<mpfr_t>::set_prec(prec);

  // Allocates space for vectors and matrices in constructors
  FP_NR<mpfr_t> max_dist;
  Z_NR<mpz_t> int_max_dist;
  Z_NR<mpz_t> itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.update_gso();
  gen_zero_vect(sol_coord, d);

  // If the last b_i* are too large, removes them to avoid an underflow
  int new_d = last_useful_index(gso.get_r_matrix());
  if (new_d < d)
  {
    // FPLLL_TRACE("Ignoring the last " << d - new_d << " vector(s)");
    d = new_d;
  }

  if (flags & SVP_DUAL)
  {
    max_dist = 1.0 / gso.get_r_exp(d - 1, d - 1);
    if (flags & SVP_VERBOSE)
    {
      cout << "max_dist = " << max_dist << endl;
    }
  }
  else
  {
    // Computes a bound for the enumeration. This bound would work for an
    //   exact algorithm, but we will increase it later to ensure that the fp
    //   algorithm finds a solution

    // Use the GSO version of get_basis_min
    get_basis_min(int_max_dist, gso, 0, d);
    max_dist.set_z(int_max_dist, GMP_RNDU);
  }

  // Initializes the evaluator of solutions
  ErrorBoundedEvaluator *evaluator;
  if (method == SVPM_FAST)
  {
    evaluator =
        new FastErrorBoundedEvaluator(d, gso.get_mu_matrix(), gso.get_r_matrix(), eval_mode,
                                      max_aux_sols + 1, EVALSTRATEGY_BEST_N_SOLUTIONS, findsubsols);
  }
  else if (method == SVPM_PROVED)
  {
    ExactErrorBoundedEvaluator *p = new ExactErrorBoundedEvaluator(
        gso, eval_mode, max_aux_sols + 1, EVALSTRATEGY_BEST_N_SOLUTIONS, findsubsols);
    p->int_max_dist = int_max_dist;
    evaluator       = p;
  }
  else
  {
    FPLLL_ABORT("shortestVector: invalid evaluator type");
  }
  evaluator->init_delta_def(prec, rho, true);

  if (!(flags & SVP_OVERRIDE_BND) && (eval_mode == EVALMODE_SV || method == SVPM_PROVED))
  {
    FP_NR<mpfr_t> ftmp1;
    bool result = evaluator->get_max_error_aux(max_dist, true, ftmp1);
    FPLLL_CHECK(result, "shortestVector: cannot compute an initial bound");
    max_dist.add(max_dist, ftmp1, GMP_RNDU);
  }

  // Main loop of the enumeration
  enumerate_svp(d, gso, max_dist, *evaluator, pruning, flags);  // Only uses r and mu

  int result = RED_ENUM_FAILURE;
  if (eval_mode != EVALMODE_SV)
  {
    result    = RED_SUCCESS;
    sol_count = evaluator->sol_count * 2;
  }
  else if (!evaluator->empty())
  {
    // FP_NR<mpfr_t> fMaxError;
    // validMaxError = evaluator->get_max_error(fMaxError);
    // max_error = fMaxError.get_d(GMP_RNDU);
    for (int i = 0; i < d; i++)
    {
      itmp1.set_f(evaluator->begin()->second[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
    }
    result = RED_SUCCESS;
  }

  if (findsubsols)
  {
    subsol_coord->clear();
    subsol_dist->clear();
    subsol_dist->resize(evaluator->sub_solutions.size());
    for (size_t i = 0; i < evaluator->sub_solutions.size(); ++i)
    {
      (*subsol_dist)[i] = evaluator->sub_solutions[i].first.get_d();

      vector<Z_NR<mpz_t>> ss_c;
      for (size_t j = 0; j < evaluator->sub_solutions[i].second.size(); ++j)
      {
        itmp1.set_f(evaluator->sub_solutions[i].second[j]);
        ss_c.emplace_back(itmp1);
      }
      subsol_coord->emplace_back(std::move(ss_c));
    }
  }
  if (findauxsols)
  {
    auxsol_coord->clear();
    auxsol_dist->clear();
    // iterators over all solutions
    auto it = evaluator->begin(), itend = evaluator->end();
    // skip shortest solution
    if (!merge_sol_in_aux)
      ++it;
    for (; it != itend; ++it)
    {
      auxsol_dist->push_back(it->first.get_d());

      vector<Z_NR<mpz_t>> as_c;
      for (size_t j = 0; j < it->second.size(); ++j)
      {
        itmp1.set_f(it->second[j]);
        as_c.emplace_back(itmp1);
      }
      auxsol_coord->emplace_back(std::move(as_c));
    }
  }

  delete evaluator;
  FP_NR<mpfr_t>::set_prec(old_prec);
  return result;
}

int shortest_vector(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso,
                    vector<Z_NR<mpz_t>> &sol_coord, SVPMethod method, int flags)
{
  long long tmp;
  return shortest_vector_ex(gso, sol_coord, method, vector<double>(), flags, EVALMODE_SV, tmp);
}

int shortest_vectors(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso,
                     vector<vector<Z_NR<mpz_t>>> &sol_coord, vector<enumf> &sol_dist,
                     const int max_sols, SVPMethod method, int flags)
{
  long long tmp;
  vector<Z_NR<mpz_t>> sol_coord_tmp;
  return shortest_vector_ex(gso, sol_coord_tmp, method, vector<double>(), flags, EVALMODE_SV, tmp,
                            nullptr, nullptr, &sol_coord, &sol_dist, max_sols - 1, true);
}

int shortest_vector_pruning(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso,
                            vector<Z_NR<mpz_t>> &sol_coord, const vector<double> &pruning,
                            int flags)
{
  long long tmp;
  return shortest_vector_ex(gso, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp);
}

int shortest_vector_pruning(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso,
                            vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &subsol_coord, vector<enumf> &subsol_dist,
                            const vector<double> &pruning, int flags)
{
  long long tmp;
  return shortest_vector_ex(gso, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp,
                            &subsol_coord, &subsol_dist);
}

int shortest_vector_pruning(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso,
                            vector<Z_NR<mpz_t>> &sol_coord,
                            vector<vector<Z_NR<mpz_t>>> &auxsol_coord, vector<enumf> &auxsol_dist,
                            const int max_aux_sols, const vector<double> &pruning, int flags)
{
  long long tmp;
  return shortest_vector_ex(gso, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp, nullptr,
                            nullptr, &auxsol_coord, &auxsol_dist, max_aux_sols);
}

///////////////////////////////////
//////END SVP FOR GSO OBJECT //////
///////////////////////////////////

/* Closest vector problem
   ====================== */

static void get_gscoords(const Matrix<FP_NR<mpfr_t>> &matrix, const Matrix<FP_NR<mpfr_t>> &mu,
                         const Matrix<FP_NR<mpfr_t>> &r, const vector<FP_NR<mpfr_t>> &v,
                         vector<FP_NR<mpfr_t>> &vcoord)
{

  int n = matrix.get_rows(), m = matrix.get_cols();

  if (static_cast<int>(vcoord.size()) != n)
    vcoord.resize(n);
  FPLLL_DEBUG_CHECK(mu.get_rows() == n && mu.get_cols() == n && r.get_rows() == n &&
                    r.get_cols() == n && static_cast<int>(v.size()) == m);

  for (int i = 0; i < n; i++)
  {
    vcoord[i] = 0.0;
    for (int j = 0; j < m; j++)
      vcoord[i].addmul(v[j], matrix(i, j));
    for (int j = 0; j < i; j++)
      vcoord[i].submul(mu(i, j), vcoord[j]);
  }
  for (int i = 0; i < n; i++)
  {
    vcoord[i].div(vcoord[i], r(i, i));
  }
}

static void babai(const FP_mat<mpfr_t> &matrix, const Matrix<FP_NR<mpfr_t>> &mu,
                  const Matrix<FP_NR<mpfr_t>> &r, const vector<FP_NR<mpfr_t>> &target,
                  vector<FP_NR<mpfr_t>> &target_coord)
{

  int d = matrix.get_rows();
  get_gscoords(matrix, mu, r, target, target_coord);
  for (int i = d - 1; i >= 0; i--)
  {
    target_coord[i].rnd(target_coord[i]);
    for (int j = 0; j < i; j++)
      target_coord[j].submul(mu(i, j), target_coord[i]);
  }
}

int closest_vector(ZZ_mat<mpz_t> &b, const vector<Z_NR<mpz_t>> &int_target,
                   vector<Z_NR<mpz_t>> &sol_coord, int method, int flags)
{
  // d = lattice dimension (note that it might decrease during preprocessing)
  int d = b.get_rows();
  // n = dimension of the space
  int n = b.get_cols();

  FPLLL_CHECK(d > 0 && n > 0, "closestVector: empty matrix");
  FPLLL_CHECK(d <= n, "closestVector: number of vectors > size of the vectors");

  // Sets the floating-point precision
  // Error bounds on GSO are valid if prec >= minprec
  double rho;
  int min_prec = gso_min_prec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec     = max(53, min_prec + 10);
  int old_prec = FP_NR<mpfr_t>::set_prec(prec);

  // Allocates space for vectors and matrices in constructors
  ZZ_mat<mpz_t> empty_mat;
  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> gso(b, empty_mat, empty_mat, GSO_INT_GRAM);
  vector<FP_NR<mpfr_t>> target_coord;
  FP_NR<mpfr_t> max_dist;
  Z_NR<mpz_t> itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.update_gso();
  gen_zero_vect(sol_coord, d);

  /* Applies Babai's algorithm. Because we use fp, it might be necessary to
      do it several times (if ||target|| >> ||b_i||) */
  FP_mat<mpfr_t> float_matrix(d, n);
  vector<FP_NR<mpfr_t>> target(n), babai_sol;
  vector<Z_NR<mpz_t>> int_new_target = int_target;

  for (int i = 0; i < d; i++)
    for (int j = 0; j < n; j++)
      float_matrix(i, j).set_z(b(i, j));

  for (int loop_idx = 0;; loop_idx++)
  {
    if (loop_idx >= 0x100 && ((loop_idx & (loop_idx - 1)) == 0))
      FPLLL_INFO("warning: possible infinite loop in Babai's algorithm");

    for (int i = 0; i < n; i++)
    {
      target[i].set_z(int_new_target[i]);
    }
    babai(float_matrix, gso.get_mu_matrix(), gso.get_r_matrix(), target, babai_sol);
    int idx;
    for (idx = 0; idx < d && babai_sol[idx] >= -1 && babai_sol[idx] <= 1; idx++)
    {
    }
    if (idx == d)
      break;

    for (int i = 0; i < d; i++)
    {
      itmp1.set_f(babai_sol[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
      for (int j = 0; j < n; j++)
        int_new_target[j].submul(itmp1, b(i, j));
    }
  }
  // FPLLL_TRACE("BabaiSol=" << sol_coord);
  get_gscoords(float_matrix, gso.get_mu_matrix(), gso.get_r_matrix(), target, target_coord);

  /* Computes a very large bound to make the algorithm work
      until the first solution is found */
  max_dist = 0.0;
  for (int i = 1; i < d; i++)
  {
    // get_r_exp(i, i) = r(i, i) because gso is initialized without GSO_ROW_EXPO
    max_dist.add(max_dist, gso.get_r_exp(i, i));
  }

  vector<int> max_indices;
  if (method & CVPM_PROVED)
  {
    // For Exact CVP, we need to reset enum below depth with maximal r_i
    max_indices = vector<int>(d);
    int cur, max_index, previous_max_index;
    previous_max_index = max_index = d - 1;
    FP_NR<mpfr_t> max_val;

    while (max_index > 0)
    {
      max_val = gso.get_r_exp(max_index, max_index);
      for (cur = previous_max_index - 1; cur >= 0; --cur)
      {
        if (max_val <= gso.get_r_exp(cur, cur))
        {
          max_val   = gso.get_r_exp(cur, cur);
          max_index = cur;
        }
      }
      for (cur = max_index; cur < previous_max_index; ++cur)
        max_indices[cur] = max_index;
      max_indices[previous_max_index] = previous_max_index;
      previous_max_index              = max_index;
      --max_index;
    }
  }
  FPLLL_TRACE("max_indices " << max_indices);

  FastErrorBoundedEvaluator evaluator(n, gso.get_mu_matrix(), gso.get_r_matrix(), EVALMODE_CV);

  // Main loop of the enumeration
  Enumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>> enumobj(gso, evaluator, max_indices);
  enumobj.enumerate(0, d, max_dist, 0, target_coord);

  int result = RED_ENUM_FAILURE;
  if (!evaluator.empty())
  {
    FPLLL_TRACE("evaluator.bestsol_coord=" << evaluator.begin()->second);
    if (flags & CVP_VERBOSE)
      FPLLL_INFO("max_dist=" << max_dist);
    for (int i = 0; i < d; i++)
    {
      itmp1.set_f(evaluator.begin()->second[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
    }
    result = RED_SUCCESS;
  }

  FP_NR<mpfr_t>::set_prec(old_prec);
  return result;
}

FPLLL_END_NAMESPACE
