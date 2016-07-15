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
static int last_useful_index(const Matrix<Float> &r)
{
  int i;
  Float rdiag_min_value;
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
static void get_basis_min(Integer &basis_min, const IntMatrix &b, int first, int last)
{
  Integer sq_norm;
  int n = b.get_cols();
  sqr_norm(basis_min, b[first], n);

  for (int i = first + 1; i < last; i++)
  {
    sqr_norm(sq_norm, b[i], n);
    if (sq_norm < basis_min)
      basis_min = sq_norm;
  }
}

static bool enumerate_svp(int d, MatGSO<Integer, Float> &gso, Float &max_dist,
                          Evaluator<Float> &evaluator, const vector<enumf> &pruning, int flags)
{
  Enumeration<Float> enumobj(gso, evaluator);
  bool dual = (flags & SVP_DUAL);
  if (d == 1 || !pruning.empty() || dual)
  {
    enumobj.enumerate(0, d, max_dist, 0, vector<Float>(), vector<enumxt>(), pruning, dual);
  }
  else
  {
    Enumerator enumerator(d, gso.get_mu_matrix(), gso.get_r_matrix());
    while (enumerator.enumNext(max_dist))
    {
      if (flags & SVP_VERBOSE)
      {
        evaluator.newSolFlag = false;
        cerr << enumerator.getSubTree();
        if (evaluator.evalMode != EVALMODE_SV)
          cerr << " (count=2*" << evaluator.solCount << ")";
      }

      /* Enumerates short vectors only in enumerator.getSubTree()
        (about maxVolume iterations or less) */
      enumobj.enumerate(0, d, max_dist, 0, FloatVect(), enumerator.getSubTree(), pruning);

      if (flags & SVP_VERBOSE)
      {
        cerr << "\r" << (char)27 << "[K";
        if (evaluator.evalMode == EVALMODE_SV && evaluator.newSolFlag)
          cerr << "Solution norm^2=" << evaluator.lastPartialDist
               << " value=" << evaluator.sol_coord << endl;
      }
    }
  }
  return !evaluator.sol_coord.empty();
}

static int shortest_vector_ex(IntMatrix &b, IntVect &sol_coord, SVPMethod method,
                              const vector<double> &pruning, int flags, EvaluatorMode evalMode,
                              long long &solCount, vector<IntVect> *subsol_coord = nullptr,
                              vector<enumf> *subsolDist = nullptr)
{
  bool findsubsols = (subsol_coord != nullptr) && (subsolDist != nullptr);

  // d = lattice dimension (note that it might decrease during preprocessing)
  int d = b.get_rows();
  // n = dimension of the space
  int n = b.get_cols();

  FPLLL_CHECK(d > 0 && n > 0, "shortestVector: empty matrix");
  FPLLL_CHECK(d <= n, "shortestVector: number of vectors > size of the vectors");

  // Sets the floating-point precision
  // Error bounds on GSO are valid if prec >= minprec
  double rho;
  int minPrec = gso_min_prec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec    = max(53, minPrec + 10);
  int oldprec = Float::set_prec(prec);

  // Allocates space for vectors and matrices in constructors
  IntMatrix empty_mat;
  MatGSO<Integer, Float> gso(b, empty_mat, empty_mat, GSO_INT_GRAM);
  Float max_dist;
  Integer int_max_dist;
  Integer itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.update_gso();
  gen_zero_vect(sol_coord, d);

  // If the last b_i* are too large, removes them to avoid an underflow
  int newD = last_useful_index(gso.get_r_matrix());
  if (newD < d)
  {
    // FPLLL_TRACE("Ignoring the last " << d - newD << " vector(s)");
    d = newD;
  }

  if (flags & SVP_DUAL)
  {
    max_dist = gso.get_r_exp(d - 1, d - 1);
    Float one;
    one = 1.0;
    max_dist.div(one, max_dist);
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
  Evaluator<Float> *evaluator;
  if (method == SVPM_FAST)
  {
    evaluator = new FastEvaluator<Float>(d, gso.get_mu_matrix(), gso.get_r_matrix(), evalMode, 0,
                                         findsubsols);
  }
  else if (method == SVPM_PROVED)
  {
    ExactEvaluator *p =
        new ExactEvaluator(d, b, gso.get_mu_matrix(), gso.get_r_matrix(), evalMode, 0, findsubsols);
    p->int_max_dist = int_max_dist;
    evaluator     = p;
  }
  else
  {
    FPLLL_ABORT("shortestVector: invalid evaluator type");
  }
  evaluator->init_delta_def(prec, rho, true);

  if (!(flags & SVP_OVERRIDE_BND) && (evalMode == EVALMODE_SV || method == SVPM_PROVED))
  {
    Float ftmp1;
    bool result = evaluator->get_max_error_aux(max_dist, true, ftmp1);
    FPLLL_CHECK(result, "shortestVector: cannot compute an initial bound");
    max_dist.add(max_dist, ftmp1, GMP_RNDU);
  }

  // Main loop of the enumeration
  enumerate_svp(d, gso, max_dist, *evaluator, pruning, flags);

  int result = RED_ENUM_FAILURE;
  if (evalMode != EVALMODE_SV)
  {
    result   = RED_SUCCESS;
    solCount = evaluator->solCount * 2;
  }
  else if (!evaluator->sol_coord.empty())
  {
    /*Float fMaxError;
    validMaxError = evaluator->get_max_error(fMaxError);
    maxError = fMaxError.get_d(GMP_RNDU);*/
    for (int i = 0; i < d; i++)
    {
      itmp1.set_f(evaluator->sol_coord[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
    }
    result = RED_SUCCESS;
  }

  if (findsubsols)
  {
    subsol_coord->clear();
    subsolDist->clear();
    subsolDist->resize(evaluator->sub_sol_coord.size());
    for (size_t i = 0; i < evaluator->sub_sol_coord.size(); ++i)
    {
      (*subsolDist)[i] = evaluator->sub_solDist[i];

      IntVect ssC;
      for (size_t j = 0; j < evaluator->sub_sol_coord[i].size(); ++j)
      {
        itmp1.set_f(evaluator->sub_sol_coord[i][j]);
        ssC.emplace_back(itmp1);
      }
      subsol_coord->emplace_back(std::move(ssC));
    }
  }

  delete evaluator;
  Float::set_prec(oldprec);
  return result;
}

int shortest_vector(IntMatrix &b, IntVect &sol_coord, SVPMethod method, int flags)
{
  long long tmp;
  return shortest_vector_ex(b, sol_coord, method, vector<double>(), flags, EVALMODE_SV, tmp);
}

int shortest_vector_pruning(IntMatrix &b, IntVect &sol_coord, const vector<double> &pruning,
                            int flags)
{
  long long tmp;
  return shortest_vector_ex(b, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp);
}

int shortest_vector_pruning(IntMatrix &b, IntVect &sol_coord, vector<IntVect> &subsol_coord,
                            vector<enumf> &subsolDist, const vector<double> &pruning, int flags)
{
  long long tmp;
  return shortest_vector_ex(b, sol_coord, SVPM_FAST, pruning, flags, EVALMODE_SV, tmp,
                            &subsol_coord, &subsolDist);
}

/* Closest vector problem
   ====================== */

static void get_gscoords(const Matrix<Float> &matrix, const Matrix<Float> &mu,
                         const Matrix<Float> &r, const FloatVect &v, FloatVect &vcoord)
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

static void babai(const FloatMatrix &matrix, const Matrix<Float> &mu, const Matrix<Float> &r,
                  const FloatVect &target, FloatVect &target_coord)
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

int closest_vector(IntMatrix &b, const IntVect &intTarget, IntVect &sol_coord, int flags)
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
  int minPrec = gso_min_prec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec    = max(53, minPrec + 10);
  int oldprec = Float::set_prec(prec);

  // Allocates space for vectors and matrices in constructors
  IntMatrix empty_mat;
  MatGSO<Integer, Float> gso(b, empty_mat, empty_mat, GSO_INT_GRAM);
  FloatVect target_coord;
  Float max_dist;
  Integer itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.update_gso();
  gen_zero_vect(sol_coord, d);

  /* Applies Babai's algorithm. Because we use fp, it might be necessary to
      do it several times (if ||target|| >> ||b_i||) */
  FloatMatrix floatMatrix(d, n);
  FloatVect target(n), babaiSol;
  IntVect intNewTarget = intTarget;

  for (int i = 0; i < d; i++)
    for (int j = 0; j < n; j++)
      floatMatrix(i, j).set_z(b(i, j));

  for (int loopIdx = 0;; loopIdx++)
  {
    if (loopIdx >= 0x100 && ((loopIdx & (loopIdx - 1)) == 0))
      FPLLL_INFO("warning: possible infinite loop in Babai's algorithm");

    for (int i = 0; i < n; i++)
    {
      target[i].set_z(intNewTarget[i]);
    }
    babai(floatMatrix, gso.get_mu_matrix(), gso.get_r_matrix(), target, babaiSol);
    int idx;
    for (idx = 0; idx < d && babaiSol[idx] >= -1 && babaiSol[idx] <= 1; idx++)
    {
    }
    if (idx == d)
      break;

    for (int i = 0; i < d; i++)
    {
      itmp1.set_f(babaiSol[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
      for (int j = 0; j < n; j++)
        intNewTarget[j].submul(itmp1, b(i, j));
    }
  }
  // FPLLL_TRACE("BabaiSol=" << sol_coord);
  get_gscoords(floatMatrix, gso.get_mu_matrix(), gso.get_r_matrix(), target, target_coord);

  /* Computes a very large bound to make the algorithm work
      until the first solution is found */
  max_dist = 0.0;
  for (int i = 1; i < d; i++)
  {
    // get_r_exp(i, i) = r(i, i) because gso is initialized without GSO_ROW_EXPO
    max_dist.add(max_dist, gso.get_r_exp(i, i));
  }

  FastEvaluator<Float> evaluator(n, gso.get_mu_matrix(), gso.get_r_matrix(), EVALMODE_CV);

  // Main loop of the enumeration
  Enumeration<Float> enumobj(gso, evaluator);
  enumobj.enumerate(0, d, max_dist, 0, target_coord, vector<enumxt>(), vector<enumf>());

  int result = RED_ENUM_FAILURE;
  if (!evaluator.sol_coord.empty())
  {
    // FPLLL_TRACE("evaluator.sol_coord=" << evaluator.sol_coord);
    if (flags & CVP_VERBOSE)
      FPLLL_INFO("max_dist=" << max_dist);
    for (int i = 0; i < d; i++)
    {
      itmp1.set_f(evaluator.sol_coord[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
    }
    result = RED_SUCCESS;
  }

  Float::set_prec(oldprec);
  return result;
}

FPLLL_END_NAMESPACE
