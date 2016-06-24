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
static int lastUsefulIndex(const Matrix<Float>& r) {
  int i;
  Float rdiagMinValue;
  rdiagMinValue.mul_2si(r(0, 0), 1);
  for (i = r.getRows() - 1; i > 0; i--) {
    if (r(i, i) <= rdiagMinValue) break;
  }
  return i + 1;
}

/* Finds the shortest vector of the basis b and returns its squared norm in
   basisMin */
static void getBasisMin(Integer& basisMin, const IntMatrix& b,
                        int first, int last) {
  Integer sqNorm;
  int n = b.getCols();
  sqrNorm(basisMin, b[first], n);

  for (int i = first + 1; i < last; i++) {
    sqrNorm(sqNorm, b[i], n);
    if (sqNorm < basisMin)
      basisMin = sqNorm;
  }
}

static bool enumerateSVP(int d, MatGSO<Integer, Float>& gso, Float& maxDist,
        Evaluator<Float>& evaluator, const vector<enumf>& pruning,
        int flags) 
{
  Enumeration<Float> enumobj(gso, evaluator);
  bool dual = (flags & SVP_DUAL);
  if (d == 1 || !pruning.empty() || dual) {
    enumobj.enumerate(0, d, maxDist, 0, vector<Float>(), vector<enumxt>(), pruning, dual);
  }
  else {
    Enumerator enumerator(d, gso.getMuMatrix(), gso.getRMatrix());
    while (enumerator.enumNext(maxDist)) {
      if (flags & SVP_VERBOSE) {
        evaluator.newSolFlag = false;
        cerr << enumerator.getSubTree();
        if (evaluator.evalMode != EVALMODE_SV)
          cerr << " (count=2*" << evaluator.solCount << ")";
      }

      /* Enumerates short vectors only in enumerator.getSubTree()
        (about maxVolume iterations or less) */
      enumobj.enumerate(0, d, maxDist, 0, FloatVect(), enumerator.getSubTree(), pruning);

      if (flags & SVP_VERBOSE) {
        cerr << "\r" << (char) 27 << "[K";
        if (evaluator.evalMode == EVALMODE_SV && evaluator.newSolFlag)
          cerr << "Solution norm^2=" << evaluator.lastPartialDist
                << " value=" << evaluator.solCoord << endl;
      }
    }
  }
  return !evaluator.solCoord.empty();
}

static int shortestVectorEx(IntMatrix& b, IntVect& solCoord,
        SVPMethod method, const vector<double>& pruning, int flags,
        EvaluatorMode evalMode,
        long long& solCount, 
        vector<IntVect>* subsolCoord = nullptr, vector<enumf>* subsolDist = nullptr) 
{
  bool findsubsols = (subsolCoord != nullptr) && (subsolDist != nullptr);
  
  // d = lattice dimension (note that it might decrease during preprocessing)
  int d = b.getRows();
  // n = dimension of the space
  int n = b.getCols();

  FPLLL_CHECK(d > 0 && n > 0, "shortestVector: empty matrix");
  FPLLL_CHECK(d <= n, "shortestVector: number of vectors > size of the vectors");

  // Sets the floating-point precision
  // Error bounds on GSO are valid if prec >= minprec
  double rho;
  int minPrec = gsoMinPrec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec = max(53, minPrec + 10);
  int oldprec = Float::setprec(prec);

  // Allocates space for vectors and matrices in constructors
  IntMatrix emptyMat;
  MatGSO<Integer, Float> gso(b, emptyMat, emptyMat, GSO_INT_GRAM);
  Float maxDist;
  Integer intMaxDist;
  Integer itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.updateGSO();
  genZeroVect(solCoord, d);

  // If the last b_i* are too large, removes them to avoid an underflow
  int newD = lastUsefulIndex(gso.getRMatrix());
  if (newD < d) {
    //FPLLL_TRACE("Ignoring the last " << d - newD << " vector(s)");
    d = newD;
  }
  
  if (flags & SVP_DUAL) {
    maxDist = gso.getRExp(d - 1, d - 1);
    Float one; one = 1.0;
    maxDist.div(one, maxDist);
    if (flags & SVP_VERBOSE) {
      cout << "maxDist = " << maxDist << endl;
    }
  } else {
      /* Computes a bound for the enumeration. This bound would work for an
         exact algorithm, but we will increase it later to ensure that the fp
         algorithm finds a solution */
      getBasisMin(intMaxDist, b, 0, d);
    maxDist.set_z(intMaxDist, GMP_RNDU);
  }

  // Initializes the evaluator of solutions
  Evaluator<Float>* evaluator;
  if (method == SVPM_FAST) {
    evaluator = new FastEvaluator<Float>(d, gso.getMuMatrix(),
            gso.getRMatrix(), evalMode, 0, findsubsols);
  }
  else if (method == SVPM_PROVED) {
    ExactEvaluator* p = new ExactEvaluator(d, b, gso.getMuMatrix(),
            gso.getRMatrix(), evalMode, 0, findsubsols);
    p->intMaxDist = intMaxDist;
    evaluator = p;
  }
  else {
    FPLLL_ABORT("shortestVector: invalid evaluator type");
  }
  evaluator->initDeltaDef(prec, rho, true);

  if (!(flags & SVP_OVERRIDE_BND) && (evalMode == EVALMODE_SV || method == SVPM_PROVED)) {
    Float ftmp1;
    bool result = evaluator->getMaxErrorAux(maxDist, true, ftmp1);
    FPLLL_CHECK(result, "shortestVector: cannot compute an initial bound");
    maxDist.add(maxDist, ftmp1, GMP_RNDU);
  }
  
  // Main loop of the enumeration
  enumerateSVP(d, gso, maxDist, *evaluator, pruning, flags);

  int result = RED_ENUM_FAILURE;
  if (evalMode != EVALMODE_SV) {
    result = RED_SUCCESS;
    solCount = evaluator->solCount * 2;
  }
  else if (!evaluator->solCoord.empty()) {
    /*Float fMaxError;
    validMaxError = evaluator->getMaxError(fMaxError);
    maxError = fMaxError.get_d(GMP_RNDU);*/
    for (int i = 0; i < d; i++) {
      itmp1.set_f(evaluator->solCoord[i]);
      solCoord[i].add(solCoord[i], itmp1);
    }
    result = RED_SUCCESS;
  }
  
  if (findsubsols)
  {
    subsolCoord->clear();
    subsolDist->clear();
    subsolDist->resize(evaluator->sub_solCoord.size());
    for (size_t i = 0; i < evaluator->sub_solCoord.size(); ++i)
    {
      (*subsolDist)[i] = evaluator->sub_solDist[i];

      IntVect ssC;
      for (size_t j = 0; j < evaluator->sub_solCoord[i].size(); ++j)
      {
        itmp1.set_f(evaluator->sub_solCoord[i][j]);
        ssC.emplace_back(itmp1);
      }
      subsolCoord->emplace_back( std::move(ssC) );
    }
  }

  delete evaluator;
  Float::setprec(oldprec);
  return result;
}

int shortestVector(IntMatrix& b, IntVect& solCoord,
                   SVPMethod method, int flags) {
  long long tmp;
  return shortestVectorEx(b, solCoord, method, vector<double>(), flags,
          EVALMODE_SV, tmp);
}

int shortestVectorPruning(IntMatrix& b, IntVect& solCoord,
                   const vector<double>& pruning, int flags) {
  long long tmp;
  return shortestVectorEx(b, solCoord, SVPM_FAST, pruning, flags,
          EVALMODE_SV, tmp);
}

int shortestVectorPruning(IntMatrix& b, IntVect& solCoord, vector<IntVect>& subsolCoord, vector<enumf>& subsolDist,
        const vector<double>& pruning, int flags)
{
  long long tmp;
  return shortestVectorEx(b, solCoord, SVPM_FAST, pruning, flags,
          EVALMODE_SV, tmp, &subsolCoord, &subsolDist);
}

/* Closest vector problem
   ====================== */

static void getGSCoords(const Matrix<Float>& matrix, const Matrix<Float>& mu,
  const Matrix<Float>& r, const FloatVect& v, FloatVect& vcoord) {

  int n = matrix.getRows(), m = matrix.getCols();

  if (static_cast<int>(vcoord.size()) != n) vcoord.resize(n);
  FPLLL_DEBUG_CHECK(mu.getRows() == n && mu.getCols() == n &&
    r.getRows() == n && r.getCols() == n &&
    static_cast<int>(v.size()) == m);

  for (int i = 0; i < n; i++) {
    vcoord[i] = 0.0;
    for (int j = 0; j < m; j++)
      vcoord[i].addmul(v[j], matrix(i, j));
    for (int j = 0; j < i; j++)
      vcoord[i].submul(mu(i, j), vcoord[j]);
  }
  for (int i = 0; i < n; i++) {
    vcoord[i].div(vcoord[i], r(i, i));
  }
}

static void babai(const FloatMatrix& matrix, const Matrix<Float>& mu,
  const Matrix<Float>& r, const FloatVect& target, FloatVect& targetcoord) {

  int d = matrix.getRows();
  getGSCoords(matrix, mu, r, target, targetcoord);
  for (int i = d - 1; i >= 0; i--) {
    targetcoord[i].rnd(targetcoord[i]);
    for (int j = 0; j < i; j++)
      targetcoord[j].submul(mu(i, j), targetcoord[i]);
  }
}

int closestVector(IntMatrix& b, const IntVect& intTarget,
                  IntVect& solCoord, int flags) {
  // d = lattice dimension (note that it might decrease during preprocessing)
  int d = b.getRows();
  // n = dimension of the space
  int n = b.getCols();

  FPLLL_CHECK(d > 0 && n > 0, "closestVector: empty matrix");
  FPLLL_CHECK(d <= n, "closestVector: number of vectors > size of the vectors");

  // Sets the floating-point precision
  // Error bounds on GSO are valid if prec >= minprec
  double rho;
  int minPrec = gsoMinPrec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec = max(53, minPrec + 10);
  int oldprec = Float::setprec(prec);

  // Allocates space for vectors and matrices in constructors
  IntMatrix emptyMat;
  MatGSO<Integer, Float> gso(b, emptyMat, emptyMat, GSO_INT_GRAM);
  FloatVect targetCoord;
  Float maxDist;
  Integer itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.updateGSO();
  genZeroVect(solCoord, d);

  /* Applies Babai's algorithm. Because we use fp, it might be necessary to
      do it several times (if ||target|| >> ||b_i||) */
  FloatMatrix floatMatrix(d, n);
  FloatVect target(n), babaiSol;
  IntVect intNewTarget = intTarget;

  for (int i = 0; i < d; i++)
    for (int j = 0; j < n; j++)
      floatMatrix(i, j).set_z(b(i, j));

  for (int loopIdx = 0;; loopIdx++) {
    if (loopIdx >= 0x100 && ((loopIdx & (loopIdx - 1)) == 0))
      FPLLL_INFO("warning: possible infinite loop in Babai's algorithm");

    for (int i = 0; i < n; i++) {
      target[i].set_z(intNewTarget[i]);
    }
    babai(floatMatrix, gso.getMuMatrix(), gso.getRMatrix(), target, babaiSol);
    int idx;
    for (idx = 0; idx < d && babaiSol[idx] >= -1 && babaiSol[idx] <= 1;
          idx++) {}
    if (idx == d) break;

    for (int i = 0; i < d; i++) {
      itmp1.set_f(babaiSol[i]);
      solCoord[i].add(solCoord[i], itmp1);
      for (int j = 0; j < n; j++)
        intNewTarget[j].submul(itmp1, b(i, j));
    }
  }
  //FPLLL_TRACE("BabaiSol=" << solCoord);
  getGSCoords(floatMatrix, gso.getMuMatrix(), gso.getRMatrix(), target,
          targetCoord);

  /* Computes a very large bound to make the algorithm work
      until the first solution is found */
  maxDist = 0;
  for (int i = 1; i < d; i++) {
    // getRExp(i, i) = r(i, i) because gso is initialized without GSO_ROW_EXPO
    maxDist.add(maxDist, gso.getRExp(i, i));
  }

  FastEvaluator<Float> evaluator(n, gso.getMuMatrix(), gso.getRMatrix(),
          EVALMODE_CV);

  // Main loop of the enumeration
  Enumeration<Float> enumobj(gso, evaluator);
  enumobj.enumerate(0, d, maxDist, 0, targetCoord, vector<enumxt>(), vector<enumf>());

  int result = RED_ENUM_FAILURE;
  if (!evaluator.solCoord.empty()) {
    //FPLLL_TRACE("evaluator.solCoord=" << evaluator.solCoord);
    if (flags & CVP_VERBOSE)
      FPLLL_INFO("maxDist=" << maxDist);
    for (int i = 0; i < d; i++) {
      itmp1.set_f(evaluator.solCoord[i]);
      solCoord[i].add(solCoord[i], itmp1);
    }
    result = RED_SUCCESS;
  }

  Float::setprec(oldprec);
  return result;
}

FPLLL_END_NAMESPACE
