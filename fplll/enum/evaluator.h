/* Copyright (C) 2008-2011 Xavier Pujol.

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

#ifndef FPLLL_EVALUATOR_H
#define FPLLL_EVALUATOR_H

#include "../util.h"

FPLLL_BEGIN_NAMESPACE

enum EvaluatorMode {
  EVALMODE_SV = 0,
  EVALMODE_CV = 0,
  EVALMODE_COUNT = 1,
  EVALMODE_PRINT = 2
};

/**
 * Evaluator stores the best solution found by enumerate. The Float
 * specialization provides additional information about the solution accuracy.
 */
template<class FT>
class Evaluator {
public:
  Evaluator() : newSolFlag(false) {}
  virtual ~Evaluator() {}

  /** Called by enumerate when a solution is found.
     Input: newSolCoord = coordinates of the solution in Gram-Schmidt basis
     newPartialDist = estimated distance between the solution and the
     orthogonal projection of target on the lattice space
     maxDist = current bound of the algorithm
     Output: maxDist can be decreased */
  virtual void evalSol(const vector<FT>& newSolCoord,
          const enumf& newPartialDist, enumf& maxDist, long normExp) = 0;

  /** Coordinates of the solution in the lattice */
  vector<FT> solCoord;

  /** Set to true when solCoord is updated */
  bool newSolFlag;
};

/**
 * Simple solution evaluator which provides a result without error bound.
 * The same instance can be used for several calls to enumerate on different
 * problems.
 */
template<class FT>
class FastEvaluator : public Evaluator<FT> {
public:
  using Evaluator<FT>::solCoord;
  using Evaluator<FT>::newSolFlag;

  FastEvaluator() : Evaluator<FT>() {}
  virtual ~FastEvaluator() {}

  /**
   * Called by enumerate when a solution is found.
   * FastEvaluator always accepts the solution and sets the bound maxDist to
   * newPartialDist.
   *
   * @param newSolCoord    Coordinates of the solution in the lattice
   * @param newPartialDist Floating-point evaluation of the norm of the solution
   * @param maxDist        Bound of the enumeration (updated by the function)
   * @param normExp        r(i, i) is divided by 2^normExp in enumerate before
   *                       being converted to double
   */
  virtual void evalSol(const vector<FT>& newSolCoord,
                       const enumf& newPartialDist, enumf& maxDist,
                       long normExp) {
    solCoord = newSolCoord;
    maxDist = newPartialDist;
    newSolFlag = true;
  }
};

/**
 * Evaluator stores the best solution found by enumerate and provides
 * information about the accuracy of this solution.
 */
template<>
class Evaluator<Float> {
public:
  Evaluator<Float>(int d, const Matrix<Float>& mu, const Matrix<Float>& r,
                   int evalMode) :
    newSolFlag(false), evalMode(evalMode), inputErrorDefined(false),
    d(d), mu(mu), r(r)
  {
    maxDRdiag.resize(d);
    maxDMu.resize(d);
  }

  virtual ~Evaluator<Float>() {};

  void initDeltaDef(int prec, double rho, bool withRoundingToEnumf);

  /**
   * Computes maxError such that
   * normOfSolution^2 <= (1 + maxError) * lambda_1(L)^2.
   * The default implementation might fail (i.e. return false).
   */
  virtual bool getMaxError(Float& maxError) = 0;

  /**
   * Called by enumerate when a solution is found.
   * The default implementation always accepts the solution and sets the bound
   * maxDist to newPartialDist.
   *
   * @param newSolCoord    Coordinates of the solution
   * @param newPartialDist Floating-point estimation of the norm of the solution
   * @param maxDist        Bound of the enumeration (updated by the function)
   * @param normExp        It is assumed that r(i, i) is divided by 2^normExp
   *                       in enumerate
   */
  virtual void evalSol(const FloatVect& newSolCoord,
          const enumf& newPartialDist, enumf& maxDist, long normExp) = 0;

  // Internal use
  bool getMaxErrorAux(const Float& maxDist, bool boundOnExactVal, Float& maxDE);

  /** Coordinates of the solution in the lattice */
  FloatVect solCoord;

  /** Set to true when solCoord is updated */
  bool newSolFlag;
  /** Incremented when solCoord is updated */
  long long solCount;
  int evalMode;

  /* To enable error estimation, the caller must set
     inputErrorDefined=true and fill maxDRdiag and maxDMu */
  bool inputErrorDefined;
  FloatVect maxDRdiag, maxDMu;  // Error bounds on input parameters
  Float lastPartialDist;        // Approx. squared norm of the last solution

  int d;
  const Matrix<Float>& mu;
  const Matrix<Float>& r;
};

/**
 * Simple solution evaluator which provides a non-certified result, but can
 * give an error bound.
 * The same object can be used for several calls to enumerate on different
 * instances.
 */
template<>
class FastEvaluator<Float> : public Evaluator<Float> {
public:
  FastEvaluator(int d = 0, const Matrix<Float>& mu = Matrix<Float>(),
                const Matrix<Float>& r = Matrix<Float>(),
                int evalMode = EVALMODE_SV) :
    Evaluator<Float>(d, mu, r, evalMode) {}
  virtual ~FastEvaluator() {}

  virtual bool getMaxError(Float& maxError);
  virtual void evalSol(const FloatVect& newSolCoord,
          const enumf& newPartialDist, enumf& maxDist, long normExp);
};

/**
 * ExactEvaluator stores the best solution found by enumerate.
 * The result is guaranteed, but the the evaluation of new solutions is longer.
 */
class ExactEvaluator : public Evaluator<Float> {
public:
  ExactEvaluator(int d, const IntMatrix& matrix, const Matrix<Float>& mu,
                 const Matrix<Float>& r, int evalMode) :
    Evaluator<Float>(d, mu, r, evalMode), matrix(matrix)
  {
    intMaxDist = -1;
  }

  /**
   * Sets maxError to 0: the result is guaranteed.
   */
  virtual bool getMaxError(Float& maxError);

  virtual void evalSol(const FloatVect& newSolCoord,
          const enumf& newPartialDist, enumf& maxDist, long normExp);

  Integer intMaxDist;       // Exact norm of the last vector

private:
  void updateMaxDist(enumf& maxDist, long normExp);

  const IntMatrix& matrix;  // matrix of the lattice
};

FPLLL_END_NAMESPACE

#endif
