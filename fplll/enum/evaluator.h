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

#include <deque>
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
  Evaluator(size_t max_aux_solutions = 0, bool findsubsolutions = false) : max_aux_sols(max_aux_solutions), findsubsols(findsubsolutions), newSolFlag(false) {}
  virtual ~Evaluator() {}

  /** Called by enumerate when a solution is found.
     Input: newSolCoord = coordinates of the solution in Gram-Schmidt basis
     newPartialDist = estimated distance between the solution and the
     orthogonal projection of target on the lattice space
     maxDist = current bound of the algorithm
     Output: maxDist can be decreased */
  virtual void evalSol(const vector<FT>& newSolCoord,
          const enumf& newPartialDist, enumf& maxDist) = 0;
          
  virtual void evalSubSol(int offset, const vector<FT>& newSubSolCoord,
          const enumf& subDist) = 0;
          
  virtual void set_normexp(long normExp) {}

  /** Coordinates of the solution in the lattice */
  vector<FT> solCoord;
  enumf solDist;
  
  /** Other solutions found in the lattice */
  size_t max_aux_sols;
  std::deque< vector<FT> > aux_solCoord;
  std::deque< enumf > aux_solDist;
  
  /** Subsolutions found in the lattice */
  bool findsubsols;
  vector< vector<FT> > sub_solCoord;
  vector< enumf > sub_solDist;

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
  using Evaluator<FT>::solDist;
  using Evaluator<FT>::newSolFlag;
  using Evaluator<FT>::aux_solCoord;
  using Evaluator<FT>::aux_solDist;
  using Evaluator<FT>::sub_solCoord;
  using Evaluator<FT>::sub_solDist;
  using Evaluator<FT>::max_aux_sols;

  FastEvaluator(size_t max_aux_solutions = 0, bool findsubsolutions = false) 
    : Evaluator<FT>(max_aux_solutions, findsubsolutions) 
  {
  }
  
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
                       const enumf& newPartialDist, enumf& maxDist) 
  {
    if (max_aux_sols != 0 && !solCoord.empty())
    {
      aux_solCoord.emplace_front( std::move(solCoord) );
      aux_solDist.emplace_front( solDist );
      if (aux_solCoord.size() > max_aux_sols)
      {
        aux_solCoord.pop_back();
        aux_solDist.pop_back();
      }
    }
    solCoord = newSolCoord;
    maxDist = solDist = newPartialDist;
    newSolFlag = true;
  }
  
  virtual void evalSubSol(int offset, const vector<FT>& newSubSolCoord, const enumf& subDist)
  {
    sub_solCoord.resize( std::max(sub_solCoord.size(), std::size_t(offset+1)) );
    sub_solDist.resize( sub_solCoord.size(), -1.0 );
    if (sub_solDist[offset] == -1.0 || subDist < sub_solDist[offset])
    {
      sub_solCoord[offset] = newSubSolCoord;
      for (int i = 0; i < offset; ++i)
        sub_solCoord[offset][i] = 0.0;
      sub_solDist[offset] = subDist;
    }
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
                   int evalMode, size_t max_aux_solutions = 0, bool findsubsolutions = false) :
    max_aux_sols(max_aux_solutions), findsubsols(findsubsolutions), newSolFlag(false), evalMode(evalMode), inputErrorDefined(false),
    d(d), mu(mu), r(r)
  {
    maxDRdiag.resize(d);
    maxDMu.resize(d);
  }

  virtual ~Evaluator<Float>() {}

  virtual void set_normexp(long norm_exp) 
  {
    normExp = norm_exp;
  }
  long normExp;

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
          const enumf& newPartialDist, enumf& maxDist) = 0;
  virtual void evalSubSol(int offset, const FloatVect& newSubSolCoord,
          const enumf& subDist) = 0;

  // Internal use
  bool getMaxErrorAux(const Float& maxDist, bool boundOnExactVal, Float& maxDE);

  /** Coordinates of the solution in the lattice */
  FloatVect solCoord;
  enumf solDist;

  /** Other solutions found in the lattice */
  size_t max_aux_sols;
  std::deque< FloatVect > aux_solCoord;
  std::deque< enumf > aux_solDist;
  
  /** Subsolutions found in the lattice */
  bool findsubsols;
  vector< FloatVect > sub_solCoord;
  vector< enumf > sub_solDist;

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
                int evalMode = EVALMODE_SV, size_t max_aux_solutions = 0, bool findsubsolutions = false) :
    Evaluator<Float>(d, mu, r, evalMode, max_aux_solutions, findsubsolutions) {}
  virtual ~FastEvaluator() {}

  virtual bool getMaxError(Float& maxError);
  virtual void evalSol(const FloatVect& newSolCoord,
          const enumf& newPartialDist, enumf& maxDist);
  virtual void evalSubSol(int offset, const FloatVect& newSubSolCoord,
          const enumf& subDist);
};

/**
 * ExactEvaluator stores the best solution found by enumerate.
 * The result is guaranteed, but the the evaluation of new solutions is longer.
 */
class ExactEvaluator : public Evaluator<Float> {
public:
  ExactEvaluator(int d, const IntMatrix& matrix, const Matrix<Float>& mu,
                 const Matrix<Float>& r, int evalMode, size_t max_aux_solutions = 0, bool findsubsolutions = false) :
    Evaluator<Float>(d, mu, r, evalMode, max_aux_solutions, findsubsolutions), matrix(matrix)
  {
    intMaxDist = -1;
  }

  /**
   * Sets maxError to 0: the result is guaranteed.
   */
  virtual bool getMaxError(Float& maxError);

  virtual void evalSol(const FloatVect& newSolCoord,
          const enumf& newPartialDist, enumf& maxDist);

  virtual void evalSubSol(int offset, const FloatVect& newSubSolCoord,
          const enumf& subDist);

  Integer intMaxDist;       // Exact norm of the last vector

  std::deque< Integer > aux_solintDist; // Exact norm of aux vectors
  vector< Integer> sub_solintDist; // Exact norm of sub vectors

private:
  void updateMaxDist(enumf& maxDist);

  const IntMatrix& matrix;  // matrix of the lattice
};

FPLLL_END_NAMESPACE

#endif
