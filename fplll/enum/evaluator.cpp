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

#include "evaluator.h"

FPLLL_BEGIN_NAMESPACE

void Evaluator<Float>::initDeltaDef(int prec, double rho, bool withRoundingToEnumf) {
  /* Computes error bounds on GSO
      For all 0 <= i < d and 0 <= j <= i we have:
      |r~_i - r_i| / r_i <= maxRelDR = d * rho ^ (i + 1) * 2 ^ (2 - prec)
      |mu~_(i,j) - mu_(i,j)| <= d * rho ^ (i + 1) * 2 ^ (4 - prec)
      The following formula is used to bound absolute error on r_i:
      |r~_i - r_i| <= r_i * maxRelDR <= (r~_i + |r~_i - r_i| * maxRelDR)
      ==> |r~_i - r_i| <= r~_i * maxRelDR / (1 - maxRelDR)  */
  Float rhoPow, maxRelDR, ftmp1;
  inputErrorDefined = true;
  for (int j = 0; j < d; j++) {
    rhoPow = rho;
    rhoPow.pow_si(rhoPow, j + 1, GMP_RNDU);      // >= rho ^ (j + 1)
    ftmp1 = d;
    ftmp1.mul_2si(ftmp1, 2 - prec);              // = d * 2 ^ (2 - prec)
    maxRelDR.mul(ftmp1, rhoPow, GMP_RNDU);       // >= |r~_j - r_j| / r_j
    ftmp1 = 1.0;
    ftmp1.sub(ftmp1, maxRelDR, GMP_RNDD);        // <= 1 - maxRelDR
    ftmp1.div(maxRelDR, ftmp1, GMP_RNDU);        // >= maxRelDR/(1-maxRelDR)
    ftmp1.mul(ftmp1, r(j, j));                   // >= |r~_j - r_j|
    maxDRdiag[j] = ftmp1;

    ftmp1 = d;
    ftmp1.mul_2si(ftmp1, 4 - prec);              // = d * 2 ^ (4 - prec)
    ftmp1.mul(ftmp1, rhoPow, GMP_RNDU);          // >= |mu~_(?,j) - mu_(?,j)|
    maxDMu[j] = ftmp1;
  }
  if (withRoundingToEnumf) {
    Float halfULP;
    halfULP = numeric_limits<double>::epsilon() * 0.5;
    for (int j = 0; j < d; j++) {
      maxDRdiag[j].addmul(r(j, j), halfULP, GMP_RNDU); // >= |double(r~_j) - r_j|
      maxDMu[j].add(maxDMu[j], halfULP, GMP_RNDU);
    }
  }
}

/* Main function for error evaluation
   Returns maxDE such that in each loop of the enumeration, if
     exact_squared_norm(vector) <= maxDist (when boundOnExactVal = true)
   or
     computed_squared_norm(vector) <= maxDist (when boundOnExactVal = false)
   then
     |exact_squared_norm(vector) - computed_squared_norm(vector)| <= maxDE
   (computed_squared_norm(vector) is the variable 'newDist' is enumerateLoop)

   In the following comments:
   - The symbols +~ and *~ represent rounded fp addition and multiplication
   - a,b,c,... represent exact values
   - a~,b~,c~,... are approx. values used or computed by the fp algorithm */

bool Evaluator<Float>::getMaxErrorAux(const Float& maxDist,
  bool boundOnExactVal, Float& maxDE) {

  FPLLL_CHECK(inputErrorDefined,
    "Evaluator: error evaluation failed because the input error is undefined");

  Float ulp, halfULP, K, tmp1, tmp2;
  Float rdiagTilde, minRDiag, maxRDiag, muTilde, maxMu, maxMuTildeX;
  Float maxC, maxCTilde, maxY, maxYTilde, maxY2Tilde, maxRY2Tilde;
  Float maxDC, maxDY, maxDY2, maxDRY2;
  FloatVect maxX(d);

  ulp = numeric_limits<double>::epsilon();
  halfULP.mul_2si(ulp, -1);
  K = 1.0;
  K.add(K, halfULP, GMP_RNDU);
  maxDE = 0.0;

  for (int i = d - 1; i >= 0; i--) {
    maxC = 0.0;
    maxCTilde = 0.0;
    maxDC = 0.0;

    long rdiagExp = r(i, i).exponent();
    tmp1.mul_2si(r(i, i), -rdiagExp);
    tmp1 = tmp1.get_d();
    rdiagTilde.mul_2si(tmp1, rdiagExp);             // = r~_i

    /* Computes bounds on:
       C = mu(d-1,i) * x_(d-1) + ... + mu(j,i) * x_j
       C~ = mu~(d-1,i) *~ x_(d-1) +~ ... +~ mu~(j,i) *~ x_j
       DC = |C - C~|  */
    for (int j = d - 1; j > i; j--) {
      maxMu.abs(mu(j, i));
      maxMu.add(maxMu, maxDMu[i], GMP_RNDU);        // >= |mu(j,i)|
      maxC.addmul(maxMu, maxX[j], GMP_RNDU);
      // now maxC >= |mu(d-1,i)| * |x_(d-1)| + ... + |mu(j,i)| * |x_j|
      muTilde = fabs(mu(j, i).get_d());             // = |mu~(j,i)|
      maxMuTildeX.mul(muTilde, maxX[j], GMP_RNDU);  // >= mu~(j,i) * x_j
      maxDC.addmul(maxDMu[i], maxX[j]);             // err1: initial error on mu
      maxDC.addmul(maxMuTildeX, halfULP, GMP_RNDU); // err2: rounding after *
      maxMuTildeX.mul(maxMuTildeX, K, GMP_RNDU);    // >= mu~(j,i) *~ x_j
      maxCTilde.addmul(maxMuTildeX, K, GMP_RNDU);
      // now maxCTilde >= |mu(d-1,i)| *~ |x_(d-1)| +~ ... + |mu(j,i)| *~ |x_j|
      maxDC.addmul(maxCTilde, halfULP, GMP_RNDU);   // err3: rounding after +
      maxCTilde.mul(maxCTilde, K, GMP_RNDU);
      // now maxCTilde >= |mu(d-1,i)| *~ |x_(d-1)| +~ ... +~ |mu(j,i)| *~ |x_j|
    }

    if (boundOnExactVal) {
      // We have dist <= maxDist
      minRDiag.sub(r(i, i), maxDRdiag[i], GMP_RNDD); // <= r_i
      if (minRDiag <= 0.0) return false;
      tmp1.div(maxDist, minRDiag, GMP_RNDU);        // >= dist / r_i
      maxY.sqrt(tmp1, GMP_RNDU);                    // >= |y_i|      
      /* DY = |y_i - y~_i|
            <= DC + |x_i - C~| * halfULP
            <= DC + (|x_i - C| + DC) * halfULP
            = DC * K + y * halfULP  */
      maxDY.mul(maxY, halfULP, GMP_RNDU);
      maxDY.addmul(maxDC, K, GMP_RNDU);             // >= |y_i - y~_i|
      maxYTilde.add(maxY, maxDY, GMP_RNDU);         // >= |y~_i|
      tmp1.add(maxY, maxC, GMP_RNDD);               // >= |x_i|
      maxX[i].floor(tmp1);                          // x_i is an integer
      tmp1 = maxY;
    }
    else {
      // We have dist~ <= maxDist
      /* dist~ >= (y~_i *~ y~_i) *~ r~_i
               >= (y~_i *~ y~_i) * r~_i - dist~ * halfULP
         ==> y~_i *~ y~_i <= dist~ * K / r~_i  */
      tmp1.mul(maxDist, K, GMP_RNDU);
      tmp1.div(tmp1, rdiagTilde, GMP_RNDU);         // >= y~_i *~ y~_i
      /* tmp1 >= y~_i *~ y~_i >= y~_i * y~_i - tmp1 * halfULP
         ==> y~_i <= sqrt(tmp1 * K)   */
      tmp1.mul(tmp1, K, GMP_RNDU);
      maxYTilde.sqrt(tmp1, GMP_RNDU);               // >= y~_i
      maxDY.mul(maxYTilde, halfULP, GMP_RNDU);      // err2: rounding after +
      maxDY.add(maxDY, maxDC, GMP_RNDU);            // err1: initial error on C
      // now maxDY >= |y_i - y~_i|
      /* maxYTilde >= |C~_i -~ x_i| >= |C~_i - x_i| - maxYTilde * halfULP
         ==> |x_i| <= C~_i + maxYTilde * K  */
      tmp1 = maxCTilde;
      tmp1.addmul(maxYTilde, K, GMP_RNDD);          // >= |x_i|
      maxX[i].floor(tmp1);                          // x_i is an integer
      tmp1 = maxYTilde;
    }

    maxDY2.mul(maxDY, tmp1);
    maxDY2.mul_2si(maxDY2, 1);
    maxDY2.addmul(maxDY, maxDY, GMP_RNDU);
    /* now, we have maxDY2 <= |y~_i * y~_i - y_i * y_i|
       Case 1: tmp1 = maxY
         |y~_i * y~_i - y_i * y_i| <= (maxY + maxDY) ^ 2 - maxY ^ 2
       Case 2: tmp1 = maxY~
         |y~_i * y~_i - y_i * y_i| <= (maxY~ + maxDY) ^ 2 - maxY~ ^ 2  */
    maxY2Tilde.mul(maxYTilde, maxYTilde, GMP_RNDU); // >= y~_i * y~_i
    maxDY2.addmul(maxY2Tilde, halfULP, GMP_RNDU);
    /* now, we have maxDY2 <= |y~_i *~ y~_i - y_i * y_i| */
    maxY2Tilde.mul(maxY2Tilde, K, GMP_RNDU);        // >= y~_i *~ y~_i
    maxRDiag.add(r(i, i), maxDRdiag[i], GMP_RNDU); // >= r_i
    maxRY2Tilde.mul(rdiagTilde, maxY2Tilde, GMP_RNDU); // >= y~_i *~ y~_i * r~_i

    maxDRY2.mul(maxRDiag, maxDY2, GMP_RNDU);
    maxDRY2.addmul(maxY2Tilde, maxDRdiag[i], GMP_RNDU);
    maxDRY2.addmul(maxRY2Tilde, halfULP, GMP_RNDU);
    /* now maxDRY2 >= |y~_i *~ y~_i - y_i * y_i| * r_i +
                      |y~_i *~ y~_i| * |r~_i - r_i| +
                      y~_i *~ y~_i * r~_i * halfULP
                   >= |y~_i *~ y~_i * r~_i - y_i * y_i * r_i| +
                      y~_i *~ y~_i * r~_i * halfULP  */

    /* Let u = dist_(i+1), u~ = dist~_(i+1), v = y_i * y_i * r_i and
       v~ = y~_i *~ y~_i *~ r~_i. We have:
         |dist_i - dist~_i| = |(u~ +~ v~) - (u + v)|
                            <= |u~ +~ v~) - (u~ + v~)| + |(u~ + v~) - (u + v)|
       Case 1: dist <= maxDist
         |dist_i - dist~_i| <= (u~ + v~) * halfULP + |(u~ + v~) - (u + v)|
                            <= (u + v) * halfULP + K * |(u~ + v~) - (u + v)|
                            <= maxDist * halfULP + K * |(u~ + v~) - (u + v)|
       Case 2: dist~ <= maxDist
         |dist_i - dist~_i| <= (u~ +~ v~) * halfULP + |(u~ + v~) - (u + v)|
                            <= (u~ +~ v~) * halfULP + K * |(u~ + v~) - (u + v)|
                            <= maxDist * halfULP + K * |(u~ + v~) - (u + v)|  */
    maxDE.add(maxDE, maxDRY2, GMP_RNDU);
    maxDE.mul(maxDE, K, GMP_RNDU);
    maxDE.addmul(maxDist, halfULP, GMP_RNDU);
  }

  return true;
}

void FastEvaluator<Float>::evalSol(const FloatVect& newSolCoord,
        const enumf& newPartialDist, enumf& maxDist) 
{
  // Assumes that the solution is valid
  if (evalMode == EVALMODE_SV) 
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
    lastPartialDist = newPartialDist; // Exact conversion
    lastPartialDist.mul_2si(lastPartialDist, normExp);
  }
  else if (evalMode == EVALMODE_PRINT) 
  {
    cout << newSolCoord << "\n";
  }
  newSolFlag = true;
  solCount++;
}

void FastEvaluator<Float>::evalSubSol(int offset, const FloatVect& newSubSolCoord,
        const enumf& subDist) 
{
  sub_solCoord.resize( std::max(sub_solCoord.size(), std::size_t(offset+1)) );
  sub_solDist.resize( sub_solCoord.size(), -1.0);
  if (sub_solDist[offset] == -1.0 || subDist < sub_solDist[offset])
  {
    sub_solCoord[offset] = newSubSolCoord;
    for (int i = 0; i < offset; ++i)
      sub_solCoord[offset][i] = 0.0;
    sub_solDist[offset] = subDist;
  }
}

bool FastEvaluator<Float>::getMaxError(Float& maxError) {
  Float maxE, maxDE, maxOptDE, minOptE, one;

  if (solCoord.empty() || !inputErrorDefined) return false;
  if (!getMaxErrorAux(lastPartialDist, false, maxDE)) return false;

  // Exact norm of an optimal solution <= Exact norm of the result <= maxE
  maxE.add(lastPartialDist, maxDE, GMP_RNDU);
  // Error on the norm of an optimal solution <= maxOptDE
  if (!getMaxErrorAux(maxE, true, maxOptDE)) return false;
  // Approx. norm of an optimal solution >= minOptE
  minOptE.sub(lastPartialDist, maxOptDE, GMP_RNDD);

  one = 1.0;
  maxError.div(maxE, minOptE, GMP_RNDU);
  maxError.sub(maxError, one, GMP_RNDU);
  return true;
}

bool ExactEvaluator::getMaxError(Float& maxError) {
  maxError = 0.0;
  return true;
}

void ExactEvaluator::evalSol(const FloatVect& newSolCoord,
        const enumf& newPartialDist, enumf& maxDist) 
{
  int n = matrix.get_cols();
  Integer newSolDist, coord;
  IntVect newSol;

  genZeroVect(newSol, n);
  newSolDist = 0;

  // Computes the distance between x and zero
  for (int i = 0; i < d; i++) {
    coord.set_f(newSolCoord[i]);
    for (int j = 0; j < n; j++)
      newSol[j].addmul(coord, matrix(i, j));
  }
  for (int i = 0; i < n; i++) {
    coord = newSol[i];
    newSolDist.addmul(coord, coord);
  }

  if (intMaxDist < 0 || newSolDist <= intMaxDist) {
    if (evalMode == EVALMODE_SV) {
      if (max_aux_sols != 0 && !solCoord.empty())
      {
        aux_solCoord.emplace_front( std::move(solCoord) );
        aux_solDist.emplace_front( solDist );
        aux_solintDist.emplace_front( intMaxDist );
        if (aux_solCoord.size() > max_aux_sols)
        {
          aux_solCoord.pop_back();
          aux_solDist.pop_back();
          aux_solintDist.pop_back();
        }
      }
      // Updates the stored solution
      lastPartialDist = newPartialDist;
      lastPartialDist.mul_2si(lastPartialDist, normExp);
      solCoord = newSolCoord;
      intMaxDist = newSolDist;
      updateMaxDist(maxDist);
      solDist = maxDist;
    }
    else if (evalMode == EVALMODE_PRINT) {
      cout << newSolCoord << "\n";
    }
    newSolFlag = true;
    solCount++;
  }
}

void ExactEvaluator::evalSubSol(int offset, const FloatVect& newSubSolCoord,
        const enumf& subDist) 
{
  sub_solCoord.resize( std::max(sub_solCoord.size(), std::size_t(offset+1)) );
  sub_solDist.resize( sub_solCoord.size(), -1.0 );
  Integer minusone;
  minusone = -1;
  sub_solintDist.resize( sub_solCoord.size(), minusone);

  int n = matrix.get_cols();
  Integer newSolDist, coord;
  IntVect newSol;

  genZeroVect(newSol, n);
  newSolDist = 0;

  // Computes the distance between x[[offset,d)] and zero
  for (int i = offset; i < d; i++) 
  {
    coord.set_f(newSubSolCoord[i]);
    for (int j = 0; j < n; j++)
      newSol[j].addmul(coord, matrix(i, j));
  }
  for (int i = 0; i < n; i++) 
  {
    coord = newSol[i];
    newSolDist.addmul(coord, coord);
  }

  if (sub_solintDist[offset] < 0 || newSolDist <= sub_solintDist[offset]) 
  {
    sub_solCoord[offset] = newSubSolCoord;
    for (int i = 0; i < offset; ++i)
      sub_solCoord[offset][i] = 0.0;
    sub_solDist[offset] = subDist;
    sub_solintDist[offset] = newSolDist;
  }
}



// Decreases the bound of the algorithm when a solution is found
void ExactEvaluator::updateMaxDist(enumf& maxDist) {
  Float fMaxDist, maxDE;
  fMaxDist.set_z(intMaxDist, GMP_RNDU);
  bool result = getMaxErrorAux(fMaxDist, true, maxDE);
  FPLLL_CHECK(result, "ExactEvaluator: error cannot be bounded");
  FPLLL_CHECK(maxDE <= r(0, 0), "ExactEvaluator: max error is too large");
  fMaxDist.add(fMaxDist, maxDE);
  fMaxDist.mul_2si(fMaxDist, -normExp);
  maxDist = fMaxDist.get_d();
}

FPLLL_END_NAMESPACE
