/* Copyright (C) 2011 Xavier Pujol.

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

#include "util.h"

#ifdef DEBUG
int debugDepth = 0;
#endif

FPLLL_BEGIN_NAMESPACE

enum MinPrecAlgo {
  MINPREC_GSO,
  MINPREC_L2
};

/* State of LDConvHelper (declared in nr.h, must be defined in exactly one
   source file) */
#ifdef FPLLL_WITH_LONG_DOUBLE
mpfr_t LDConvHelper::temp;
bool LDConvHelper::tempInitialized = false;
#endif

/* State of the random generator (declared in nr.h, must be defined in exactly
   one source file) */
bool RandGen::initialized = false;
gmp_randstate_t RandGen::gmpState;

static int computeMinPrec(double& rho, int d, double delta, double eta,
        double epsilon, MinPrecAlgo algo) {
  int oldprec = Float::setprec(53);
  Float fMinprec, fRho, fD, fEta, fDelta, fEpsilon, tmp1, tmp2;

  // These four conversions are exact
  fD = static_cast<double>(d);
  fEta = eta;
  fDelta = delta;
  fEpsilon = epsilon;
  if (algo == MINPREC_L2) {
    // eta - 0.5 is an exact fp operation
    if (fEpsilon > eta - 0.5) fEpsilon = eta - 0.5;
    tmp1 = 1.0;
    tmp1.sub(tmp1, fDelta, GMP_RNDD);
    if (fEpsilon > tmp1) fEpsilon = tmp1;
    // now fEpsilon <= min(epsilon, eta - 0.5, 1 - delta);
  }
  // Computes tmp1 >= (1 + eta) ^ 2 + epsilon
  tmp1 = 1.0;                          // exact
  tmp1.add(fEta, tmp1, GMP_RNDU);      // >= 1 + eta
  tmp1.mul(tmp1, tmp1, GMP_RNDU);      // >= (1 + eta) ^ 2
  tmp1.add(tmp1, fEpsilon, GMP_RNDU);
  // Computes tmp2 <= delta - eta ^ 2
  tmp2.mul(fEta, fEta, GMP_RNDU);
  tmp2.sub(fDelta, tmp2, GMP_RNDD);
  FPLLL_CHECK(tmp2 > 0, "invalid LLL parameters, eta must be < sqrt(delta)");
  // Computes rho >= ((1 + eta) ^ 2 + epsilon) / (delta - eta ^ 2)
  fRho.div(tmp1, tmp2, GMP_RNDU);
  rho = fRho.get_d(GMP_RNDU);

  /* Computes minprec >= constant + 2 * log2(d) - log2(epsilon) + d * log2(rho)
     (constant = 5 for GSO, 10 for LLL) */
  tmp1.log(fD, GMP_RNDU);              // >= log(d)
  tmp1.mul_2si(tmp1, 1);               // >= 2 * log(d)
  tmp2.log(fEpsilon, GMP_RNDD);        // <= log(epsilon) (<= 0)
  tmp1.sub(tmp1, tmp2, GMP_RNDU);      // >= 2 * log(d) - log(epsilon)
  tmp2.log(fRho, GMP_RNDU);            // >= log(rho)
  tmp2.mul(fD, tmp2, GMP_RNDU);        // >= d * log(rho)
  tmp1.add(tmp1, tmp2, GMP_RNDU);      // >= 2*log(d)-log(epsilon)+d*log(rho)
  tmp2 = 2.0;                          // exact
  tmp2.log(tmp2, GMP_RNDD);            // <= log(2)
  tmp1.div(tmp1, tmp2, GMP_RNDU);      // >= 2*log2(d)-log2(epsilon)+d*log2(rho)
  tmp2 = (algo == MINPREC_L2) ? 10.0 : 5.0;
  fMinprec.add(tmp1, tmp2, GMP_RNDU);
  int minprec = static_cast<int>(ceil(fMinprec.get_d(GMP_RNDU)));
  mpfr_free_cache();
  Float::setprec(oldprec);
  return minprec;
}

int gsoMinPrec(double& rho, int d, double delta, double eta, double epsilon) {
  return computeMinPrec(rho, d, delta, eta, epsilon, MINPREC_GSO);
}

int l2MinPrec(int d, double delta, double eta, double epsilon) {
  double rho;
  return computeMinPrec(rho, d, delta, eta, epsilon, MINPREC_L2);
}

/**
 * Computes the volume of a d-dimensional hypersphere of radius 1.
 */
void sphereVolume(Float& volume, int d) {
  Float rtmp1;
  volume = pow(M_PI, (double)(d / 2));

  if (d % 2 == 0)
    for (int i = 1; i <= d / 2; i++) {
      rtmp1 = (double) i;
      volume.div(volume, rtmp1);
    }
  else
    for (int i = 0; i <= d / 2; i++) {
      rtmp1 = 2.0 / (double)(2 * i + 1);
      volume.mul(volume, rtmp1);
    }
}

/**
 * Estimates the cost of the enumeration for SVP.
 */
void costEstimate(Float& cost, const Float& bound,
                  const Matrix<Float>& r, int dimMax) {
  Float det, levelCost, tmp1;
  det = 1.0;
  cost = 0.0;

  for (int i = dimMax - 1; i >= 0; i--) {
    tmp1.div(bound, r(i, i));
    det.mul(det, tmp1);

    levelCost.sqrt(det);
    sphereVolume(tmp1, dimMax - i);
    levelCost.mul(levelCost, tmp1);

    cost.add(cost, levelCost);
  }
}

#ifdef FPLLL_V3_COMPAT

void gramSchmidt(const IntMatrix& b, Matrix<Float>& mu, FloatVect& rdiag) {
  int d = b.getRows();
  int n = b.getCols();
  Matrix<Float> r(d, d);
  Integer dotProd;
  Float coeff;

  FPLLL_DEBUG_CHECK(mu.getRows() == d && mu.getCols() == d);
  if (static_cast<int>(rdiag.size()) != d) rdiag.resize(d);

  for (int i = 0; i < d; i++) {
    for (int j = 0; j <= i; j++) {
      dotProd = 0;
      for (int k = 0; k < n; k++) {
        dotProd.addmul(b(i, k), b(j, k));
      }
      coeff.set_z(dotProd);
      for (int k = 0; k < j; k++) {
        coeff.submul(mu(j, k), r(i, k));
      }
      r(i, j) = coeff;
      mu(i, j).div(coeff, r(j, j));
    }
    rdiag[i].set(r(i, i));
  }
}

#endif // #ifdef FPLLL_V3_COMPAT

FPLLL_END_NAMESPACE
