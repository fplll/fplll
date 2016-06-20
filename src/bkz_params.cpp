#include "bkz_params.h"

FPLLL_BEGIN_NAMESPACE

Pruning::Pruning(int blockSize, int level) {
  int startDescent = blockSize - level;

  if (startDescent > blockSize)
    startDescent = blockSize;

  if (startDescent < 1)
    startDescent = 1;

  coefficients.resize(blockSize);
  for (int k = 0; k < startDescent; k++) {
    coefficients[k] = 1.0;
  }
  for (int k = 0; k < blockSize - startDescent; k++) {
    coefficients[startDescent + k] = ((double)(blockSize - k - 1)) / blockSize;
  }
  // TODO: need to adapt probability
  radiusFactor = 1.0;
  probability = 1.0;
}

const Pruning &BKZParam::getPruning(double radius, double gh) const {
  double ghFactor = radius/gh;
  double closestDist = pow(2, 80);
  auto best = pruning.begin();

  for(auto it = pruning.begin(); it != pruning.end(); ++it) {
    if (abs(it->radiusFactor - ghFactor) < closestDist) {
      closestDist = abs(it->radiusFactor - ghFactor);
      best = it;
    }
  }

  return *best;
}

FPLLL_END_NAMESPACE
