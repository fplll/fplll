#include "bkz_params.h"

FPLLL_BEGIN_NAMESPACE

const Pruning &Strategy::getPruning(double radius, double gh) const {
  double ghFactor = radius/gh;
  double closestDist = pow(2, 80);
  auto best = pruning_parameters.begin();

  for(auto it = pruning_parameters.begin(); it != pruning_parameters.end(); ++it) {
    if (abs(it->radiusFactor - ghFactor) < closestDist) {
      closestDist = abs(it->radiusFactor - ghFactor);
      best = it;
    }
  }

  return *best;
}

FPLLL_END_NAMESPACE
