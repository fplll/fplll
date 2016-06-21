#include "bkz_params.h"

FPLLL_BEGIN_NAMESPACE

const Pruning &Strategy::get_pruning(double radius, double gh) const
{
  double gh_factor    = radius / gh;
  double closest_dist = pow(2, 80);
  auto best           = pruning_parameters.begin();

  for (auto it = pruning_parameters.begin(); it != pruning_parameters.end(); ++it)
  {
    if (abs(it->radius_factor - gh_factor) < closest_dist)
    {
      closest_dist = abs(it->radius_factor - gh_factor);
      best         = it;
    }
  }

  return *best;
}

FPLLL_END_NAMESPACE
