#include <cstdio>
#include "bkz_params.h"
#include "io/json.hpp"
using json = nlohmann::json;

FPLLL_BEGIN_NAMESPACE

const char *default_strategy_directory()
{
  return FPLLL_DEFAULT_STRATEGY_DIR;
}

const char *default_strategy()
{
  return FPLLL_DEFAULT_STRATEGY;
}

const std::string strategy_full_path(const char *strategy_path)
{
  FILE *fh = fopen(strategy_path, "r");
  string path;

  if (fh)
  {
    fclose(fh);
    path = strategy_path;
    return path;
  } else if (strcmp(FPLLL_DEFAULT_STRATEGY_DIR, "") > 0)
  {
    path += FPLLL_DEFAULT_STRATEGY_DIR;
    path += "/";
    path += strategy_path;
    fh = fopen(path.c_str(), "r");
    if (fh) {
      fclose(fh);
    }
    return path;
  }
  path = "";
  return path;
}

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

vector<Strategy> load_strategies_json(const char* filename)
{
  std::fstream fs;
  fs.open(filename, std::fstream::in);
  json js;
  js << fs;
  fs.close();

  vector<Strategy> strategies;

  for(auto it = js.begin(); it != js.end(); ++it)
  {
    const json &j_strat = *it;
    const size_t block_size = j_strat["block_size"];
    FPLLL_DEBUG_CHECK(block_size < 4096); // some arbitrary upper limit

    // ensure to add this strategy in the right place
    while(strategies.size() <= block_size)
    {
      strategies.emplace_back();
    }

    Strategy strategy;

    if (j_strat.find("preprocessing") != j_strat.end())
    {
      for (auto p_it = j_strat["preprocessing"].begin(); p_it != j_strat["preprocessing"].end();
           ++p_it)
      {
        if ((*p_it).is_number())
        {
          strategy.preprocessing_blocksizes.emplace_back((*p_it).get<int>());
        }
        else
        {
          strategy.preprocessing_blocksizes.emplace_back((*p_it)["block_size"]);
        }
      }
    }

    if (j_strat.find("pruning_coefficients") != j_strat.end())
    {
      for (auto p_it = j_strat["pruning_coefficients"].begin();
           p_it != j_strat["pruning_coefficients"].end(); ++p_it)
      {
        const json &j_prun = *p_it;
        Pruning pruning;
        pruning.radius_factor = j_prun[0];

        // fplll enforces that the first pruning coefficient is 1.0
        FPLLL_DEBUG_CHECK(j_prun[1][0] == 1.0);

        for (auto c_it = j_prun[1].begin(); c_it != j_prun[1].end(); ++c_it)
        {
          double c = (*c_it);
          pruning.coefficients.emplace_back(c);
        }
        pruning.probability = j_prun[2];
        FPLLL_DEBUG_CHECK(pruning.probability > 0.0 && pruning.probability <= 1.0);

        strategy.pruning_parameters.emplace_back(pruning);
      }
    }

    strategies[block_size] = std::move(strategy);
  }

  // finally, we make sure all strategies are sound
  for(auto it = strategies.begin(); it != strategies.end(); ++it) {
    if(it->pruning_parameters.size() == 0)
      it->pruning_parameters.emplace_back(Pruning());
  }

  return strategies;
}

FPLLL_END_NAMESPACE

