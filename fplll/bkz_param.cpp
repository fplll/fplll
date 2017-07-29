
#include "bkz_param.h"
#include "io/json.hpp"
#include "nr/nr.h"
#include <cstdio>

using json = nlohmann::json;

FPLLL_BEGIN_NAMESPACE

PruningParams PruningParams::LinearPruningParams(int block_size, int level)
{

  PruningParams pruning = PruningParams();
  int start_descent     = block_size - level;

  if (start_descent > block_size)
    start_descent = block_size;

  if (start_descent < 1)
    start_descent = 1;

  pruning.coefficients.resize(block_size);
  for (int k = 0; k < start_descent; k++)
  {
    pruning.coefficients[k] = 1.0;
  }
  for (int k = 0; k < block_size - start_descent; k++)
  {
    pruning.coefficients[start_descent + k] = ((double)(block_size - k - 1)) / block_size;
  }
  pruning.gh_factor   = 1.0;
  pruning.metric      = PRUNER_METRIC_PROBABILITY_OF_SHORTEST;
  pruning.expectation = fplll::svp_probability<FP_NR<double>>(pruning.coefficients).get_d();

  return pruning;
}

const std::string &default_strategy_path()
{
  static const std::string ret(FPLLL_DEFAULT_STRATEGY_PATH);
  return ret;
}

const std::string strategy_full_path(const std::string &strategy_path)
{
  if (std::ifstream(strategy_path).good())
    return strategy_path;
  std::string path = default_strategy_path() + "/" + strategy_path;
  if (std::ifstream(path).good())
    return path;
  path.clear();
  return path;
}

const std::string &default_strategy()
{
  static const std::string ret(FPLLL_DEFAULT_STRATEGY);
  return ret;
}

const PruningParams &Strategy::get_pruning(double radius, double gh) const
{
  double gh_factor    = radius / gh;
  double closest_dist = pow(2, 80);
  auto best           = pruning_parameters.begin();

  for (auto it = pruning_parameters.begin(); it != pruning_parameters.end(); ++it)
  {
    if (fabs(it->gh_factor - gh_factor) < closest_dist)
    {
      closest_dist = fabs(it->gh_factor - gh_factor);
      best         = it;
    }
  }

  return *best;
}

vector<Strategy> load_strategies_json(const std::string &filename)
{
  json js;
  {
    std::ifstream fs(filename);
    if (fs.fail())
      throw std::runtime_error("Cannot open strategies file.");
    fs >> js;
  }

  vector<Strategy> strategies;

  for (auto it = js.begin(); it != js.end(); ++it)
  {
    const json &j_strat     = *it;
    const size_t block_size = j_strat["block_size"];
    FPLLL_DEBUG_CHECK(block_size < 4096);  // some arbitrary upper limit

    // ensure to add this strategy in the right place
    while (strategies.size() <= block_size)
    {
      strategies.emplace_back();
    }

    Strategy strategy;
    strategy.block_size = block_size;

    if (j_strat.find("preprocessing_block_sizes") != j_strat.end())
    {
      for (auto p_it = j_strat["preprocessing_block_sizes"].begin();
           p_it != j_strat["preprocessing_block_sizes"].end(); ++p_it)
      {
        if ((*p_it).is_number())
        {
          strategy.preprocessing_block_sizes.emplace_back((*p_it).get<int>());
        }
        else
        {
          strategy.preprocessing_block_sizes.emplace_back((*p_it)["block_size"]);
        }
      }
    }

    if (j_strat.find("pruning_parameters") != j_strat.end())
    {
      for (auto p_it = j_strat["pruning_parameters"].begin();
           p_it != j_strat["pruning_parameters"].end(); ++p_it)
      {
        const json &j_prun = *p_it;
        PruningParams pruning;
        pruning.gh_factor = j_prun[0];

        for (auto c_it = j_prun[1].begin(); c_it != j_prun[1].end(); ++c_it)
        {
          double c = (*c_it);
          pruning.coefficients.emplace_back(c);
        }
        pruning.expectation = j_prun[2];
        // TODO: don't hardcode success probability as metric
        pruning.metric = PRUNER_METRIC_PROBABILITY_OF_SHORTEST;
        FPLLL_DEBUG_CHECK(pruning.expectation > 0.0 && pruning.expectation <= 1.0);

        strategy.pruning_parameters.emplace_back(pruning);
      }
    }

    strategies[block_size] = std::move(strategy);
  }

  // finally, we make sure all strategies are sound
  for (auto it = strategies.begin(); it != strategies.end(); ++it)
  {
    if (it->pruning_parameters.size() == 0)
      it->pruning_parameters.emplace_back(PruningParams());
  }

  return strategies;
}

FPLLL_END_NAMESPACE
