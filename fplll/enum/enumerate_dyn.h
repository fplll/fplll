#ifndef FPLLL_ENUMERATE_DYN_H
#define FPLLL_ENUMERATE_DYN_H

#include <array>
#include <fplll/enum/enumerate_base.h>
#include <fplll/enum/evaluator.h>
#include <fplll/gso_interface.h>
#include <memory>

FPLLL_BEGIN_NAMESPACE

template <typename ZT, typename FT> class EnumerationDyn : public EnumerationBase
{
public:
  EnumerationDyn(MatGSOInterface<ZT, FT> &gso, Evaluator<FT> &evaluator,
                 const vector<int> &max_indices = vector<int>())
      : _gso(gso), _evaluator(evaluator)
  {
    _max_indices = max_indices;
    std::fill(nodes.begin(), nodes.end(), 0);
  }

  void enumerate(int first, int last, FT &fmaxdist, long fmaxdistexpo,
                 const vector<FT> &target_coord = vector<FT>(),
                 const vector<enumxt> &subtree  = vector<enumxt>(),
                 const vector<enumf> &pruning = vector<enumf>(), bool dual = false,
                 bool subtree_reset = false);

  inline uint64_t get_nodes(const int level = -1) const
  {
    if (level == -1)
    {
      return std::accumulate(nodes.cbegin(), nodes.cend(), 0);
    }
    return nodes[level];
  }

  inline array<uint64_t, FPLLL_MAX_ENUM_DIM> get_nodes_array() { return nodes; }

private:
  MatGSOInterface<ZT, FT> &_gso;
  Evaluator<FT> &_evaluator;
  vector<FT> target;

  vector<enumf> pruning_bounds;
  enumf maxdist;
  vector<FT> fx;

  void prepare_enumeration(const vector<enumxt> &subtree, bool solvingsvp, bool subtree_reset);

  void do_enumerate();

  void set_bounds();
  void reset(enumf cur_dist, int cur_depth);
  virtual void process_solution(enumf newmaxdist);
  virtual void process_subsolution(int offset, enumf newdist);
};

FPLLL_END_NAMESPACE

#endif