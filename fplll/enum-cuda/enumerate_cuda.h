#ifndef FPLLL_ENUMERATE_CUDA_H
#define FPLLL_ENUMERATE_CUDA_H

#include <array>
#include <fplll/enum/enumerate_base.h>
#include <fplll/enum/enumerate_ext.h>
#include <fplll/enum/evaluator.h>
#include <fplll/gso_interface.h>
#include <functional>
#include <memory>

FPLLL_BEGIN_NAMESPACE

template<typename FT>
struct AllPointsEvaluator;

template <typename ZT, typename FT> class CudaEnumeration : public ExternalEnumeration<ZT, FT>
{
private:
  using ExternalEnumeration<ZT, FT>::_gso;
  using ExternalEnumeration<ZT, FT>::_evaluator;
  using ExternalEnumeration<ZT, FT>::_pruning;
  using ExternalEnumeration<ZT, FT>::_normexp;

  using ExternalEnumeration<ZT, FT>::_nodes;
  using ExternalEnumeration<ZT, FT>::_dual;
  using ExternalEnumeration<ZT, FT>::_d;
  using ExternalEnumeration<ZT, FT>::_first;
  using ExternalEnumeration<ZT, FT>::_maxdist;
  using ExternalEnumeration<ZT, FT>::_fx;

  using ExternalEnumeration<ZT, FT>::callback_set_config;

public:
  CudaEnumeration(MatGSOInterface<ZT, FT> &gso, Evaluator<FT> &evaluator)
      : ExternalEnumeration<ZT, FT>(gso, evaluator)
  {
  }

  virtual ~CudaEnumeration() = default;

  virtual bool enumerate(int first, int last, FT &fmaxdist, long fmaxdistexpo,
                         const vector<enumf> &pruning = vector<enumf>(), bool dual = false);

private:
  std::unique_ptr<AllPointsEvaluator<FT>> enumerate_start_points(size_t start_dims, FT &fmaxdist, long fmaxdistexpo);
};

FPLLL_END_NAMESPACE

#endif