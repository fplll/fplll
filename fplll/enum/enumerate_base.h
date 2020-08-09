/* Copyright (C) 2008-2011 Xavier Pujol
   (C) 2015 Michael Walter.
   (C) 2016 Marc Stevens. (generic improvements, auxiliary solutions, subsolutions)
   (C) 2016 Guillaume Bonnoron. (CVP improvements)

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

#ifndef FPLLL_ENUMERATE_BASE_H
#define FPLLL_ENUMERATE_BASE_H

#include "fplll/fplll_config.h"
#include "fplll/nr/nr.h"
#include <array>
#include <cfenv>
#include <cmath>
#include <vector>

FPLLL_BEGIN_NAMESPACE

inline void roundto(int &dest, const double &src) { dest = (int)round(src); }
inline void roundto(double &dest, const double &src) { dest = round(src); }

/* config */
#define MAXTEMPLATEDDIMENSION 80  // unused
//#define FORCE_ENUM_INLINE // not recommended
/* end config */

#ifndef __has_attribute
#define __has_attribute(x) 0  // Compatibility with non - GCC/clang compilers.
#endif
#if __has_attribute(always_inline)
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define ALWAYS_INLINE
#endif

#ifndef FORCE_ENUM_INLINE
#define ENUM_ALWAYS_INLINE
#else
#define ENUM_ALWAYS_INLINE ALWAYS_INLINE
#endif

/***
 * WholeTreeCounter. This class counts every single node visited in enumeration and stores it in the
 * _nodes field. This class can be viewed as a direct replacement for the old method of naively
 * counting every element in the enumeration tree.
 *
 * Note: in order to retain backwards compatibility with the previous mechanisms in fplll, this is
 * the default template parameter for all enumeration classes. This may not be the best case for
 * your application; instead, it may be better to use the LevelTreecounter.
 */
class WholeTreeCounter
{
public:
  // This is just here so we can replace the type later if necessary.
  // For almost all applications, uint64_t is likely to be fine.
  using UnderlyingIndividualType = uint64_t;
  using UnderlyingCounterType    = UnderlyingIndividualType;

  WholeTreeCounter(UnderlyingCounterType starting_count = 0) : _nodes{starting_count} {}

  // This is a getter method that is typically called from the outside world.
  inline UnderlyingCounterType get_nodes() const { return _nodes; }
  // In this instance, this method does exactly the same as get_nodes().
  inline UnderlyingCounterType get_total_nodes() const { return _nodes; }

  // This adds amount to _nodes.
  inline void update_nodes_count(const unsigned int index, const uint64_t amount)
  {
    _nodes += amount;
  }

  // We provide an equals operator here for copying.
  // We accept a copy: we're only copying 64 bits at a time, and so doing so is cheap.
  WholeTreeCounter &operator=(const UnderlyingCounterType value)
  {
    _nodes = value;
    return *this;
  }

  // This is an override for the += operator: this method just operates on the underlying type
  // directly.
  WholeTreeCounter &operator+=(const UnderlyingCounterType &value)
  {
    _nodes += value;
    return *this;
  }

  // This is an override for the += operator: here we accept another WholeTreecounter, which we just
  // add to our own counter.
  WholeTreeCounter &operator+=(const WholeTreeCounter &value)
  {
    _nodes += value._nodes;
    return *this;
  }

  // This resets the number of nodes in the tree to 0.
  inline void reset() { _nodes = 0; }

  // This method tells us if the value in the counter is valid.
  inline bool is_valid() const { return _nodes != ~uint64_t(0); }

private:
  // This is the underlying counter for the number of nodes in the tree.
  UnderlyingCounterType _nodes;
};

/***
 * LevelTreeCounter. This counter is used to count the number of nodes that are visited on each
 * level on the tree. The structure of this class is the following: we keep an array in this object
 * that is indexed by the level. For example, if the current level k is k, then the number of nodes
 * visited at level k can be found in _nodes[k].
 *
 * We also keep a track of a memoised quantity, _total_nodes. This keeps track of the sum of the
 * elements in the array in total. This is done so that we don't need to add all of the elements in
 * the array every time we want to ascertain the total number of nodes that were visited during
 * enumeration.
 */
class LevelTreeCounter
{
public:
  // These using declarations are simply for neater return types & easier changing.
  // Similarly to WholeTreeCounter, it's likely that uint64_t is the perfect type for your
  // application.
  using UnderlyingIndividualType = uint64_t;
  using UnderlyingCounterType    = std::array<UnderlyingIndividualType, FPLLL_MAX_ENUM_DIM>;

  // Default, zero initialisation for the constructors
  LevelTreeCounter() : _nodes{} {}
  LevelTreeCounter(const UnderlyingIndividualType total_nodes) : _nodes{} {}
  // Copy constructor for an array input
  LevelTreeCounter(const UnderlyingCounterType &node_set) : _nodes{node_set} {}
  // Move constructor for an array input
  LevelTreeCounter(const UnderlyingCounterType &&node_set)
      : _nodes{std::move(node_set)} {}
  // This returns the array containing the nodes on each level.
  inline UnderlyingCounterType get_nodes() const { return _nodes; }
  // This returns the total number of nodes in the tree.
  inline UnderlyingIndividualType get_all_nodes() const { 
      UnderlyingIndividualType total = 0;
      for(unsigned i = 0; i < _nodes.size();i++) {
          total+= _nodes[i];
      }
      return total;
  }

  // This updates the number of nodes visited in the tree.
  inline void update_nodes_count(const unsigned int index, const uint64_t amount)
  {
    _nodes[index] += amount;
  }

  // This is a copy assignment operator for the number of nodes in the tree.
  // In this instance, this is an array copy.
  LevelTreeCounter &operator=(const UnderlyingCounterType &value)
  {
    _nodes = value;
    return *this;
  }

  // These operators add an array and another LevelTreeCounter to this tree counter.
  // Usage: in this instance we'd have the following.
  // LevelTreeCounter a;
  // ....
  // a += enum_obj.get_nodes();
  LevelTreeCounter &operator+=(const UnderlyingCounterType &value)
  {
    // Note: we use FPLLL_MAX_ENUM_DIM. This is safe because the counter type is defined with width
    // FPLLL_MAX_ENUM. A sensible compiler should unroll this loop.
    for (unsigned i = 0; i < FPLLL_MAX_ENUM_DIM; i++)
    {
      _nodes[i] += value[i];
    }
    return *this;
  }
  // Usage:
  // LevelTreeCounter a = .....
  // LevelTreeCounter b = .....
  // a += b;
  LevelTreeCounter &operator+=(const LevelTreeCounter &value)
  {
    for (unsigned i = 0; i < FPLLL_MAX_ENUM_DIM; i++)
    {
      _nodes[i] += value._nodes[i];
    }
    return *this;
  }

  // Reset. We simply clear the _nodes array and set the _total_nodes to 0.
  inline void reset()
  {
    std::fill(std::begin(_nodes), std::end(_nodes), 0);
  }

  // Similarly to WholeTreeCounter, we denote failure as ~uint64_t(0).
  inline bool is_valid() const { return _nodes[0] != ~uint64_t(0); }

private:
  UnderlyingCounterType _nodes;
};

/***
 * CounterClassWrapper.
 * This class is a wrapper class for the underlying counter that is used. In particular, this is an
 * instance of the policy pattern from Alexandrescu's "Modern C++ design". The idea is to provide a
 * unified interface to the outside world, regardless of the underlying type of the counter. This is
 * useful for two reasons. Firstly, this essentially mimmicks C++20's concepts. Rather than
 * explicitly constraining the template parameter via concepts, the compiler will instead error if
 * one supplies a CounterClass that doesn't implement the underlying methods. This is useful because
 * it ensures that this class cannot be instantiated incorrectly. The second advantage is that we
 * don't need to resort to virtual functions: since the type is known at compile-time, the compiler
 * can correctly link to the correct function, and so there's limited performance differences. It
 * also means that new counters can be dropped into place easily, simply by writing similar
 * functions to those above.
 */

template <class CounterClass> class CounterClassWrapper
{
  // This is the counter that's being wrapped.
  CounterClass base_counter;

public:
  // We pass this type up so that the EnumerationBase class can access the counter's store.
  // This is so that we can get neater return types.
  using UnderlyingCounterType    = typename CounterClass::UnderlyingCounterType;
  using UnderlyingIndividualType = typename CounterClass::UnderlyingIndividualType;

  CounterClassWrapper() : base_counter{} {}
  CounterClassWrapper(const UnderlyingIndividualType starting_value) : base_counter{starting_value}
  {
  }
  // Every method in this class just delegates to the underlying counter's methods.
  // As every method here is inline, this should be done without paying the overhead of a function
  // call.
  inline UnderlyingCounterType get_nodes() const { return base_counter.get_nodes(); }

  inline void update_nodes_count(const unsigned int index, const uint64_t amount = 1)
  {
    base_counter.update_nodes_count(index, amount);
  }

  CounterClassWrapper<CounterClass> &operator=(const UnderlyingCounterType &value)
  {
    base_counter = value;
    return *this;
  }

  CounterClassWrapper<CounterClass> &operator+=(const UnderlyingCounterType &value)
  {
    base_counter += value;
    return *this;
  }

  CounterClassWrapper<CounterClass> &operator+=(const CounterClass &value)
  {
    base_counter += value;
    return *this;
  }

  inline void reset() { base_counter.reset(); }

  inline bool is_valid() { return base_counter.is_valid(); }

  inline UnderlyingIndividualType get_total_nodes() const { return base_counter.get_total_nodes(); }
};

template <class CounterClass> class EnumerationBase
{
public:
  static const int maxdim = FPLLL_MAX_ENUM_DIM;

  // These methods all delegate to the interface from the CounterClassWrapper (which, in turn,
  // delegate)
  using UnderlyingCounterType = typename CounterClassWrapper<CounterClass>::UnderlyingCounterType;
  using UnderlyingIndividualType =
      typename CounterClassWrapper<CounterClass>::UnderlyingIndividualType;

  inline UnderlyingCounterType get_nodes() const { return nodes_counter.get_nodes(); }
  inline UnderlyingIndividualType get_total_nodes() const
  {
    return nodes_counter.get_total_nodes();
  }

  inline void update_nodes_count(const unsigned int index = 0, const long amount = 1)
  {
    nodes_counter.update_nodes_count(index, amount);
  }

  inline void reset_nodes_count() { nodes_counter.reset(); }

  virtual ~EnumerationBase() {}

protected:
  /* configuration */
  bool dual;
  bool is_svp;
  bool resetflag;

  /* enumeration input */
  enumf mut[maxdim][maxdim];
  array<enumf, maxdim> rdiag, partdistbounds;
  int d, k_end;  // dimension, subtreelevel

  /* partial sum cache */
  enumf center_partsums[maxdim][maxdim];
  array<enumf, maxdim> center_partsum;
  array<int, maxdim> center_partsum_begin;

  /* enumeration data for each level */
  array<enumf, maxdim> partdist, center, alpha;
  array<enumxt, maxdim> x, dx, ddx;
  array<enumf, maxdim> subsoldists;

  /* CVP reset informations */
  vector<int> _max_indices;
  int reset_depth;

  int k, k_max;
  bool finished;
  /* Node counter wrapper */
  CounterClassWrapper<CounterClass> nodes_counter;

  template <int kk, int kk_start, bool dualenum, bool findsubsols, bool enable_reset> struct opts
  {
  };

  /* need templated function argument for support of integer specialization for kk==-1 */
  template <int kk, int kk_start, bool dualenum, bool findsubsols, bool enable_reset>
  inline void enumerate_recursive(opts<kk, kk_start, dualenum, findsubsols, enable_reset>)
      ENUM_ALWAYS_INLINE;
  template <int kk_start, bool dualenum, bool findsubsols, bool enable_reset>
  inline void enumerate_recursive(opts<-1, kk_start, dualenum, findsubsols, enable_reset>)
  {
  }

  /* simple wrapper with no function argument as helper for dispatcher */
  template <int kk, bool dualenum, bool findsubsols, bool enable_reset>
  void enumerate_recursive_wrapper()
  {
    // kk < maxdim-1:
    // kk < kk_end                         (see enumerate_loop(), enumerate_base.cpp)
    // kk_end = d - subtree.size() <= d    (see prepare_enumeration(), enumerate.cpp)
    // d < maxdim                          (see enumerate(), enumerate.cpp)
    enumerate_recursive(
        opts<(kk < (maxdim - 1) ? kk : -1), 0, dualenum, findsubsols, enable_reset>());
  }

  template <bool dualenum, bool findsubsols, bool enable_reset>
  inline void enumerate_recursive_dispatch(int kk);

  template <bool dualenum, bool findsubsols, bool enable_reset> void enumerate_loop();

  virtual void reset(enumf, int)                              = 0;
  virtual void process_solution(enumf newmaxdist)             = 0;
  virtual void process_subsolution(int offset, enumf newdist) = 0;

  int rounding_backup;
  void save_rounding()
  {
    rounding_backup = std::fegetround();
    std::fesetround(FE_TONEAREST);
  }
  void restore_rounding() { std::fesetround(rounding_backup); }

  inline bool next_pos_up()
  {
    ++k;
    if (partdist[k] != 0.0)
    {
      x[k] += dx[k];
      ddx[k] = -ddx[k];
      dx[k]  = ddx[k] - dx[k];
    }
    else
    {
      if (k >= k_end)
        return false;
      k_max = k;
      ++x[k];
    }
    return true;
  }
};

FPLLL_END_NAMESPACE

#endif
