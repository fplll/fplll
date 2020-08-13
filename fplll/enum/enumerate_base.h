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
 * LevelTreeCounter. This counter is used to count the number of nodes that are visited on each
 * level on the tree. The structure of this class is the following: we keep an array in this object
 * that is indexed by the level. For example, if the current level k is k, then the number of nodes
 * visited at level k can be found in _nodes[k].
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
  inline LevelTreeCounter() : _nodes{} {}
  // This returns the array containing the nodes on each level.
  inline UnderlyingCounterType get_nodes() const { return _nodes; }
  // This returns the total number of nodes in the tree.
  inline UnderlyingIndividualType get_total_nodes() const
  {
    UnderlyingIndividualType total = 0;
    for (unsigned i = 0; i < _nodes.size(); i++)
    {
      total += _nodes[i];
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
  inline LevelTreeCounter &operator=(const UnderlyingCounterType &value)
  {
    _nodes = value;
    return *this;
  }

  inline bool operator==(const UnderlyingCounterType value) { return _nodes == value; }

  inline bool operator==(const UnderlyingIndividualType value)
  {
    return is_valid() && get_total_nodes() == value;
  }

  inline bool operator!=(const UnderlyingIndividualType value) { return !(*this == value); }

  // These operators add an array and another LevelTreeCounter to this tree counter.
  // Usage: in this instance we'd have the following.
  // LevelTreeCounter a;
  // ....
  // a += enum_obj.get_nodes();
  inline LevelTreeCounter &operator+=(const UnderlyingCounterType &value)
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
  inline LevelTreeCounter &operator+=(const LevelTreeCounter &value)
  {
    for (unsigned i = 0; i < FPLLL_MAX_ENUM_DIM; i++)
    {
      _nodes[i] += value._nodes[i];
    }
    return *this;
  }

  // Reset. We simply clear the _nodes array and set the _total_nodes to 0.
  inline void reset() { std::fill(std::begin(_nodes), std::end(_nodes), 0); }

  // Similarly to WholeTreeCounter, we denote failure as ~uint64_t(0).
  inline bool is_valid() const { return _nodes[0] != ~uint64_t(0); }
  inline void invalidate() { _nodes[0] = ~uint64_t(0); }
  inline constexpr static UnderlyingCounterType produce_invalid_entry()
  {
    return UnderlyingCounterType{{~uint64_t(0)}};
  }

private:
  UnderlyingCounterType _nodes;
};

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

  inline WholeTreeCounter(UnderlyingCounterType starting_count = 0) : _nodes{starting_count} {}

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
  inline WholeTreeCounter &operator=(const UnderlyingCounterType value)
  {
    _nodes = value;
    return *this;
  }

  // This is an override for the += operator: this method just operates on the underlying type
  // directly.
  inline WholeTreeCounter &operator+=(const UnderlyingCounterType &value)
  {
    _nodes += value;
    return *this;
  }

  inline WholeTreeCounter &operator+=(const LevelTreeCounter::UnderlyingCounterType &value)
  {
    // If valid, add up the entries
    for (unsigned i = 0; i < value.size(); i++)
    {
      _nodes += value[i];
    }

    return *this;
  }

  inline WholeTreeCounter &operator=(const LevelTreeCounter::UnderlyingCounterType &value)
  {
    _nodes = 0;
    for (unsigned i = 0; i < value.size(); i++)
    {
      _nodes += value[i];
    }
    return *this;
  }

  // This is an override for the += operator: here we accept another WholeTreecounter, which we just
  // add to our own counter.
  inline WholeTreeCounter &operator+=(const WholeTreeCounter &value)
  {
    _nodes += value._nodes;
    return *this;
  }

  inline bool operator==(const UnderlyingCounterType value) { return _nodes == value; }

  inline bool operator!=(const UnderlyingCounterType value) { return !(*this == value); }

  // This resets the number of nodes in the tree to 0.
  inline void reset() { _nodes = 0; }

  // This sets the counter to the invalid counter.
  inline void invalidate() { _nodes = ~uint64_t(0); }

  // This method tells us if the value in the counter is valid.
  inline bool is_valid() const { return _nodes != ~uint64_t(0); }
  inline constexpr static UnderlyingCounterType produce_invalid_entry() { return ~uint64_t(0); }

private:
  // This is the underlying counter for the number of nodes in the tree.
  UnderlyingCounterType _nodes;
};

/**
 * InvalidCounterFactory.
 * To separate producing invalid entries from the particular instance of a counter, this namespace
 * uses template resolution rules to delegate to the correct counter class, at compile time. This
 * then lets the compiler return the right definition of "invalid" depending on the counter that's
 * used.
 *
 * Note: we leave the generic cases empty. If this function is called with a type that isn't
 * explicitly specified, then the linker will fail. This is necessary for correctness.
 */
namespace InvalidCounterFactory
{
template <typename T> constexpr T produce_invalid_entry();

template <typename extenum_return_type, typename in_data_type,
          std::size_t in_width = FPLLL_MAX_PARALLEL_ENUM_DIM>
extenum_return_type
convert_array_counter_to_extenum_rt(const std::array<in_data_type, in_width> &input);

template <> constexpr LevelTreeCounter::UnderlyingCounterType produce_invalid_entry()
{
  return LevelTreeCounter::produce_invalid_entry();
}

template <> constexpr WholeTreeCounter::UnderlyingCounterType produce_invalid_entry()
{
  return WholeTreeCounter::produce_invalid_entry();
}

}  // namespace InvalidCounterFactory

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

  CounterClassWrapper<CounterClass> &operator=(const LevelTreeCounter::UnderlyingCounterType &value)
  {
    base_counter = value;
    return *this;
  }

  CounterClassWrapper<CounterClass> &operator=(const WholeTreeCounter::UnderlyingCounterType &value)
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

  bool operator==(const UnderlyingCounterType &value) { return base_counter == value; }

  template <typename = std::enable_if<
                std::is_same<UnderlyingCounterType, UnderlyingIndividualType>::value>>
  bool operator==(const UnderlyingIndividualType &value)
  {
    return base_counter == value;
  }
  inline void reset() { base_counter.reset(); }

  inline bool is_valid() { return base_counter.is_valid(); }
  inline void invalidate() { return base_counter.invalidate(); }
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
