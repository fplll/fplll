#ifndef FPLLL_COUNTERS_HEADER_H
#define FPLLL_COUNTERS_HEADER_H

#include "fplll/fplll_config.h"
#include "fplll/nr/nr.h"
#include <array>
#include <cmath>

/***
 * This file contains the implementations of all of the various counters for the enumeration
 * process.
 */

FPLLL_BEGIN_NAMESPACE

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
  LevelTreeCounter() : _nodes{} {}
  // This returns the array containing the nodes on each level.
  inline UnderlyingCounterType get_nodes() const { return _nodes; }
  // This returns the total number of nodes in the tree.
  inline UnderlyingIndividualType get_all_nodes() const
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
  inline void reset() { std::fill(std::begin(_nodes), std::end(_nodes), 0); }

  // Similarly to WholeTreeCounter, we denote failure as ~uint64_t(0).
  inline bool is_valid() const { return _nodes[0] != ~uint64_t(0); }
  inline void invalidate() { _nodes[0] = ~uint64_t(0); }
  constexpr static UnderlyingCounterType produce_invalid_entry()
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

  WholeTreeCounter &operator+=(const LevelTreeCounter::UnderlyingCounterType &value)
  {
    // If valid, add up the entries
    for (unsigned i = 0; i < value.size(); i++)
    {
      _nodes += value[i];
    }

    return *this;
  }

  WholeTreeCounter &operator=(const LevelTreeCounter::UnderlyingCounterType &value)
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
  WholeTreeCounter &operator+=(const WholeTreeCounter &value)
  {
    _nodes += value._nodes;
    return *this;
  }

  // This resets the number of nodes in the tree to 0.
  inline void reset() { _nodes = 0; }

  // This sets the counter to the invalid counter.
  inline void invalidate() { _nodes = ~uint64_t(0); }

  // This method tells us if the value in the counter is valid.
  inline bool is_valid() const { return _nodes != ~uint64_t(0); }
  constexpr static UnderlyingCounterType produce_invalid_entry() { return ~uint64_t(0); }

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

  inline void reset() { base_counter.reset(); }

  inline bool is_valid() { return base_counter.is_valid(); }
  inline void invalidate() { return base_counter.invalidate(); }
  inline UnderlyingIndividualType get_total_nodes() const { return base_counter.get_total_nodes(); }
};

FPLLL_END_NAMESPACE
#endif
