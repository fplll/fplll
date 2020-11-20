/*
MIT License

Copyright (c) 2016 Marc Stevens

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <fplll/fplll_config.h>

#ifdef ENUMDIMENSION
#if ENUMDIMENSION <= FPLLL_MAX_PARALLEL_ENUM_DIM

#include "enumeration.h"
#include "enumlib.h"

FPLLL_BEGIN_NAMESPACE

namespace enumlib
{

template <int dimension> struct enumerate_traits
{
  static const int SWIRLY          = 1 + (dimension / 20);
  static const int SWIRLY2BUF      = 1 << 10;
  static const int SWIRLY1FRACTION = 4;
};

template <int dimension, bool findsubsols>
array<uint64_t, FPLLL_MAX_ENUM_DIM>
enumerate_dim_detail(int dim, float_type maxdist,
                     std::function<extenum_cb_set_config> cb_set_config,
                     std::function<extenum_cb_process_sol> cb_process_sol,
                     std::function<extenum_cb_process_subsol> cb_process_subsol, bool dual)
{
  static const int SWIRLY          = enumerate_traits<dimension>::SWIRLY;
  static const int SWIRLY2BUF      = enumerate_traits<dimension>::SWIRLY2BUF;
  static const int SWIRLY1FRACTION = enumerate_traits<dimension>::SWIRLY1FRACTION;
  typedef lattice_enum_t<dimension, SWIRLY, SWIRLY2BUF, SWIRLY1FRACTION, findsubsols> lat_t;

  globals_t<dimension> globals;
  globals.A              = maxdist;
  globals.process_sol    = cb_process_sol;
  globals.process_subsol = cb_process_subsol;

  lat_t lat(globals);

  cb_set_config(&lat.muT[0][0], dimension, true, &lat.risq[0], &lat.pr[0]);
  lat.pr2 = lat.pr;

  lat.activeswirly = false;

  lat.enumerate_recursive();

  if (findsubsols)
  {
    for (int j = 0; j < dimension; ++j)
    {
      if (lat._subsolL[j] < lat.risq[j])
      {
        cb_process_subsol(lat._subsolL[j], &lat._subsol[j][0], j);
      }
    }
  }
  static_assert(FPLLL_MAX_ENUM_DIM > dimension,
                "The maximum enumeration dimension needs to be bigger than, or equal to, the "
                "parallel enumeration dim");
  array<uint64_t, FPLLL_MAX_ENUM_DIM> arr{};
  std::copy(lat._counts.cbegin(), lat._counts.cend(), arr.begin());
  return arr;
}

template <int dimension>
array<uint64_t, FPLLL_MAX_ENUM_DIM>
enumerate_dim(int dim, float_type maxdist, std::function<extenum_cb_set_config> cb_set_config,
              std::function<extenum_cb_process_sol> cb_process_sol,
              std::function<extenum_cb_process_subsol> cb_process_subsol, bool dual,
              bool findsubsols)
{
  if (findsubsols)
    return enumerate_dim_detail<dimension, true>(dim, maxdist, cb_set_config, cb_process_sol,
                                                 cb_process_subsol, dual);
  else
    return enumerate_dim_detail<dimension, false>(dim, maxdist, cb_set_config, cb_process_sol,
                                                  cb_process_subsol, dual);
}

#ifndef ENUMDIMENSION
#error "ENUMDIMENSION not defined"
#endif
#define DIMFUNCNAME(DIM) enumerate##DIM
#define DIMFUNC(DIM) DIMFUNCNAME(DIM)
#define GENENUM(d)                                                                                 \
  case (d):                                                                                        \
    return enumerate_dim<(d)>(dim, maxdist, cb_set_config, cb_process_sol, cb_process_subsol,      \
                              dual, findsubsols);                                                  \
    break;

array<uint64_t, FPLLL_MAX_ENUM_DIM> DIMFUNC(ENUMDIMENSION)(
    int dim, float_type maxdist, std::function<extenum_cb_set_config> cb_set_config,
    std::function<extenum_cb_process_sol> cb_process_sol,
    std::function<extenum_cb_process_subsol> cb_process_subsol, bool dual, bool findsubsols)
{
  static const int d = ENUMDIMENSION;
  switch (dim)
  {
    GENENUM(d - 9);
    GENENUM(d - 8);
    GENENUM(d - 7);
    GENENUM(d - 6);
    GENENUM(d - 5);
    GENENUM(d - 4);
    GENENUM(d - 3);
    GENENUM(d - 2);
    GENENUM(d - 1);
    GENENUM(d);
  }

  cout << "[enumlib] " << ENUMDIMENSION << ":" << dim << " wrong dimension!" << endl;
  array<uint64_t, FPLLL_MAX_ENUM_DIM> arr{};
  arr[0] = ~uint64_t(0);
  return arr;
}

}  // namespace enumlib

FPLLL_END_NAMESPACE

#endif
#endif
