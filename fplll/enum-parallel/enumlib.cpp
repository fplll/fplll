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

#include "enumlib.h"
#include "enumeration.h"
#include "fplll_types.h"
#include <iostream>

FPLLL_BEGIN_NAMESPACE

namespace enumlib
{

using namespace std;

NodeCountType enumlib_enumerate(int dim, fplll::enumf maxdist,
                           std::function<fplll::extenum_cb_set_config> cbfunc,
                           std::function<fplll::extenum_cb_process_sol> cbsol,
                           std::function<fplll::extenum_cb_process_subsol> cbsubsol, bool dual,
                           bool findsubsols);

#define ENUMFUNCNAME(DIM)                                                                          \
    NodeCountType enumerate##DIM(int, float_type, std::function<extenum_cb_set_config>,                   \
                          std::function<extenum_cb_process_sol>,                                   \
                          std::function<extenum_cb_process_subsol>, bool, bool);

#if 20 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(20);
#endif
#if 30 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(30);
#endif
#if 40 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(40);
#endif
#if 50 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(50);
#endif
#if 60 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(60);
#endif
#if 70 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(70);
#endif
#if 80 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(80);
#endif
#if 90 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(90);
#endif
#if 100 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(100);
#endif
#if 110 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(110);
#endif
#if 120 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(120);
#endif
#if 130 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(130);
#endif
#if 140 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(140);
#endif
#if 150 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(150);
#endif
#if 160 <= FPLLL_MAX_PARALLEL_ENUM_DIM
ENUMFUNCNAME(160);
#endif

NodeCountType enumlib_enumerate(int dim, fplll_float maxdist,
                           std::function<extenum_cb_set_config> cbfunc,
                           std::function<extenum_cb_process_sol> cbsol,
                           std::function<extenum_cb_process_subsol> cbsubsol, bool dual,
                           bool findsubsols)
{
  // dual svp enumeration not supported yet
  if (dim <= 10 || dual)
    return ~uint64_t(0);

#if 20 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 20)
    return enumerate20(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 30 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 30)
    return enumerate30(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 40 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 40)
    return enumerate40(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 50 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 50)
    return enumerate50(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 60 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 60)
    return enumerate60(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 70 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 70)
    return enumerate70(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 80 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 80)
    return enumerate80(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 90 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 90)
    return enumerate90(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 100 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 100)
    return enumerate100(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 110 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 110)
    return enumerate110(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 120 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 120)
    return enumerate120(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 130 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 130)
    return enumerate130(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 140 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 140)
    return enumerate140(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 150 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 150)
    return enumerate150(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
#if 160 <= FPLLL_MAX_PARALLEL_ENUM_DIM
  if (dim <= 160)
    return enumerate160(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
#endif
  return ~uint64_t(0);
}

}  // namespace enumlib

FPLLL_END_NAMESPACE
