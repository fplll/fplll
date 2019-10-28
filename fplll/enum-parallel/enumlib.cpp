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

thread_pool::thread_pool threadpool(std::thread::hardware_concurrency());

int enumlib_nrthreads = std::thread::hardware_concurrency();
int enumlib_loglevel  = 0;

uint64_t enumlib_enumerate(int dim, fplll::enumf maxdist,
                           std::function<fplll::extenum_cb_set_config> cbfunc,
                           std::function<fplll::extenum_cb_process_sol> cbsol,
                           std::function<fplll::extenum_cb_process_subsol> cbsubsol, bool dual,
                           bool findsubsols);

extern "C"
{

  void fplll_enumlib_set_loglevel(int level)
  {
    if (level < 0)
      level = -1;
    if (level > 2)
      level = 2;
    enumlib_loglevel               = level;
    static const char *levelstrs[] = {"quiet", "normal", "verbose", "very verbose"};
    cout << "[enumlib] setting verbose level to " << levelstrs[level + 1] << "." << endl;
  }

  void fplll_enumlib_set_numthreads(int th)
  {
    if (th <= 0 || th > (int)(std::thread::hardware_concurrency()))
      th = std::thread::hardware_concurrency();
    cout << "[enumlib] setting number of threads to " << th << "." << endl;
    enumlib_nrthreads = th;
    threadpool.resize(th);
  }

  void fplll_register_enumlib() { fplll::set_external_enumerator(enumlib_enumerate); }

}  // extern "C"

#define ENUMFUNCNAME(DIM)                                                                          \
  uint64_t enumerate##DIM(int, float_type, std::function<extenum_cb_set_config>,                   \
                          std::function<extenum_cb_process_sol>,                                   \
                          std::function<extenum_cb_process_subsol>, bool, bool);

ENUMFUNCNAME(20);
ENUMFUNCNAME(30);
ENUMFUNCNAME(40);
ENUMFUNCNAME(50);
ENUMFUNCNAME(60);
ENUMFUNCNAME(70);
ENUMFUNCNAME(80);
ENUMFUNCNAME(90);
ENUMFUNCNAME(100);
ENUMFUNCNAME(110);
ENUMFUNCNAME(120);
ENUMFUNCNAME(130);
ENUMFUNCNAME(140);
ENUMFUNCNAME(150);
ENUMFUNCNAME(160);

uint64_t enumlib_enumerate(int dim, fplll_float maxdist,
                           std::function<extenum_cb_set_config> cbfunc,
                           std::function<extenum_cb_process_sol> cbsol,
                           std::function<extenum_cb_process_subsol> cbsubsol, bool dual,
                           bool findsubsols)
{
  // dual svp enumeration not supported yet
  if (dim <= 10 || dual)
    return ~uint64_t(0);

  if (dim <= 20)
    return enumerate20(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 30)
    return enumerate30(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 40)
    return enumerate40(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 50)
    return enumerate50(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 60)
    return enumerate60(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 70)
    return enumerate70(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 80)
    return enumerate80(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 90)
    return enumerate90(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 100)
    return enumerate100(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 110)
    return enumerate110(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 120)
    return enumerate120(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 130)
    return enumerate130(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 140)
    return enumerate140(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 150)
    return enumerate150(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  if (dim <= 160)
    return enumerate160(dim, maxdist, cbfunc, cbsol, cbsubsol, dual, findsubsols);
  return ~uint64_t(0);
}

}  // namespace enumlib

FPLLL_END_NAMESPACE
