/* Copyright (C) 2019 Marc Stevens.

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

#include <fplll/threadpool.h>

FPLLL_BEGIN_NAMESPACE

thread_pool::thread_pool threadpool;

/* get and set number of threads in threadpool, both return the (new) number of threads */
int get_threads() { return threadpool.size() + 1; }

int set_threads(int th)
{
  if (th > int(std::thread::hardware_concurrency()) || th == -1)
    th = std::thread::hardware_concurrency();
  if (th < 1)
    th = 1;
  threadpool.resize(th - 1);
  return get_threads();
}

FPLLL_END_NAMESPACE
