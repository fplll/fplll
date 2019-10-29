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

#ifndef FPLLL_THREADPOOL_H
#define FPLLL_THREADPOOL_H

#include <fplll/defs.h>
#include <fplll/io/thread_pool.hpp>

FPLLL_BEGIN_NAMESPACE

/* import C++11 standard mutex, lock_guard<mutex>, unique_lock<mutex> for convenience */
typedef std::mutex mutex;
typedef std::lock_guard<std::mutex> lock_guard;
typedef std::unique_lock<std::mutex> unique_lock;

/* import thread_pool::barrier
	class barrier {
	public:
		barrier(size_t N);
		void wait(); // wait until N threads have called wait()
	}
*/
typedef thread_pool::barrier barrier;

/* fplll's threadpool

	Note that for N threads, the total threadpool consists of the main thread and N-1 pooled threads.
	Default use is to submit N jobs and call wait_work in the main thread, 
	which will cause the main thread to also start processing jobs.

	class threadpool {
	public:
	
		// simple way to add a job
		void push(const std::function<void()>& f);
		void push(std::function<void()>&& f);
	
		// advanced way to add job
		template<typename F, typename... Args>
		auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>;
		
		// process tasks with main thread & then wait until all threads are idle
		// DO NOT CALL FROM A JOB FUNCTION: WILL CAUSE DEADLOCK
		void wait_work();

	}
*/

extern thread_pool::thread_pool threadpool;

/* get and set number of threads in threadpool, both return the (new) number of threads */
int get_threads();
int set_threads(int th = -1); // -1 defaults number of threads to machine's number of cores

FPLLL_END_NAMESPACE

#endif
