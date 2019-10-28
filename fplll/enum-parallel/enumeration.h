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

#ifndef ENUMLIB_ENUMERATION_HPP
#define ENUMLIB_ENUMERATION_HPP

#include "fplll_types.h"
#include "thread_pool.hpp"

#include <iostream>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include <memory>
#include <functional>

#include <fplll/defs.h>

FPLLL_BEGIN_NAMESPACE

namespace enumlib
{

//#define SINGLE_THREADED
//#define NOCOUNTS
#define NOLOCALUPDATE

using namespace std;

typedef fplll_float float_type;
typedef std::mutex mutex_type;
typedef std::lock_guard<std::mutex> lock_type;

typedef atomic<float_type> global_A_t;
typedef array<atomic_int_fast8_t, 256> global_signal_t;

template<int N>
struct globals_t {
	typedef array<int, N>        introw_t;
	typedef pair<introw_t, pair<float_type, float_type> > swirl_item_t;

	mutex_type      mutex;
	global_A_t      A;
	global_signal_t signal;

	std::function<extenum_cb_process_sol> process_sol;
	std::function<extenum_cb_process_subsol> process_subsol;

	vector< vector< swirl_item_t > > swirlys;
};


extern mutex_type global_mutex;

extern int enumlib_nrthreads;
extern int enumlib_loglevel;

extern thread_pool::thread_pool threadpool;


template <int N, int SWIRLY, int SWIRLY2BUF, int SWIRLY1FRACTION, bool findsubsols = false>
struct lattice_enum_t
{
	typedef array<float_type, N> fltrow_t;
	typedef array<int, N>        introw_t;
	typedef pair<introw_t, pair<float_type, float_type> > swirl_item_t;

	/* inputs */
	float_type muT[N][N];
	fltrow_t risq;
	fltrow_t pr;
	fltrow_t pr2;

	/* config */
	bool activeswirly;

	/* internals */
	int threadid;
	globals_t<N>& globals;
	
	float_type _A; // overall enumeration bound
	fltrow_t _AA, _AA2; // enumeration pruning bounds
	introw_t _x, _Dx, _D2x; 
	fltrow_t _sol; // to pass to fplll
	fltrow_t _c;
	introw_t _r;
	array<float_type, N + 1> _l;
	array<uint64_t, N + 1> _counts;

	float_type _sigT[N][N];
	
	fltrow_t _subsolL;
	array<fltrow_t, N> _subsol;
	
	std::chrono::system_clock::time_point starttime;

	lattice_enum_t(globals_t<N>& globals_)
		: globals(globals_)
		, starttime(std::chrono::system_clock::now())
		, activeswirly(false)
	{
	}

	inline int myround(double a)
	{
		return (int)(round(a));
	}
	inline int myround(float a)
	{
		return (int)(roundf(a));
	}
	inline int myround(long double a)
	{
		return (int)(roundl(a));
	}

	inline void _thread_local_update()
	{
		if (globals.signal[threadid])
		{
			globals.signal[threadid] = 0;
			_A = globals.A;
			_update_AA();
		}
	}

	inline void _update_AA()
	{
		for (int k = 0; k < N; ++k)
			_AA[k] = _A * pr[k];
		for (int k = 0; k < N; ++k)
			_AA2[k] = _A * pr2[k];
	}

	// compile time parameters for enumerate_recur (without ANY runtime overhead)
	// allows specialization for certain specific cases, e.g., i=0, or i=swirl
	template<int i, bool svp, int swirl, int swirlid> struct i_tag {};

	template<int i, bool svp, int swirl, int swirlid>
	inline void enumerate_recur(i_tag<i, svp, swirl, swirlid>)
	{
		if (_r[i] > _r[i - 1])
			_r[i - 1] = _r[i];
		float_type ci = _sigT[i][i];
		float_type yi = round(ci);
		int xi = (int)(yi);
		yi = ci - yi;
		float_type li = _l[i + 1] + (yi * yi * risq[i]);
#ifndef NOCOUNTS
		++_counts[i];
#endif

		if (findsubsols && li < _subsolL[i] && li != 0.0)
		{
			_subsolL[i] = li;
			_subsol[i][i] = xi;
			for (int j = i + 1; j < N; ++j)
				_subsol[i][j] = _x[j];
		}
		if (li > _AA[i])
			return;

		_Dx[i] = _D2x[i] = (((int)(yi >= 0) & 1) << 1) - 1;
		_c[i] = ci;
		_x[i] = xi;
		_l[i] = li;


		for (int j = _r[i - 1]; j > i - 1; --j)
			_sigT[i - 1][j - 1] = _sigT[i - 1][j] - _x[j] * muT[i - 1][j];

#ifndef NOLOCALUPDATE
		_thread_local_update();
#endif

		while (true)
		{
			enumerate_recur(i_tag<i - 1, svp, swirl, swirlid>());

			if (_l[i + 1] == 0.0)
			{
				++_x[i];
				_r[i - 1] = i;
				float_type yi2 = _c[i] - _x[i];
				float_type li2 = _l[i + 1] + (yi2 * yi2 * risq[i]);
				if (li2 > _AA2[i])
					return;
				_l[i] = li2;
				_sigT[i - 1][i - 1] = _sigT[i - 1][i] - _x[i] * muT[i - 1][i];
			}
			else
			{
				_x[i] += _Dx[i]; _D2x[i] = -_D2x[i]; _Dx[i] = _D2x[i] - _Dx[i];
				_r[i - 1] = i;
				float_type yi2 = _c[i] - _x[i];
				float_type li2 = _l[i + 1] + (yi2 * yi2 * risq[i]);
				if (li2 > _AA2[i])
					return;
				_l[i] = li2;
				_sigT[i - 1][i - 1] = _sigT[i - 1][i] - _x[i] * muT[i - 1][i];
			}
		}
	}

	template<bool svp, int swirl, int swirlid>
	inline void enumerate_recur(i_tag<0, svp, swirl, swirlid>)
	{
		static const int i = 0;
		float_type ci = _sigT[i][i];
		float_type yi = round(ci);
		int xi = (int)(yi);
		yi = ci - yi;
		float_type li = _l[i + 1] + (yi * yi * risq[i]);
#ifndef NOCOUNTS
		++_counts[i];
#endif

		if (findsubsols && li < _subsolL[i] && li != 0.0)
		{
			_subsolL[i] = li;
			_subsol[i][i] = xi;
			for (int j = i + 1; j < N; ++j)
				_subsol[i][j] = _x[j];
		}
		if (li > _AA[i])
			return;

		_Dx[i] = _D2x[i] = (((int)(yi >= 0) & 1) << 1) - 1;
		_c[i] = ci;
		_x[i] = xi;
		_l[i] = li;


#ifndef NOLOCALUPDATE
		_thread_local_update();
#endif

		while (true)
		{
			enumerate_recur(i_tag<i - 1, svp, swirl, swirlid>());

			if (_l[i + 1] == 0.0)
			{
				++_x[i];
				float_type yi2 = _c[i] - _x[i];
				float_type li2 = _l[i + 1] + (yi2 * yi2 * risq[i]);
				if (li2 > _AA2[i])
					return;
				_l[i] = li2;
			}
			else
			{
				_x[i] += _Dx[i]; _D2x[i] = -_D2x[i]; _Dx[i] = _D2x[i] - _Dx[i];
				float_type yi2 = _c[i] - _x[i];
				float_type li2 = _l[i + 1] + (yi2 * yi2 * risq[i]);
				if (li2 > _AA2[i])
					return;
				_l[i] = li2;
			}
		}
	}


	template<bool svp, int swirl, int swirlid>
	inline void enumerate_recur(i_tag<-1, svp, swirl, swirlid>)
	{
		if (_l[0] > _AA[0] || _l[0] == 0.0)
			return;

		lock_type lock(globals.mutex);
		
		for (int j = 0; j < N; ++j)
			_sol[j] = _x[j];
		globals.A = globals.process_sol(_l[0], &_sol[0]);

		// if it has changed then signal all threads to update and update ourselves
		if (globals.A != _A)
		{
			for (int j = 0; j < globals.signal.size(); ++j)
				globals.signal[j] = 1;
		
			_thread_local_update();
		}
	}

	template<bool svp, int swirl, int swirlid>
	inline void enumerate_recur(i_tag<N+1, svp, swirl, swirlid>)
	{
	}
	template<bool svp, int swirl, int swirlid>
	inline void enumerate_recur(i_tag<N+2, svp, swirl, swirlid>)
	{
	}

	template<int i, bool svp, int swirlid>
	inline void enumerate_recur(i_tag<i, svp, i, swirlid>)
	{
		if (_r[i] > _r[i - 1])
			_r[i - 1] = _r[i];

		float_type ci = _sigT[i][i];
		float_type yi = round(ci);
		int xi = (int)(yi);
		yi = ci - yi;
		float_type li = _l[i + 1] + (yi * yi * risq[i]);
#ifndef NOCOUNTS
		++_counts[i];
#endif

		if (findsubsols && li < _subsolL[i] && li != 0.0)
		{
			_subsolL[i] = li;
			_subsol[i][i] = xi;
			for (int j = i + 1; j < N; ++j)
				_subsol[i][j] = _x[j];
		}
		if (li > _AA[i])
			return;
		_c[i] = ci;
		_x[i] = xi;
		_l[i] = li;
		_Dx[i] = _D2x[i] = (((int)(yi >= 0) & 1) << 1) - 1;


		for (int j = _r[i - 1]; j > i - 1; --j)
			_sigT[i - 1][j - 1] = _sigT[i - 1][j] - _x[j] * muT[i - 1][j];

		while (true)
		{
			float_type ci2 = _sigT[i - 1][i - 1];
			int xi2 = myround(ci2);
			float_type yi2 = ci2 - xi2;
			float_type li2 = _l[i] + (yi2 * yi2 * risq[i - 1]);

			globals.swirlys[swirlid].emplace_back();
			for (int j = i; j < N; ++j)
				globals.swirlys[swirlid].back().first[j] = _x[j];
			globals.swirlys[swirlid].back().second.first = _l[i];
			globals.swirlys[swirlid].back().second.second = li2;

			if (_l[i + 1] == 0.0)
			{
				++_x[i];
				_r[i - 1] = i;
				float_type yi2 = _c[i] - _x[i];
				float_type li = _l[i + 1] + (yi2 * yi2 * risq[i]);
				if (li > _AA2[i])
					return;
				_l[i] = li;
				_sigT[i - 1][i - 1] = _sigT[i - 1][i] - _x[i] * muT[i - 1][i];
			}
			else
			{
				_x[i] += _Dx[i]; _D2x[i] = -_D2x[i]; _Dx[i] = _D2x[i] - _Dx[i];
				_r[i - 1] = i;
				float_type yi2 = _c[i] - _x[i];
				float_type li = _l[i + 1] + (yi2 * yi2 * risq[i]);
				if (li > _AA2[i])
					return;
				_l[i] = li;
				_sigT[i - 1][i - 1] = _sigT[i - 1][i] - _x[i] * muT[i - 1][i];
			}
		}

	}


	template<bool svp = true>
	void enumerate_recursive()
	{
		for (int i = 0; i < globals.signal.size(); ++i)
			globals.signal[i] = 0;
		threadid = enumlib_nrthreads;

		_A = globals.A;
		_update_AA();

		for (int j = 0; j < N; ++j)
		{
			_x[j] = _Dx[j] = 0; _D2x[j] = -1;
			_sol[j] = 0;
			_c[j] = _l[j] = 0.0;
			_subsolL[j] = risq[j];
			for (int k = 0; k < N; ++k)
			{
				_sigT[j][k] = 0.0;
				_subsol[j][k] = 0;
			}
			_r[j] = N - 1;
			_counts[j] = 0;
		}
		_l[N] = 0.0;
		_counts[N] = 0;

#ifdef SINGLE_THREADED
		enumerate_recur(i_tag<N-1, svp, -2, 0>());
#else
		auto& swirlys = globals.swirlys;
		swirlys.resize(2);
		swirlys[0].clear();
		enumerate_recur(i_tag<N - 1, svp, N - SWIRLY, 0>());
		if (enumlib_loglevel >= 1) cout << "[enumlib] Swirly0: #=" << swirlys[0].size() << endl;

		const auto swirl_less = [](const swirl_item_t& l, const swirl_item_t& r) { return l.second.second < r.second.second; };
		if (activeswirly)
		{
			sort(swirlys[0].begin(), swirlys[0].end(), swirl_less);
		}

		int swirly0idx = 0;
		swirlys[1].clear();
		while (swirly0idx < swirlys[0].size())
		{
			int swirly1newstart = (int)(swirlys[1].size());
			while (swirly0idx < swirlys[0].size() && swirlys[1].size() < SWIRLY2BUF)
			{
				const int i = N - SWIRLY;
				_x = swirlys[0][swirly0idx].first;
				_l[i] = swirlys[0][swirly0idx].second.first;
				for (int j = 0; j < N; ++j)
					_r[j] = N - 1;
				for (int j = N - 1; j > i - 1; --j)
					_sigT[i - 1][j - 1] = _sigT[i - 1][j] - _x[j] * muT[i - 1][j];

				enumerate_recur(i_tag<N - SWIRLY - 1, svp, N - 2 * SWIRLY, 1>());

				++swirly0idx;
			}
			if (enumlib_loglevel >= 1) cout << "[enumlib] Swirly1: #=" << swirlys[1].size() << " (" << swirly0idx << "/" << swirlys[0].size() << ")" << endl;

			int swirly1end = (int)(swirlys[1].size());
			if (activeswirly)
			{
				// sort the new additions to swirly1
				sort(swirlys[1].begin() + swirly1newstart, swirlys[1].end(), swirl_less);
				// merge with previous elms in swirly1
				inplace_merge(swirlys[1].begin(), swirlys[1].begin() + swirly1newstart, swirlys[1].end(), swirl_less);

				// process portion of swirly[1] and then add more
				swirly1end = (SWIRLY2BUF >> SWIRLY1FRACTION);
				if (swirly1end > (int)(swirlys[1].size()))
					swirly1end = (int)(swirlys[1].size());
			}

			auto& swirly_ref = swirlys[1];
			std::atomic<std::size_t> swirly_i(0);
			unsigned threadid = 0;
			auto f = [this, &swirly_ref, &swirly_i, swirly1end, &threadid]()
			{
				auto mylat = *this;
				{
					lock_type lock(globals.mutex);
					mylat.threadid = threadid++;
				}
				if (enumlib_loglevel >= 2) cout << "[enumlib] Thread " << mylat.threadid << " started." << endl;
				for (int j = 0; j < N - SWIRLY; ++j)
					mylat._counts[j] = 0;
				while (true)
				{
					std::size_t idx = swirly_i++;
					if (idx >= swirly1end)
						break;

					const int i = N - 2 * SWIRLY;
					mylat._x = swirly_ref[idx].first;
					mylat._l[i] = swirly_ref[idx].second.first;
					for (int j = 0; j < N; ++j)
						mylat._r[j] = N - 1;
					for (int j = N - 1; j > i - 1; --j)
						mylat._sigT[i - 1][j - 1] = mylat._sigT[i - 1][j] - mylat._x[j] * mylat.muT[i - 1][j];

					mylat._thread_local_update();

					mylat.enumerate_recur(i_tag<N - 2 * SWIRLY - 1, svp, -2, -1>());
				}

				lock_type lock(globals.mutex);
				for (int j = 0; j < N - SWIRLY; ++j)
					this->_counts[j] += mylat._counts[j];
				for (int j = 0; j < N; ++j)
					if (mylat._subsolL[j] < this->_subsolL[j])
					{
						this->_subsolL[j] = mylat._subsolL[j];
						this->_subsol[j] = mylat._subsol[j];
					}
			};
			if (threadpool.size() != enumlib_nrthreads)
				cout << "[enumlib] threadpool size mismatch!" << enumlib_nrthreads << "!=" << threadpool.size() << endl;
			for (int i = 0; i < enumlib_nrthreads; ++i)
				threadpool.push(f);
			threadpool.wait_work();

			swirlys[1].erase(swirlys[1].begin(), swirlys[1].begin() + swirly1end);
		}
#ifndef NOCOUNTS
//		if (enumlib_loglevel >= 1) cout << "[enumlib] counts: " << _counts << endl;
#endif
#endif
	}

};

} // namespace enumblib

FPLLL_END_NAMESPACE

#endif // ENUMLIB_ENUMERATION_HPP
