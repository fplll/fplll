#include "sieve_Gauss.h"

//#define DEBUG_CHECK_2RED

/**
 * reduction function: reduces recuirsively the point with all the
 * points with smaller norm. Adds the point to the list if we 
 * don't have colission and puts to the queue all the points 
 * with bigger norm that can be reduced with it.
 */
template<class ZT, class F>
Z_NR<ZT> Gauss_sieve<ZT, F>::update_p_2reduce (ListPoint<ZT>* p)
{

#if 0
  struct timeval time;
  gettimeofday(&time, 0);
  long startt = 1000000 * time.tv_sec + time.tv_usec;
#endif

  typename list<ListPoint<ZT> *>::iterator lp_it, tmp_lp_it;
  ListPoint<ZT>* v;
  bool loop = true;

  /* 1. this loop should be stopping eventually hopefully */
  int count = 0;
  while(loop) {

    count ++;
    loop = false;
    
    /* if |p| >= |v_i| for any v_i in L, reduce p */
    for (lp_it = List.begin(); lp_it != List.end(); ++lp_it) {
      v = *lp_it;
      if (p->norm < v->norm)
        break;
      
      /* if there is one reduction the vector should re-pass the list */
      if (half_2reduce(p, v)) {
        loop = true;
      }
    }
  }


#if 0
  gettimeofday(&time, 0);
  long endt = 1000000 * time.tv_sec + time.tv_usec;
  cout << "[info] loop times " << count << endl;
  cout << "[info] updatelist1 time " << endt-startt << endl;
#endif

  /* 2. if collision, remove point */
  if (p->norm == 0) {
    del_listpoint(p);
    Z_NR<ZT> t;
    t = 0;
    return t;
  }

  /* 3. lp_it shows to the first point with bigger norm
     this is where we will insert the new point */
  List.insert(lp_it, p);

  /* 4. reduce List by p */
  while (lp_it != List.end()) {
    tmp_lp_it = lp_it;
    v = *lp_it;
    ++lp_it;
    if (half_2reduce(v, p)) {
      List.erase(tmp_lp_it);
      Queue.push(v);
    }
  }

#if 0
  gettimeofday(&time, 0);
  long endt2 = 1000000 * time.tv_sec + time.tv_usec;
  cout << "[info] updatelist2 time " << endt2-endt << endl;
#endif

  return p->norm;
}


/**
 * 2-sieve
 */
template<class ZT, class F>
bool Gauss_sieve<ZT, F>::run_2sieve () {

  ListPoint<ZT>* current_point;
  NumVect<Z_NR<ZT> > vec(nc);
  Z_NR<ZT> current_norm;

#ifdef REDUCE_TIMING
  struct timeval time, time2;
  gettimeofday(&time, 0);
  long startt = 1000000 * time.tv_sec + time.tv_usec;
  double sec, sec2 = 0.0;
  long nred = 0;
#endif
  
  /*
    Loop till you find a short enough vector,
    or enough collisions.
    The 0.1 * max_list_size_ + 200 is much more
    than needed in general so we are almost sure
    we have found the shortest vector. */
  while ( (best_sqr_norm > goal_sqr_norm) &&
          (collisions < mult * max_list_size + add) )
  {
    /* update stats */
    iterations++;
    max_list_size = max(max_list_size, long(List.size()));

    /* sample new or fetch from queue */
    if (Queue.empty()) {
      vec = Sampler->sample();
      current_point = NumVectToListPoint (vec, nc);
      samples++;
    }
    else {
      current_point = Queue.front();
      Queue.pop();
    }

#ifdef REDUCE_TIMING
    gettimeofday(&time2, 0);
    long startt2 = 1000000 * time2.tv_sec + time2.tv_usec;
    nred ++;
#endif

    /* sieve current_point */
    current_norm = update_p_2reduce(current_point);

#ifdef REDUCE_TIMING
    gettimeofday(&time2, 0);
    long endt2 = 1000000 * time2.tv_sec + time2.tv_usec;
    sec2 += (endt2-startt2) / 1000000.0;
#endif

    /* record collisions */
    if (current_norm == 0)
      collisions++;
    if (current_norm > 0 && current_norm < best_sqr_norm) {
      best_sqr_norm = current_norm;
    }


#if 1
    print_curr_info();
    /*if (samples+nr-List.size()-Queue.size() != collisions)
      exit(1);
      */
#endif

    /* tuples of (iters, max_list_size) */
    iters_norm.push_back (best_sqr_norm);
    iters_ls.push_back (max_list_size);
  }

#ifdef REDUCE_TIMING
  gettimeofday(&time, 0);
  long endt = 1000000 * time.tv_sec + time.tv_usec;
  sec = (endt-startt) / 1000000.0;
  cout << "# [info] total-reducuction time " << sec2;
  cout << ", average-reducuction time " << sec2/nred << " (" << nred << ")" << endl;
  cout << "# [info] total time " << sec << endl;
#endif

  /* finished main procedure, output some information */
  print_final_info();

#ifdef DEBUG_CHECK_2RED
  if (check_2reduce_order_list<ZT>(List))
    cout << "# [info] check 2-reduced OK" << endl;
  else
    cout << "# Error: check 2-reduced not OK" << endl;
#endif

  if (best_sqr_norm > goal_sqr_norm)
    return false;
  return true;
}
