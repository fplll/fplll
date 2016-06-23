#include "sieve_Gauss.h"
#include <iterator>

//#define DEBUG_CHECK_4RED


/**
 * auxiliary function, return the indicator of larger norms
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::update_p_4reduce_aux
(ListPoint<ZT>* p, typename list<ListPoint<ZT> *>::iterator &lp_it_k)
{
  typename list<ListPoint<ZT> *>::iterator lp_it;
  ListPoint<ZT>* v;
  for (lp_it = List.begin(); lp_it != List.end(); ++lp_it) {
    v = *lp_it;
    if (p->norm < v->norm)
      break;
  }
  lp_it_k = lp_it;
}


/**
 * 3-reduce the point p w.r.t to the list
 * Note p could be deleted !
 */
template<class ZT, class F>
Z_NR<ZT> Gauss_sieve<ZT, F>::update_p_4reduce_3reduce (ListPoint<ZT>* p)
{
  typename list<ListPoint<ZT> *>::iterator lp_it1, lp_it2,
    lp_it_k, tmp_lp_it;
  ListPoint<ZT> *v1, *v2;
  ListPoint<ZT> *vnew = new_listpoint<ZT>(nc);
  ListPoint<ZT> *vnew2;
  int count = 0;
  Z_NR<ZT> current_norm;
  bool loop = true;

  /*
    check_0_list(List);
    cout << "#  -- before inserting s1 " << p->v << endl;
    cout << "#  -- before inserting s1 " << p->norm << endl;
    print_list(List);
  */
  
  while(loop) {
    count ++;
    loop = false;
    lp_it_k = List.begin();

    /* now p and L are 2-reduced and lp_it_k is the larger-norm-borderline */
    current_norm = update_p_3reduce_2reduce (p, lp_it_k);

    /* p has been deleted in another function, need only to delete vnew */
    if (current_norm == 0) {
      del_listpoint(vnew);
      return current_norm;
    }

    /* ordered (v1, v2, p), 3-reduce p */
    for (lp_it1 = List.begin(); lp_it1 != lp_it_k && 
           lp_it1 != List.end(); ++lp_it1) {
      if (loop)
        break;
      v1 = *lp_it1;
      lp_it2 = lp_it1;
      for (; lp_it2 != lp_it_k && lp_it2 != List.end(); ++lp_it2) {
        v2 = *lp_it2;
        if (v1->norm >= v2->norm || v2->norm >= p->norm ||
            v1->norm >= p->norm)
          continue;
        if (check_3reduce(v1, v2, p, vnew)!=1) {
          clone_listpoint(vnew, p);
          loop = true;
          break;
        }
      }
    }
  }

  /*
  cout << "#  -- before inserting s2 " << p->v << endl;
  cout << "#  -- before inserting s2 " << p->norm << endl;
  check_0_list(List);
  print_list(List);
  */
  
  del_listpoint(vnew);
  if (p->norm == 0) {
    del_listpoint(p);
    Z_NR<ZT> t;
    t = 0;
    return t;
  }

  /* 3-reduce (v1, p, v2) or (p, v1, v2) */
  lp_it1 = List.begin();
  while (lp_it1 != List.end()) {
    v1 = *lp_it1;
    if (v1->norm == p->norm) {
      ++lp_it1; continue;
    }
    lp_it2 = lp_it_k;
    while (lp_it2 != List.end()) {
      tmp_lp_it = lp_it2;
      v2 = *lp_it2;
      if (v2->norm <= v1->norm || v2->norm <= p->norm) {
        ++lp_it2;
        if (lp_it2 == List.end())
          break;
        continue;
      }
      ++lp_it2;
      if (v1->norm < p->norm) {
        /*cout << "#   --- here 1 " << endl;
        cout << v1->norm << endl;
        cout << p->norm << endl;
        cout << v2->norm << endl;
        check_0_list(List);
        */
        vnew2 = new_listpoint<ZT>(nc);
        if (check_3reduce(v1, p, v2, vnew2)!=1) {
          if (tmp_lp_it == lp_it_k) lp_it_k ++;
          List.erase(tmp_lp_it);
          Queue.push(vnew2);
          del_listpoint(v2);
        }
        else {
          del_listpoint(vnew2);
        }
      }
      else {
        vnew2 = new_listpoint<ZT>(nc);
        if (check_3reduce(p, v1, v2, vnew2)!=1) {
          if (tmp_lp_it == lp_it_k) lp_it_k ++;
          List.erase(tmp_lp_it);
          del_listpoint(v2);
          Queue.push(vnew2);
        }
        else
          del_listpoint(vnew2);
      }
    }
    ++lp_it1;
  }
  
  return p->norm;
}


/**
 * 4-reduce the point p w.r.t to the list
 * after this, the list should be 4-reduced tuple-wisely.
 * note p could be deleted
 */
template<class ZT, class F>
Z_NR<ZT> Gauss_sieve<ZT, F>::update_p_4reduce (ListPoint<ZT>* p)
{
  typename list<ListPoint<ZT> *>::iterator lp_it1, lp_it2, lp_it3,
    lp_it_k, tmp_lp_it;
  ListPoint<ZT> *v1, *v2, *v3;
  ListPoint<ZT> *vnew = new_listpoint<ZT>(nc);
  ListPoint<ZT> *vnew2;
  int count = 0;
  Z_NR<ZT> current_norm;
  bool loop = true;

  while (loop) {
    count ++;
    loop = false;

    /* 3-reduce p w.r.t list */
    current_norm = update_p_4reduce_3reduce (p);

   /* Here p has been deleted in another function, need only to 
      delete vnew */
    if (current_norm == 0) {
      del_listpoint(vnew);
      return current_norm;
    }
    
    /* find indicator for larger norms */
    update_p_4reduce_aux (p, lp_it_k);

#if 0
    if (!check_3reduce_order_list<ZT>(List)) {
      cout << "!!! Failed [check 3-red]" << endl;
      exit(1);
    }
    else {
      cout << "[check 3-red] OK " << endl;
    }
#endif

    /* case (v1, v2, v3, p) when p has largest norm  */
    for (lp_it1 = List.begin(); lp_it1 != lp_it_k && 
           lp_it1 != List.end(); ++lp_it1) {
      if (loop) break;
      v1 = *lp_it1;
      lp_it2 = lp_it1;
      for (; lp_it2 != lp_it_k && lp_it2 != List.end(); ++lp_it2) {
        if (loop) break;
        v2 = *lp_it2;
        lp_it3 = lp_it2;
        if (v1->norm >= v2->norm || v2->norm >= p->norm ||
            v1->norm >= p->norm)
          continue;
        for (; lp_it3 != lp_it_k && lp_it3 != List.end(); ++lp_it3) {
          v3 = *lp_it3;
          if (v1->norm >= v3->norm || v2->norm >= v3->norm ||
              v3->norm >= p->norm )
            continue;
          if (check_4reduce(v1, v2, v3, p, vnew)!=1) {
            clone_listpoint(vnew, p);
            loop = true;
            break;
          }
        }
      }
    }
  }
  
  del_listpoint(vnew);
  List.insert(lp_it_k, p);
  lp_it_k --;
  
  /* 4-reduce (p, v1, v2, v3) or (v1, p, v2, v3) or (v1, v2, p, v3) */
  lp_it1 = List.begin();
  while (lp_it1 != List.end()) {
    v1 = *lp_it1;
    if (v1->norm == p->norm) {
      ++lp_it1; continue;
    }
    lp_it2 = List.begin();
    while (lp_it2 != List.end()) {
      v2 = *lp_it2;
      if ((v2->norm == p->norm) || (v2->norm == v1->norm)) {
        ++lp_it2; continue;
      }
      lp_it3 = lp_it_k;
      while (lp_it3 != List.end()) {
        tmp_lp_it = lp_it3;
        v3 = *lp_it3;
        if (v3->norm <= v1->norm || v3->norm <= v2->norm
            || v3->norm <= p->norm) {
          ++lp_it3;
          if (lp_it3 == List.end())
            break;
          continue;
        }
        ++lp_it3;
        /* (v1, p, v2, v3) or (v1, v2, p, v3) */
        if (v1->norm < p->norm) {
          /* (v1, p, v2, v3) */
          if (v2->norm > p->norm) {
            vnew2 = new_listpoint<ZT>(nc);
            if (check_4reduce(v1, p, v2, v3, vnew2)!=1) {
              List.erase(tmp_lp_it);
              Queue.push(vnew2);
              del_listpoint(v3);
            }
            else {
              del_listpoint(vnew2);
            }
          }
          /* (v1, v2, p, v3) */
          else {
            vnew2 = new_listpoint<ZT>(nc);
            if (check_4reduce(v1, v2, p, v3, vnew2)!=1) {
              List.erase(tmp_lp_it);
              Queue.push(vnew2);
              del_listpoint(v3);
            }
            else {
              del_listpoint(vnew2);
            }
          }
        }
        /* (p, v1, v2, v3) */
        else {
          vnew2 = new_listpoint<ZT>(nc);
          if (check_4reduce(p, v1, v2, v3, vnew2)!=1) {
            List.erase(tmp_lp_it);
            Queue.push(vnew2);
            del_listpoint(v3);
          }
          else {
            del_listpoint(vnew2);
          }
        }
      }
      ++lp_it2;
    }
    ++lp_it1;
  }
  
  return p->norm;
}


/**
 * 4-sieve
 */
template<class ZT, class F>
bool Gauss_sieve<ZT, F>::run_4sieve () {

  ListPoint<ZT>* current_point;
  NumVect<Z_NR<ZT> > vec(nc);
  Z_NR<ZT> current_norm;

  /* main iteration */
  while ( (best_sqr_norm > goal_sqr_norm) &&
          (collisions < mult * max_list_size + 200) )
  {
    iterations++;
    max_list_size = max(max_list_size, long(List.size()));

    if (Queue.empty()) {
      vec = Sampler->sample();
      current_point = NumVectToListPoint (vec, nc);
      samples++;
    }
    else {
      current_point = Queue.front();
      Queue.pop();
    }

    current_norm = update_p_4reduce(current_point);

#if 0
    if (!check_4reduce_order_list<ZT>(List)) {
      cout << "!!! Failed " << endl;
      //print_list(List);
      exit(1);
    }
    else {
      cout << "[check 4-red] OK " << endl;
      //print_list(List);
    }
#endif

    if (current_norm == 0)
      collisions++;
    if (current_norm > 0 && current_norm < best_sqr_norm) {
      best_sqr_norm = current_norm;
    }

#if 1
    print_curr_info();
#endif

    /* tuples of (iters, max_list_size) */
    iters_norm.push_back (best_sqr_norm);
    iters_ls.push_back (max_list_size);
  }

  /* finished main procedure, output some information */
  print_final_info();

#ifdef DEBUG_CHECK_4RED
  if (check_4reduce_order_list<ZT>(List))
    cout << "# [info] check 4-reduced OK" << endl;
  else
    cout << "# Error: check 4-reduced not OK" << endl;
#endif
  //print_list(List);

  if (best_sqr_norm > goal_sqr_norm)
    return false;
  return true;
}
