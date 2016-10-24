#include "sieve_gauss.h"
#include <iterator>

//#define DEBUG_CHECK_3RED
#define EXTENSION_FILTERING

/**
 * auxiliary functionused in 3-reduction, make sure p and L is
 * pairwisely 2-reduced and return lp_it_k as indicator of larger norms
 */
template <class ZT, class F>
Z_NR<ZT>
GaussSieve<ZT, F>::update_p_3reduce_2reduce(ListPoint<ZT> *p,
                                            typename list<ListPoint<ZT> *>::iterator &lp_it_k)
{
  typename list<ListPoint<ZT> *>::iterator lp_it, tmp_lp_it;
  ListPoint<ZT> *v;
  bool loop = true;
  int count = 0;

  /* iteratively reduce p w.r.t to v \in L for every
   * v such that v.norm() < p.norm() */
  while (loop)
  {
    count++;
    loop = false;
    for (lp_it = List.begin(); lp_it != List.end(); ++lp_it)
    {
      v = *lp_it;
      if (p->norm < v->norm)
        break;
      if (half_2reduce(p, v))
      {
        loop = true;
      }
    }
  }

  if (p->norm == 0)
  {
    del_listpoint(p);
    Z_NR<ZT> t;
    t = 0;
    return t;
  }

  /* record position after which the v would be larger than p */
  lp_it_k = lp_it;

  /* now every v after lp_it is larger than p, we should try to reduce
   * v instead */
  while (lp_it != List.end())
  {
    tmp_lp_it = lp_it;
    v         = *tmp_lp_it;
    ++lp_it;

    /* corner case, do not remove lp_it_k as it is needed to be passed
     * back to the other function */
    if (half_2reduce(v, p))
    {
      if (tmp_lp_it == lp_it_k)
      {
        lp_it_k++;
      }
      List.erase(tmp_lp_it);
      Queue.push(v);
    }
  }

  return p->norm;
}

/**
 * 3-reduce the point p w.r.t to the list
 */
template <class ZT, class F> Z_NR<ZT> GaussSieve<ZT, F>::update_p_3reduce(ListPoint<ZT> *p)
{
  typename list<ListPoint<ZT> *>::iterator lp_it1, lp_it2, lp_it_k, tmp_lp_it;
  ListPoint<ZT> *v1, *v2;
  ListPoint<ZT> *vnew = new_listpoint<ZT>(nc);
  ListPoint<ZT> *vnew2;
  int count = 0;
  Z_NR<ZT> current_norm;
  bool loop = true;

  while (loop)
  {
    count++;
    loop    = false;
    lp_it_k = List.begin();

    /* now p and L are 2-reduced and lp_it_k is the larger-norm-borderline */
    current_norm = update_p_3reduce_2reduce(p, lp_it_k);

    /* p has been deleted in another function, need only to delete vnew */
    if (current_norm == 0)
    {
      del_listpoint(vnew);
      return current_norm;
    }

    /* ordered (v1, v2, p), 3-reduce p */
    for (lp_it1 = List.begin(); lp_it1 != lp_it_k && lp_it1 != List.end(); ++lp_it1)
    {
      if (loop)
        break;
      v1     = *lp_it1;
      lp_it2 = lp_it1;

#ifdef EXTENSION_FILTERING
      if (apply_filtering(p, v1))
        continue;
#endif

      for (; lp_it2 != lp_it_k && lp_it2 != List.end(); ++lp_it2)
      {
        v2 = *lp_it2;
        if (v1->norm >= v2->norm || v2->norm >= p->norm || v1->norm >= p->norm)
          continue;
        if (check_3reduce(v1, v2, p, vnew) != 1)
        {
          clone_listpoint(vnew, p);
          loop = true;
          break;
        }
      }
    }
  }
  del_listpoint(vnew);
  List.insert(lp_it_k, p);
  /* this is important to prevent from it to be deleted in
     the later operation List.erase(tmp_lp_it) */
  lp_it_k--;

  /* 3-reduce (v1, p, v2) or (p, v1, v2) */
  lp_it1 = List.begin();

  while (lp_it1 != List.end())
  {
    v1 = *lp_it1;
    if (v1->norm == p->norm)
    {
      ++lp_it1;
      continue;
    }

#ifdef EXTENSION_FILTERING
    if (apply_filtering(p, v1))
    {
      ++lp_it1;
      continue;
    }
#endif

    lp_it2 = lp_it_k;
    while (lp_it2 != List.end())
    {
      tmp_lp_it = lp_it2;
      v2        = *lp_it2;
      if (v2->norm <= v1->norm || v2->norm <= p->norm)
      {
        ++lp_it2;
        if (lp_it2 == List.end())
          break;
        continue;
      }
      ++lp_it2;
      if (v1->norm < p->norm)
      {
        vnew2 = new_listpoint<ZT>(nc);
        if (check_3reduce(v1, p, v2, vnew2) != 1)
        {
          List.erase(tmp_lp_it);
          Queue.push(vnew2);
          del_listpoint(v2);
        }
        else
        {
          del_listpoint(vnew2);
        }
      }
      else
      {
        vnew2 = new_listpoint<ZT>(nc);
        if (check_3reduce(p, v1, v2, vnew2) != 1)
        {
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
 * 3-sieve
 */
template <class ZT, class F> bool GaussSieve<ZT, F>::run_3sieve()
{

  ListPoint<ZT> *current_point;
  NumVect<Z_NR<ZT>> vec(nc);
  Z_NR<ZT> current_norm;

  /* main iteration */
  while ((best_sqr_norm > goal_sqr_norm) && (collisions < mult * max_list_size + add))
  {
    iterations++;
    max_list_size = max(max_list_size, long(List.size()));

    if (Queue.empty())
    {
      vec           = Sampler->sample();
      current_point = num_vec_to_list_point(vec, nc);
      samples++;
    }
    else
    {
      current_point = Queue.front();
      Queue.pop();
    }

    current_norm = update_p_3reduce(current_point);

#if 0
    if (!check_3reduce_order_list<ZT>(List)) {
      cout << "!!! Failed " << endl;
      print_list(List);
      exit(1);
    }
    else {
      cout << "[check] OK " << endl;
      cout << endl;
      print_list(List);
      
    }
#endif

    if (current_norm == 0)
      collisions++;
    if (current_norm > 0 && current_norm < best_sqr_norm)
    {
      best_sqr_norm = current_norm;
    }

#if 1
    print_curr_info();
#endif

    /* tuples of (iters, max_list_size) */
    iters_norm.push_back(best_sqr_norm);
    iters_ls.push_back(max_list_size);
  }

  /* finished main procedure, output some information */
  print_final_info();

#ifdef DEBUG_CHECK_3RED
  if (check_3reduce_order_list<ZT>(List))
    cout << "# [info] check 3-reduced OK" << endl;
  else
    cout << "# Error: check 3-reduced not OK" << endl;
#endif
  // print_list(List);

  if (best_sqr_norm > goal_sqr_norm)
    return false;
  return true;
}
