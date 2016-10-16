#ifndef FPLLL_SIEVE_GAUSS_STR_H
#define FPLLL_SIEVE_GAUSS_STR_H

#include "sieve_common.h"

#if 0
extern long dot_time;
extern long dot_num;
extern long count_bad;
#endif

/**
 * list point struct
 */
template <class ZT> class ListPoint
{

public:
  /* vector */
  NumVect<Z_NR<ZT>> v;

  /* square L2 norm of the vector */
  Z_NR<ZT> norm;
};

template <class ZT> inline ListPoint<ZT> *new_listpoint(int n)
{
  ListPoint<ZT> *p = new ListPoint<ZT>;
  p->norm          = 0;
  p->v.resize(n);
  for (int i = 0; i < n; ++i)
  {
    p->v[i] = 0;
  }
  return p;
}

/**
 * copy p1 to p2
 */
template <class ZT> void clone_listpoint(ListPoint<ZT> *p1, ListPoint<ZT> *p2)
{
  int n = p1->v.size();
  if (p2->v.size() != n)
    p2->v.resize(n);
  p2->norm = p1->norm;
  p2->v    = p1->v;
}

/**
 * set a of norm2 to p
 */
template <class ZT>
void set_listpoint_numvect(NumVect<Z_NR<ZT>> a, Z_NR<ZT> norm2, ListPoint<ZT> *p)
{
  int n = a.size();
  if (p->v.size() != n)
    p->v.resize(n);
  p->v    = a;
  p->norm = norm2;
}

template <class ZT> inline void del_listpoint(ListPoint<ZT> *p) { delete p; }

/**
 * reduce p1 w.r.t to p2
 * (TBA: optimize this function)
 */
template <class ZT> inline bool half_2reduce(ListPoint<ZT> *p1, const ListPoint<ZT> *p2)
{

#if 0
  struct timeval time;
  gettimeofday(&time, 0);
  long startt = 1000000 * time.tv_sec + time.tv_usec;
#endif
  /* debug info
  int debug = 0;
  if (debug) {
    cout << "[debug] start reduction " << endl;
    cout << "[debug]  p1.norm() = " << p1->norm << endl;
    cout << "[debug]  p2.norm() = " << p2->norm << endl;
    cout << "[debug]  p1 = " << p1->v << endl;
    cout << "[debug]  p2 = " << p2->v << endl;
  }
  */
  /* if 2*<p1,p2> <= |p2|^2, do nothing */
  Z_NR<ZT> dot, t, s;

#if 0
  struct timeval time;
  gettimeofday(&time, 0);
  long startt = 1000000 * time.tv_sec + time.tv_usec;
#endif

  dot_product(dot, p1->v, p2->v);

#if 0
  gettimeofday(&time, 0);
  long endt = 1000000 * time.tv_sec + time.tv_usec;
  dot_time += endt-startt;
  dot_num ++;
#endif

  t.mul_ui(dot, 2);
  t.abs(t);
  if (t <= p2->norm)
    return false;

  /* reduce p1 w.r.t to p2 (use dpe_t for now) */
  // FP_NR<dpe_t> t1, t2;
  FP_NR<double> t1, t2;
  t1.set_z(dot);
  t2.set_z(p2->norm);
  t1.div(t1, t2, GMP_RNDN);
  t1.rnd(t1);
  t.set_f(t1);
  NumVect<Z_NR<ZT>> vt((p1->v).size());
  vt.mul(p2->v, t);
  (p1->v).sub(vt);

  /* update new norm of |p1| */
  /* p1->norm = p1->norm + t * t * p2->norm - 2 * t * dot */
  s.mul(t, t);
  s.mul(s, p2->norm);
  s.add(s, p1->norm);
  t.mul_ui(t, 2);
  t.mul(t, dot);
  (p1->norm).sub(s, t);
/*
if (debug) {
  cout << "[debug]  p1new = " << p1->v << endl;
  cout << "[debug]  |p1new| = " << p1->norm << endl;
  cout << "[debug]  |v| = " << p2->norm << endl;
     if (p1->norm <= p2->norm && p1->norm != 0)
     exit(1);
}*/

#if 0
  gettimeofday(&time, 0);
  long endt = 1000000 * time.tv_sec + time.tv_usec;
  cout << "[info] red time " << endt-startt << endl;
#endif

  return true;
}

/***************************************************************/
/***************************************************************/

/**
 * check if p1 and p2 are 2-reduced for ordered p1 and p2
 */
template <class ZT>
inline bool check_2reduce_order(const ListPoint<ZT> *p1, const ListPoint<ZT> *p2)
{
  Z_NR<ZT> dot, t;
  dot_product(dot, p1->v, p2->v);
  t.mul_ui(dot, 2);
  t.abs(t);
  if (t > p2->norm)
    return 0;
  return 1;
}

/**
 * check if p1 and p2 are 2-reduced
 */
template <class ZT> inline bool check_2reduce(const ListPoint<ZT> *p1, const ListPoint<ZT> *p2)
{
  if (p1->norm <= p2->norm)
    return check_2reduce_order(p1, p2);
  else
    return check_2reduce_order(p2, p1);
}

/**
 * check if the list is pairwisely 2-reduced
 * the list needs to be ordered by norm which is ok in our context
 */
template <class ZT> inline bool check_2reduce_order_list(const list<ListPoint<ZT> *> List)
{
  typename list<ListPoint<ZT> *>::const_iterator i, j;
  ListPoint<ZT> *v1, *v2;
  for (i = List.begin(); i != List.end(); ++i)
  {
    v1 = *i;
    j  = i;
    advance(j, 1);
    for (; j != List.end(); ++j)
    {
      v2 = *j;
      if (!check_2reduce(v1, v2))
      {
        cout << v1->v << endl;
        cout << v2->v << endl;
        return 0;
      }
    }
  }
  return 1;
}

/***************************************************************/
/***************************************************************/

/**
 * check if p1, p2 and p3 are 3-reduced for ordered p1, p2, p3
 */
template <class ZT>
inline int check_3reduce_order(const ListPoint<ZT> *p1, const ListPoint<ZT> *p2,
                               const ListPoint<ZT> *p3, ListPoint<ZT> *pnew)
{
  /* check 2-reduced condition */
  if (!check_2reduce(p1, p2))
    return 0;
  if (!check_2reduce(p2, p3))
    return 0;
  if (!check_2reduce(p1, p3))
    return 0;

  /*
  cout << " checking p1->v " << p1->v << endl;
  cout << " checking p2->v " << p2->v << endl;
  cout << " checking p3->v " << p3->v << endl;
  cout << endl;
  */

  /* check 3-reduced condition */
  Z_NR<ZT> dot12, dot13, dot23;
  dot_product(dot12, p1->v, p2->v);
  dot_product(dot13, p1->v, p3->v);
  dot_product(dot23, p2->v, p3->v);

  if (dot12.sgn() * dot13.sgn() * dot23.sgn() != -1)
    return 1;
  else
  {
    NumVect<Z_NR<ZT>> a, b, c;
    Z_NR<ZT> t;
    a = p1->v;
    b = p2->v;
    c = p3->v;
    a.addmul_si(b, -dot12.sgn());
    a.addmul_si(c, -dot13.sgn());
    dot_product(t, a, a);
    if (t < p3->norm)
    {
      set_listpoint_numvect(a, t, pnew);

      /*
      cout << p3->v << endl;
      cout << p3->norm << endl;
      cout << " new is " << endl;
      cout << a << endl;
      cout << t << endl;
      */
      return -1;
    }
  }
  return 1;
}

/**
 * check if v1, v2 and v3 are 3-reduced
 */
template <class ZT>
inline int check_3reduce(const ListPoint<ZT> *p1, const ListPoint<ZT> *p2, const ListPoint<ZT> *p3,
                         ListPoint<ZT> *pnew)
{
  if (p1->norm <= p2->norm)
  {
    if (p2->norm <= p3->norm)
      return check_3reduce_order(p1, p2, p3, pnew);
    else if (p1->norm <= p3->norm)
      return check_3reduce_order(p1, p3, p2, pnew);
    else
      return check_3reduce_order(p3, p1, p2, pnew);
  }
  else
  {
    if (p1->norm <= p3->norm)
      return check_3reduce_order(p2, p1, p3, pnew);
    else if (p2->norm <= p3->norm)
      return check_3reduce_order(p2, p3, p1, pnew);
    else
      return check_3reduce_order(p3, p2, p1, pnew);
  }
}

/**
 * check if the list is 3-reduced
 * the list needs to be ordered by norm which is ok in our context
 */
template <class ZT> inline bool check_3reduce_order_list(const list<ListPoint<ZT> *> List)
{
  if (List.size() < 3)
    return 1;
  typename list<ListPoint<ZT> *>::const_iterator i, j, k;
  ListPoint<ZT> *v1, *v2, *v3;
  i                   = List.begin();
  v1                  = *i;
  ListPoint<ZT> *vnew = new_listpoint<ZT>(v1->v.size());
  for (i = List.begin(); i != List.end(); ++i)
  {
    v1 = *i;
    j  = i;
    for (; j != List.end(); ++j)
    {
      v2 = *j;
      k  = j;
      if (v1->norm >= v2->norm)
        continue;
      for (; k != List.end(); ++k)
      {
        v3 = *k;
        if (v2->norm >= v3->norm || v1->norm >= v3->norm)
          continue;
        if (check_3reduce(v1, v2, v3, vnew) != 1)
        {
          cout << " Error on vector: " << endl;
          cout << v1->v << endl;
          cout << v2->v << endl;
          cout << v3->v << endl;
#if 0
          count_bad ++;
          cout << " # bad pairs " << count_bad << endl;
#endif
          // del_listpoint(vnew);
          // return 0;
        }
      }
    }
  }
  del_listpoint(vnew);
  return 1;
}

/***************************************************************/
/***************************************************************/

/**
 * check if p1, p2, p3 and p4 are 4-reduced for ordered p1, p2, p3, p4
 */
template <class ZT>
inline int check_4reduce_order(const ListPoint<ZT> *p1, const ListPoint<ZT> *p2,
                               const ListPoint<ZT> *p3, const ListPoint<ZT> *p4,
                               ListPoint<ZT> *p4new)
{

  if ((p1->norm >= p2->norm) || (p2->norm >= p3->norm) || (p3->norm >= p4->norm))
  {
    cout << "error , please debug, this function can only used in order" << endl;
    exit(1);
  }

  /* check 2-reduced condition */
  if (!check_2reduce(p1, p2))
    return 0;
  if (!check_2reduce(p1, p3))
    return 0;
  if (!check_2reduce(p1, p4))
    return 0;
  if (!check_2reduce(p2, p3))
    return 0;
  if (!check_2reduce(p2, p4))
    return 0;
  if (!check_2reduce(p3, p4))
    return 0;

  /* check 3-reduced condition */
  ListPoint<ZT> *pnew = new_listpoint<ZT>(p1->v.size());
  if (!check_3reduce(p1, p2, p3, pnew))
  {
    del_listpoint(pnew);
    return 0;
  }
  if (!check_3reduce(p1, p2, p4, pnew))
  {
    del_listpoint(pnew);
    return 0;
  }
  if (!check_3reduce(p1, p3, p4, pnew))
  {
    del_listpoint(pnew);
    return 0;
  }
  if (!check_3reduce(p2, p3, p4, pnew))
  {
    del_listpoint(pnew);
    return 0;
  }
  del_listpoint(pnew);

  /*
  cout << " checking p1->v " << p1->v << endl;
  cout << " checking p2->v " << p2->v << endl;
  cout << " checking p3->v " << p3->v << endl;
  cout << " checking p4->v " << p4->v << endl;
  cout << endl;
  */

  /* check 4-reduced condition */
  ListPoint<ZT> *p4_update = new_listpoint<ZT>(p4->v.size());
  set_listpoint_numvect(p4->v, p4->norm, p4_update);
  int flag = 1;
  for (int i = -1; i <= 1; i += 2)
  {
    for (int j = -1; j <= 1; j += 2)
    {
      for (int k = -1; k <= 1; k += 2)
      {
        NumVect<Z_NR<ZT>> a, b, c, sum;
        Z_NR<ZT> t;
        a   = p1->v;
        b   = p2->v;
        c   = p3->v;
        sum = p4_update->v;
        sum.addmul_si(a, i);
        sum.addmul_si(b, j);
        sum.addmul_si(c, k);
        dot_product(t, sum, sum);
        if (t < p4_update->norm)
        {
          flag = 0;
          /*
            cout << "t:" << t << endl;
            cout << "sum:" << sum << endl;
            cout << "p4->norm:" << p4_update->norm << endl;
            cout << "p4->v:" << p4_update->v << endl;
          */
          set_listpoint_numvect(sum, t, p4_update);
          /*
            cout << "ne p4->v:" << p4_update->v << endl;
            cout << p4_update->v << endl;
            cout << p4_update->norm << endl;
            cout << " new is " << endl;
            cout << endl;
          */
        }
      }
    }
  }

  if (flag == 0)
  {
    set_listpoint_numvect(p4_update->v, p4_update->norm, p4new);
    del_listpoint(p4_update);
    return -1;
  }
  else
  {
    del_listpoint(p4_update);
    return 1;
  }
}

/**
 * check if v1, v2, v3 and v4 are 4-reduced
 */
template <class ZT>
inline int check_4reduce(ListPoint<ZT> *p1, ListPoint<ZT> *p2, ListPoint<ZT> *p3, ListPoint<ZT> *p4,
                         ListPoint<ZT> *p4new)
{
  NumVect<Z_NR<ZT>> index;
  index.push_back(p1->norm);
  index.push_back(p2->norm);
  index.push_back(p3->norm);
  index.push_back(p4->norm);
  std::sort(index.begin(), index.end());
  ListPoint<ZT> *s1 = NULL;
  ListPoint<ZT> *s2 = NULL;
  ListPoint<ZT> *s3 = NULL;
  ListPoint<ZT> *s4 = NULL;

  if (p1->norm == index[0])
    s1 = p1;
  else if (p2->norm == index[0])
    s1 = p2;
  else if (p3->norm == index[0])
    s1 = p3;
  else if (p4->norm == index[0])
    s1 = p4;

  if (p1->norm == index[1])
    s2 = p1;
  else if (p2->norm == index[1])
    s2 = p2;
  else if (p3->norm == index[1])
    s2 = p3;
  else if (p4->norm == index[1])
    s2 = p4;

  if (p1->norm == index[2])
    s3 = p1;
  else if (p2->norm == index[2])
    s3 = p2;
  else if (p3->norm == index[2])
    s3 = p3;
  else if (p4->norm == index[2])
    s3 = p4;

  if (p1->norm == index[3])
    s4 = p1;
  else if (p2->norm == index[3])
    s4 = p2;
  else if (p3->norm == index[3])
    s4 = p3;
  else if (p4->norm == index[3])
    s4 = p4;

  return check_4reduce_order(s1, s2, s3, s4, p4new);
}

/**
 * check if the list is 4-reduced
 * the list needs to be ordered by norm which is ok in our context
 */
template <class ZT> inline bool check_4reduce_order_list(const list<ListPoint<ZT> *> List)
{
  typename list<ListPoint<ZT> *>::const_iterator i, j, k, l;
  ListPoint<ZT> *v1, *v2, *v3, *v4;
  i                    = List.begin();
  v1                   = *i;
  ListPoint<ZT> *p4new = new_listpoint<ZT>(v1->v.size());
  for (i = List.begin(); i != List.end(); ++i)
  {
    v1 = *i;
    j  = i;
    for (; j != List.end(); ++j)
    {
      v2 = *j;
      k  = j;
      if (v1->norm >= v2->norm)
        continue;
      for (; k != List.end(); ++k)
      {
        v3 = *k;
        l  = k;
        if (v2->norm >= v3->norm || v1->norm >= v3->norm)
          continue;
        for (; l != List.end(); ++l)
        {
          v4 = *l;
          if (v3->norm >= v4->norm || v2->norm >= v4->norm || v1->norm >= v4->norm)
            continue;
          if (check_4reduce(v1, v2, v3, v4, p4new) != 1)
          {
            cout << " Error on vector: " << endl;
            cout << v1->v << endl;
            cout << v2->v << endl;
            cout << v3->v << endl;
            cout << v4->v << endl;
            del_listpoint(p4new);
            return 0;
          }
        }
      }
    }
  }
  del_listpoint(p4new);
  return 1;
}

/**
 * print current list
 */
template <class ZT> inline void print_list(const list<ListPoint<ZT> *> List)

{
  typename list<ListPoint<ZT> *>::const_iterator lp_it;
  for (lp_it = List.begin(); lp_it != List.end(); ++lp_it)
  {
    ///*
    cout << (*lp_it)->v << ", ";
    cout << (*lp_it)->norm << endl;
    //*/
    // cout << (*lp_it)->v << endl;
  }
}

/**
 * print current list
 */
template <class ZT> inline void check_0_list(const list<ListPoint<ZT> *> List)

{
  typename list<ListPoint<ZT> *>::const_iterator lp_it;
  for (lp_it = List.begin(); lp_it != List.end(); ++lp_it)
  {
    if ((*lp_it)->norm == 0)
    {
      cout << " error, list containing zero vector " << endl;
      exit(1);
    }
  }
}

/**
 * Use to convert MatrixRow to ListPoint
 */
template <class ZT>
inline void matrix_row_to_list_point(const MatrixRow<Z_NR<ZT>> &row, ListPoint<ZT> *p)
{
  int dims = row.size();
  //(p->v) (dims);
  p->v.resize(dims);
  p->norm = 0;
  Z_NR<ZT> t;
  for (int i = 0; i < dims; i++)
  {
    p->v[i] = row[i];
    t.mul(p->v[i], p->v[i]);
    p->norm.add(p->norm, t);
  }
}

/**
 * Use to convert sample() results NumVect to ListPoint
 */
template <class ZT> inline ListPoint<ZT> *num_vec_to_list_point(const NumVect<Z_NR<ZT>> &vec, int n)
{
  ListPoint<ZT> *p = new_listpoint<ZT>(n);
  int dims         = vec.size();
  p->v.resize(dims);
  p->norm = 0;
  Z_NR<ZT> t;
  for (int i = 0; i < dims; i++)
  {
    p->v[i] = vec[i];
    t.mul(p->v[i], p->v[i]);
    p->norm.add(p->norm, t);
  }
  return p;
}

template <class ZT> bool apply_filtering(const ListPoint<ZT> *p1, const ListPoint<ZT> *p2)
{
  Z_NR<ZT> dot;
  dot_product(dot, p1->v, p2->v);
  // cout << " dot is " << dot << endl;
  double t, t1, t2;
  t  = dot.get_d();
  t  = t * t;
  t1 = (p1->norm).get_d();
  t2 = (p2->norm).get_d();
  t  = t / t1;
  t  = t / t2;
  t  = abs(sqrt(t));
  if (t >= 1 / 3.0)
    return false;
  else
    return true;
}

#endif
