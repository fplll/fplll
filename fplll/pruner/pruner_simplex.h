/**
 *  Inline functions for computing the volume of even simplex (used in
 *  pruner_prob.cpp and pruner_cost.cpp.
 */
template <class FT> inline FT Pruner<FT>::eval_poly(const int ld, /*i*/ const poly &p, const FT x)
{
  FT acc;
  acc = 0.0;
  for (int i = ld; i >= 0; --i)
  {
    acc = acc * x;
    acc = acc + p[i];
  }
  return acc;
}

template <class FT> inline void Pruner<FT>::integrate_poly(const int ld, /*io*/ poly &p)
{
  for (int i = ld; i >= 0; --i)
  {
    FT tmp;
    tmp      = i + 1.;
    p[i + 1] = p[i] / tmp;
  }
  p[0] = 0.0;
}

/**
 * volume of even simplex
 */
template <class FT>
inline FT Pruner<FT>::relative_volume(const int rd,
                                      /*i*/ const evec &b)
{
  poly P(rd + 1);
  P[0]   = 1;
  int ld = 0;
  for (int i = rd - 1; i >= 0; --i)
  {
    integrate_poly(ld, P);
    ld++;
    P[0] = -1.0 * eval_poly(ld, P, b[i] / b[rd - 1]);
  }
  FT res = P[0] * tabulated_factorial[rd];
  return (rd % 2) ? -res : res;
}
