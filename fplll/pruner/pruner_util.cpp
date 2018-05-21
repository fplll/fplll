#include "fplll.h"

FPLLL_BEGIN_NAMESPACE

template <class FT> void Pruner<FT>::print_coefficients(const vector<double> &b)
{
  cout << "# b = ";
  for (auto i = b.begin(); i != b.end(); ++i)
    std::cout << *i << ' ';
  cout << endl;
}

template <class FT> void Pruner<FT>::print_coefficients(const vec &b)
{
  cout << "# b = ";
  for (auto i = b.begin(); i != b.end(); ++i)
    std::cout << *i << ' ';
  cout << endl;
}

/**
    load gso
 */
template <class FT>
void Pruner<FT>::load_basis_shape(const vector<double> &gso_r, bool reset_normalization)
{
  shape_loaded = true;
  FT tmp;
  logvol = 0.0;
  r.resize(n);
  ipv.resize(n);
  r_old.resize(n);
  for (int i = 0; i < n; ++i)
  {
    r[i]     = gso_r[n - 1 - i];
    r_old[i] = gso_r[i];
    logvol += log(r[i]);
  }

  // by default reset
  if (reset_normalization)
  {
    /* just (1 / det^n)^2, the squared normalization factor */
    normalization_factor = exp(logvol / ((float)(-n)));

    /* re-normalise R = radius_d * (1 / sqr_det^n) */
    normalized_radius = sqrt(enumeration_radius * normalization_factor);
  }

  // r[i] is normalized by (1 / det^n)^2
  for (int i = 0; i < n; ++i)
  {
    r[i] *= normalization_factor;
  }

  tmp = 1.;
  // hence ipv is also normalized where ipv[n-1] = 1
  for (int i = 0; i < 2 * d; ++i)
  {
    tmp *= sqrt(r[i]);
    ipv[i] = 1.0 / tmp;
  }
}

template <class FT> void Pruner<FT>::load_basis_shapes(const vector<vector<double>> &gso_rs)
{
  n = gso_rs[0].size();
  vec sum_ipv(n);

  for (int i = 0; i < n; ++i)
  {
    sum_ipv[i] = 0.;
  }
  int count = gso_rs.size();
  for (int k = 0; k < count; ++k)
  {
    if (gso_rs[k].size() != static_cast<size_t>(n))
    {
      throw std::runtime_error("loading several bases with different dimensions");
    }
    load_basis_shape(gso_rs[k], (k == 0));
    for (int i = 0; i < n; ++i)
    {
      sum_ipv[i] += ipv[i];
    }
  }
  for (int i = 0; i < n; ++i)
  {
    ipv[i] = sum_ipv[i] / (1.0 * count);
  }
}

template <class FT> FT Pruner<FT>::gaussian_heuristic()
{
  return exp(2. * log(tabulated_ball_vol[n]) / ((float)-n)) / normalization_factor;
}

template <class FT>
void Pruner<FT>::save_coefficients(/*o*/ vector<double> &pr, /*i*/ const evec &b)
{
  pr.resize(n);
  int dn = b.size();
  if (dn == d)
  {
    for (int i = 0; i < d; ++i)
    {
      pr[n - 1 - 2 * i] = b[i].get_d();
      pr[n - 2 - 2 * i] = b[i].get_d();
    }
  }
  else
  {
    for (int i = 0; i < n; ++i)
    {
      pr[n - 1 - i] = b[i].get_d();
    }
  }
  pr[0] = 1.;
}

/**
 *   load pruning coefficients
 */
template <class FT> void Pruner<FT>::load_coefficients(/*o*/ vec &b, /*i*/ const vector<double> &pr)
{
  int dn = b.size();
  int c  = (dn == d) ? 2 : 1;
  for (int i = 0; i < dn; ++i)
  {
    b[i] = pr[n - c * i - 1];
  }
  if (enforce(b))
  {
    throw std::runtime_error(
        "Ill formed pruning coefficients (must be decreasing, starting with two 1.0)");
  }
}

FPLLL_END_NAMESPACE
