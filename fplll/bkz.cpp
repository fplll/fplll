/* Copyright (C) 2011 Xavier Pujol
   (C) 2014-2016 Martin R. Albrecht.

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

#include <iomanip>

#include "bkz.h"
#include "bkz_params.h"
#include "enum/enumerate.h"
#include <iomanip>

FPLLL_BEGIN_NAMESPACE

template <class FT>
BKZReduction<FT>::BKZReduction(MatGSO<Integer, FT> &m, LLLReduction<Integer, FT> &lll_obj,
                               const BKZParam &param)
    : status(RED_SUCCESS), param(param), m(m), lll_obj(lll_obj)
{
  for (num_rows = m.d; num_rows > 0 && m.b[num_rows - 1].is_zero(); num_rows--)
  {
  }
  this->delta = param.delta;
}

template <class FT> BKZReduction<FT>::~BKZReduction() {}

template <class FT> double get_current_slope(MatGSO<Integer, FT> &m, int start_row, int stop_row)
{
  FT f, log_F;
  long expo;
  vector<double> x;
  x.resize(stop_row);
  for (int i = start_row; i < stop_row; i++)
  {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    log_F.log(f, GMP_RNDU);
    x[i] = log_F.get_d() + expo * std::log(2.0);
  }
  int n         = stop_row - start_row;
  double i_mean = (n - 1) * 0.5 + start_row, x_mean = 0, v1 = 0, v2 = 0;
  for (int i = start_row; i < stop_row; i++)
  {
    x_mean += x[i];
  }
  x_mean /= n;
  for (int i = start_row; i < stop_row; i++)
  {
    v1 += (i - i_mean) * (x[i] - x_mean);
    v2 += (i - i_mean) * (i - i_mean);
  }
  return v1 / v2;
}

template<class FT>
FT get_root_det(MatGSO<Integer, FT>& m, int start, int end) {
  FT root_det = 0;
  start = max(0, start);
  end = min (m.d, end);
  FT f,h;
  for(int i = start; i < end; ++i)
  {
      m.getR(h,i,i);
      h.log(h);
      root_det += h;
  }
  h = (double)(end - start);
  root_det /= h;
  root_det.exponential(root_det);
  return root_det;
}

template <class FT>
void compute_gauss_heur_dist(FT &root_det, FT &max_dist, long max_dist_expo, int kappa,
                             int block_size, double gh_factor)
{
  double t = (double)block_size / 2.0 + 1;
  t        = tgamma(t);
  t        = pow(t, 2.0 / (double)block_size);
  t        = t / M_PI;
  FT f = t;
  f = f * root_det;
  f.mul_2si(f, -max_dist_expo);
  f = f * gh_factor;
  if (f < max_dist)
  {
    max_dist = f;
  }
}

/** Randomize basis between from ``min_row`` and ``max_row`` (exclusive)

    1. permute rows
    2. apply lower triangular matrix with coefficients in -1,0,1
    3. LLL reduce result

    @param min_row start in this row

    @param max_row stop at this row (exclusive)

    @param density number of non-zero coefficients in lower triangular
    transformation matrix
**/

template <class FT> void BKZReduction<FT>::rerandomize_block(int min_row, int max_row, int density)
{
  if (max_row - min_row < 2)
    return;

  // 1. permute rows
  size_t niter = 4 * (max_row - min_row);  // some guestimate

  for (size_t i = 0; i < niter; ++i)
  {
    size_t a = gmp_urandomm_ui(RandGen::getGMPState(), max_row - min_row - 1) + min_row;
    size_t b = a;
    while (b == a)
    {
      b = gmp_urandomm_ui(RandGen::getGMPState(), max_row - min_row - 1) + min_row;
    }
    m.moveRow(b, a);
  }

  // 2. triangular transformation matrix with coefficients in -1,0,1
  m.rowOpBegin(min_row, max_row);
  FT x;
  for (long a = min_row; a < max_row - 2; ++a)
  {
    for (long i = 0; i < density; i++)
    {
      size_t b = gmp_urandomm_ui(RandGen::getGMPState(), max_row - (a + 1) - 1) + a + 1;
      if (gmp_urandomm_ui(RandGen::getGMPState(), 2))
        m.row_add(a, b);
      else
        m.row_sub(a, b);
    }
  }
  m.rowOpEnd(min_row, max_row);

  // 3. LLL reduce
  if (!lll_obj.lll(0, min_row, max_row))
    throw lll_obj.status;
  return;
}

template <class FT>
const Pruning &BKZReduction<FT>::get_pruning(int kappa, int block_size, const BKZParam &par) const
{

  FPLLL_DEBUG_CHECK(param.strategies.size() > blockSize);

  Strategy &strat = par.strategies[block_size];

  long max_dist_expo;
  FT max_dist = m.getRExp(kappa, kappa, max_dist_expo);

  FT gh_max_dist;
  FT root_det = get_root_det(m, kappa, kappa + block_size);
  compute_gauss_heur_dist(root_det, gh_max_dist, max_dist_expo, kappa, block_size, 1.0);
  return strat.get_pruning(max_dist.get_d() * pow(2, max_dist_expo),
                           gh_max_dist.get_d() * pow(2, max_dist_expo));
}

template <class FT>
bool BKZReduction<FT>::svp_postprocessing(int kappa, int block_size, const vector<FT> &solution)
{
  // Is it already in the basis ?
  int nz_vectors = 0, i_vector = -1;
  for (int i = 0; i < block_size; i++)
  {
    if (!solution[i].is_zero())
    {
      nz_vectors++;
      if (i_vector == -1 && fabs(solution[i].get_d()) == 1)
        i_vector = i;
    }
  }

  FPLLL_DEBUG_CHECK(nzVectors > 0);

  if (nz_vectors == 1)
  {
    // Yes, it is another vector
    FPLLL_DEBUG_CHECK(i_vector != -1 && i_vector != 0);
    m.moveRow(kappa + i_vector, kappa);
    if (!lll_obj.sizeReduction(kappa, kappa + 1))
      throw lll_obj.status;
  }
  else
  {
    // No, general case
    int d = m.d;
    m.createRow();
    m.rowOpBegin(d, d + 1);
    for (int i = 0; i < block_size; i++)
    {
      m.row_addmul(d, kappa + i, solution[i]);
    }
    m.rowOpEnd(d, d + 1);
    m.moveRow(d, kappa);
    if (!lll_obj.lll(kappa, kappa, kappa + block_size + 1))
      throw lll_obj.status;
    FPLLL_DEBUG_CHECK(m.b[kappa + block_size].is_zero());
    m.moveRow(kappa + block_size, d);
    m.removeLastRow();
  }
  return false;
}

template <class FT>
bool BKZReduction<FT>::svp_reduction(int kappa, int block_size, const BKZParam &par)
{
  bool clean = true;

  int lll_start = (par.flags & BKZ_BOUNDED_LLL) ? kappa : 0;

  if (!lll_obj.lll(lll_start, kappa, kappa + block_size))
  {
    throw std::runtime_error(RED_STATUS_STR[lll_obj.status]);
  }

  clean &= (lll_obj.nSwaps == 0);

  size_t trial                 = 0;
  double remaining_probability = 1.0;

  while (remaining_probability > 1. - par.min_success_probability)
  {
    if (trial > 0)
    {
      rerandomize_block(kappa + 1, kappa + block_size, par.rerandomization_density);
    }

    clean &= svp_preprocessing(kappa, block_size, par);

    long max_dist_expo;
    FT max_dist = m.getRExp(kappa, kappa, max_dist_expo);
    FT delta_max_dist;
    delta_max_dist.mul(delta, max_dist);

    if ((par.flags & BKZ_GH_BND) && block_size > 30)
    {
      FT root_det = get_root_det(m, kappa, kappa + block_size);
      compute_gauss_heur_dist(root_det, max_dist, max_dist_expo, kappa, block_size, par.gh_factor);
    }

    const Pruning &pruning = get_pruning(kappa, block_size, par);

    vector<FT> &solCoord = evaluator.solCoord;
    solCoord.clear();
    Enumeration::enumerate(m, max_dist, max_dist_expo, evaluator, empty_sub_tree, empty_sub_tree,
                           kappa, kappa + block_size, pruning.coefficients);

    nodes += Enumeration::getNodes();

    if (solCoord.empty())
    {
      if (pruning.coefficients[0] == 1 && !(par.flags & BKZ_GH_BND))
      {
        throw std::runtime_error(RED_STATUS_STR[RED_ENUM_FAILURE]);
      }
    }

    if (max_dist < delta_max_dist)
    {
      clean &= svp_postprocessing(kappa, block_size, solCoord);
    }
    remaining_probability *= (1 - pruning.probability);
    trial += 1;
  }
  return clean;
}

template <class FT>
bool BKZReduction<FT>::tour(const int loop, int &kappa_max, const BKZParam &par, int min_row,
                            int max_row)
{
  bool clean = true;
  clean &= trunc_tour(kappa_max, par,  min_row, max_row);
  clean &= hkz(kappa_max, par, max_row - par.block_size, max_row);
  
  if (par.flags & BKZ_VERBOSE)
  {
    print_tour(loop, min_row, max_row);
  }

  if (par.flags & BKZ_DUMP_GSO)
  {
    std::ostringstream prefix;
    prefix << "End of BKZ loop " << std::setw(4) << loop;
    prefix << " (" << std::fixed << std::setw(9) << std::setprecision(3)
           << (cputime() - cputime_start) * 0.001 << "s)";
    dump_gso(par.dump_gso_filename, prefix.str());
  }
  
  return clean;
}

template <class FT>
bool BKZReduction<FT>::trunc_tour(int &kappa_max, const BKZParam &par, int min_row,
                            int max_row)
{
  bool clean = true;
  int block_size = par.block_size;
  for (int kappa = min_row; kappa < max_row - block_size; ++kappa)
  {
    clean &= svp_reduction(kappa, block_size, par);
    if ((par.flags & BKZ_VERBOSE) && kappa_max < kappa && clean)
    {
      cerr << "Block [1-" << setw(4) << kappa + 1 << "] BKZ-" << setw(0) << par.block_size
           << " reduced for the first time" << endl;
      kappa_max = kappa;
    }
  }

  return clean;
}

template <class FT>
bool BKZReduction<FT>::trunc_dtour(const BKZParam &par, int min_row,
                            int max_row)
{
  bool clean = true;
  int block_size = par.block_size;
  
  for (int kappa = max_row - block_size; kappa > 0; --kappa)
  {
    clean &= dsvp_reduction(kappa, block_size, par);
  }

  return clean;
}

template <class FT>
bool BKZReduction<FT>::hkz(int &kappa_max, const BKZParam &param, int min_row, int max_row)
{
  bool clean = true;
  for (int kappa = min_row; kappa < max_row - 1; ++kappa) {
    int block_size = max_row - kappa;
    clean &= svp_reduction(kappa, block_size, param);
    if ((param.flags & BKZ_VERBOSE) && kappa_max < kappa && clean)
    {
      cerr << "Block [1-" << setw(4) << kappa + 1 << "] BKZ-" << setw(0) << param.block_size
           << " reduced for the first time" << endl;
      kappa_max = kappa;
    }
  }
  cerr << endl;
  
  return clean;
}

template <class FT>
bool BKZReduction<FT>::sd_tour(const int loop, const BKZParam &par, int min_row,
                            int max_row)
{
  int dummy_kappa_max = num_rows;
  bool clean = true;
  clean &= trunc_dtour(par, min_row, max_row);
  clean &= trunc_tour(dummy_kappa_max, par, min_row, max_row);
  
  if (par.flags & BKZ_VERBOSE)
  {
    print_tour(loop, min_row, max_row);
  }

  if (par.flags & BKZ_DUMP_GSO)
  {
    std::ostringstream prefix;
    prefix << "End of SD-BKZ loop " << std::setw(4) << loop;
    prefix << " (" << std::fixed << std::setw(9) << std::setprecision(3)
           << (cputime() - cputime_start) * 0.001 << "s)";
    dump_gso(par.dump_gso_filename, prefix.str());
  }

  return clean;
}

template <class FT> bool BKZReduction<FT>::bkz()
{
  int flags        = param.flags;
  int final_status = RED_SUCCESS;
  nodes            = 0;

  if (flags & BKZ_DUMP_GSO)
  {
    std::ostringstream prefix;
    prefix << "Input";
    dump_gso(param.dump_gso_filename, prefix.str(), false);
  }

  if (param.block_size < 2)
    return set_status(RED_SUCCESS);

  int i = 0;

  BKZAutoAbort<FT> auto_abort(m, num_rows);
  
  if ((flags & BKZ_SD_VARIANT) && 
      !(flags & (BKZ_MAX_LOOPS | BKZ_MAX_TIME | BKZ_AUTO_ABORT))) 
  {
    cerr << "Warning: SD Variant of BKZ requires explicit termination condition. Turning auto abort on!" << endl;
    flags |= BKZ_AUTO_ABORT;
  }

  if (flags & BKZ_VERBOSE)
  {
    cerr << "Entering BKZ:" << endl;
    print_params(param, cerr);
    cerr << endl;
  }
  cputime_start = cputime();

  m.discoverAllRows();

  int kappaMax;
  bool clean = true;
  for (i = 0;; i++)
  {
    if ((flags & BKZ_MAX_LOOPS) && i >= param.max_loops)
    {
      final_status = RED_BKZ_LOOPS_LIMIT;
      break;
    }
    if ((flags & BKZ_MAX_TIME) && (cputime() - cputime_start) * 0.001 >= param.max_time)
    {
      final_status = RED_BKZ_TIME_LIMIT;
      break;
    }
    if ((flags & BKZ_AUTO_ABORT) &&
        autoAbort.test_abort(param.auto_abort_scale, param.auto_abort_max_no_dec))
    {
      break;
    }

    try
    {
      if (flags & BKZ_SD_VARIANT)
      {
        clean = sd_tour(i, param, 0, num_rows);
      } 
      else
      {
        clean = tour(i, kappa_max, param, 0, num_rows);
      }
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }

    if (clean || param.block_size >= num_rows)
      break;
  }
  
  if (flags & BKZ_SD_VARIANT) {
    int dummy_kappa_max = num_rows;
    try
    {
      hkz(dummy_kappa_max, param, num_rows - param.block_size, num_rows);
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }
  
  if (flags & BKZ_DUMP_GSO)
  {
    std::ostringstream prefix;
    prefix << "Output ";
    prefix << " (" << std::fixed << std::setw(9) << std::setprecision(3)
           << (cputime() - cputime_start) * 0.001 << "s)";
    dump_gso(param.dump_gso_filename, prefix.str());
  }
  return set_status(final_status);
}

template <class FT> void BKZReduction<FT>::print_tour(const int loop, int min_row, int max_row)
{
  FT r0;
  Float fr0;
  long expo;
  r0  = m.getRExp(min_row, min_row, expo);
  fr0 = r0.get_d();
  fr0.mul_2si(fr0, expo);
  cerr << "End of BKZ loop " << std::setw(4) << loop << ", time = " << std::fixed << std::setw(9)
       << std::setprecision(3) << (cputime() - cputime_start) * 0.001 << "s";
  cerr << ", r_" << min_row << " = " << fr0;
  cerr << ", slope = " << std::setw(9) << std::setprecision(6)
       << get_current_slope(m, min_row, max_row);
  cerr << ", log2(nodes) = " << std::setw(9) << std::setprecision(6) << log2(nodes) << endl;
}

template <class FT> void BKZReduction<FT>::print_params(const BKZParam &param, ostream &out)
{
  out << "block size: " << std::setw(3) << param.block_size << ", ";
  out << "flags: 0x" << std::setw(4) << setfill('0') << std::hex << param.flags << ", " << std::dec
      << std::setfill(' ');
  out << "max_loops: " << std::setw(3) << param.max_loops << ", ";
  out << "max_time: " << std::setw(0) << std::fixed << std::setprecision(1) << param.max_time
      << ", ";
  if (param.flags & BKZ_AUTO_ABORT)
  {
    out << "autoAbort: (" << std::setw(0) << std::fixed << std::setprecision(4)
        << param.auto_abort_scale;
    out << ", " << std::setw(2) << param.auto_abort_max_no_dec << "), ";
  }
  else
  {
    out << "autoAbort: (     -,  -), ";
  }
  out << endl;
}

template <class FT> bool BKZReduction<FT>::set_status(int newStatus)
{
  status = newStatus;
  if (param.flags & BKZ_VERBOSE)
  {
    if (status == RED_SUCCESS)
      cerr << "End of BKZ: success" << endl;
    else
      cerr << "End of BKZ: failure: " << RED_STATUS_STR[status] << endl;
  }
  return status == RED_SUCCESS;
}

template <class FT>
void BKZReduction<FT>::dump_gso(const std::string filename, const std::string prefix, bool append)
{
  ofstream dump;
  if (append)
    dump.open(filename.c_str(), std::ios_base::app);
  else
    dump.open(filename.c_str());
  dump << std::setw(4) << prefix << ": ";
  FT f, logF;
  long expo;
  for (int i = 0; i < num_rows; i++)
  {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    dump << std::setprecision(8) << logF.get_d() + expo * std::log(2.0) << " ";
  }
  dump << std::endl;
  dump.close();
}

template <class FT> bool BKZAutoAbort<FT>::test_abort(double scale, int maxNoDec)
{
  double new_slope = -get_current_slope(m, start_row, num_rows);
  if (no_dec == -1 || new_slope < scale * old_slope)
    no_dec = 0;
  else
    no_dec++;
  old_slope = min(old_slope, new_slope);
  return no_dec >= maxNoDec;
}

/**
   Force instantiation of templates
*/

template class BKZReduction<FP_NR<double> >;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class BKZReduction<FP_NR<long double> >;
#endif

#ifdef FPLLL_WITH_QD
template class BKZReduction<FP_NR<dd_real> >;

template class BKZReduction<FP_NR<qd_real> >;
#endif

#ifdef FPLLL_WITH_DPE
template class BKZReduction<FP_NR<dpe_t> >;
#endif

template class BKZReduction<FP_NR<mpfr_t> >;

FPLLL_END_NAMESPACE
