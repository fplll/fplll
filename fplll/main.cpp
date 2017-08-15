/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2008-2011 Xavier Pujol.

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

#include "main.h"
#include <config.h>

template <class ZT> int lll(Options &o, ZZ_mat<ZT> &b)
{
  // Stupid intialization of u and u_inv to be not empty.
  ZZ_mat<ZT> u(1, 1), u_inv(1, 1);
  const char *format = o.output_format ? o.output_format : "b";
  int status, flags = 0;
  if (o.verbose)
    flags |= LLL_VERBOSE;
  if (o.early_red)
    flags |= LLL_EARLY_RED;
  if (o.siegel)
    flags |= LLL_SIEGEL;

  if (strchr(format, 'v') != NULL)
  {
    // LLL-reduction with transform and inverse transform
    status = lll_reduction(b, u, u_inv, o.delta, o.eta, o.method, o.float_type, o.precision, flags);
  }
  else if (strchr(format, 'u') != NULL)
  {
    // LLL-reduction with transform
    status = lll_reduction(b, u, o.delta, o.eta, o.method, o.float_type, o.precision, flags);
  }
  else
  {
    status = lll_reduction(b, o.delta, o.eta, o.method, o.float_type, o.precision, flags);
  }

  for (int i = 0; format[i]; i++)
  {
    switch (format[i])
    {
    case 'b':
      cout << b << endl;
      break;
    case 'u':
      cout << u << endl;
      break;
    case 'v':
      cout << u_inv << endl;
      break;
    case 't':
      cout << status << endl;
      break;
    case ' ':
      cout << endl;
      break;
    }
  }
  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
  }
  return status;
}

/* BKZ reduction */

void read_pruning_vector(const char *file_name, PruningParams &pr, int n)
{
  double x;
  FILE *file = fopen(file_name, "r");
  CHECK(file, "Cannot open '" << file_name << "'");

  pr.coefficients.clear();
  for (int i = 0; i <= n && fscanf(file, "%lf", &x) == 1; i++)
  {
    pr.coefficients.push_back(x);
    CHECK(x > 0 && x <= 1, "Number " << x << " in file '" << file_name
                                     << "' is not in the interval (0,1]");
    if (i == 0)
    {
      CHECK(x == 1, "The first number in file '" << file_name << "' should be 1");
    }
    else
    {
      CHECK(pr.coefficients[i] <= pr.coefficients[i - 1],
            "File '" << file_name << "' should contain a non-increasing sequence of numbers");
    }
  }
  CHECK(static_cast<int>(pr.coefficients.size()) == n,
        "File '" << file_name << "' should contain exactly " << n << " numbers");
}

template <class ZT> int bkz(Options &o, ZZ_mat<ZT> &b) { ABORT_MSG("mpz required for BKZ"); }

template <> int bkz(Options &o, ZZ_mat<mpz_t> &b)
{
  CHECK(o.block_size > 0, "Option -b is missing");
  vector<Strategy> strategies;
  if (!o.bkz_strategy_file.empty())
  {
    strategies = load_strategies_json(strategy_full_path(o.bkz_strategy_file));
  }

  BKZParam param(o.block_size, strategies);
  // Stupid intialization of u to be not empty.
  ZZ_mat<mpz_t> u(1, 1);
  const char *format = o.output_format ? o.output_format : "b";
  int status;

  param.delta = o.delta;
  param.flags = o.bkz_flags;

  if (o.bkz_flags & BKZ_DUMP_GSO)
    param.dump_gso_filename = o.bkz_dump_gso_filename;
  if (o.bkz_flags & BKZ_GH_BND)
    param.gh_factor = o.bkz_gh_factor;
  if (o.bkz_flags & BKZ_MAX_LOOPS)
    param.max_loops = o.bkz_max_loops;
  if (o.verbose)
    param.flags |= BKZ_VERBOSE;
  if (o.no_lll)
    param.flags |= BKZ_NO_LLL;

  status = bkz_reduction(&b, strchr(format, 'u') ? &u : NULL, param, o.float_type, o.precision);

  for (int i = 0; format[i]; i++)
  {
    switch (format[i])
    {
    case 'b':
      cout << b << endl;
      break;
    case 'u':
      cout << u << endl;
      break;
    case 't':
      cout << status << endl;
      break;
    case ' ':
      cout << endl;
      break;
    }
  }
  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
  }
  return status;
}

/* HKZ reduction
   Note: since we only force |mu_i,j| <= eta with eta > 0.5, the solution
   is not unique even for a generic matrix */

template <class ZT> int hkz(Options &o, ZZ_mat<ZT> &b) { ABORT_MSG("mpz required for HKZ"); }

template <> int hkz(Options &o, ZZ_mat<mpz_t> &b)
{
  int flags = 0;
  if (o.verbose)
    flags |= HKZ_VERBOSE;
  int status = hkz_reduction(b, flags, o.float_type, o.precision);
  cout << b << endl;
  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
  }
  return status;
}

/* Shortest vector problem and closest vector problem */

template <class ZT> int svpcvp(Options &o, ZZ_mat<ZT> &b, const vector<Z_NR<ZT>> &target)
{
  if (target.empty())
  {
    ABORT_MSG("mpz required for SVP");
  }
  else
  {
    ABORT_MSG("mpz required for CVP");
  }
}

template <> int svpcvp(Options &o, ZZ_mat<mpz_t> &b, const vector<Z_NR<mpz_t>> &target)
{
  const char *format = o.output_format ? o.output_format : "s";
  vector<Z_NR<mpz_t>> sol_coord;    // In the LLL-reduced basis
  vector<Z_NR<mpz_t>> sol_coord_2;  // In the initial basis
  vector<Z_NR<mpz_t>> solution;
  ZZ_mat<mpz_t> u;
  bool with_coord     = strchr(format, 'c') != NULL;
  bool with_coord_std = strchr(format, 's') != NULL;
  int flags           = SVP_DEFAULT | (o.verbose ? SVP_VERBOSE : 0);
  int flagsLLL        = LLL_DEFAULT | (o.verbose ? LLL_VERBOSE : 0);
  int status;

  if (!o.no_lll)
  {
    if (with_coord)
    {
      status = lll_reduction(b, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, flagsLLL);
    }
    else
    {
      status = lll_reduction(b, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, flagsLLL);
    }
    if (status != RED_SUCCESS)
    {
      cerr << "LLL reduction failed: " << get_red_status_str(status) << endl;
      return status;
    }
  }

  if (target.empty())
    status = shortest_vector(b, sol_coord, SVPM_PROVED, flags);
  else
    status = closest_vector(b, target, sol_coord, flags);

  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
    return status;
  }
  if (with_coord)
  {
    if (o.no_lll)
      sol_coord_2 = sol_coord;
    else
      vector_matrix_product(sol_coord_2, sol_coord, u);
  }
  if (with_coord_std)
  {
    vector_matrix_product(solution, sol_coord, b);
  }

  for (int i = 0; format[i]; i++)
  {
    switch (format[i])
    {
    case 'c':
      cout << sol_coord_2 << endl;
      break;
    case 's':
      cout << solution << endl;
      break;
    case 't':
      cout << status << endl;
      break;
    case ' ':
      cout << endl;
      break;
    }
  }
  return status;
}

template <class ZT> int run_action(Options &o)
{
  istream *is;
  ZZ_mat<ZT> m;
  vector<Z_NR<ZT>> target;

  if (o.input_file)
    is = new ifstream(o.input_file);
  else
    is = &cin;

  *is >> m;
  if (o.action == ACTION_CVP)
  {
    *is >> target;
  }
  if (!*is)
    ABORT_MSG("invalid input");
  if (o.input_file)
    delete is;

  int result = 0;
  switch (o.action)
  {
  case ACTION_LLL:
    result = lll(o, m);
    break;
  case ACTION_SVP:
    result = svpcvp(o, m, target);
    break;
  case ACTION_CVP:
    result = svpcvp(o, m, target);
    break;
  case ACTION_HKZ:
    result = hkz(o, m);
    break;
  case ACTION_BKZ:
    result = bkz(o, m);
    break;
  default:
    ABORT_MSG("unimplemented action");
    break;
  }
  return result;
}

/* Command line parsing */

void read_options(int argc, char **argv, Options &o)
{
  for (int ac = 1; ac < argc; ac++)
  {
    if (strcmp(argv[ac], "-a") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -a switch");
      if (strcmp(argv[ac], "lll") == 0)
        o.action = ACTION_LLL;
      else if (strcmp(argv[ac], "hkz") == 0)
        o.action = ACTION_HKZ;
      else if (strcmp(argv[ac], "bkz") == 0)
        o.action = ACTION_BKZ;
      else if (strcmp(argv[ac], "svp") == 0)
        o.action = ACTION_SVP;
      else if (strcmp(argv[ac], "cvp") == 0)
        o.action = ACTION_CVP;
      else if (strcmp(argv[ac], "sdb") == 0)
      {
        o.action = ACTION_BKZ;
        o.bkz_flags |= BKZ_SD_VARIANT;
      }
      else if (strcmp(argv[ac], "sld") == 0)
      {
        o.action = ACTION_BKZ;
        o.bkz_flags |= BKZ_SLD_RED;
      }
      else
        ABORT_MSG("parse error in -a switch: lll or svp expected");
    }
    else if (strcmp(argv[ac], "-b") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -b switch");
      o.block_size = atoi(argv[ac]);
    }
    else if (strcmp(argv[ac], "-bkzboundedlll") == 0)
    {
      o.bkz_flags |= BKZ_BOUNDED_LLL;
    }
    else if (strcmp(argv[ac], "-bkzmaxloops") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after '-bkzmaxloops'");
      o.bkz_max_loops = atoi(argv[ac]);
      o.bkz_flags |= BKZ_MAX_LOOPS;
    }
    else if (strcmp(argv[ac], "-bkzmaxtime") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after '-bkzmaxtime'");
      o.bkz_max_time = atof(argv[ac]);
      o.bkz_flags |= BKZ_MAX_TIME;
    }
    else if (strcmp(argv[ac], "-bkzautoabort") == 0)
    {
      o.bkz_flags |= BKZ_AUTO_ABORT;
    }
    else if (strcmp(argv[ac], "-bkzdumpgso") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing filename after -bkzdumpgso switch");
      o.bkz_dump_gso_filename = argv[ac];
      o.bkz_flags |= BKZ_DUMP_GSO;
    }
    else if (strcmp(argv[ac], "-c") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -c switch");
      // o.c=atoi(argv[ac]); // ignored (was the number of columns)
    }
    else if (strcmp(argv[ac], "-bkzghbound") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after '-bkzghbound'");
      o.bkz_gh_factor = atof(argv[ac]);
      o.bkz_flags |= BKZ_GH_BND;
    }
    else if (strcmp(argv[ac], "-d") == 0 || strcmp(argv[ac], "-delta") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -d switch");
      o.delta = atof(argv[ac]);
    }
    else if (strcmp(argv[ac], "-e") == 0 || strcmp(argv[ac], "-eta") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -e switch");
      o.eta = atof(argv[ac]);
    }
    else if (strcmp(argv[ac], "-f") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -f switch");
      if (strcmp("mpfr", argv[ac]) == 0)
        o.float_type = FT_MPFR;
      else if (strcmp("dpe", argv[ac]) == 0)
        o.float_type = FT_DPE;
      else if (strcmp("dd", argv[ac]) == 0)
        o.float_type = FT_DD;
      else if (strcmp("qd", argv[ac]) == 0)
        o.float_type = FT_QD;
      else if (strcmp("double", argv[ac]) == 0)
        o.float_type = FT_DOUBLE;
      else if (strcmp("longdouble", argv[ac]) == 0)
        o.float_type = FT_LONG_DOUBLE;
      else
        ABORT_MSG("parse error in -f switch : mpfr, qd, dd, dpe or double expected");
    }
    else if (strcmp(argv[ac], "-s") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -s switch");
      o.bkz_strategy_file = argv[ac];
    }
    else if (strcmp(argv[ac], "-l") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -l switch");
      o.siegel = !atoi(argv[ac]);
    }
    else if (strcmp(argv[ac], "-m") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -m switch");
      if (strcmp("wrapper", argv[ac]) == 0)
        o.method = LM_WRAPPER;
      else if (strcmp("proved", argv[ac]) == 0)
        o.method = LM_PROVED;
      else if (strcmp("heuristic", argv[ac]) == 0)
        o.method = LM_HEURISTIC;
      else if (strcmp("fast", argv[ac]) == 0)
        o.method = LM_FAST;
      else if (strcmp("fastearly", argv[ac]) == 0)
      {
        o.method    = LM_FAST;
        o.early_red = true;
      }
      else if (strcmp("heuristicearly", argv[ac]) == 0)
      {
        o.method    = LM_HEURISTIC;
        o.early_red = true;
      }
      else
        ABORT_MSG("parse error in -m switch : proved, heuristic, fast, "
                  << "or wrapper expected");
    }
    else if (strcmp(argv[ac], "-nolll") == 0)
    {
      o.no_lll = true;
    }
    else if (strcmp(argv[ac], "-of") == 0)
    {
      ac++;
      CHECK(ac < argc, "missing value after -of switch");
      o.output_format = argv[ac];
    }
    else if (strcmp(argv[ac], "-p") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -p switch");
      o.precision = atoi(argv[ac]);
    }
    else if (strcmp(argv[ac], "-r") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -r switch");
      // o.r = atoi(argv[ac]); // ignored (was the number of rows)
    }
    else if (strcmp(argv[ac], "-v") == 0)
    {
      o.verbose = true;
    }
    else if (strcmp(argv[ac], "-y") == 0)
    {
      o.early_red = true;
    }
    else if (strcmp(argv[ac], "-z") == 0)
    {
      ++ac;
      CHECK(ac < argc, "missing value after -z switch");
      if (strcmp("mpz", argv[ac]) == 0)
        o.int_type = ZT_MPZ;
      else if (strcmp("long", argv[ac]) == 0 || strcmp("int", argv[ac]) == 0)
        o.int_type = ZT_LONG;
      else if (strcmp("double", argv[ac]) == 0)
        o.int_type = ZT_DOUBLE;
      else
        ABORT_MSG("parse error in -z switch : int, double or mpz expected");
    }
    else if ((strcmp(argv[ac], "-h") == 0) || (strcmp(argv[ac], "--help") == 0))
    {
      cout << "Usage: " << argv[0] << " [options] [file]\n"

           << "List of options:\n"
           << "  -a [lll|bkz|hkz|svp|sdb|sld|cvp]\n"
           << "       lll = LLL-reduce the input matrix (default)\n"
           << "       bkz = BKZ-reduce the input matrix\n"
           << "       hkz = HKZ-reduce the input matrix\n"
           << "       svp = compute a shortest non-zero vector of the lattice\n"
           << "       sdb = reduce the input matrix using the self dual BKZ variant\n"
           << "       sld = slide reduce the input matrix\n"
           << "       cvp = compute the vector in the input lattice closest to an input vector\n"
           << "  -v\n"
           << "       Enable verbose mode\n"
           << "  -nolll\n"
           << "       Does not apply intial LLL-reduction (for bkz, hkz and svp)\n"
           << "  -c <size>\n"
           << "       Was the number of columns (ignored)\n"
           << "  -r <size>\n"
           << "       Was the number of rows (ignored)\n"

           << "  -d <delta> (default=0.99; alias to -delta <delta>)\n"
           << "  -e <eta> (default=0.51; alias to -eta <eta>)\n"
           << "  -l <lovasz>\n"
           << "       If <lovasz> != 0, Lovasz's condition, otherwise, Siegel's condition\n"
           << "  -f [mpfr|dd|qd|dpe|double|longdouble]\n"
           << "       Floating-point type in LLL\n"
           << "  -p <precision>\n"
           << "       Floating-point precision (only with -f mpfr)\n"
           << "  -z [mpz|int|long|double]\n"
           << "       Integer type in LLL (default: mpz; long is an alias to int)\n"
           << "  -m [wrapper|fast|fastearly|heuristic|heuristicearly|proved]\n"
           << "       LLL version (default: wrapper)\n"
           << "  -y\n"
           << "       Enable early reduction\n"

           << "  -b <block_size>\n"
           << "       Size of BKZ blocks\n"
           << "  -bkzmaxloops <loops>\n"
           << "       Maximum number of full loop iterations\n"
           << "  -bkzmaxtime <time>\n"
           << "        Stops after <time> seconds\n"
           << "  -bkzautoabort\n"
           << "        Stops when the average slope does not decrease fast enough\n"
           << "  -s <filename.json>\n"
           << "        Load BKZ strategies from filename\n"
           << "  -bkzghbound <factor>\n"
           << "        Multiplies the Gaussian heuristic by <factor> (of float type)\n"
           << "  -bkzboundedlll\n"
           << "        Restricts the LLL call\n"
           << "  -bkzdumpgso <file_name>\n"
           << "        Dumps the log of the Gram-Schmidt vectors in specified file\n"

           << "  -of [bcstuv]\n"
           << "        Output formats.\n"

           << "Please refer to https://github.com/fplll/fplll/README.md for more informations.\n";
      exit(0);
    }
    else if (strcmp(argv[ac], "--version") == 0)
    {
      cout << "fplll " << VERSION << endl
           << "Copyright 2005-2012 Damien Stehle, David Cade, Xavier Pujol." << endl
           << "fplll is free software. You can redistribute it and/or modify" << endl
           << "it under the terms of the GNU Lesser General Public License as published by" << endl
           << "the Free Software Foundation, either version 2.1 of the License, or" << endl
           << "(at your option) any later version." << endl;
      exit(0);
    }
    else if (argv[ac][0] == '-')
    {
      ABORT_MSG("invalid option '" << argv[ac] << "'.\n"
                                                  "Try '"
                                   << argv[0] << " --help' for more information.");
    }
    else
    {
      CHECK(!o.input_file, "too many input files");
      o.input_file = argv[ac];
    }
  }
}

int main(int argc, char **argv)
{
  int result;
  Options o;
  read_options(argc, argv, o);
  ZZ_mat<mpz_t>::set_print_mode(MAT_PRINT_REGULAR);
  switch (o.int_type)
  {
  case ZT_MPZ:
    result = run_action<mpz_t>(o);
    break;
#ifdef FPLLL_WITH_ZLONG
  case ZT_LONG:
    result = run_action<long int>(o);
    break;
#endif
#ifdef FPLLL_WITH_ZDOUBLE
  case ZT_DOUBLE:
    result = run_action<double>(o);
    break;
#endif
  default:
    ABORT_MSG("compiled without support for this integer type");
  }
  return result;
}
