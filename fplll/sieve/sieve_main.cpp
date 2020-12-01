/*
  This provides an implementation of Gauss sieving, including using
  tuples of vectors in fplll framework. The Gauss Sieve code is
  based on Panagiotis Voulgaris's implementation of the Gauss sieve.
  In addition a CVP(P) functionality is added based on using sieving
  as the preprocessing.
*/
#include "sieve_main.h"
#include "fplll.h"

/**
 * help function
 */
static void main_usage(char *myself)
{
  cout << "Usage: " << myself << " [options]\n"
       << "List of options:\n"
       << "  -a [2|3|4]\n"
       << "     2- or 3- or 4-sieve;\n"
       << "  -f filename\n"
       << "     Input filename\n"
       << "  -r nnn\n"
       << "     Generate a random instance of dimension nnn\n"
       << "  -t nnn\n"
       << "     Targeted norm^2=nnn\n"
       << "  -s nnn\n"
       << "     Using seed=nnn\n"
       << "  -b nnn\n"
       << "     BKZ preprocessing of blocksize=nnn\n"
       << "  -p filename\n"
       << "     Output filename for the list\n"
       << "  -L filename\n"
       << "     Input filename used for CVPP (includes the list)\n"
       << "  -C filename\n"
       << "     Input filename for CVPP (includes the target(s))\n"
       << "  -v\n"
       << "     Verbose mode\n";
}

/**
 * run sieve
 */
template <class ZT>
int main_run_sieve(ZZ_mat<ZT> B, Z_NR<ZT> target_norm, int alg, bool ver, int seed,
                   char *file_name)  //----
{
  clock_t stime = clock(), etime = clock();
  double secs;

  GaussSieve<ZT, FP_NR<double>> gsieve(B, alg, ver, seed);
  gsieve.sieve(target_norm);

  etime = clock();
  secs  = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (ver)
    cout << "# [info] sieve took time " << secs << " s" << endl;

  gsieve.print_list_to_file(file_name);
  return 0;
}

/**
 * run cvpp
 */
template <class ZT>
int main_run_cvpp(ZZ_mat<ZT> B, Z_NR<ZT> target_norm, int alg, bool ver, int seed,
                  const char *input_list_file, const char *targets_file_name,
                  const char *output_list_file)
{
  GaussSieve<ZT, FP_NR<double>> gsieve(B, alg, ver, seed);
  gsieve.approx_voronoi_cvpp(input_list_file, targets_file_name, output_list_file);
  return 0;
}

/**
 * main function
 */
int main(int argc, char **argv)
{
  char *input_file_name  = NULL;
  char *target_norm_s    = NULL;
  char *output_list_file = NULL, *input_list_file = NULL, *targets_file_name = NULL;
  bool flag_verbose = true, flag_file = false, flag_cvpp = false;
  int option, alg, dim = 10, seed = 0, bs = 0;

#if 0
  dot_time = 0;
  dot_num = 0;
  count_bad = 0;
#endif

  alg = 2;

  /* parse */
  if (argc == 1)
  {
    main_usage(argv[0]);
    return -1;
  }
  while ((option = getopt(argc, argv, "a:f:r:t:s:b:p:L:C:v")) != -1)
  {
    switch (option)
    {
    case 'a':
      alg = atoi(optarg);
      if (alg != 2 && alg != 3 && alg != 4)
        throw std::invalid_argument("only support 2-, 3- and 4-sieve");
      break;
    case 'f':
      input_file_name = optarg;
      flag_file       = true;
      break;
    case 'r':
      dim       = atoi(optarg);
      flag_file = false;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'b':
      bs = atoi(optarg);
      break;
    case 'v':
      flag_verbose = true;
      break;
    case 't':
      // ntarget_norm = atol(optarg);
      // cout << optarg << endl;   ---- not necessary ----
      target_norm_s = optarg;
      break;
    case 'p':
      output_list_file = optarg;
      break;
    case 'L':
      input_list_file = optarg;
      break;
    case 'C':
      targets_file_name = optarg;
      flag_cvpp         = true;
      break;
    case 'h':
      main_usage(argv[0]);
      return -1;
    case '?':
      main_usage(argv[0]);
      return -1;
    case ':':
      main_usage(argv[0]);
      return -1;
    }
  }

  /* set lattice */
  ZZ_mat<mpz_t> B;
  if (flag_file)
  {
    ifstream input_file(input_file_name);
    if (input_file.is_open())
    {
      input_file >> B;
      input_file.close();
    }
    else
    {
      cin >> B;
    }
    if (flag_verbose)
    {
      cout << "# [info] reading lattice of dimension " << B.get_rows() << "x" << B.get_cols()
           << endl;
    }
  }
  else
  {
    if (flag_verbose)
    {
      cout << "# [info] generating random lattice of dimension " << dim << endl;
    }
    srand(time(NULL));
    B.resize(dim, dim);
    B.gen_trg(1.1);
  }

  /* set targeted norm */
  Z_NR<mpz_t> target_norm, max;
  if (target_norm_s != NULL)
  {
    target_norm.set_str(target_norm_s);
  }
  if (target_norm < 0)
    target_norm = 0;
  if (flag_verbose)
    cout << "# [info] target norm^2 is " << target_norm << endl;

  /* preprocessing of basis */
  clock_t stime = clock();
  if (bs > 0)
    bkz_reduction(B, bs, BKZ_DEFAULT, FT_DEFAULT, 0);
  else
    lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

  clock_t etime = clock();
  double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (flag_verbose)
  {
    if (bs > 0)
      cout << "# [info] BKZ took time " << secs << " s" << endl;
    else
      cout << "# [info] LLL took time " << secs << " s" << endl;
  }
  // cout << B << endl;

  /* decide integer type */
  max = B.get_max();

#if 1
  if (max < std::numeric_limits<int>::max())
  {
    long target_norm_l = abs(target_norm.get_si());
    Z_NR<long> target_norm_lt;
    target_norm_lt = target_norm_l;
    ZZ_mat<long> B2(B.get_rows(), B.get_cols());
    for (int i = 0; i < B.get_rows(); i++)
      for (int j = 0; j < B.get_cols(); j++)
        B2(i, j) = B(i, j).get_si();
    /* Choose if we solve cvpp or svp */
    if (flag_cvpp)
      main_run_cvpp<long>(B2, target_norm_lt, alg, flag_verbose, seed, input_list_file,
                          targets_file_name, output_list_file);
    else
      main_run_sieve<long>(B2, target_norm_lt, alg, flag_verbose, seed, output_list_file);
  }
  else
  {
#endif
    /* Choose if we solve cvpp or svp */
    if (flag_cvpp)
      main_run_cvpp<mpz_t>(B, target_norm, alg, flag_verbose, seed, input_list_file,
                           targets_file_name, output_list_file);
    else
      main_run_sieve<mpz_t>(B, target_norm, alg, flag_verbose, seed, output_list_file);
  }

  if (flag_verbose)
  {
/* dot product time */
#if 0
    cout << "# [info] dot_time " << dot_time << endl;
    cout << "# [info] dot_num " << dot_num << endl;
    cout << "# [info] dot_time/dot_number " << (double) dot_time/dot_num << endl;
#endif
  }
  return 1;
}
