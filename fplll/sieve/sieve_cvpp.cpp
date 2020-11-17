#include "sieve_gauss.h"

/**
 * sets the filename for the output
 */
string solutions_filename(const char *targets_file_name)
{
  string file_name = string(targets_file_name);
  string old_name;
  old_name = file_name.substr(0, file_name.find_last_of('.'));
  return old_name + string("-cvpp_solutions.txt");
}

/**
 * imports target vectors from file
 * Note: the targets are imported without their norm included.
 */
template <class ZT> list<ListPoint<ZT> *> target_import(const char *targets_file_name, const int nc)
{
  list<ListPoint<ZT> *> target_list;
  NumVect<Z_NR<ZT>> t_vec(nc);
  ListPoint<ZT> *p;

  ifstream input_file(targets_file_name);
  if (!(input_file.is_open()))
  {
    cout << "# [ERROR] The input file " << targets_file_name;
    cout << " was not opened properly. The target(s) were not imported." << endl;
  }
  else
  {
    p = new_listpoint<ZT>(nc);
    while (input_file >> t_vec)
    {
      p = num_vec_to_list_point(t_vec, nc);
      if (p->norm > 0)
        target_list.push_back(p);
    }
    input_file.close();
  }
  return target_list;
}

/**
 * imports List from a file
 * Note: the list of vectors is imported without their norms
 */
template <class ZT, class F>
void GaussSieve<ZT, F>::import_list_from_file(const char *list_file_name)
{
  NumVect<Z_NR<ZT>> t_vec;
  ListPoint<ZT> *p;

  ifstream input_file(list_file_name);
  if (!(input_file.is_open()))
  {
    cout << "[ERROR] The input file " << list_file_name;
    cout << " was not opened properly. The list was not imported." << endl;
  }
  else
  {
    free_list_queue();
    p = new_listpoint<ZT>(nc);
    while (input_file >> t_vec)
    {
      p = num_vec_to_list_point(t_vec, nc);
      if (p->norm > 0)
        List.push_back(p);
    }
    input_file.close();

    /* Ensure that the imported list is sorted */
    List.sort(compare_listpoints<ZT>);

    if (verbose)
    {
      cout << "# [info] Imported a list of " << List.size();
      cout << " lattice vectors from file " << list_file_name << endl;
    }
  }
}

/**
 * implements the iterative slicer
 */
template <class ZT, class F> void GaussSieve<ZT, F>::iterative_slicer(ListPoint<ZT> *target)
{
  typename list<ListPoint<ZT> *>::iterator lp_it;
  ListPoint<ZT> *vc;
  bool loop = true;

  /* Loop while the target could get reduced */
  int count = 0;
  while (loop)
  {
    count++;
    loop = false;

    /* Loop through the list of vectors */
    for (lp_it = List.begin(); lp_it != List.end(); ++lp_it)
    {
      vc = *lp_it;
      // if ((v->norm) >= 4*(target->norm))
      //  break;

      /* If there is a reduction the vector should re-pass the list */
      if (half_2reduce(target, vc))
      {
        loop = true;
        break;
      }
    }
  }
}

/**
 * computes the inverse of the base-probability of the
 * randomised slicer as given in DLW20
 */
double dlw20_base(double alpha)
{
  int i;
  double beta, num, base, sq, u, v, x_i, x_i1;

  beta = (alpha * alpha) / (4.0 * (alpha - 1.0));
  num =
      ceil(-0.5 + sqrt(pow(4.0 * beta - alpha, 2.0) - 8.0 * (2.0 * beta - alpha)) / (2.0 * alpha));

  if (num == 1)
  {
    base = 1.0 / (alpha - pow(alpha + beta - 1.0, 2.0) / (4.0 * beta));
  }
  else
  {
    sq   = sqrt(pow(alpha * num * num - beta - 1.0, 2.0) + 4.0 * beta * (num * num - 1.0));
    u    = ((beta + 1.0 - alpha) * num - sq) / (pow(num, 3.0) - num);
    v    = ((alpha - 2.0 * beta) * pow(num, 2.0) + beta - 1.0 + sq * num) / (pow(num, 3.0) - num);
    base = 1.0;
    x_i1 = beta;
    for (i = 1; i < num + 1; i++)
    {
      x_i = u * i * i + v * i + beta;
      base *= (1.0 / (alpha - pow(alpha + x_i1 - x_i, 2.0) / (4.0 * x_i1)));
      x_i1 = x_i;
    }
  }
  base = sqrt(base);
  return base;
}

/**
 * implements the randomised slicer
 */
template <class ZT, class F>
void GaussSieve<ZT, F>::randomised_slicer(ListPoint<ZT> *target, const int bound)
{
  Z_NR<ZT> dot;
  NumVect<Z_NR<ZT>> v_rand(nc);
  ListPoint<ZT> *start_t, *t_red;
  int i;

  start_t = new_listpoint<ZT>(nc);
  t_red   = new_listpoint<ZT>(nc);

  clone_listpoint(target, start_t);
  clone_listpoint(target, t_red);

  /* Call the iterative slicer #bound times */
  for (i = 0; i < bound; i++)
  {
    iterative_slicer(t_red);

    if (t_red->norm < target->norm)
      clone_listpoint(t_red, target);

    /* Add a random lattice point to the target. */
    do
    {
      v_rand = Sampler->sample();
      squared_norm(dot, v_rand);
    } while (dot == 0);
    (v_rand).add(start_t->v);
    squared_norm(dot, v_rand);
    set_listpoint_numvect(v_rand, dot, t_red);
  }

  del_listpoint(start_t);
  del_listpoint(t_red);
}

/**
 * runs the randomised slicer for a given list of targets
 * target vector type will interact with the lattice vector type (for now the same)
 */
template <class ZT, class F>
void GaussSieve<ZT, F>::approx_voronoi_cvpp(const char *list_file_name,
                                            const char *targets_file_name,
                                            const char *output_list_file)
{
  Z_NR<ZT> target_norm;
  NumVect<Z_NR<ZT>> cvp(nc), t(nc);
  ListPoint<ZT> *t_red = new_listpoint<ZT>(nc);
  list<ListPoint<ZT> *> target_list;
  typename list<ListPoint<ZT> *>::iterator lp_it;

  int bound, targets_num = 0, targets_count;
  double alpha, expb, base;
  clock_t stime, etime;
  double secs = 0;

  /* Import a list or create a list and store it */
  if (verbose)
    cout << endl << "# [****] Preprocessing starts..." << endl;

  if (list_file_name != NULL)
    import_list_from_file(list_file_name);
  else
  {
    target_norm = 0;
    sieve(target_norm);
    if (output_list_file != NULL)
    {
      print_list_to_file(output_list_file);
      if (verbose)
      {
        cout << "# [info] Created a list of " << List.size();
        cout << " lattice vectors and stored it in the file " << output_list_file << endl;
      }
    }
  }

  /* Import targets from file */
  target_list = target_import<ZT>(targets_file_name, nc);
  targets_num = target_list.size();
  if (targets_num == 0)
  {
    cout << "# [ERROR] No targets given in the file " << targets_file_name << " ." << endl;
    return;
  }
  if (verbose)
  {
    cout << endl << "# [****] Query phase starts..." << endl;
    cout << "# [info] Imported a list of " << targets_num;
    cout << " target vector(s) from file " << targets_file_name << endl;
  }

  // debug
  /*cout << return_first() << endl;
  for (lp_it = target_list.begin(); lp_it != target_list.end(); ++lp_it)
  {
    cout << (*lp_it)->v <<endl;
  }*/

  /* Obtain a number of rerandomisations based on DLW20 */
  alpha = pow(List.size(), 2.0 / (double)nr);
  base  = dlw20_base(alpha);
  expb  = pow(base, nr);
  expb *= 20;
  expb  = fmax(1.0, expb);
  bound = round(expb);

  /* Run randomised slicer for target(s) */
  string solution_file = solutions_filename(targets_file_name);
  ofstream output_file(solution_file);

  targets_count = 1;
  stime         = clock();
  for (lp_it = target_list.begin(); lp_it != target_list.end(); ++lp_it)
  {
    t_red = *lp_it;
    t     = t_red->v;

    randomised_slicer(t_red, bound);

    etime = clock();
    secs  = (etime - stime) / (double)CLOCKS_PER_SEC;

    if (verbose)
    {
      cout << "\33[2K\r# [info] Processing targets: " << targets_count << "/" << targets_num;
      cout << " in " << secs << " s" << flush;
      targets_count++;
    }

    /* Print the answer(s) to file */
    cvp = t;
    cvp.sub(t_red->v);
    output_file << cvp << endl;

    // debug
    /*output_file << "Target vector: " << endl;
    output_file << t << endl;
    output_file << "Error vector: " << endl;
    output_file << (t_red)->v << " ";
    output_file << (t_red)->norm << endl;
    output_file << "Closest lattice vector: " << endl;
    cvp = t;
    cvp.sub(t_red->v);
    output_file << cvp << endl << endl;*/
  }
  output_file.close();

  if (verbose)
  {
    cout << endl;
    cout << "# [****] done!" << endl;
    cout << "# [info] Listsize = " << List.size() << endl;
    cout << "# [info] Number of targets = " << targets_num << endl;
    cout << "# [info] Number of rerandomisations = 20*2^(" << log2(base) << "*" << nr
         << ") = " << bound << endl;
    secs = secs / (double)targets_num;
    cout << "# [info] Average query time = " << secs << " s" << endl;
  }

  del_listpoint(t_red);
}
