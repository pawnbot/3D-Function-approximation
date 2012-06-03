#include "approximator.h"
#include "precond.h"
#include "iter_alg.h"

#define FREE_ARRAY(A) {if (A) delete[] (A); (A) = 0;}

static inline double
length (double x1, double y1, double x2, double y2)
{
  return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

approximator::approximator(QObject *parent) :
  QThread(parent)
{
  threads_count = 1;

  old_answer.res = 0;

  res = 0;
  workspace = 0;
  matrix_rnz = 0;
  matrix = 0;
  rhs = 0;

  n1 = 50;
  n2 = 50;

  /// Jacobian = MINOR / MAJOR * r
  jacobian = polynom (MINOR / MAJOR, 0., 0.);
  /// Polynom class checking here
  point_t pij = point_t (5.58273, 4.20386);
  point_t p2 = point_t (5.64029, 4.24906);
  point_t p = point_t (5.64029, 4.20386);
  polynom matrix_pol = jacobian * kurant_pol (pij, p2, p) * kurant_pol (pij, p2, p);
  printf ("Polynom: ");
  matrix_pol.print ();
  printf ("Points: (%g, %g), (%g, %g), (%g, %g)\n", pij.x, pij.y, p2.x, p2.y, p.x, p.y);
  printf ("Integral: %.15lf\n", matrix_pol.integrate (p, pij, p2));
}

void approximator::clean_up ()
{
  FREE_ARRAY (res)
  FREE_ARRAY (workspace)
  FREE_ARRAY (matrix_rnz)
  FREE_ARRAY (matrix)
  FREE_ARRAY (rhs)
}

approximator::~approximator ()
{
  clean_up ();
  FREE_ARRAY (old_answer.res)
}

void approximator::update_answer ()
{
  if (solver_done)
    {
      old_answer.n1 = n1;
      old_answer.n2 = n2;
      FREE_ARRAY (old_answer.res)
      old_answer.res = res;
      res = 0;
    }
}

void approximator::update_values (
  int n1,
  int n2,
  functions function_index
  )
{
  if (this->function_index != function_index)
    {
      this->function_index = function_index;
      FREE_ARRAY (old_answer.res)
    }
  this->n1 = n1;
  this->n2 = n2;
}

void approximator::run ()
{
  solver_done = false;
  clean_up ();
  /// GENERATE MATRIX FOR THREADS
  printf ("Generating matrix...\n");
  h1 = MAJOR / (double) (n1 - 1);
  h2 = 2 * M_PI / (double) (n2 - 1);

  size_t n_rows = n1 * n2;
  size_t max_non_zeros = 6;
  matrix_rnz = new int [n_rows * max_non_zeros];
  matrix = new double [n_rows * (max_non_zeros + 1)];
  rhs = new double [n_rows];
  /// Generating matrix and rhs
#pragma omp parallel for
  for (int i = 0; i < n1; ++i)
    for (int j = 0; j < n2; ++j)
      {
        int u1, u2;
        int index = get_index (i, j, n1, n2);
        /// Generating Matrix
        u1 = 0;
        u2 = 0;
        /// Diagonal element
        matrix[index * (max_non_zeros + 1) + 0]     = l2_prod (i, j, i + u1, j + u2);
        if (matrix[index * (max_non_zeros + 1) + 0] < 1.e-12)
          {
            printf ("Trash at diagonal elemnt: %g\n", matrix[index * (max_non_zeros + 1) + 0]);
          }
        u1 = 1;
        u2 = 0;
        matrix_rnz[index * max_non_zeros]           = get_index (i + u1, j + u2, n1, n2);
        matrix[index * (max_non_zeros + 1) + 1]     = l2_prod (i, j, i + u1, j + u2);
        u1 = 1;
        u2 = 1;
        matrix_rnz[index * max_non_zeros + 1]       = get_index (i + u1, j + u2, n1, n2);
        matrix[index * (max_non_zeros + 1) + 2]     = l2_prod (i, j, i + u1, j + u2);
        u1 = 0;
        u2 = 1;
        matrix_rnz[index * max_non_zeros + 2]       = get_index (i + u1, j + u2, n1, n2);
        matrix[index * (max_non_zeros + 1) + 3]     = l2_prod (i, j, i + u1, j + u2);
        u1 = -1;
        u2 = 0;
        matrix_rnz[index * max_non_zeros + 3]       = get_index (i + u1, j + u2, n1, n2);
        matrix[index * (max_non_zeros + 1) + 4]     = l2_prod (i, j, i + u1, j + u2);
        u1 = -1;
        u2 = -1;
        matrix_rnz[index * max_non_zeros + 4]       = get_index (i + u1, j + u2, n1, n2);
        matrix[index * (max_non_zeros + 1) + 5]     = l2_prod (i, j, i + u1, j + u2);
        u1 = 0;
        u2 = -1;
        matrix_rnz[index * max_non_zeros + 5]       = get_index (i + u1, j + u2, n1, n2);
        matrix[index * (max_non_zeros + 1) + 6]     = l2_prod (i, j, i + u1, j + u2);
        /// Generating RHS
        rhs[index] = rhs_integrate (i, j);
      }
  double rhs_norm = 0;
  res = new double [n_rows];
  workspace = new double [7 * n_rows];

  memset (res, 0, n_rows * sizeof (double));

  /// Solve system
  printf ("Solving system...\n");
  omp_set_num_threads (threads_count);
  PRECISION = 1.e-14;
  double *precond = 0;
  if (0)
    {
      precond = new double [n_rows];
      solve_linear_system_bcgs (n_rows, max_non_zeros, matrix_rnz, PRECONDITIONER_JACOBI, 1000, PRECISION,
                                matrix, rhs, res, res, precond, 0, workspace, &rhs_norm);
    }
  else
    {
      precond = new double [n_rows * (2 * max_non_zeros + 1)];
      int iters =  solve_linear_system_bcgs (n_rows, max_non_zeros, matrix_rnz, PRECONDITIONER_ILU,
                                             1000, PRECISION, matrix, rhs, res, res, precond,
                                             0, workspace, &rhs_norm);
      printf ("Solved in %d iterations!\n", iters);
      if (PRECISION < 1.)
        solver_done = true;
    }
  FREE_ARRAY (precond)
}

void approximator::set_threads_count (int count)
{
  threads_count = count;
}

void approximator::stop ()
{
#pragma omp atomic
  PRECISION += 1.e80;
#pragma omp atomic
  n1 -= 10000;
#pragma omp atomic
  n2 -= 10000;
}

double approximator::evaluate (double x, double y)
{
  if (!old_answer.res)
    return 0;

  /// Elipse to Circle
  x  = x;
  y *= MAJOR / MINOR;
  /// Circle to Rectangle
  double radius = sqrt (x * x + y * y);
  double angle = 0;
  if (radius > 1.e-14)
    {
      angle = atan2 (y, x);
      if (angle <= 0)
        angle += 2 * M_PI;
    }

  if (radius >= MAJOR)
    radius -= MAJOR - 1.e-15;
  if (angle >= 2 * M_PI)
    angle = 2 * M_PI - 1.e-15;

  x = radius;
  y = angle;

  double h1 = MAJOR / (double) (old_answer.n1 - 1);
  double h2 = 2 * M_PI / (double) (old_answer.n2 - 1);
  int i = x / h1;
  int j = y / h2;

  double result = evaluate_phi (i, j, h1, h2, x, y);
  result += evaluate_phi (i, j + 1, h1, h2, x, y);
  result += evaluate_phi (i + 1, j, h1, h2, x, y);
  result += evaluate_phi (i + 1, j + 1, h1, h2, x, y);
  return result;
}
