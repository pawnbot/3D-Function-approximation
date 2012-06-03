#include <math.h>
#include <stdio.h>
#include <string.h>

#include "array_op.h"
#include "precond.h"
#include "iter_alg.h"

////////////////////////////////////////////////////////////////////////////////
// Make array linear combination: RES = X + alpha Y
static __inline__ void
linear_combination_1 (int n,    // array size
                      const double *x,  // input array X
                      const double *y,  // input array Y
                      double alpha,     // multiplier
                      double *res       // output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    res[i] = x[i] + alpha * y[i];

  // Unrolled loop
  for (i = m; i < n; i += 8)
    {
      res[i] = x[i] + alpha * y[i];
      res[i + 1] = x[i + 1] + alpha * y[i + 1];
      res[i + 2] = x[i + 2] + alpha * y[i + 2];
      res[i + 3] = x[i + 3] + alpha * y[i + 3];
      res[i + 4] = x[i + 4] + alpha * y[i + 4];
      res[i + 5] = x[i + 5] + alpha * y[i + 5];
      res[i + 6] = x[i + 6] + alpha * y[i + 6];
      res[i + 7] = x[i + 7] + alpha * y[i + 7];
    }
}

////////////////////////////////////////////////////////////////////////////////
// Make array linear combination: RES = X + beta (RES + alpha Y)
static __inline__ void
linear_combination_2 (int n,    // array size
                      const double *x,  // input array X
                      const double *y,  // input array Y
                      double alpha,     // multiplier
                      double beta,      // multiplier
                      double *res       // input/output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    res[i] = x[i] + beta * (res[i] + alpha * y[i]);

  // Unrolled loop
  for (i = m; i < n; i += 8)
    {
      res[i] = x[i] + beta * (res[i] + alpha * y[i]);
      res[i + 1] = x[i + 1] + beta * (res[i + 1] + alpha * y[i + 1]);
      res[i + 2] = x[i + 2] + beta * (res[i + 2] + alpha * y[i + 2]);
      res[i + 3] = x[i + 3] + beta * (res[i + 3] + alpha * y[i + 3]);
      res[i + 4] = x[i + 4] + beta * (res[i + 4] + alpha * y[i + 4]);
      res[i + 5] = x[i + 5] + beta * (res[i + 5] + alpha * y[i + 5]);
      res[i + 6] = x[i + 6] + beta * (res[i + 6] + alpha * y[i + 6]);
      res[i + 7] = x[i + 7] + beta * (res[i + 7] + alpha * y[i + 7]);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Make array linear combination: RES = RES + alpha X + beta Y
static __inline__ void
linear_combination_5 (int n,    // array size
                      const double *x,  // input array X
                      const double *y,  // input array Y
                      double alpha,     // multiplier
                      double beta,      // multiplier
                      double *res       // input/output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    res[i] += alpha * x[i] + beta * y[i];

  // Unrolled loop
  for (i = m; i < n; i += 8)
    {
      res[i] += alpha * x[i] + beta * y[i];
      res[i + 1] += alpha * x[i + 1] + beta * y[i + 1];
      res[i + 2] += alpha * x[i + 2] + beta * y[i + 2];
      res[i + 3] += alpha * x[i + 3] + beta * y[i + 3];
      res[i + 4] += alpha * x[i + 4] + beta * y[i + 4];
      res[i + 5] += alpha * x[i + 5] + beta * y[i + 5];
      res[i + 6] += alpha * x[i + 6] + beta * y[i + 6];
      res[i + 7] += alpha * x[i + 7] + beta * y[i + 7];
    }
}

////////////////////////////////////////////////////////////////////////////////
/*
 * Solve system A u = b for non-symmetric matrix A
 * Return negative value if error, number of iterations on success
 */
int
solve_linear_system_bcgs (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  int preconditioner_type,      // type of preconditioner
  unsigned int itmax,           // maximal number of iterations
  const double& itprec,          // stopping precision: |Ap - b| / |b| < prec
  const double *a,              // (n x (max_non_zeros+1)) matrix
  const double *rhs,            // (n) right hand side
  const double *u_ini,          // (n) first approximation for solution
  // modified values:
  double *u,                    // (n) on output: solution
  double *a_precond,            // (n x (max_non_zeros+1)) workspace for preconditioner for matrix
                                // (n x (2 x max_non_zeros+1)) workspace for ILU1 preconditioner for matrix
  int *precond_rnz,             // (n x 2 x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in ILU1 preconditioner matrix row
  double *workspace,            // workspace vector of length (7 * n)
  double *p_rhs_norm            // pointer to right hand side norm to be stored
  )
{
  // workspace arrays
  double *r;                    // residual
  double *p;                    // temporary vector
  double *s;                    // temporary vector
  double *t;                    // temporary vector
  double *v;                    // temporary vector
  double *w;                    // temporary vector
  double *z;                    // temporary vector

  double u_ini_norm;            // initial approximation norm
  double rhs_norm;              // initial norm of right hand side
  double residual_norm;         // norm of current residual
  double stop_test;             // stopping criterion
  unsigned int iter;            // Output: number of iterations proceeded
  double masheps = 1.e-16;      // ANSI/IEEE 754-1985 floating point precision
  double min_for_division = 1.e-64;     // Minimal for division

  double c_1, alpha, beta, rho, omega;  // Algorithm variables
  double scalp;                 // Temporary variable

  // initialize workspace arrays
  r = workspace;
  p = r + n;
  s = p + n;
  t = s + n;
  v = t + n;
  w = v + n;
  z = w + n;

  // Save original value of u_ini in u: u = u_ini
  copy_array (n, u_ini, u);

  // Compute initial approximation norm
  u_ini_norm = e2_norm_array (n, u_ini);

  if (u_ini_norm > masheps)
    {
      // Non-zero initial approximation
      // Compute initial residual
      // r = A u
      mult_sparse_matrix_vector (n, max_non_zeros, row_non_zeros, a, u, r);
      // make residual: r = rhs - A u
      sub_arrays (n, rhs, r);
      // Compute norm of residual
      residual_norm = e2_norm_array (n, r);

      // Compute norm of right hand side
      rhs_norm = e2_norm_array (n, rhs);
    }
  else
    {
      // Zero initial approximation
      // r = rhs
      copy_array (n, rhs, r);
      // Compute norm of residual
      residual_norm = e2_norm_array (n, r);
      // Compute norm of right hand side
      rhs_norm = residual_norm;
    }

  if (rhs_norm < masheps)
    rhs_norm = masheps;

  // Store value if pointer non-zero
  if (p_rhs_norm)
    *p_rhs_norm = rhs_norm;
  stop_test = residual_norm / rhs_norm;
  printf ("\tit=%3d, residual = %le, rhs = %le\n", 0, residual_norm, rhs_norm);

  if (stop_test < itprec)
    return 0;   // solution found (in u copy of u_ini)

  // Build preconditioner matrix
  if (make_preconditioner_matrix (n, max_non_zeros, row_non_zeros,
                                  preconditioner_type, a,
                                  a_precond, precond_rnz))
    {
      // Preconditioning method cannot be applied
      printf ("Fatal error: preconditioning failed\n");
      return -1;  // error
    }
  iter = 1;

  // Compute residual z = a_precond^{-1} r for preconditioned system B^{-1} Au = B^{-1} b
  apply_preconditioner_matrix (n, max_non_zeros, row_non_zeros,
                               preconditioner_type, a_precond,
                               precond_rnz, r, z);
  // Initialize r: r = z
  copy_array (n, z, r);

  // Initialize p: p = z
  copy_array (n, z, p);

  // Compute rho = (z, r)
  rho = e2_scalar_product_array (n, z, r);
  if (fabs (rho) < min_for_division)
    {
      // Method cannot be applied with this initial guesss
      printf ("Fatal error: acceleration algorithm breaks down\n");
      return -1;
    }

  // Iterative loop
  for (; iter <= itmax; iter++)
    {
      // Compute v = B^{-1}A(p)
      // First w = A(p)
      mult_sparse_matrix_vector (n, max_non_zeros, row_non_zeros, a, p, w);
      // Next v = B^{-1} w
     apply_preconditioner_matrix (n, max_non_zeros, row_non_zeros,
                                   preconditioner_type, a_precond,
                                   precond_rnz, w, v);
      // Compute alpha = rho / (z, v)
      // Compute first c_1 = (z, v)
      c_1 = e2_scalar_product_array (n, z, v);
      if (fabs (c_1) < min_for_division)
        {
          // Method cannot be applied with this initial guesss
          printf ("Fatal error: acceleration algorithm breaks down\n");
          return -2;
        }
      // Compute alpha = rho / (z, v)
      alpha = rho / c_1;

      // s = r - alpha v
      linear_combination_1 (n, r, v, -alpha, s);

      // Compute t = B^{-1}A(s)
      // First w = A(s)
      mult_sparse_matrix_vector (n, max_non_zeros, row_non_zeros, a, s, w);
      // Next t = B^{-1} w
      apply_preconditioner_matrix (n, max_non_zeros, row_non_zeros,
                                   preconditioner_type, a_precond,
                                   precond_rnz, w, t);

      // Compute omega = (t, s) / (t, t)
      scalp = e2_scalar_product_array (n, t, t);
      if (fabs (scalp) < min_for_division)
        {
          // t = B^{-1}A(s) = 0 => s = 0 that is r = alpha v = alpha B^{-1}A(p)
          // so u = u + alpha is precise solution
          // u = u + alpha p + omega s
          linear_combination_5 (n, p, s, alpha, 0., u);
          break;
        }
      // Compute omega = (t, s) / (t, t)
      omega = e2_scalar_product_array (n, t, s) / scalp;
      if (fabs (omega) < min_for_division)
        {
          // Method cannot be applied with this initial guesss
          printf ("Fatal error: acceleration algorithm breaks down\n");
          return -4;
        }

      // u = u + alpha p + omega s
      linear_combination_5 (n, p, s, alpha, omega, u);

      // r = s - omega t
      linear_combination_1 (n, s, t, -omega, r);

      // Compute residual norm
      residual_norm = e2_norm_array (n, r);

      // Check stopping criterion
      stop_test = residual_norm / rhs_norm;

      if (iter % 100 == 0)
        printf ("\tit=%3d, residual = %le\n", iter, residual_norm);

      if (stop_test < itprec)
        break;

      // Compute rho_next = c_1 = (z, r)
      c_1 = e2_scalar_product_array (n, z, r);
      if (fabs (c_1) < min_for_division)
        {
          // Method cannot be applied with this initial guesss
          printf ("Fatal error: acceleration algorithm breaks down\n");
          return -5;
        }
      // beta = (rho_next / rho) * (alpha / omega)
      beta = c_1 / rho;
      beta *= alpha / omega;
      // rho = rho_next
      rho = c_1;

      // p = r + beta (p - omega v)
      linear_combination_2 (n, r, v, - omega, beta, p);
    }

  if (iter >= itmax)
    {
      // Convergence too slow or failed
      printf ("Fatal error: failure to converge in max %d iterations\n",
                  itmax);

      return -10;
    }

  return iter;
}

