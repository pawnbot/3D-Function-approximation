#ifndef ITER_ALG_H
#define ITER_ALG_H

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
  const double *a,              // (n x (max_non_zeros + 1)) matrix
  const double *rhs,            // (n) right hand side
  const double *u_ini,          // (n) first approximation for solution
  // modified values:
  double *u,                    // (n) on output: solution
  double *a_precond,            // (n x (max_non_zeros+1)) workspace for preconditioner for matrix
                                // (n x (2 x max_non_zeros+1)) workspace for ILU1 preconditioner for matrix
  int *precond_rnz,             // (n x 2 x max_non_zeros) array with column, for ILU1 preconditioner
                                // numbers of non-zero off-diagonal elements in preconditioner matrix row
  double *workspace,            // workspace vector of length (7 * n)
  double *p_rhs_norm            // pointer to right hand side norm to be stored
  );

#endif /* ITER_ALG_H */
