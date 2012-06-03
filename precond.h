#ifndef PRECOND_H
#define PRECOND_H

// Preconditioner types
#define PRECONDITIONER_NO               0
#define PRECONDITIONER_ILU             1
#define PRECONDITIONER_JACOBI           2
#define PRECONDITIONER_LUMP_JACOBI      3
#define PRECONDITIONER_ILU1             4

/*
 * Build ILU1 Factorization
 * Return non-zero if error
 */
int
make_sparse_matrix_ilu1_factorisation (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  const double *a,              // (n x (max_non_zeros+1)) matrix A
  double *r,                    // (n x (2 x max_non_zeros+1)) result matrices L,U
  int *r_rnz                    // (n x 2 x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in preconditioner matrix row              
  );

/*
 * Build ILU Factorization
 * Return non-zero if error
 */
int
make_sparse_matrix_ilu_factorisation (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  const double *a,              // (n x (max_non_zeros+1)) matrix A
  double *r                     // (n x (max_non_zeros+1)) result matrices L,U
  );

/*
 * Solve system U*X=Y for Incomplete LU Factorization sparse
 * upper triangle matrix U
 * X may coinside with Y
 */
void
solve_system_sparse_up_matrix_ilu_factorized (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  const double *a,              //!< (n x (max_non_zeros+1)) matrix A
  const double *y,              //!< (n) right hand size
  double *x                     //!< (n) result vector
  );

/*
 * Solve system L*X=Y for Incomplete LU Factorization sparse
 * lower triangle matrix L
 */
void
solve_system_sparse_down_matrix_ilu_factorized (
  unsigned int n,                       // order of system
  unsigned int max_non_zeros,           // maximal number of non-zero off-diagonal
                                        // elements in matrix row
  const int *row_non_zeros,             // (n x max_non_zeros) array with column
                                        // numbers of non-zero off-diagonal elements in matrix row
  const double *a,                      // (n x (max_non_zeros+1)) matrix A
  const double *y,                      // (n) right hand size
  double *x                             // (n) result vector
  );

int
make_jacobian_preconditioner_matrix (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const double *a,              // (n x (max_non_zeros+1)) matrix A
  double *r                     // (n x (max_non_zeros+1)) result preconditioner matrix
  );

int
make_lumped_jacobian_preconditioner_matrix (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  const double *a,              // (n x (max_non_zeros+1)) matrix A
  double *r                     // (n x (max_non_zeros+1)) result preconditioner matrix
  );

/*
 * Prepare preconditioner for linear solver
 * Return non-zero if error
 */
int
make_preconditioner_matrix (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  int preconditioner_type,      // type of preconditioner
  const double *a,              // (n x (max_non_zeros+1)) matrix A
  double *r,                    // (n x (max_non_zeros+1)) workspace for preconditioner for matrix
                                // (n x (2 x max_non_zeros+1)) workspace for ILU1 preconditioner for matrix
  int *r_rnz                    // (n x 2 x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in preconditioner matrix row
  );

/*
 * Apply preconditioner for linear solver
 * Return non-zero if error
 */
void
apply_preconditioner_matrix (
  unsigned int n,                       // order of system
  unsigned int max_non_zeros,           // maximal number of non-zero off-diagonal
                                        // elements in matrix row
  const int *row_non_zeros,             // (n x max_non_zeros) array with column
                                        // numbers of non-zero off-diagonal elements in matrix row
  int preconditioner_type,              // type of preconditioner
  const double *a_precond,              // (n x (max_non_zeros+1)) workspace for preconditioner for matrix
                                        // (n x (2 x max_non_zeros+1)) workspace for ILU1 preconditioner for matrix
  const int *precond_rnz,               // (n x 2 x max_non_zeros) array with column
                                        // numbers of non-zero off-diagonal elements in preconditioner matrix row
  const double *y,                      // (n) right hand size
  double *x                             // (n) result vector
  );

#endif /* PRECOND_H */
