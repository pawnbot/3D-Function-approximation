#ifndef ARRAY_OP_H
#define ARRAY_OP_H

/*
 * Вычисляет произведение разреженной матрицы на вектор
 */
void
mult_sparse_matrix_vector (
  unsigned int n_rows,          // number of rows in matrix
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal 
                                // elements in matrix row
  const int *row_non_zeros,     // n_rows x max_non_zeros array with
                                // column numbers of non-zero off-diagonal
                                // elements in matrix row
  const double *matrix,         // n_rows x (max_non_zeros + 1) array
  const double *u,              // funcion u
  double *destin                // result - matrix * u
  );

/*
 * Вычисляет отличие разреженной матрицы от симметричной
 */
double
sparse_matrix_symmetric_p (
  unsigned int n_rows,          // number of rows in matrix
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal 
                                // elements in matrix row
  const unsigned int *row_non_zeros,    // n_rows x max_non_zeros array with
                                // column numbers of non-zero off-diagonal
                                // elements in matrix row
  const double *matrix          // n_rows x (max_non_zeros + 1) array
  );

/* Copy array x to array y. */
void
copy_array (int n,              // array size
            const double *src,  // input array
            double *dst         // output array
  );

/* Compute euclidean norm of array */
double
e2_norm_array (int n,           //!< array size
               const double *a  //!< array
  );

/* Compute euclidean scalar product of arrays */
double
e2_scalar_product_array (int n, //!< array size
                         const double *a,       //!< array 1
                         const double *b        //!< array 2
  );

/* Subtract arrays: dst = src - dst */
void
sub_arrays (int n,              // array size
            const double *src,  // input array
            double *dst         // output array
  );

/* Scale array x */
void
scale_array (int n,              // array size
            double *a,          // input/output array
            double mult
  );

/*
 * Write sparse matrix A into binary file filename.
 * Return negative value if error, number of bytes written on success
 */
int
write_sparse_matrix_file (
  const char *filename,         // file name
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  const double *a               // (n x (max_non_zeros+1)) matrix
  );

/*
 * Read sparse matrix A from binary file filename.
 * Return negative value if error, number of bytes read on success
 */
int
read_sparse_matrix_file (
  const char *filename,         // file name
  unsigned int *p_n,            // order of system
  unsigned int *p_max_non_zeros,// maximal number of non-zero off-diagonal
                                // elements in matrix row
  int **p_row_non_zeros,        // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  double **p_a                  // (n x (max_non_zeros+1)) matrix
  );

/*
 * Write vector V into binary file filename.
 * Return negative value if error, number of bytes written on success
 */
int
write_vector_file (
  const char *filename,         // file name
  unsigned int n,               // order of system
  const double *v               // (n) vector
  );

/*
 * Read vector from binary file filename.
 * Return negative value if error, number of bytes read on success
 */
int
read_vector_file (
  const char *filename,         // file name
  unsigned int *p_n,            // vector length
  double **p_v                  // (n) matrix
  );

/*
 * Распечатать разреженную матрицу
 */
void
print_sparse_matrix (
  unsigned int n_rows,          // number of rows in matrix
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal 
                                // elements in matrix row
  const int *row_non_zeros,    // n_rows x max_non_zeros array with
                                // column numbers of non-zero off-diagonal
                                // elements in matrix row
  const double *matrix          // n_rows x (max_non_zeros + 1) array
  );

/*
 * Распечатать вектор
 */
void
print_vector (
  unsigned int n_rows,          // number of elements in vector
  const double *vector          // n_rows array
  );

#endif /* ARRAY_OP_H */
