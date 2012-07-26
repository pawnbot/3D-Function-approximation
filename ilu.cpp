#include <math.h>
#include <stdio.h>
#include <string.h>
#include "precond.h"

/******************************************************************************/
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
  )
{
  unsigned int i, k, l, m, ni;
  int j, nk;
  const double *pa, *pa2;
  double *pr, *pr2, *pr3;
  const int *pn, *pn3;
  int *prn, *pn2;
  double s;
  //new variables
  int added;
  int count;

  // Minimal value alloved on diagonal of factorized matrix
  const double min_on_diagonal = 1.e-12;

  // Zero result
  memset (r, 0, (2 * max_non_zeros + 1) * n * sizeof (double));
  count = 0;
  // Copy A zero pattern
  for (i = 0; i < n; ++i) {
    pn = row_non_zeros + i * max_non_zeros;  // pointer to column numbers of non-zero elements of i-th row
    prn = r_rnz + i * (2 * max_non_zeros);   // pointer to column numbers of non-zero elements of precond

    for (l = 0; l < max_non_zeros; ++l)
      prn[l] = pn[l];
    for (l = max_non_zeros; l < 2 * max_non_zeros; ++l)
      prn[l] = -1;
  }

  // Generating ILU1 zero pattern
  for (i = 0; i < n; ++i) {
    prn = r_rnz + i * (2 * max_non_zeros);   // pointer to column numbers of non-zero elements of precond

    for (l = 0; l < max_non_zeros; ++l) {
      j = prn[l];
      if (j >= 0 && j < (int) i)
        for (k = l+1; k < max_non_zeros; ++k) {
          nk = prn[k];
          if (nk >= 0 && nk < (int) i) {
            // Looking for a[j][nk] and a[nk][j]
            // Pattern is symmetrical
            added = 0;
            pn2 = r_rnz + j * (2 * max_non_zeros);  // pointer to column numbers of non-zero elements of j-th row

            for (m = 0; m < 2*max_non_zeros; ++m)
              if (pn2[m] == nk) {
                added = 1;
                break;
              }
            if (added != 0) continue;
            if ((r_rnz[(j+1)*2*max_non_zeros - 1] >= 0) || (r_rnz[(nk+1)*2*max_non_zeros - 1] >= 0)) {
              printf("NOT ENOUGH SPACE IN ROWS: %d | %d |\n", j , nk);
              return -1;
            }
            ++count;
            // Adding a[j][nk]
            for (m = max_non_zeros; m < 2*max_non_zeros; ++m)
              if  (pn2[m] < 0) {
                pn2[m] = nk;
                break;
              }
            // Adding a[nk][j]
            pn2 = r_rnz + nk * (2 * max_non_zeros);  // pointer to column numbers of non-zero elements of nk-th row

            for (m = max_non_zeros; m < 2*max_non_zeros; ++m)
              if  (pn2[m] < 0) {
                pn2[m] = j;
                break;
              }
          }
        }
    }
  }
  printf("Number of elements added to pattern during ILU1: %d\n", 2*count);
  // Loop over matrix rows
  for (i = 0; i < n; i++)
    {
      pa = a + i * (max_non_zeros + 1);     // pointer to i-th row a
      pr = r + i * (2*max_non_zeros + 1);   // pointer to i-th row r
      pn = r_rnz + i * 2 * max_non_zeros;   // pointer to column numbers of non-zero elements

      // Compute diagonal element of L
      // Sum: s = \sum_{j=1}^{i-1} l_{ij}u_{ji}
      s = 0.;
      for (k = 0; k < 2 * max_non_zeros; k++)
        {
          j = pn[k];
          if (j >= 0 && j < (int) i)
            {
              // l_{ij} = pr[k + 1] - element of L in row i and column j
              // Look for element of U in row j and column i
              pr2 = r + j * (2*max_non_zeros + 1);  // pointer to j-th row
              pn2 = r_rnz + j * 2 * max_non_zeros;  // pointer to column numbers of non-zero elements of j-th row

              // Look for element in row j and column i
              for (m = 0; m < 2 * max_non_zeros; m++)
                {
                  if (pn2[m] == (int) i)
                    {
                      // there is element
                      // s += l_{ij} u_{ji}
                      s += pr[k + 1] * pr2[m + 1];
                      break;
                    }
                }
            }
        }
      // l_{ii} = a_{ii} - \sum_{j=1}^{i-1} l_{ij} u_{ji}
      // Correct matrix to finish process
      if (fabs (pa[0] - s) >= min_on_diagonal)
        pr[0] = pa[0] - s;
      else {
        if (pa[0] - s >= 0.)
         pr[0] = min_on_diagonal;
        else
          pr[0] = - min_on_diagonal;
      }
      // off-diagonal elements of L
      for (l = 0; l < 2 * max_non_zeros; l++)
        {
          j = pn[l];
          if (j >= 0 && j > (int) i)
            {
              pa2 = a + j * (max_non_zeros + 1);   // pointer to j-th row A
              pr2 = r + j * (2*max_non_zeros + 1); // pointer to j-th row L
              pn2 = r_rnz + j * 2 * max_non_zeros; // pointer to column numbers of non-zero elements of j-th row

              // l_{ji} = a_{ji}-\sum_{k=1}^{i-1} l_{jk}u_{ki}
              // Summ : s = \sum_{k=1}^{i-1} l_{jk}u_{ki}
              s = 0.;
              ni = (unsigned int) -1;
              for (k = 0; k < 2 * max_non_zeros; k++)
                {
                  nk = pn2[k];
                  if (nk >= 0 && nk < (int) i)
                    {
                      // l_{jk} = pr2[k + 1]
                      pr3 = r + nk * (2*max_non_zeros + 1); // pointer to nk-th row U
                      pn3 = r_rnz + nk * 2 * max_non_zeros; // pointer to column numbers of non-zero elements of nk-th row

                      // Look for element in row nk and column i
                      for (m = 0; m < 2 * max_non_zeros; m++)
                        {
                          if (pn3[m] == (int) i)
                            {
                              // there is element
                              // u_{ki} = pr3[m + 1]
                              // s += l_{jk}u_{ki}
                              s += pr2[k + 1] * pr3[m + 1];
                              break;
                            }
                        }
                    }
                  else if (nk == (int) i)
                    ni = k; // column number of l_{ji} = l_{j,ni}
                }
              // Check arrays consistence
              if (ni == (unsigned int) (-1))
                return -1;      // wrong off-diagonal element numbers
              // r_{ji} = a_{ji} - s;
              if (ni < max_non_zeros)
                pr2[ni + 1] = pa2[ni + 1] - s;
              else
                pr2[ni + 1] = -s;
            }
        }

      // off-diagonal elements of U
      for (l = 0; l < 2 * max_non_zeros; l++)
        {
          j = pn[l];
          if (j > (int) i)
            {
              // u_{ij}=(a_{ij}-\sum_{k=1}^{i-1} l_{ik}u_{kj})/l_{ii}
              // Summ : s = \sum_{k=1}^{i-1} l_{ik}u_{kj}
              s = 0.;
              for (k = 0; k < 2 * max_non_zeros; k++)
                {
                  nk = pn[k];
                  if (nk >= 0 && nk < (int) i)
                    {
                      // l_{ik} = pr[k + 1]
                      pr2 = r + nk * (2*max_non_zeros + 1); // pointer to nk-th row U
                      pn2 = r_rnz + nk * 2 * max_non_zeros; // pointer to column numbers of non-zero elements of nk-th row

                      // Look for element in row nk and column j
                      for (m = 0; m < 2 * max_non_zeros; m++)
                        {
                          if (pn2[m] == j)
                            {
                              // there is element
                              // u_{kj} = pr2[m + 1]
                              // s += l_{ik}u_{kj}
                              s += pr[k + 1] * pr2[m + 1];
                              break;
                            }
                        }
                    }
                }
              // u_{ij} = (a_{ij} - s) / l_{ii}
              if (l < max_non_zeros)
                pr[l + 1] = (pa[l + 1] - s) / pr[0];
              else
                pr[l + 1] = -s / pr[0];
            }
        }
    }
  return 0;
}

/******************************************************************************/

/******************************************************************************/
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
  )
{
  unsigned int i, k, l, m, ni;
  int j, nk;
  const double *pa, *pa2;
  double *pr, *pr2, *pr3;
  const int *pn, *pn2, *pn3;
  double s;

  // Minimal value alloved on diagonal of factorized matrix
  const double min_on_diagonal = 1.e-12;

  // Zero result
  memset (r, 0, (max_non_zeros + 1) * n * sizeof (double));

  // Loop over matrix rows
  for (i = 0; i < n; i++)
    {
      pa = a + i * (max_non_zeros + 1); // pointer to i-th row a
      pr = r + i * (max_non_zeros + 1); // pointer to i-th row r
      pn = row_non_zeros + i * max_non_zeros;   // pointer to column numbers of non-zero elements

      // Compute diagonal element of L
      // Sum: s = \sum_{j=1}^{i-1} l_{ij}u_{ji}
      s = 0.;
      for (k = 0; k < max_non_zeros; k++)
        {
          j = pn[k];
          if (j >= 0 && j < (int) i)
            {
              // l_{ij} = pr[k + 1] - element of L in row i and column j
              // Look for element of U in row j and column i
              pr2 = r + j * (max_non_zeros + 1);    // pointer to j-th row
              pn2 = row_non_zeros + j * max_non_zeros; // pointer to column numbers of non-zero elements of j-th row

              // Look for element in row j and column i
              for (m = 0; m < max_non_zeros; m++)
                {
                  if (pn2[m] == (int) i)
                    {
                      // there is element
                      // s += l_{ij} u_{ji}
                      s += pr[k + 1] * pr2[m + 1];
                      break;
                    }
                }
            }
        }
      // l_{ii} = a_{ii} - \sum_{j=1}^{i-1} l_{ij} u_{ji}
      // Correct matrix to finish process
      if (fabs (pa[0] - s) >= min_on_diagonal)
        pr[0] = pa[0] - s;
      else if (pa[0] - s >= 0.)
        pr[0] = min_on_diagonal;
      else
        pr[0] = - min_on_diagonal;

      // off-diagonal elements of L
      for (l = 0; l < max_non_zeros; l++)
        {
          j = pn[l];
          if (j >= 0 && j > (int) i)
            {
              pa2 = a + j * (max_non_zeros + 1); // pointer to j-th row A
              pr2 = r + j * (max_non_zeros + 1); // pointer to j-th row L
              pn2 = row_non_zeros + j * max_non_zeros; // pointer to column numbers of non-zero elements of j-th row

              // l_{ji} = a_{ji}-\sum_{k=1}^{i-1} l_{jk}u_{ki}
              // Summ : s = \sum_{k=1}^{i-1} l_{jk}u_{ki}
              s = 0.;
              ni = (unsigned int) -1;
              for (k = 0; k < max_non_zeros; k++)
                {
                  nk = pn2[k];
                  if (nk >= 0 && nk < (int) i)
                    {
                      // l_{jk} = pr2[k + 1]
                      pr3 = r + nk * (max_non_zeros + 1); // pointer to nk-th row U
                      pn3 = row_non_zeros + nk * max_non_zeros; // pointer to column numbers of non-zero elements of nk-th row

                      // Look for element in row nk and column i
                      for (m = 0; m < max_non_zeros; m++)
                        {
                          if (pn3[m] == (int) i)
                            {
                              // there is element
                              // u_{ki} = pr3[m + 1]
                              // s += l_{jk}u_{ki}
                              s += pr2[k + 1] * pr3[m + 1];
                              break;
                            }
                        }
                    }
                  else if (nk == (int) i)
                    ni = k; // column number of l_{ji} = l_{j,ni}
                }
              // Check arrays consistence
              if (ni == (unsigned int) (-1))
                {
                  printf ("Cannot find %d\n", ni);
                  return -1;      // wrong off-diagonal element numbers
                }
              // r_{ji} = a_{ji} - s;
              pr2[ni + 1] = pa2[ni + 1] - s;
            }
        }

      // off-diagonal elements of U
      for (l = 0; l < max_non_zeros; l++)
        {
          j = pn[l];
          if (j > (int) i)
            {
              // u_{ij}=(a_{ij}-\sum_{k=1}^{i-1} l_{ik}u_{kj})/l_{ii}
              // Summ : s = \sum_{k=1}^{i-1} l_{ik}u_{kj}
              s = 0.;
              for (k = 0; k < max_non_zeros; k++)
                {
                  nk = pn[k];
                  if (nk >= 0 && nk < (int) i)
                    {
                      // l_{ik} = pr[k + 1]
                      pr2 = r + nk * (max_non_zeros + 1); // pointer to nk-th row U
                      pn2 = row_non_zeros + nk * max_non_zeros; // pointer to column numbers of non-zero elements of nk-th row

                      // Look for element in row nk and column j
                      for (m = 0; m < max_non_zeros; m++)
                        {
                          if (pn2[m] == j)
                            {
                              // there is element
                              // u_{kj} = pr2[m + 1]
                              // s += l_{ik}u_{kj}
                              s += pr[k + 1] * pr2[m + 1];
                              break;
                            }
                        }
                    }
                }
              // u_{ij} = (a_{ij} - s) / l_{ii}
              pr[l + 1] = (pa[l + 1] - s) / pr[0];
            }
        }
    }
  return 0;
}

/******************************************************************************/

/******************************************************************************/

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
  )
{
  int i, j;
  unsigned int i_row;           // index of row data
  unsigned int i_col;           // index of column numbers
  int l;
  double res;

  // Sequential exclusion starting from the last row
  for (i = n - 1, i_row = (i + 1) * (max_non_zeros + 1) - 1,
       i_col = (i + 1) * max_non_zeros - 1; i >= 0; i--)
    {
      // x_i = (y_i - \sum_{j=i+1}^{n}u_{i,j}x_j)
      res = y[i];
      // off-diagonal elements in upper triangle
      for (j = max_non_zeros - 1; j >= 0; j--, i_row--, i_col--)
        {
          if ((l = row_non_zeros[i_col]) >= i)
            res -= a[i_row] * x[l];
        }
      // Diagonal element of U not stored
      i_row--;
      x[i] = res;
    }
}

/******************************************************************************/
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
  )
{
  unsigned int i, j;
  unsigned int i_row;           // index of row data
  unsigned int i_col;           // index of column numbers
  int l;
  double res, v;

  // Sequential exclusion starting from the first row
  for (i = 0, i_row = 0, i_col = 0; i < n; i++)
    {
      // x_i = (y_i - \sum_{j=1}^{i-1}l_{i,j}x_j) / r_{ii}
      // Save diagonal element
      v = a[i_row++];
      res = y[i];
      for (j = 0; j < max_non_zeros; j++, i_row++, i_col++)
        {
          if ((l = row_non_zeros[i_col]) >= 0 && l < (int) i) {
            res = res - a[i_row] * x[l];
          }
        }
      // Diagonal element tested for zero during factorisation process
      x[i] = res / v;
    }
}

/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
// Copy array x of double to array y.
static inline void
copy_array (int n,              // array size
            const double *src,  // input array
            double *dst         // output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  for (i = 0; i < m; i++)
    dst[i] = src[i];

  // Unrolled loop
  for (i = m; i < n; i += 8)
    {
      dst[i] = src[i];
      dst[i + 1] = src[i + 1];
      dst[i + 2] = src[i + 2];
      dst[i + 3] = src[i + 3];
      dst[i + 4] = src[i + 4];
      dst[i + 5] = src[i + 5];
      dst[i + 6] = src[i + 6];
      dst[i + 7] = src[i + 7];
    }
}

////////////////////////////////////////////////////////////////////////////////
int
make_jacobian_preconditioner_matrix (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const double *a,              // (n x (max_non_zeros+1)) matrix A
  double *r                     // (n x (max_non_zeros+1)) result preconditioner matrix
  )
{
  // Minimal value alloved on diagonal of factorized matrix
  const double min_on_diagonal = 1.e-12;

  // Loop over matrix rows
  for (unsigned int i = 0; i < n; i++)
    {
      const double *pa = a + i * (max_non_zeros + 1); // pointer to i-th row a
      // Get diagonal value
      double val = pa[0];
      if (fabs (val) < min_on_diagonal)
        {
          if (val > 0)
            val = min_on_diagonal;
          else
            val = -min_on_diagonal;
        }

      // Build preconditioner diagonal matrix
      r[i] = 1. / val;
    }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
int
make_lumped_jacobian_preconditioner_matrix (
  unsigned int n,               // order of system
  unsigned int max_non_zeros,   // maximal number of non-zero off-diagonal
                                // elements in matrix row
  const int *row_non_zeros,     // (n x max_non_zeros) array with column
                                // numbers of non-zero off-diagonal elements in matrix row
  const double *a,              // (n x (max_non_zeros+1)) matrix A
  double *r                     // (n x (max_non_zeros+1)) result preconditioner matrix
  )
{
  // Minimal value alloved on diagonal of factorized matrix
  const double min_on_diagonal = 1.e-12;

  // Loop over matrix rows
  for (unsigned int i = 0; i < n; i++)
    {
      const double *pa = a + i * (max_non_zeros + 1); // pointer to i-th row a
      const int *pn = row_non_zeros + i * max_non_zeros;   // pointer to column numbers of non-zero elements
      // Get diagonal value
      double val = pa[0];

      // Sum off-diagonal values
      for (unsigned int k = 0; k < max_non_zeros; k++)
        {
          if (pn[k] >= 0)
            val += pa[k + 1];
        }

      if (fabs (val) < min_on_diagonal)
        {
          if (val > 0)
            val = min_on_diagonal;
          else
            val = -min_on_diagonal;
        }

      if (fabs (val) < 1.e-2)
        r[i] = 1.;
      else
        // Build preconditioner diagonal matrix
        r[i] = 1. / val;
    }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
/*
 * Scale vector x by vector mult and put result into y
 */
static inline void
scale_vector_coordinates2 (
  unsigned int n,               // vector length
  const double *x,              // (n) vector x
  const double *mults,          // (n) multipliers for backward substitution
  double *y                     // (n) vector y
  )
{
  unsigned int m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (unsigned int i = 0; i < m; i++)
    y[i] = x[i] * mults[i];

  // Unrolled loop
  for (unsigned int i = m; i < n; i += 8)
    {
      y[i] = x[i] * mults[i];
      y[i + 1] = x[i + 1] * mults[i + 1];
      y[i + 2] = x[i + 2] * mults[i + 2];
      y[i + 3] = x[i + 3] * mults[i + 3];
      y[i + 4] = x[i + 4] * mults[i + 4];
      y[i + 5] = x[i + 5] * mults[i + 5];
      y[i + 6] = x[i + 6] * mults[i + 6];
      y[i + 7] = x[i + 7] * mults[i + 7];
    }
}

////////////////////////////////////////////////////////////////////////////////
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
  )
{
  switch (preconditioner_type)
    {
    case PRECONDITIONER_NO:
      return 0;

    case PRECONDITIONER_ILU:
      return make_sparse_matrix_ilu_factorisation
               (n, max_non_zeros, row_non_zeros, a, r);

    case PRECONDITIONER_JACOBI:
      return make_jacobian_preconditioner_matrix
               (n, max_non_zeros, a, r);

    case PRECONDITIONER_LUMP_JACOBI:
      return make_lumped_jacobian_preconditioner_matrix
               (n, max_non_zeros, row_non_zeros, a, r);

    case PRECONDITIONER_ILU1:
      return make_sparse_matrix_ilu1_factorisation
               (n, max_non_zeros,
               row_non_zeros, a, r, r_rnz);
    default:
      return -1;
    }
}

////////////////////////////////////////////////////////////////////////////////
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
  )
{
  switch (preconditioner_type)
    {
    default:
    case PRECONDITIONER_NO:
      // x = y
      copy_array (n, y, x);
      return;

    case PRECONDITIONER_ILU:
      // Solve system with right hand side y and matrix = L: L x = y
      solve_system_sparse_down_matrix_ilu_factorized
        (n, max_non_zeros, row_non_zeros, a_precond, y, x);
      // solve system with right hand side x and matrix = U: U x = x
      solve_system_sparse_up_matrix_ilu_factorized
        (n, max_non_zeros, row_non_zeros, a_precond, x, x);
      return;

    case PRECONDITIONER_JACOBI:
    case PRECONDITIONER_LUMP_JACOBI:
      scale_vector_coordinates2 (n, y, a_precond, x);
      return;

    case PRECONDITIONER_ILU1:
      // Solve system with right hand side y and matrix = L: L x = y
      solve_system_sparse_down_matrix_ilu_factorized
        (n, 2 * max_non_zeros, precond_rnz, a_precond, y, x);
      // solve system with right hand side x and matrix = U: U x = x
      solve_system_sparse_up_matrix_ilu_factorized
        (n, 2 * max_non_zeros, precond_rnz, a_precond, x, x);
      return;
    }
}

