#include <math.h>
#include <stdio.h>
#include <string.h>
#include "array_op.h"

/******************************************************************************/
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
  )
{
  /* Обходим строки матрицы */
  for (unsigned int i = 0; i < n_rows; i++)
    {
      const double *pr;
      const int *pn;
      double res;
      // указатель на строку, соответствующую узлу
      pr = matrix + i * (max_non_zeros  + 1);
      // указатель на строку соседей для узла
      pn = row_non_zeros + i * max_non_zeros;

      // диагональный элемент
      res = pr[0] * u[i];

      // внедиагональные элементы
      for (unsigned int j = 0; j < max_non_zeros; j++)
        {
          if (pn[j] != (-1))
            res += pr[j + 1] * u[pn[j]];
        }
      destin[i] = res;
    }
}

/******************************************************************************/
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
  )
{
  unsigned int i;
  const double *pr, *prk;
  const unsigned int *pn, *pnk;
  unsigned int j, k, l;
  double res = 0., v1, v2;

  /* Обходим строки матрицы */
  for (i = 0; i < n_rows; i++)
    {
      // Указатель на строку матрицы
      pr = matrix + i * (max_non_zeros + 1);
      // Указатель на строку номеров внедиагональных элементов
      pn = row_non_zeros + i * max_non_zeros;

      // диагональный элемент pr[0] не расматривается

      // внедиагональные элементы
      for (j = 0; j < max_non_zeros; j++)
        {
          // k - номер столбца, где стоит элемент pr[j + 1]
          if ((k = pn[j]) == (unsigned int) (-1))
            continue;
          v1 = pr[j + 1];
          // берем в строке k элемент в i-м столбце
          // указатель на k-ю строку
          prk = matrix + k * (max_non_zeros + 1);
          // указатель на номера столбцов k-й строки
          pnk = row_non_zeros + k * max_non_zeros;
          v2 = 0.;
          for (l = 0; l < max_non_zeros; l++)
            {
              if (pnk[l] != (unsigned int) (-1)  && pnk[l] == i)
                {
                  v2 = prk[l + 1];
                  break;
                }
            }
          res += fabs (v1 - v2);
        }
    }
  return res;
}

/******************************************************************************/

// Copy array x to array y.
void
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
// Compute euclidean norm of array
double
e2_norm_array (int n,           //!< array size
               const double *a  //!< array
  )
{
  int i, m;
  double res = 0.;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    res += a[i] * a[i];

  // Unrolled loop
  for (i = m; i < n; i += 8)
    res += a[i] * a[i] + a[i + 1] * a[i + 1]
      + a[i + 2] * a[i + 2] + a[i + 3] * a[i + 3]
      + a[i + 4] * a[i + 4] + a[i + 5] * a[i + 5]
      + a[i + 6] * a[i + 6] + a[i + 7] * a[i + 7];

  return sqrt (res);
}

////////////////////////////////////////////////////////////////////////////////
// Compute euclidean scalar product of arrays
double
e2_scalar_product_array (int n, //!< array size
                         const double *a,       //!< array 1
                         const double *b        //!< array 2
  )
{
  int i, m;
  double res = 0.;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    res += a[i] * b[i];

  // Unrolled loop
  for (i = m; i < n; i += 8)
    res += a[i] * b[i] + a[i + 1] * b[i + 1]
      + a[i + 2] * b[i + 2] + a[i + 3] * b[i + 3]
      + a[i + 4] * b[i + 4] + a[i + 5] * b[i + 5]
      + a[i + 6] * b[i + 6] + a[i + 7] * b[i + 7];

  return res;
}

////////////////////////////////////////////////////////////////////////////////
// Subtract arrays: dst = src - dst
void
sub_arrays (int n,              // array size
            const double *src,  // input array
            double *dst         // output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    dst[i] = src[i] - dst[i];

  // Unrolled loop
  for (i = m; i < n; i += 8)
    {
      dst[i] = src[i] - dst[i];
      dst[i + 1] = src[i + 1] - dst[i + 1];
      dst[i + 2] = src[i + 2] - dst[i + 2];
      dst[i + 3] = src[i + 3] - dst[i + 3];
      dst[i + 4] = src[i + 4] - dst[i + 4];
      dst[i + 5] = src[i + 5] - dst[i + 5];
      dst[i + 6] = src[i + 6] - dst[i + 6];
      dst[i + 7] = src[i + 7] - dst[i + 7];
    }
}

/******************************************************************************/

// Scale array x
void
scale_array (int n,              // array size
             double *a,          // input/output array
             double mult
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  for (i = 0; i < m; i++)
    a[i] *= mult;

  // Unrolled loop
  for (i = m; i < n; i += 8)
    {
      a[i] *= mult;
      a[i + 1] *= mult;
      a[i + 2] *= mult;
      a[i + 3] *= mult;
      a[i + 4] *= mult;
      a[i + 5] *= mult;
      a[i + 6] *= mult;
      a[i + 7] *= mult;
    }
}

