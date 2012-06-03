#ifndef APPROXIMATOR_H
#define APPROXIMATOR_H

#include <QThread>
#include <omp.h>
#include <cmath>
#include <cstdio>

#include "polynom.h"

//#define DEBUG_OUTPUT

enum functions
{
  APPRPOX_FUNC = -1,
  SIN_FUNC = 0,
  x2y2_FUNC = 1,
  LINEAR_FUNC = 2,
  X2X = 3
};

struct answer
{
  double* res;
  double  n1, n2;
};

const double MAJOR = 8;
const double MINOR = 7;

static inline bool
is_inside (point_t u)
{
  if (-1.e-15 < u.x && u.x < MAJOR + 1.e-15 && -1.e-15 < u.y && u.y < 2 * M_PI + 1.e-15)
    return true;
  else
    return false;
}

inline double
scal (const point_t& a, const point_t& b)
{
  return a.x * b.x + a.y * b.y;
}

class approximator : public QThread
{
  Q_OBJECT
public:
  explicit approximator(QObject *parent = 0);
  ~approximator ();

  void run ();

  void update_values (
    int n1,
    int n2,
    functions function_index
    );

  inline double function (double radius, double angle)
  {
    double x = radius * cos (angle);
    double y = radius * sin (angle) * MINOR / MAJOR;
    switch (function_index)
      {
        case (SIN_FUNC):
          return sin (x * x + y * y);
        default:
        case (x2y2_FUNC):
          return (x * x + y * y);
        case (LINEAR_FUNC):
          return x + y;
        case (X2X):
          return x * x + x;
      }
  }

  inline polynom
  L (point_t a, point_t b)
  {
    /// L = (x - u1) * (v2 - v1) - (y - v1) * (u2 - u1) = {v2 - v1, u1 - u2, u1(v1-v2) + v1(u2-u1)}
    return polynom (b.y - a.y, a.x - b.x, a.x * (a.y - b.y) + a.y * (b.x - a.x));
  }

  inline polynom
  psi (int index, point_t p1, point_t p2, point_t p3)
  {
    switch (index)
    {
      default:
#ifdef DEBUG_OUTPUT
        printf ("PSI FUNCTION FAIL!\n");
#endif
        return polynom (1.e100, 1.e100, 1.e100);
      case (1):
        return L (p2, p3);
      case (2):
        return L (p1, p3);
      case (3):
        return L (p1, p2);
    }
  }

  inline polynom
  phi (int index, point_t p1, point_t p2, point_t p3)
  {
    polynom psi_pol = psi (index, p1, p2, p3);
    switch (index)
    {
      default:
#ifdef DEBUG_OUTPUT
        printf ("PHI FUNCTION FAIL!\n");
#endif
        psi_pol.scale (1.e100);
        break;
      case (1):
        psi_pol.div (psi_pol.eval (p1.x, p1.y));
        break;
      case (2):
        psi_pol.div (psi_pol.eval (p2.x, p2.y));
        break;
      case (3):
        psi_pol.div (psi_pol.eval (p3.x, p3.y));
        break;
    }
#ifdef DEBUG_OUTPUT
    switch (index)
    {
      default:
      case (1):
        if (fabs (psi_pol.eval (p1.x, p1.y) - 1.) > 1.e-12 ||
            fabs (psi_pol.eval (p2.x, p2.y) - 0.) > 1.e-12 ||
            fabs (psi_pol.eval (p3.x, p3.y) - 0.) > 1.e-12)
          printf ("PSI_POL FAIL: 1\n");
        break;
      case (2):
        if (fabs (psi_pol.eval (p1.x, p1.y) - 0.) > 1.e-12 ||
            fabs (psi_pol.eval (p2.x, p2.y) - 1.) > 1.e-12 ||
            fabs (psi_pol.eval (p3.x, p3.y) - 0.) > 1.e-12)
          printf ("PSI_POL FAIL: 2\n");
        break;
      case (3):
        if (fabs (psi_pol.eval (p1.x, p1.y) - 0.) > 1.e-12 ||
            fabs (psi_pol.eval (p2.x, p2.y) - 0.) > 1.e-12 ||
            fabs (psi_pol.eval (p3.x, p3.y) - 1.) > 1.e-12)
          printf ("PSI_POL FAIL: 3\n");
        break;
    }
#endif
    return psi_pol;
  }

  inline polynom
  func_pol (point_t p1, point_t p2, point_t p3)
  {
    polynom result;
    polynom tmp;
    /// 1
    tmp = phi (1, p1, p2, p3);
    tmp.scale (function (p1.x, p1.y));
    result = tmp;
    /// 2
    tmp = phi (2, p1, p2, p3);
    tmp.scale (function (p2.x, p2.y));
    result = result + tmp;
    /// 3
    tmp = phi (3, p1, p2, p3);
    tmp.scale (function (p3.x, p3.y));
    result = result + tmp;
#ifdef DEBUG_OUTPUT
    if (fabs (result.eval (p1.x, p1.y) - function (p1.x, p1.y)) > 1.e-10)
      {
        printf ("FUNC_POL FAIL\n");
        printf ("%.16lf\n", phi (1, p1, p2, p3).eval (p1));
        phi (1, p1, p2, p3).print ();
        printf ("%.16lf\n", phi (2, p1, p2, p3).eval (p2));
        phi (2, p1, p2, p3).print ();
        printf ("%.16lf\n", phi (3, p1, p2, p3).eval (p3));
        phi (3, p1, p2, p3).print ();
      }
    if (fabs (result.eval (p2.x, p2.y) - function (p2.x, p2.y)) > 1.e-10)
      printf ("FUNC_POL FAIL\n");
    if (fabs (result.eval (p3.x, p3.y) - function (p3.x, p3.y)) > 1.e-10)
      printf ("FUNC_POL FAIL\n");
#endif
    return result;
  }

  inline polynom
  func_pol (double x1, double y1, double x2, double y2, double x3, double y3)
  {
    return func_pol (point_t (x1, y1), point_t (x2, y2), point_t (x3, y3));
  }

  inline polynom
  kurant_pol (point_t p, point_t p1, point_t p2)
  {
    return phi (1, p, p1, p2);
  }

  inline double
  rhs_integrate (int i, int j)
  {
    double result = 0.;
    polynom rhs_pol;
    point_t p (i * h1, j * h2);
    point_t p1, p2;
    /// 1
    p1 = point_t (i * h1, (j + 1) * h2);
    p2 = point_t ((i + 1) * h1, (j + 1) * h2);
    rhs_pol = jacobian * kurant_pol (p, p1, p2) * func_pol (p, p1, p2);
    if (is_inside (p1) && is_inside (p2))
      result += rhs_pol.integrate (p, p1, p2);
    /// 2
    p1 = p2;
    p2 = point_t ((i + 1) * h1, j * h2);
    rhs_pol = jacobian * kurant_pol (p, p1, p2) * func_pol (p, p1, p2);
    if (is_inside (p1) && is_inside (p2))
      result += rhs_pol.integrate (p, p1, p2);
    /// 3
    p1 = p2;
    p2 = point_t (i * h1, (j - 1) * h2);
    rhs_pol = jacobian * kurant_pol (p, p1, p2) * func_pol (p, p1, p2);
    if (is_inside (p1) && is_inside (p2))
      result += rhs_pol.integrate (p, p1, p2);
    /// 4
    p1 = p2;
    p2 = point_t ((i - 1) * h1, (j - 1) * h2);
    rhs_pol = jacobian * kurant_pol (p, p1, p2) * func_pol (p, p1, p2);
    if (is_inside (p1) && is_inside (p2))
      result += rhs_pol.integrate (p, p1, p2);
    /// 5
    p1 = p2;
    p2 = point_t ((i - 1) * h1, j * h2);
    rhs_pol = jacobian * kurant_pol (p, p1, p2) * func_pol (p, p1, p2);
    if (is_inside (p1) && is_inside (p2))
      result += rhs_pol.integrate (p, p1, p2);
    /// 6
    p1 = p2;
    p2 = point_t (i * h1, (j + 1) * h2);
    rhs_pol = jacobian * kurant_pol (p, p1, p2) * func_pol (p, p1, p2);
    if (is_inside (p1) && is_inside (p2))
      result += rhs_pol.integrate (p, p1, p2);
#ifdef DEBUG_OUTPUT
    if (fabs (result) > 1.e3)
      printf ("RHS INTEGRAL FAIL!\n");
    //printf ("%g | %g | %lf\n", i * h1, j * h2, result);
#endif
    return result;
  }


  inline double
  l2_prod (int i, int j, int i2, int j2)
  {
    double result = 0.;
    if (get_index (i, j, n1, n2) < 0 || get_index (i2, j2, n1, n2) < 0)
      return result;
    point_t pij (i * h1, j * h2);
    point_t p2 (i2 * h1, j2 * h2);
    point_t p;
    polynom matrix_pol;
    /// Diagonal
    if (i == i2 && j == j2)
      {
        /// 1
        p2 = point_t ((i + 0) * h1, (j + 1) * h2);
        p = point_t ((i + 1) * h1, (j + 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p2, p) * kurant_pol (pij, p2, p);
        if (is_inside (p) && is_inside (p2))
          result += matrix_pol.integrate (p, pij, p2);
        /// 2
        p2 = p;
        p = point_t ((i + 1) * h1, (j + 0) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p2, p) * kurant_pol (pij, p2, p);
        if (is_inside (p) && is_inside (p2))
          result += matrix_pol.integrate (p, pij, p2);
        /// 3
        p2 = p;
        p = point_t ((i + 0) * h1, (j - 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p2, p) * kurant_pol (pij, p2, p);
        if (is_inside (p) && is_inside (p2))
          result += matrix_pol.integrate (p, pij, p2);
        /// 4
        p2 = p;
        p = point_t ((i - 1) * h1, (j - 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p2, p) * kurant_pol (pij, p2, p);
        if (is_inside (p) && is_inside (p2))
          result += matrix_pol.integrate (p, pij, p2);
        /// 5
        p2 = p;
        p = point_t ((i - 1) * h1, (j + 0) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p2, p) * kurant_pol (pij, p2, p);
        if (is_inside (p) && is_inside (p2))
          result += matrix_pol.integrate (p, pij, p2);
        /// 6
        p2 = p;
        p = point_t ((i + 0) * h1, (j + 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p2, p) * kurant_pol (pij, p2, p);
        if (is_inside (p) && is_inside (p2))
          result += matrix_pol.integrate (p, pij, p2);
      }
    /// Top
    if (i == i2 && j + 1 == j2)
      {
        /// 6
        p = point_t ((i - 1) * h1, (j + 0) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
        /// 1
        p = point_t ((i + 1) * h1, (j + 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
      }
    /// Top Right
    if (i + 1 == i2 && j + 1 == j2)
      {
        /// 1
        p = point_t ((i + 0) * h1, (j + 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
        /// 2
        p = point_t ((i + 1) * h1, (j + 0) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
      }
    /// Right
    if (i + 1 == i2 && j == j2)
      {
        /// 2
        p = point_t ((i + 1) * h1, (j + 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
        /// 3
        p = point_t ((i + 0) * h1, (j - 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
      }
    /// Bottom
    if (i == i2 && j - 1 == j2)
      {
        /// 3
        p = point_t ((i + 1) * h1, (j + 0) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
        /// 4
        p = point_t ((i - 1) * h1, (j - 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
      }
    /// Bottom Left
    if (i - 1 == i2 && j - 1 == j2)
      {
        /// 4
        p = point_t ((i + 0) * h1, (j - 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
        /// 5
        p = point_t ((i - 1) * h1, (j + 0) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
      }
    /// Left
    if (i - 1 == i2 && j == j2)
      {
        /// 5
        p = point_t ((i - 1) * h1, (j - 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
        /// 6
        p = point_t ((i + 0) * h1, (j + 1) * h2);
        matrix_pol = jacobian * kurant_pol (pij, p, p2) * kurant_pol (p2, p, pij);
        if (is_inside (p))
          result += matrix_pol.integrate (p, pij, p2);
      }
#ifdef DEBUG_OUTPUT
    if (fabs (result) > 1.e3)
      {
        printf ("L2 INTEGRAL FAIL in points pair: %d, %d | %d, %d\n", i, j, i2, j2);
      }
#endif
    return result;
  }
  void clean_up ();
  inline int get_index (int i, int j, int n1, int n2)
  {
    if (0 <= i && i < n1 && 0 <= j && j < n2)
      return i * n2 + j;
    else
      return -1;
  }
  inline double evaluate_phi (int i, int j, double h1, double h2, double x, double y)
  {
    int index = get_index (i, j, old_answer.n1, old_answer.n2);
    if (index < 0)
      {
#ifdef DEBUG_OUTPUT
        printf ("TRASH INDEX!\n");
#endif
        return 1.e100;
      }
    const double phi_x = i * h1;
    const double phi_y = j * h2;
    double vec_x;
    /// 1
    if (phi_x <= x && phi_y <= y && (x - phi_x) * h2 <= (y - phi_y) * h1)
      {
        if (y - phi_y < h2)
          return (1. - (y - phi_y) / h2) * old_answer.res[index];
        else
          return 0.;
      }
    /// 4
    if (phi_x >= x && phi_y >= y && (phi_x - x) * h2 <= (phi_y - y) * h1)
      {
        if (phi_y - y < h2)
          return (1. - (phi_y - y) / h2) * old_answer.res[index];
        else
          return 0.;
      }
    /// 2
    if (phi_x <= x && phi_y <= y && (x - phi_x) * h2 >= (y - phi_y) * h1)
      {
        if (x - phi_x < h1)
          return (1. - (x - phi_x) / h1) * old_answer.res[index];
        else
          return 0.;
      }
    /// 5
    if (phi_x >= x && phi_y >= y && (phi_x - x) * h2 >= (phi_y - y) * h1)
      {
        if (phi_x - x < h1)
          return (1. - (phi_x - x) / h1) * old_answer.res[index];
        else
          return 0.;
      }
    /// 3
    if (phi_x <= x && phi_y >= y)
      {
        vec_x = (phi_y - y) / h2;
        vec_x = x + vec_x * h1;
        if (vec_x - phi_x < h1)
          return (1. - (vec_x - phi_x) / h1) * old_answer.res[index];
        else
          return 0.;
      }
    /// 6
    if (phi_x >= x && phi_y <= y)
      {
        vec_x = (phi_y - y) / h2;
        vec_x = x + vec_x * h1;
        if (phi_x - vec_x  < h1)
          return (1. - (phi_x - vec_x) / h1) * old_answer.res[index];
        else
          return 0.;
      }
#ifdef DEBUG_OUTPUT
    printf ("PHI EVALUATE FAIL!\n");
#endif
    return 1.e100;
  }
  double evaluate (double x, double y);

signals:

public slots:
  void set_threads_count (int count);
  void stop ();
  void update_answer ();

private:
  polynom jacobian;
  int* matrix_rnz;
  double* matrix;
  double* rhs;
  double* res;
  double* workspace;
  answer old_answer;

  functions function_index;
  double PRECISION;
  int threads_count;
  int n1, n2;
  double h1, h2;
  bool solver_done;
};

#endif // APPROXIMATOR_H
