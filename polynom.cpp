#include <cstdio>
#include <cmath>

#include "polynom.h"

polynom::polynom (size_t deg_x, size_t deg_y)
{
  resize (deg_x, deg_y);
}

inline double po (double a, size_t k)
{
  double result = 1.;
  for (size_t i = 0; i < k; ++i)
    result *= a;
  return result;
}

polynom::polynom (double a, double b, double c)
{
  resize (1, 1);
  coef[0][0] = c;
  coef[0][1] = b;
  coef[1][0] = a;
  coef[1][1] = 0.;
}

void polynom::resize (size_t deg_x, size_t deg_y)
{
  if (deg_y == (size_t) -1)
    deg_y = deg_x;
  if (deg_x > 0 && deg_y > 0)
    {
      coef.resize (deg_x + 1);
      for (size_t i = 0; i <= deg_x; ++i)
        coef[i] = vector<double> (deg_y + 1, 0.);
    }
}

double polynom::eval (double x, double y)
{
  double result = 0;
  double dx = 1;
  for (size_t i = 0; i <= deg_x (); ++i)
    {
      double dy = 1;
      for (size_t j = 0; j <= deg_y (); ++j)
        {
          if (fabs (coef[i][j]) > 1.e-32)
            result += coef[i][j] * dx * dy;
          dy *= y;
        }
      dx *= x;
    }
  return result;
}

size_t polynom::deg_x () const
{
  return coef.size () - 1;
}

size_t polynom::deg_y () const
{
  if (deg_x () > 0)
    return coef[0].size () - 1;
  else
    return 0;
}

void polynom::print ()
{
  bool printed = false;
  for (int i = deg_x (); i >= 0; --i)
    for (int j = deg_y (); j >= 0; --j)
      {
        if (fabs (coef[i][j]) > 1.e-32)
          {
            if (coef[i][j] >= 0. && printed)
              printf (" + ");
            if (coef[i][j] < 0.)
              printf (" - ");
            printf ("%.3lf", fabs (coef[i][j]));
            printed = true;
            if (i > 0)
              printf (" x^%d", i);
            if (j > 0)
              printf (" y^%d", j);
          }
      }
  printf ("\n");
}

double polynom::integrate (double x_begin, double x_end, double y_begin, double y_end)
{
  double result = 0;
  double dx_begin = x_begin;
  double dx_end = x_end;
  for (size_t i = 0; i <= deg_x (); ++i)
    {
      double dy_begin = y_begin;
      double dy_end = y_end;
      for (size_t j = 0; j <= deg_y (); ++j)
        {
          if (fabs (coef[i][j]) > 1.e-32)
            result += coef[i][j] * (dx_end - dx_begin) * (dy_end - dy_begin) / ((i + 1) * (j + 1));
          dy_begin *= y_begin;
          dy_end *= y_end;
        }
      dx_begin *= x_begin;
      dx_end *= x_end;
    }
  return result;
}

void polynom::scale (double mult)
{
  for (size_t i = 0; i <= deg_x (); ++i)
    {
      for (size_t j = 0; j <= deg_y (); ++j)
        {
          if (fabs (coef[i][j]) > 1.e-32)
            coef[i][j] *= mult;
        }
    }
}

void polynom::div (double mult)
{
  for (size_t i = 0; i <= deg_x (); ++i)
    {
      for (size_t j = 0; j <= deg_y (); ++j)
        {
          if (fabs (coef[i][j]) > 1.e-32)
            coef[i][j] /= mult;
        }
    }
}

double polynom::lower_triangle_int (double x1, double x2, double y1, double y2)
{
  const double m = (y2 - y1) / (x2 - x1);
  double result = 0.;
  for (size_t i = 0; i <= deg_x (); ++i)
    {
      for (size_t j = 0; j <= deg_y (); ++j)
        {
          if (fabs (coef[i][j]) > 1.e-32)
            {
              if (i == 1 && j == 0)
                {
                  result += (coef[i][j]/6)*m*po(x1-x2, 2)*(x1+2*x2);
                }
              else
              if (i == 1 && j == 1)
                {
                  result += (coef[i][j]/24)*m*po(x1-x2, 2)*(-m*(x1-x2)*(x1+3*x2)+4*(x1+2*x2)*y1);
                }
              else
              if (i == 2 && j == 0)
                {
                  result += (coef[i][j]/12)*m*(po(x1, 4) - 4*x1*po(x2, 3)+3*po(x2, 4));
                }
              else
              if (i == 1 && j == 2)
                {
                  result += (coef[i][j]/60)*m*po(x1-x2, 2) * (po(m, 2)*po(x1-x2, 2)*(x1+4*x2)-5*m*(x1-x2)*(x1+3*x2)*y1+10*(x1+2*x2)*po(y1,2));
                }
              else
              if (i == 2 && j == 1)
                {
                  result += (coef[i][j]/60)*m*(-m*po(x1-x2, 3)*(po(x1,2)+3*x1*x2+6*po(x2,2))+5*(po(x1,4)-4*x1*po(x2,3)+3*po(x2,4))*y1);
                }
              else
              if (i == 3 && j == 0)
                {
                  result += (coef[i][j]/20)*m*(po(x1,5)-5*x1*po(x2,4)+4*po(x2,5));
                }
              else
                {
                  print ();
                  printf ("Integration of this polynom not supported!\n");
                  return 1.e100;
                }
            }
        }
    }
  return result;
}

double polynom::upper_triangle_int (double x1, double x2, double y1, double y2)
{
  return integrate (x1, x2, y1, y2) - lower_triangle_int (x1, x2, y1, y2);
}

double polynom::upper_triangle_int2 (double x1, double x2, double y1, double y2)
{
  return integrate (x1, x2, y1, y2) - lower_triangle_int2 (x1, x2, y1, y2);
}

double polynom::lower_triangle_int2 (double x1, double x2, double y1, double y2)
{
  printf ("Triangle (%g, %g), (%g, %g), (%g, %g) not supported!\n",
          x1, y2, x1, y1, x2, y1);
  return 1.e100;
}

double polynom::integrate (point_t a, point_t b, point_t c)
{
  double x1, x2;
  double y1, y2;
  const double EPS = 1.e-14;
  if (fabs (a.y - b.y) < EPS)
    {
      x1 = min (a.x, b.x);
      x2 = max (a.x, b.x);
      y1 = min (a.y, c.y);
      y2 = max (a.y, c.y);
      if (c.y < a.y)
        {
          if (fabs (c.x - x1) < EPS)
            return upper_triangle_int (x1, x2, y1, y2);
          if (fabs (c.x - x2) < EPS)
            return upper_triangle_int2 (x1, x2, y1, y2);
        }
      else
        {
          if (fabs (c.x - x2) < EPS)
            return lower_triangle_int (x1, x2, y1, y2);
          if (fabs (c.x - x1) < EPS)
            return lower_triangle_int2 (x1, x2, y1, y2);
        }
    }
  if (fabs (a.y - c.y) < EPS)
    {
      x1 = min (a.x, c.x);
      x2 = max (a.x, c.x);
      y1 = min (a.y, b.y);
      y2 = max (a.y, b.y);
      if (b.y < a.y)
        {
          if (fabs (b.x - x1) < EPS)
            return upper_triangle_int (x1, x2, y1, y2);
          if (fabs (b.x - x2) < EPS)
            return upper_triangle_int2 (x1, x2, y1, y2);
        }
      else
        {
          if (fabs (b.x - x2) < EPS)
            return lower_triangle_int (x1, x2, y1, y2);
          if (fabs (b.x - x1) < EPS)
            return lower_triangle_int2 (x1, x2, y1, y2);
        }
    }
  if (fabs (b.y - c.y) < EPS)
    {
      x1 = min (b.x, c.x);
      x2 = max (b.x, c.x);
      y1 = min (b.y, a.y);
      y2 = max (b.y, a.y);
      if (a.y < b.y)
        {
          if (fabs (a.x - x1) < EPS)
            return upper_triangle_int (x1, x2, y1, y2);
          if (fabs (a.x - x2) < EPS)
            return upper_triangle_int2 (x1, x2, y1, y2);
        }
      else
        {
          if (fabs (a.x - x2) < EPS)
            return lower_triangle_int (x1, x2, y1, y2);
          if (fabs (a.x - x1) < EPS)
            return lower_triangle_int2 (x1, x2, y1, y2);
        }
    }
  printf ("Cannot classify triangle!\n");
  return 1.e100;
}
