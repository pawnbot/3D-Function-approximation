#ifndef POLYNOM_H
#define POLYNOM_H

#include <vector>

class point_t
{
public:
  point_t (double x = 0, double y = 0) : x(x), y(y) {}
  double x, y;
};

using std::vector;
using std::min;
using std::max;

class polynom
{
public:
  polynom (size_t deg_x = 0, size_t deg_y = (size_t) -1);
  /// P = a * x + b * y + c
  polynom (double a, double b, double c);
  void resize (size_t deg_x, size_t deg_y = (size_t) -1);
  double eval (double x, double y);
  double eval (point_t u)
  {
    return eval (u.x, u.y);
  }

  size_t deg_x () const;
  size_t deg_y () const;
  friend inline polynom operator* (const polynom& a, const polynom& b)
  {
    polynom c (a.deg_x () + b.deg_x (), a.deg_y () + b.deg_y ());
    for (size_t ax = 0; ax <= a.deg_x (); ++ax)
      for (size_t ay = 0; ay <= a.deg_y (); ++ay)
        {
          for (size_t bx = 0; bx <= b.deg_x (); ++bx)
            for (size_t by = 0; by <= b.deg_y (); ++by)
              {
                c.coef[ax + bx][ay + by] += a.coef[ax][ay] * b.coef[bx][by];
              }
        }
    return c;
  }
  friend inline polynom operator+ (const polynom& a, const polynom& b)
  {
    polynom c (max (a.deg_x (), b.deg_x ()), max (a.deg_y (), b.deg_y ()));
    for (size_t ax = 0; ax <= a.deg_x (); ++ax)
      for (size_t ay = 0; ay <= a.deg_y (); ++ay)
        {
          c.coef[ax][ay] += a.coef[ax][ay];
        }
    for (size_t bx = 0; bx <= b.deg_x (); ++bx)
      for (size_t by = 0; by <= b.deg_y (); ++by)
        {
          c.coef[bx][by] += b.coef[bx][by];
        }
    return c;
  }
  double integrate (point_t a, point_t b, point_t c);
  double integrate (double x_begin, double x_end, double y_begin, double y_end);
  void scale (double mult);
  void div (double mult);

  void print ();

private:
  double lower_triangle_int (double x1, double x2, double y1, double y2);
  double upper_triangle_int (double x1, double x2, double y1, double y2);
  double lower_triangle_int2 (double x1, double x2, double y1, double y2);
  double upper_triangle_int2 (double x1, double x2, double y1, double y2);
  vector<vector<double> > coef;
};

#endif // POLYNOM_H
