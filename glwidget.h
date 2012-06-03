#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QtGui>
#include <QGLWidget>

#include "approximator.h"

enum draw_options
{
  ONLY_FUNCTION = 0,
  ONLY_APPROXIMATION = 1,
  DRAW_BOTH = 2
};

static inline bool
is_outside (double x, double y)
{
  if (x * x / (MAJOR * MAJOR) + y * y / (MINOR * MINOR) >= 1)
    return true;
  else
    return false;
}

/* cross_product: compute the cross product of two vectors
 *
 * u - array of 3 GLfloats (GLfloat u[3])
 * v - array of 3 GLfloats (GLfloat v[3])
 * n - array of 3 GLfloats (GLfloat n[3]) to return the cross product in
 */
static inline void
cross_product (GLfloat* u, GLfloat* v, GLfloat* n)
{
  n[0] = u[1] * v[2] - u[2] * v[1];
  n[1] = u[2] * v[0] - u[0] * v[2];
  n[2] = u[0] * v[1] - u[1] * v[0];
}

class GLWidget : public QGLWidget
{
  Q_OBJECT

public:
  GLWidget (QWidget *parent = 0);
  ~GLWidget ();

  void zoom_out ();
  void zoom_in ();
  QSize minimumSizeHint () const;
  QSize sizeHint () const;
  void setClearColor (const QColor &color);
  void set_default ();
  void draw_function (functions index = SIN_FUNC, bool draw_mesh = false);
  inline double function (double x, double y)
  {
    if (!draw_approx)
      {
        switch (function_index)
          {
            case (SIN_FUNC):
              return sin (x * x + y * y);
            default:
            case (x2y2_FUNC):
              return (x * x + y * y);
            case (LINEAR_FUNC):
              return x + y;;
            case (X2X):
              return x * x + x;
          }
      }
    else
      {
        double result = least_squares->evaluate (x, y);
        switch (function_index)
          {
            case (SIN_FUNC):
              residual = qMax (residual, fabs (sin (x * x + y * y) - result));
              return result;
            default:
            case (x2y2_FUNC):
              residual = qMax (residual,fabs ((x * x + y * y) - result));
              return result;
            case (LINEAR_FUNC):
              residual = qMax (residual, fabs (x + y - result));
              return result;
            case (X2X):
              residual = qMax (residual, fabs (x * x + x - result));
              return result;
          }
      }
  }

  void update_approximation ();
  approximator *least_squares;

public slots:
  void rotate_by (int xAngle = 8, int yAngle = 8, int zAngle = 2);
  void set_sin () { function_index = SIN_FUNC; update_approximation ();}
  void set_x2y2 () { function_index = x2y2_FUNC; update_approximation ();}
  void set_linear () { function_index = LINEAR_FUNC; update_approximation ();}
  void set_x2x () { function_index = X2X; update_approximation ();}
  void mesh_draw_changed ();
  void only_func ();
  void only_approx ();
  void draw_both ();
  void update_n1 (int count);
  void update_n2 (int count);
  void approximation_done ();

signals:
  void sendmsg (QString);
  void clicked ();
  void released ();

protected:
  void initializeGL ();
  void paintGL ();
  void resizeGL (int width, int height);
  void mousePressEvent (QMouseEvent *event);
  void mouseMoveEvent (QMouseEvent *event);
  void mouseReleaseEvent (QMouseEvent *event);
  void wheelEvent (QWheelEvent *event);
  void mouseDoubleClickEvent (QMouseEvent *);
  void keyPressEvent (QKeyEvent *event);

private:
  double residual;
  bool draw_approx;
  draw_options options;
  functions function_index;
  QColor clearColor;
  QPoint last_mouse_pos;
  int xRot;
  int yRot;
  int zRot;
  double scale;
  int steps;
  int n1;
  int n2;
  double MIN_X;
  double MAX_X;
  double MIN_Y;
  double MAX_Y;
  bool draw_mesh;
};

#endif
