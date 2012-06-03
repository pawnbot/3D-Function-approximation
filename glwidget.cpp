#include <QtGui>
#include <QtOpenGL>
#include <iostream>

#include "glwidget.h"

const GLfloat LightPosition[]= { 0.0f, 0.0f, 10.0f, 1.0f };
const GLfloat LightAmbient[]= { 0.5f, 0.5f, 0.5f, 1.0f };
const GLfloat LightDiffuse[]= { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat gold_color[] = {0.85f, 0.85f, 0.10f};
const GLfloat blue_color[] = {0.137255f, 0.419608f, 0.556863f};
const GLfloat white_color[] = {1.0f, 1.0f, 1.0f};

GLWidget::GLWidget (QWidget *parent) :
  QGLWidget(parent)
{
  clearColor = Qt::black;
  function_index = SIN_FUNC;
  draw_mesh = false;
  options = ONLY_FUNCTION;
  least_squares = new approximator (this);
  set_default ();
  update_n1 (50);
  update_n2 (50);
  draw_approx = false;
  connect (least_squares, SIGNAL (finished ()), this, SLOT (approximation_done ()));
}

GLWidget::~GLWidget ()
{
}

QSize GLWidget::minimumSizeHint () const
{
  return QSize (320, 240);
}

QSize GLWidget::sizeHint() const
{
  return QSize (640, 480);
}

void GLWidget::rotate_by (int xAngle, int yAngle, int zAngle)
{
  xRot += xAngle;
  yRot += yAngle;
  zRot += zAngle;
  updateGL ();
}

void GLWidget::setClearColor (const QColor &color)
{
  clearColor = color;
  updateGL ();
}

void GLWidget::initializeGL ()
{
  glEnable (GL_DEPTH_TEST);
  glEnable (GL_CULL_FACE);
  glEnable (GL_NORMALIZE);
}

void GLWidget::paintGL ()
{
  qglClearColor (clearColor);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glLoadIdentity ();
  glTranslatef (0.0f, 0.0f, -10.0f);
  glRotatef (zRot / 16.0f, 0.0f, 0.0f, 1.0f);
  glRotatef (xRot / 16.0f, 1.0f, 0.0f, 0.0f);
  glRotatef (yRot / 16.0f, 0.0f, 1.0f, 0.0f);
  glScalef (scale, scale, scale);

  glDisable (GL_LIGHTING);
  /// Coordinate Lines
  glLineWidth(3);
  glBegin(GL_LINES);
    glColor4f (1, 0, 0, 1);
    glVertex3d(100, 0, 0);
    glVertex3d(-100, 0, 0);
  glEnd();
  qglColor(Qt::white);
  renderText(1, 0 , 0, QString::fromUtf8("X"), QFont() , 300);
  glBegin(GL_LINES);
    glColor4f (0, 1, 0, 1);
    glVertex3d(0, 100, 0);
    glVertex3d(0, -100, 0);
  glEnd();
  qglColor(Qt::white);
  renderText(0, 1 , 0, QString::fromUtf8("Y"), QFont() , 300);
  glBegin(GL_LINES);
    glColor4f (0, 0, 1, 1);
    glVertex3d(0, 0, 100);
    glVertex3d(0, 0, -100);
  glEnd();
  qglColor(Qt::white);
  renderText(0, 0 , 1, QString::fromUtf8("Z"), QFont() , 300);
  glLineWidth(1);
  /// Lights
  glEnable (GL_LIGHTING);
  glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
  glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
  glEnable (GL_LIGHT1);

  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, gold_color);
  if (options != ONLY_APPROXIMATION)
    draw_function (function_index, draw_mesh);
  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, blue_color);
  residual = 0.;
  if (options != ONLY_FUNCTION)
    {
      draw_approx = true;
      draw_function (function_index, draw_mesh);
      draw_approx = false;
    }
  emit sendmsg (trUtf8 ("Zoom %1, Residual: %3").arg (scale).arg (residual));
}

void GLWidget::draw_function (functions index, bool draw_mesh)
{
  functions old_index = function_index;
  function_index = index;
  const double step_x = (MAX_X - MIN_X) / (double) steps;
  const double step_y = (MAX_Y - MIN_Y) / (double) steps;
  GLfloat u[3], v[3], n[3];
  for (double x = MIN_X; x <= MAX_X; x += step_x)
    {
      for (double y = MIN_Y; y <= MAX_Y; y += step_y)
        {
          /// Upper triangle
          if (!(is_outside (x, y) || is_outside (x, y + step_y) ||
                is_outside (x + step_x, y + step_y)))
            {
              glBegin (GL_TRIANGLES);
              {
                u[0] = 0.0f;
                u[1] = step_y;
                u[2] = function (x + u[0], y + u[1]) - function (x, y);
                v[0] = step_x;
                v[1] = step_y;
                v[2] = function (x + v[0], y + v[1]) - function (x, y);
                cross_product (v, u, n);
                glNormal3fv (n);
                glVertex3f(x, y, function (x, y));
                glVertex3f(x, y + step_y, function (x, y + step_y));
                glVertex3f(x + step_x, y + step_y, function (x + step_x, y + step_y));
              }
              glEnd ();
            }
          /// Lower triangle
          if (!(is_outside (x, y) || is_outside (x + step_x, y) ||
                is_outside (x + step_x, y + step_y)))
            {
              glBegin (GL_TRIANGLES);
              {
                u[0] = step_x;
                u[1] = step_y;
                u[2] = function (x + u[0], y + u[1]) - function (x, y);
                v[0] = step_x;
                v[1] = 0;
                v[2] = function (x + v[0], y + v[1]) - function (x, y);
                cross_product (v, u, n);
                glNormal3fv (n);
                glVertex3f(x, y, function (x, y));
                glVertex3f(x + step_x, y + step_y, function (x + step_x, y + step_y));
                glVertex3f(x + step_x, y, function (x + step_x, y));
              }
              glEnd ();
            }
        }
    }
  if (draw_mesh)
    {
      glDisable (GL_LIGHTING);
      glColor3f (1, 1, 1);
      for (double x = MIN_X; x <= MAX_X; x += step_x)
        {
          for (double y = MIN_Y; y <= MAX_Y; y += step_y)
            {
              /// Upper triangle
              if (!(is_outside (x, y) || is_outside (x, y + step_y) ||
                    is_outside (x + step_x, y + step_y)))
                {
                  glBegin (GL_LINE_STRIP);
                  {
                    glVertex3f(x, y, function (x, y));
                    glVertex3f(x, y + step_y, function (x, y + step_y));
                    glVertex3f(x + step_x, y + step_y, function (x + step_x, y + step_y));
                    glVertex3f(x, y, function (x, y));
                  }
                  glEnd ();
                }
              /// Lower triangle
              if (!(is_outside (x, y) || is_outside (x + step_x, y) ||
                    is_outside (x + step_x, y + step_y)))
                {
                  glBegin (GL_LINE_STRIP);
                  {
                    glVertex3f(x, y, function (x, y));
                    glVertex3f(x + step_x, y + step_y, function (x + step_x, y + step_y));
                    glVertex3f(x + step_x, y, function (x + step_x, y));
                    glVertex3f(x, y, function (x, y));
                  }
                  glEnd ();
                }
            }
        }
      glEnable (GL_LIGHTING);
    }
  function_index = old_index;
}

void GLWidget::resizeGL(int width, int height)
{
  glViewport (0, 0, width, height);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  double max = qMax (MAJOR, MINOR);
  double min = qMin (-MAJOR, -MINOR);
  glOrtho (min * 1.05, max * 1.05, max * 1.05, min * 1.05, -100.0, 100.0);
  glMatrixMode (GL_MODELVIEW);
}

void GLWidget::set_default ()
{
  scale = 1.;
  xRot = 920;
  yRot = -112;
  zRot = 0;
  MIN_X = -MAJOR;
  MAX_X = MAJOR;
  MIN_Y = -MINOR;
  MAX_Y = MINOR;
  steps = 150;
  residual = 0.;
  zoom_out ();
}

void GLWidget::zoom_out ()
{
  MAX_X *= scale;
  MIN_X *= scale;
  MAX_Y *= scale;
  MIN_Y *= scale;

  scale *= 0.9;

  MAX_X /= scale;
  MIN_X /= scale;
  MAX_Y /= scale;
  MIN_Y /= scale;
}

void GLWidget::zoom_in ()
{
  MAX_X *= scale;
  MIN_X *= scale;
  MAX_Y *= scale;
  MIN_Y *= scale;

  scale *= 1.1;

  MAX_X /= scale;
  MIN_X /= scale;
  MAX_Y /= scale;
  MIN_Y /= scale;
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
  last_mouse_pos = event->pos ();
  emit clicked ();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
  last_mouse_pos = event->pos ();
  emit released ();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
  int dx = event->x () - last_mouse_pos.x ();
  int dy = event->y () - last_mouse_pos.y ();

  if (event->buttons () & Qt::LeftButton)
    {
      rotate_by (8 * dy, -8 * dx, 0);
    }
  else
    if (event->buttons () & Qt::RightButton)
      {
        rotate_by (8 * dy, 0, 8 * dx);
      }
  last_mouse_pos = event->pos ();
}

void GLWidget::wheelEvent (QWheelEvent *event)
{
  if (event->delta () > 0)
    {
      zoom_in ();
    }
  else
    {
      zoom_out ();
    }

  updateGL ();
}

void GLWidget::keyPressEvent (QKeyEvent *event)
{
  switch (event->key ())
    {
      case (Qt::Key_0):
        xRot = 0;
        yRot = 0;
        zRot = 0;
        break;
      case (Qt::Key_Plus):
        steps += 10;
        break;
      case (Qt::Key_Minus):
        steps -= 10;
        if (steps < 10)
          steps = 10;
        break;
      default:
        break;
    }
  updateGL ();
}

void GLWidget::mouseDoubleClickEvent (QMouseEvent *)
{
  set_default ();
  update_approximation ();
}

void GLWidget::mesh_draw_changed ()
{
  draw_mesh = !draw_mesh;
}

void GLWidget::only_func ()
{
  options = ONLY_FUNCTION;
}

void GLWidget::only_approx ()
{
  options = ONLY_APPROXIMATION;
}

void GLWidget::draw_both ()
{
  options = DRAW_BOTH;
}

void GLWidget::update_n1 (int count)
{
  count = qMax (count, 10);
  n1 = count;
}

void GLWidget::update_n2 (int count)
{
  count = qMax (count, 10);
  n2 = count;
}

void GLWidget::update_approximation ()
{
  least_squares->stop ();
  least_squares->wait ();
  least_squares->update_values (n1, n2, function_index);
  least_squares->start ();
}

void GLWidget::approximation_done ()
{
  least_squares->update_answer ();
  /// new results are aviable for GUI
  updateGL ();
}
