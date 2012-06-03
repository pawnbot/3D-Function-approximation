#include "main_window.h"
#include "ui_main_window.h"

main_window::main_window (QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::main_window)
{
  ui->setupUi(this);
  ui->statusBar->showMessage (trUtf8 ("Выберите функцию"));
  setWindowTitle (trUtf8 ("Метод наименьших квадратов"));
  ui->function_menu->setTitle (trUtf8 ("Функции"));
  ui->menuView->setTitle (trUtf8 ("Показать"));
  ui->menu_func->setText (trUtf8 ("Только функцию"));
  ui->menu_approx->setText (trUtf8 ("Только приближение"));
  ui->menu_both->setText (trUtf8 ("Функцию и приближение"));
  ui->menu_mesh->setText (trUtf8 ("Сетку"));
  ui->threads_label->setText (trUtf8 ("Число потоков:"));

  QColor clearColor = Qt::black;

  browser = new GLWidget(this);
  browser->setClearColor(clearColor);
  connect (ui->sin, SIGNAL (triggered ()), browser, SLOT (set_sin ()));
  connect (ui->linear, SIGNAL (triggered ()), browser, SLOT (set_linear ()));
  connect (ui->x2y2, SIGNAL (triggered ()), browser, SLOT (set_x2y2 ()));
  connect (ui->x2x, SIGNAL (triggered ()), browser, SLOT (set_x2x ()));
  connect (ui->menu_mesh, SIGNAL (triggered ()), browser, SLOT (mesh_draw_changed ()));
  connect (ui->menu_func, SIGNAL (triggered ()), browser, SLOT (only_func ()));
  connect (ui->menu_approx, SIGNAL (triggered ()), browser, SLOT (only_approx ()));
  connect (ui->menu_both, SIGNAL (triggered ()), browser, SLOT (draw_both ()));

  connect (browser, SIGNAL (sendmsg (QString)), ui->statusBar, SLOT (showMessage (QString)));
  ui->centralWidget->layout ()->addWidget (browser);
  browser->setFocus ();
  old_n1 = ui->n1_box->value ();
  old_n2 = ui->n2_box->value ();
  browser->update_approximation ();
}

main_window::~main_window ()
{
  delete ui;
}

void main_window::on_threads_box_valueChanged(int arg1)
{
  browser->least_squares->set_threads_count (arg1);
}

void main_window::on_n1_box_editingFinished()
{
  browser->update_n1 (ui->n1_box->value ());
  browser->setFocus ();
  if (old_n1 != ui->n1_box->value ())
    browser->update_approximation ();
  old_n1 = ui->n1_box->value ();
}

void main_window::on_n2_box_editingFinished()
{
  browser->update_n2 (ui->n2_box->value ());
  browser->setFocus ();
  if (old_n2 != ui->n2_box->value ())
    browser->update_approximation ();
  old_n2 = ui->n2_box->value ();
}
