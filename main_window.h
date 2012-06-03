#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>
#include "glwidget.h"

namespace Ui
{
  class main_window;
}

class main_window : public QMainWindow
{
  Q_OBJECT

public:
  explicit main_window (QWidget *parent = 0);
  ~main_window ();

private slots:
  void on_threads_box_valueChanged(int arg1);

  void on_n1_box_editingFinished();

  void on_n2_box_editingFinished();

private:
  int old_n1, old_n2;
  Ui::main_window *ui;
  GLWidget *browser;
};

#endif // MAIN_WINDOW_H
