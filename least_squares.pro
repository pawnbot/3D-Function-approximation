#-------------------------------------------------
#
# Project created by QtCreator 2012-04-06T15:12:08
#
#-------------------------------------------------

QT       += core gui
QT       += opengl

LIBS           += -fopenmp
QMAKE_CXXFLAGS += -fopenmp

TARGET   = least_squares
TEMPLATE = app


SOURCES += main.cpp\
        main_window.cpp \
    glwidget.cpp \
    approximator.cpp \
    ilu.cpp \
    bcgstab.cpp \
    array_op.cpp \
    polynom.cpp

HEADERS  += main_window.h \
    glwidget.h \
    approximator.h \
    precond.h \
    iter_alg.h \
    array_op.h \
    polynom.h

FORMS    += main_window.ui
