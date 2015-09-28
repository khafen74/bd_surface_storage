#-------------------------------------------------
#
# Project created by QtCreator 2015-09-15T16:24:00
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = PotentialStorage
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp

win32: LIBS += -L$$PWD/../../../../../../../MinGW/msys/1.0/local/lib/ -llibgdal

INCLUDEPATH += $$PWD/../../../../../../../MinGW/msys/1.0/local/include
DEPENDPATH += $$PWD/../../../../../../../MinGW/msys/1.0/local/include
