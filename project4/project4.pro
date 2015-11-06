TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp

HEADERS += \
    lib.h

# MPI Settings
#QMAKE_CXX = /usr/lib64/openmpi/bin/mpic++
QMAKE_CXX = mpic++
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

#QMAKE_CFLAGS += $$system(/usr/lib64/openmpi/bin/mpic++ --showme:compile)
#QMAKE_LFLAGS += $$system(/usr/lib64/openmpi/bin/mpic++ --showme:link)
#QMAKE_CXXFLAGS += $$system(/usr/lib64/openmpi/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(/usr/lib64/openmpi/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK
