
bin_PROGRAMS = eigen_demo1 eigen_demo2
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

eigen_demo1_SOURCES =	eigen_demo1.cpp
eigen_demo2_SOURCES =	eigen_demo2.cpp

eigen_demo1_LDADD  = $(top_builddir)/src/libofeli.a
eigen_demo2_LDADD  = $(top_builddir)/src/libofeli.a

eigen_demo1dir   = $(datadir)/ofeli/demos/solvers/EIGEN
eigen_demo2dir   = $(datadir)/ofeli/demos/solvers/EIGEN

eigen_demo1_DATA = eigen_demo1.cpp
eigen_demo2_DATA = eigen_demo2.cpp \
                   test.m \
                   eigen_demo2.dat \
                   README

dist_eigen_demo1_DATA =	eigen_demo1.cpp
dist_eigen_demo2_DATA =	eigen_demo2.cpp \
                        test.m \
                        eigen_demo2.dat \
                        README

clean-local:
	-rm -rf *.pos
