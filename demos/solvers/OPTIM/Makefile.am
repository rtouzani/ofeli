
bin_PROGRAMS = opt_demo1 opt_demo2 opt_demo3
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

opt_demo1_SOURCES =	opt_demo1.cpp
opt_demo2_SOURCES =	opt_demo2.cpp
opt_demo3_SOURCES =	opt_demo3.cpp

opt_demo1_LDADD  = $(top_builddir)/src/libofeli.a
opt_demo2_LDADD  = $(top_builddir)/src/libofeli.a
opt_demo3_LDADD  = $(top_builddir)/src/libofeli.a

opt_demo1dir   =	$(datadir)/ofeli/demos/solvers/OPTIM
opt_demo2dir   =	$(datadir)/ofeli/demos/solvers/OPTIM
opt_demo3dir   =	$(datadir)/ofeli/demos/solvers/OPTIM

opt_demo1_DATA =	opt_demo1.cpp
opt_demo2_DATA =	opt_demo2.cpp
opt_demo3_DATA =	opt_demo3.cpp Opt.h test.dat test.m Makefile.standalone README

dist_opt_demo1_DATA = opt_demo1.cpp
dist_opt_demo2_DATA = opt_demo2.cpp
dist_opt_demo3_DATA = opt_demo3.cpp Opt.h test.dat test.m Makefile.standalone README
