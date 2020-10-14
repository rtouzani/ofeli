
bin_PROGRAMS = lp_demo1
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

lp_demo1_SOURCES = lp_demo1.cpp

lp_demo1_LDADD  = $(top_builddir)/src/libofeli.a

lp_demo1dir =	$(datadir)/ofeli/demos/solvers/LP

lp_demo1_DATA = lp_demo1.cpp \
                Makefile.standalone \
                README

dist_lp_demo1_DATA = lp_demo1.cpp \
                     Makefile.standalone \
                     README
