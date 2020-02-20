
bin_PROGRAMS = ls_demo1 ls_demo2
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

ls_demo1_SOURCES = ls_demo1.cpp
ls_demo2_SOURCES = ls_demo2.cpp

ls_demo1_LDADD = $(top_builddir)/src/libofeli.a
ls_demo2_LDADD = $(top_builddir)/src/libofeli.a

ls_demo1dir = $(datadir)/ofeli/demos/solvers/LS
ls_demo2dir = $(datadir)/ofeli/demos/solvers/LS

ls_demo1_DATA = ls_demo1.cpp
ls_demo2_DATA = ls_demo2.cpp \
                mesh.m \
                Makefile.standalone \
                README

dist_ls_demo1_DATA = ls_demo1.cpp
dist_ls_demo2_DATA = ls_demo2.cpp \
                     mesh.m \
                     Makefile.standalone \
                     README

clean-local:
	-rm -rf *.pos
