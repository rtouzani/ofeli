
bin_PROGRAMS = helmholtz
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

helmholtz_SOURCES = main.cpp
helmholtz_LDADD  = $(top_builddir)/src/libofeli.a

helmholtzdir   = $(datadir)/ofeli/demos/electromagnetics/helmholtz
helmholtz_DATA = main.cpp \
                 proj.m \
                 proj.dat \
                 Makefile.standalone \
                 README
dist_helmholtz_DATA = main.cpp \
                      proj.m \
                      proj.dat \
                      Makefile.standalone \
                      README

