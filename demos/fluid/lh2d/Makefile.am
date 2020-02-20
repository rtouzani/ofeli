
bin_PROGRAMS = lh2d
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

lh2d_SOURCES = main.cpp
lh2d_HEADERS = 

lh2d_LDADD  = $(top_builddir)/src/libofeli.a

lh2ddir   =	$(datadir)/ofeli/demos/cl/lh2d

lh2d_DATA = main.cpp \
            proj.dat \
            rectangle.m \
            rectangle.v \
            Makefile.standalone \
            README

dist_lh2d_DATA = main.cpp \
                 proj.dat \
                 rectangle.m \
                 rectangle.v \
                 Makefile.standalone \
                 README

clean-local:
	-rm -rf rectangle.t

