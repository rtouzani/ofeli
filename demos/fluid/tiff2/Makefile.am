
bin_PROGRAMS = tiff2
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

tiff2_SOURCES = main.cpp

tiff2_LDADD  = $(top_builddir)/src/libofeli.a

tiff2dir   = $(datadir)/ofeli/demos/fluid/tiff2

tiff2_DATA = main.cpp \
             cavity.dat \
             cavity.m \
             Makefile.standalone \
             README

dist_tiff2_DATA = main.cpp \
                  cavity.dat \
                  cavity.m \
                  Makefile.standalone \
                  README

clean-local:
	-rm -rf cavity.v cavity.p

