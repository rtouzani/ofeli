
bin_PROGRAMS = lelas3d
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

lelas3d_SOURCES = main.cpp

lelas3d_LDADD  = $(top_builddir)/src/libofeli.a

lelas3ddir   =	$(datadir)/ofeli/demos/solid/lelas3d

lelas3d_DATA = main.cpp \
               beam.dat \
               beam.m \
               Makefile.standalone \
               README

dist_lelas3d_DATA = main.cpp \
                    beam.dat \
                    beam.m \
                    Makefile.standalone \
                    README

clean-local:
	-rm -rf *.d beam-1.m

