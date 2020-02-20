
bin_PROGRAMS = beam
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

beam_SOURCES = main.cpp

beam_LDADD  = $(top_builddir)/src/libofeli.a

beamdir = $(datadir)/ofeli/demos/solid/beam

beam_DATA = main.cpp \
            beam.dat \
            Makefile.standalone \
            README

dist_beam_DATA = main.cpp \
                 beam.dat \
                 Makefile.standalone \
	         README

clean-local:
	-rm -rf beam_tecplot.dat deformed_beam.m
