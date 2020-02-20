
bin_PROGRAMS = std3
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

std3_SOURCES = main.cpp

std3_LDADD  = $(top_builddir)/src/libofeli.a

std3dir   =	$(datadir)/ofeli/demos/thermal/std3

std3_DATA =	main.cpp \
            beam.dat \
            beam.m \
            Makefile.standalone \
            README

dist_std3_DATA = main.cpp \
                 beam.dat \
                 beam.m \
                 Makefile.standalone \
                 README

clean-local:
	-rm -rf beam.pos
