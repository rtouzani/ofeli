
bin_PROGRAMS = tins2
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

tins2_SOURCES = main.cpp

tins2_LDADD = $(top_builddir)/src/libofeli.a

tins2dir = $(datadir)/ofeli/demos/fluid/tins2

tins2_DATA = main.cpp \
             cavity.dat \
             cavity-60.dat \
             cavity-30.m \
             cavity-60.m \
             step.dom \
             step.dat \
             Makefile.standalone \
             README

dist_tins2_DATA =	main.cpp \
                  cavity.dat \
                  cavity-60.dat \
                  cavity-30.m \
                  cavity-60.m \
                  step.dom \
                  step.dat \
                  Makefile.standalone \
                  README

clean-local:
	-rm -rf *.v *.p

