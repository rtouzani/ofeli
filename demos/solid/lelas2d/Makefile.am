
bin_PROGRAMS = lelas2d
AUTOMAKE_OPTIONS = no-dependencies

lelas2d_SOURCES = main.cpp

lelas2d_LDADD  = $(top_builddir)/src/libofeli.a

lelas2ddir   =	$(datadir)/ofeli/demos/solid/lelas2d
lelas2d_DATA = main.cpp \
               beam.dat \
               beam-10x5.m \
               beam-20x10.m \
               beam.m \
               Makefile.standalone \
               README
dist_lelas2d_DATA = main.cpp \
                    beam.dat \
                    beam-10x5.m \
                    beam-20x10.m \
                    beam.m \
                    Makefile.standalone \
                    README
