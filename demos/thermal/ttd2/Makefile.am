
bin_PROGRAMS = ttd2
AUTOMAKE_OPTIONS = no-dependencies

ttd2_SOURCES = main.cpp error.cpp

ttd2_LDADD  = $(top_builddir)/src/libofeli.a

ttd2dir   =	$(datadir)/ofeli/demos/thermal/ttd2
ttd2_DATA = main.cpp \
            error.cpp \
            proj.dat \
            proj-5x10.m \
            proj-5x20.m \
            proj-5x40.m \
            proj-5x80.m \
            Makefile.standalone \
            README
dist_ttd2_DATA = main.cpp \
                 error.cpp \
                 proj.dat \
                 proj-5x10.m \
                 proj-5x20.m \
                 proj-5x40.m \
                 proj-5x80.m \
                 Makefile.standalone \
                 README
