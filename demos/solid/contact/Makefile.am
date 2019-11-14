
bin_PROGRAMS = contact
AUTOMAKE_OPTIONS = no-dependencies

contact_SOURCES = main.cpp

contact_LDADD  = $(top_builddir)/src/libofeli.a

contactdir   =	$(datadir)/ofeli/demos/solid/contact
contact_DATA = main.cpp \
               beam.dat \
               hertz-40x15.m \
               hertz.dat \
               MyMaterial.md \
               Makefile.standalone \
               README
dist_contact_DATA = main.cpp \
                    beam.dat \
                    hertz-40x15.m \
                    hertz.dat \
                    MyMaterial.md \
                    Makefile.standalone \
                    README
