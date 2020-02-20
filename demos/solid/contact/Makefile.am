
bin_PROGRAMS = contact
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

contact_SOURCES = main.cpp

contact_LDADD  = $(top_builddir)/src/libofeli.a

contactdir   =	$(datadir)/ofeli/demos/solid/contact

contact_DATA = main.cpp \
               beam.dat \
               beam.m \
               MyMaterial.md \
               Makefile.standalone \
               README

dist_contact_DATA = main.cpp \
                    beam.dat \
                    beam.m \
                    MyMaterial.md \
                    Makefile.standalone \
                    README

clean-local:
	-rm -rf *-1.m *.pos
