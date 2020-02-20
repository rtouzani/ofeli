
bin_PROGRAMS = stdc2
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

stdc2_SOURCES = main.cpp

stdc2_LDADD  = $(top_builddir)/src/libofeli.a

stdc2dir   = $(datadir)/ofeli/demos/thermal/stdc2
stdc2_DATA = main.cpp \
             proj.dat \
             proj.m \
             proj.v \
             Makefile.standalone \
             README

dist_stdc2_DATA =	main.cpp \
                  proj.dat \
                  proj.m \
                  proj.v \
                  Makefile.standalone \
                  README

clean-local:
	-rm -rf proj.pos 

