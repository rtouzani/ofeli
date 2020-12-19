
bin_PROGRAMS = ad1 ad2
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

ad1_SOURCES = adapt_demo1.cpp
ad2_SOURCES = adapt_demo2.cpp

ad1_LDADD  = $(top_builddir)/src/libofeli.a
ad2_LDADD  = $(top_builddir)/src/libofeli.a

ad1dir   =	$(datadir)/ofeli/demos/adapt
ad2dir   =	$(datadir)/ofeli/demos/adapt

ad1_DATA = adapt_demo1.cpp \
           rect.dom \
           Makefile.standalone \
           README

dist_ad1_DATA = adapt_demo1.cpp \
                rect.dom \
                Makefile.standalone \
                README

ad2_DATA = adapt_demo2.cpp \
           L.dom

dist_ad2_DATA = adapt_demo2.cpp \
                L.dom

clean-local:
	-rm -rf u.pos v.pos mesh1.m mesh2.m test.m 
