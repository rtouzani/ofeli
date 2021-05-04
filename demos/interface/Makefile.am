
bin_PROGRAMS = fmm_demo
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

fmm_demo_SOURCES = fmm_demo.cpp

fmm_demo_LDADD  = $(top_builddir)/src/libofeli.a

fmm_demodir = $(datadir)/ofeli/demos/interface

fmm_demo_DATA = fmm_demo.cpp README Makefile.standalone

dist_fmm_demo_DATA = fmm_demo.cpp README Makefile.standalone

clean-local:
	-rm -rf u.pos
