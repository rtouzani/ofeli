
bin_PROGRAMS = truss
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

truss_SOURCES = main.cpp

truss_LDADD  = $(top_builddir)/src/libofeli.a

trussdir   = $(datadir)/ofeli/demos/solid/truss

truss_DATA = main.cpp \
             truss.dat \
             truss.m \
             tower1.m \
             tower1.dat \
             tower2.m \
             tower2.dat \
             HoweBridge.m \
             PrattBridge.m \
             Makefile.standalone \
             README

dist_truss_DATA = main.cpp \
                  truss.dat \
                  truss.m \
                  tower1.m \
                  tower1.dat \
                  tower2.m \
                  tower2.dat \
                  HoweBridge.m \
                  PrattBridge.m \
                  Makefile.standalone \
                  README

clean-local:
	-rm -rf *.d truss-1.m
