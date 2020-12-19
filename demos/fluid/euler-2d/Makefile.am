
bin_PROGRAMS = euler2d
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

euler2d_SOURCES = main.cpp
euler2d_HEADERS = 

euler2d_LDADD  = $(top_builddir)/src/libofeli.a

euler2ddir   =	$(datadir)/ofeli/demos/cl/euler-2d

euler2d_DATA = main.cpp \
               proj.dat \
               rectangle.m \
               Makefile.standalone \
               README

dist_euler2d_DATA = main.cpp \
                    proj.dat \
                    rectangle.m \
                    Makefile.standalone \
                    README

clean-local:
	-rm -rf rect.c rect.p rect.r rect.v

