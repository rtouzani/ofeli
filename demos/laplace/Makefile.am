
bin_PROGRAMS = laplace_demo1 laplace_demo2 laplace_demo3 laplace_demo4
AUTOMAKE_OPTIONS = no-dependencies

laplace_demo1_SOURCES = laplace_demo1.cpp
laplace_demo2_SOURCES = laplace_demo2.cpp
laplace_demo3_SOURCES = laplace_demo3.cpp
laplace_demo4_SOURCES = laplace_demo4.cpp

laplace_demo1_LDADD  = $(top_builddir)/src/libofeli.a
laplace_demo2_LDADD  = $(top_builddir)/src/libofeli.a
laplace_demo3_LDADD  = $(top_builddir)/src/libofeli.a
laplace_demo4_LDADD  = $(top_builddir)/src/libofeli.a

laplace_demo1dir =	$(datadir)/ofeli/demos/laplace
laplace_demo2dir =	$(datadir)/ofeli/demos/laplace
laplace_demo3dir =	$(datadir)/ofeli/demos/laplace
laplace_demo4dir =	$(datadir)/ofeli/demos/laplace

laplace_demo1_DATA = laplace_demo1.cpp \
                     proj2.dat \
                     proj3.dat \
                     proj.dat \
                     square.geo \
                     square.m \
                     cube.geo \
                     cube.m \
                     annulus.geo \
                     annulus.m \
                     Makefile.standalone \
                     README
laplace_demo2_DATA =	laplace_demo2.cpp
laplace_demo3_DATA =	laplace_demo3.cpp
laplace_demo4_DATA =	laplace_demo4.cpp

dist_laplace_demo1_DATA = laplace_demo1.cpp \
                          proj2.dat \
                          proj3.dat \
                          proj.dat \
                          square.geo \
                          square.m \
                          cube.geo \
                          cube.m \
                          annulus.geo \
                          annulus.m \
                          Makefile.standalone \
                          README
dist_laplace_demo2_DATA =	laplace_demo2.cpp
dist_laplace_demo3_DATA =	laplace_demo3.cpp
dist_laplace_demo4_DATA =	laplace_demo4.cpp
