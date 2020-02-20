
bin_PROGRAMS = nl_demo1 nl_demo2 nl_demo3 nl_demo4
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

nl_demo1_SOURCES =	nl_demo1.cpp
nl_demo2_SOURCES =	nl_demo2.cpp
nl_demo3_SOURCES =	nl_demo3.cpp
nl_demo4_SOURCES =	nl_demo4.cpp

nl_demo1_LDADD  = $(top_builddir)/src/libofeli.a
nl_demo2_LDADD  = $(top_builddir)/src/libofeli.a
nl_demo3_LDADD  = $(top_builddir)/src/libofeli.a
nl_demo4_LDADD  = $(top_builddir)/src/libofeli.a

nl_demo1dir  =	$(datadir)/ofeli/demos/solvers/NLAS
nl_demo2dir  =	$(datadir)/ofeli/demos/solvers/NLAS
nl_demo3dir  =	$(datadir)/ofeli/demos/solvers/NLAS
nl_demo4dir  =	$(datadir)/ofeli/demos/solvers/NLAS

nl_demo1_DATA =	nl_demo1.cpp
nl_demo2_DATA =	nl_demo2.cpp
nl_demo3_DATA =	nl_demo3.cpp
nl_demo4_DATA =	nl_demo4.cpp

dist_nl_demo1_DATA =	nl_demo1.cpp
dist_nl_demo2_DATA =	nl_demo2.cpp
dist_nl_demo3_DATA =	nl_demo3.cpp
dist_nl_demo4_DATA =	nl_demo4.cpp
