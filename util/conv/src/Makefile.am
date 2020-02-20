
bin_PROGRAMS = cmesh cfield
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

cmesh_SOURCES = cmesh.cpp
cfield_SOURCES = cfield.cpp

cmesh_LDADD  = $(top_builddir)/src/libofeli.a
cfield_LDADD = $(top_builddir)/src/libofeli.a

convdir   =	$(datadir)/ofeli/util/conv/src
conv_DATA =	cmesh.cpp cfield.cpp
