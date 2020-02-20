AUTOMAKE_OPTIONS = no-dependencies
g2mdir   =	$(datadir)/ofeli/util/g2m
AM_CPPFLAGS = -std=c++1y

bin_PROGRAMS = g2m

g2m_SOURCES = main.cpp
g2m_LDADD = $(top_builddir)/src/libofeli.a

g2m_DATA = main.cpp test.dom README
dist_g2m_DATA = main.cpp test.dom README

