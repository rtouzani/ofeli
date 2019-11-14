
bin_PROGRAMS = elliptic heat transport
AUTOMAKE_OPTIONS = no-dependencies

elliptic_SOURCES = elliptic.cpp
heat_SOURCES = heat.cpp
transport_SOURCES = transport.cpp

elliptic_LDADD  = $(top_builddir)/src/libofeli.a
heat_LDADD  = $(top_builddir)/src/libofeli.a
transport_LDADD  = $(top_builddir)/src/libofeli.a

ellipticdir = $(datadir)/ofeli/demos/laplace
heatdir = $(datadir)/ofeli/demos/laplace
transportdir = $(datadir)/ofeli/demos/laplace

elliptic_DATA = elliptic.cpp README
heat_DATA = heat.cpp
transport_DATA = transport.cpp

dist_elliptic_DATA = elliptic.cpp README
dist_heat_DATA = heat.cpp
dist_transport_DATA = transport.cpp
