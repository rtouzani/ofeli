
bin_PROGRAMS = helmholtz
AUTOMAKE_OPTIONS = no-dependencies

helmholtz_SOURCES = main.cpp
helmholtz_LDADD  = $(top_builddir)/src/libofeli.a

helmholtzdir   =	$(datadir)/ofeli/demos/electromagnetics/helmholtz
helmholtz_DATA =	main.cpp \
						proj.m \
						proj.dat \
						Makefile.standalone \
						Helmholtz.cbp \
						README
dist_helmholtz_DATA =	main.cpp \
								proj.m \
								proj.dat \
								Makefile.standalone \
								Helmholtz.cbp \
								README

