
bin_PROGRAMS = helmholtz
AUTOMAKE_OPTIONS = no-dependencies

helmholtz_SOURCES = main.cpp error.cpp User.cpp
helmholtz_HEADERS = User.h

helmholtz_LDADD  = $(top_builddir)/src/libofeli.a

helmholtzdir   =	$(datadir)/ofeli/demos/electromagnetics/helmholtz
helmholtz_DATA =	main.cpp \
			error.cpp \
			User.cpp \
			User.h \
			proj-10.m \
			proj-20.m \
			proj-40.m \
			proj-80.m \
			proj.dat \
			Makefile.standalone \
			Helmholtz.cbp \
			README
dist_helmholtz_DATA =	main.cpp \
			error.cpp \
			User.cpp \
			User.h \
			proj-10.m \
			proj-20.m \
			proj-40.m \
			proj-80.m \
			proj.dat \
			Makefile.standalone \
			Helmholtz.cbp \
			README

