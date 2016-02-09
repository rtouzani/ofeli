
bin_PROGRAMS = tiff2
AUTOMAKE_OPTIONS = no-dependencies

tiff2_SOURCES = main.cpp
tiff2_HEADERS = User.h

tiff2_LDADD  = $(top_builddir)/src/libofeli.a

tiff2dir   =	$(datadir)/ofeli/demos/fluid/tiff2
tiff2_DATA =	main.cpp \
		User.h \
		cavity.dat \
		cavity.m \
		Makefile.standalone \
		tiff2.cbp \
		README
dist_tiff2_DATA =	main.cpp \
			User.h \
			cavity.dat \
			cavity.m \
			Makefile.standalone \
			tiff2.cbp \
			README

