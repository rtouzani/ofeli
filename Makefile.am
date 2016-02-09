## Process this file with automake to produce Makefile.in

DDIR = $(datadir)/ofeli/doc \
       $(datadir)/ofeli/tools

AUX_DIST =	config.guess \
		config.sub \
		install-sh \
		ltconfig \
		ltmain.sh \
		mdate-sh \
		missing \
		mkinstalldirs \
		ofeli.cbp

AUX_DIST_EXTRA=	readline.m4 \
		sys_errlist.m4 \
		sys_siglist.m4

EXTRA_DIST = bootstrap ChangeLog
AUTOMAKE_OPTIONS = foreign
SUBDIRS = doc include src util tools material tutorial demos
bin_SCRIPTS = ofeli-config test_ofeli.sh

MAINTAINERCLEANFILES = aclocal.m4 configure $(AUX_DIST)

DISTCLEANFILES = $(top_builddir)/include/datadir.h

ACLOCAL = aclocal

tests:
	chmod a+x test_ofeli.sh
	./test_ofeli.sh
