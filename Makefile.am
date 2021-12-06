## Process this file with automake to produce Makefile.in

DDIR = $(datadir)/ofeli/doc

AUX_DIST = config.guess \
           config.sub \
           install-sh \
           ltconfig \
           ltmain.sh \
           mdate-sh \
           missing \
           mkinstalldirs

AUX_DIST_EXTRA = readline.m4 \
                 sys_errlist.m4 \
                 sys_siglist.m4

EXTRA_DIST = bootstrap ChangeLog install-sdk.sh README-sdk
AUTOMAKE_OPTIONS = foreign
SUBDIRS = doc include src util material demos
bin_SCRIPTS = ofeli-config test_ofeli.sh install-sdk.sh

MAINTAINERCLEANFILES = aclocal.m4 configure $(AUX_DIST)

DISTCLEANFILES = $(top_builddir)/include/datadir.h

ACLOCAL = aclocal

tests:
	chmod a+x test_ofeli.sh
	./test_ofeli.sh

clean-local:
	-rm -rf autom4te.cache *.log include/stamp-h1
