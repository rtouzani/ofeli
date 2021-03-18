## Process this file with automake to produce Makefile.in

DDIR = $(datadir)/ofeli/doc

AUX_DIST = config.guess \
           config.sub \
           install-sh \
           ltconfig \
           ltmain.sh \
           mdate-sh \
           missing \
           mkinstalldirs \
           ofeli.cbp

AUX_DIST_EXTRA = readline.m4 \
                 sys_errlist.m4 \
                 sys_siglist.m4

EXTRA_DIST = bootstrap ChangeLog
AUTOMAKE_OPTIONS = foreign
SUBDIRS = doc include src util material demos
bin_SCRIPTS = ofeli-config test_ofeli.sh

MAINTAINERCLEANFILES = aclocal.m4 configure $(AUX_DIST)

DISTCLEANFILES = $(top_builddir)/include/datadir.h

ACLOCAL = aclocal

tests:
	chmod a+x test_ofeli.sh
	./test_ofeli.sh

sdk:
	mkdir -p ofeli-sdk
	mkdir -p ofeli-sdk/bin
	mkdir -p ofeli-sdk/lib
	mkdir -p ofeli-sdk/include
	mkdir -p ofeli-sdk/material
	cp src/libofeli.a ofeli-sdk/lib/.
	cp util/conv/src/cmesh ofeli-sdk/bin/.
	cp util/conv/src/cfield ofeli-sdk/bin/.
	cp util/g2m/g2m ofeli-sdk/bin/.
	cp util/vmesh/vmesh ofeli-sdk/bin/.
	cp material/* ofeli-sdk/material/.
	cp ofeli-config ofeli-sdk/bin/.
	cp -rf include/ ofeli-sdk/include/
	cp ../tools/install-sdk.sh ofeli-sdk/.
	find ofeli-sdk/ -name "Makefile.am" -exec rm {} \;
	find ofeli-sdk/ -name "Makefile.in" -exec rm {} \;
	find ofeli-sdk/ -name "Makefile" -exec rm {} \;
	tar czf ofeli-sdk.tar.gz ofeli-sdk/
	rm -rf ofeli-sdk

clean-local:
	-rm -rf autom4te.cache *.log include/stamp-h1
