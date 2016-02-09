## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = no-dependencies

SUBDIRS = equations shape_functions io linear_algebra mesh post solvers util

#includedir = $(pkgincludedir)

pkginclude_HEADERS =	OFELI.h \
			OFELI_Config.h \
			config.h.in \
			config.h \
			datadir.h \
			Laplace.h \
			Therm.h \
			Solid.h \
			Fluid.h \
			Electromagnetics.h \
			CL.h

datadir.h: Makefile
	echo '#define PATH_MATERIAL "$(datadir)/ofeli/material"' >$@
