
noinst_LIBRARIES = libequations.a
AUTOMAKE_OPTIONS = no-dependencies
SUBDIRS = cl electromagnetics fluid interface laplace porous solid therm
AM_CPPFLAGS = -std=c++1y

libequations_a_SOURCES = Equa.cpp \
                         DG.cpp

libequations_a_LIBADD =	\
			cl/libcl.a \
			electromagnetics/libelectromagnetics.a \
			fluid/libfluid.a \
			interface/libinterface.a \
			laplace/liblaplace.a \
			porous/libporous.a \
			solid/libsolid.a \
			therm/libtherm.a
