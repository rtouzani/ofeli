
bin_PROGRAMS = ode_demo1 ode_demo2 ode_demo3 ode_demo4 ts_demo1 ts_demo2
AUTOMAKE_OPTIONS = no-dependencies

ode_demo1_SOURCES =	ode_demo1.cpp
ode_demo2_SOURCES =	ode_demo2.cpp
ode_demo3_SOURCES =	ode_demo3.cpp
ode_demo4_SOURCES =	ode_demo4.cpp
ts_demo1_SOURCES =	ts_demo1.cpp
ts_demo2_SOURCES =	ts_demo2.cpp

ode_demo1_LDADD  = $(top_builddir)/src/libofeli.a
ode_demo2_LDADD  = $(top_builddir)/src/libofeli.a
ode_demo3_LDADD  = $(top_builddir)/src/libofeli.a
ode_demo4_LDADD  = $(top_builddir)/src/libofeli.a
ts_demo1_LDADD  = $(top_builddir)/src/libofeli.a
ts_demo2_LDADD  = $(top_builddir)/src/libofeli.a

ode_demo1dir   =	$(datadir)/ofeli/demos/solvers/ODE
ode_demo2dir   =	$(datadir)/ofeli/demos/solvers/ODE
ode_demo3dir   =	$(datadir)/ofeli/demos/solvers/ODE
ode_demo4dir   =	$(datadir)/ofeli/demos/solvers/ODE
ts_demo1dir   =	$(datadir)/ofeli/demos/solvers/ODE
ts_demo2dir   =	$(datadir)/ofeli/demos/solvers/ODE

ode_demo1_DATA =	ode_demo1.cpp ode_demo1.cbp
ode_demo2_DATA =	ode_demo2.cpp ode_demo2.cbp
ode_demo3_DATA =	ode_demo3.cpp ode_demo3.cbp
ode_demo4_DATA =	ode_demo4.cpp ode_demo4.cbp
ts_demo1_DATA =	ts_demo1.cpp \
		ts_demo1.cbp \
		mesh1.m \
		Makefile.standalone \
		README
ts_demo2_DATA =	ts_demo2.cpp \
		ts_demo2.cbp \
		mesh2.m

dist_ode_demo1_DATA =	ode_demo1.cpp ode_demo1.cbp
dist_ode_demo2_DATA =	ode_demo2.cpp ode_demo2.cbp
dist_ode_demo3_DATA =	ode_demo3.cpp ode_demo3.cbp
dist_ode_demo4_DATA =	ode_demo4.cpp ode_demo4.cbp
dist_ts_demo1_DATA =	ts_demo1.cpp \
			ts_demo1.cbp \
			mesh1.m \
			Makefile.standalone \
			README
dist_ts_demo2_DATA =	ts_demo2.cpp \
			ts_demo2.cbp \
			mesh2.m
