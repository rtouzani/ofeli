
bin_PROGRAMS = eigen_demo1 eigen_demo2
AUTOMAKE_OPTIONS = no-dependencies

eigen_demo1_SOURCES =	eigen_demo1.cpp
eigen_demo2_SOURCES =	eigen_demo2.cpp

eigen_demo1_LDADD  = $(top_builddir)/src/libofeli.a
eigen_demo2_LDADD  = $(top_builddir)/src/libofeli.a

eigen_demo1dir   =	$(datadir)/ofeli/demos/solvers/EIGEN
eigen_demo2dir   =	$(datadir)/ofeli/demos/solvers/EIGEN

eigen_demo1_DATA =	eigen_demo1.cpp eigen_demo1.cbp
eigen_demo2_DATA =	eigen_demo2.cpp eigen_demo2.cbp \
			test.m \
			eigen_demo2.dat \
			README

dist_eigen_demo1_DATA =	eigen_demo1.cpp eigen_demo1.cbp
dist_eigen_demo2_DATA =	eigen_demo2.cpp eigen_demo2.cbp \
			test.m \
			eigen_demo2.dat \
			README
