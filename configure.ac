dnl ------------------------------------------------------------------------------
dnl OFELI - Copyright (C) 1998-2020 Rachid Touzani

dnl This file is part of OFELI.

dnl OFELI is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU Lesser General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.

dnl OFELI is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU Lesser General Public License for more details.

dnl You should have received a copy of the GNU Lesser General Public License
dnl along with OFELI. If not, see <http://www.gnu.org/licenses/>.
dnl ------------------------------------------------------------------------------
dnl Configure script for the OFELI Library
dnl Process this file with autoconf to produce a configure script.
dnl ------------------------------------------------------------------------------

dnl Disable caching.
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl

# Define current version and release
sinclude([VERSION])
AC_SUBST(VERSION)
AC_SUBST(RELEASE)

AC_INIT([ofeli], 4.0.0)
AC_CONFIG_SRCDIR([src/mesh/Mesh.cpp])
AM_INIT_AUTOMAKE

dnl Report what was found
AC_MSG_CHECKING(Version)
AC_MSG_RESULT($VERSION)
AC_MSG_CHECKING(Release)
AC_MSG_RESULT($RELEASE)

AM_CONFIG_HEADER(include/config.h)
AM_MAINTAINER_MODE

AC_CYGWIN()

dnl Prevents setting flags.
CPPFLAGS="-Wall "
CXXFLAGS=" "
CFLAGS="-Wall "

if test "$OFELI_CFLAGS" == "" ; then
   OFELI_CFLAGS="-I$includedir/ofeli "
fi
if test "$OFELI_LIBS" == "" ; then
   OFELI_LIBS="-lofeli "
fi
if test "$OFELI_LDFLAGS" == "" ; then
   OFELI_LDFLAGS="-L$libdir "
fi
if test "$PETSC_CFLAGS" == "" ; then
   PETSC_CFLAGS="-I$includedir/petsc "
fi
if test "$PETSC_LIBS" == "" ; then
   PETSC_LIBS="-lpetsc "
fi
if test "$MPI_CFLAGS" == "" ; then
   MPI_CFLAGS="-I$includedir/mpi "
fi
if test "$MPI_LIBS" == "" ; then
   MPI_LIBS="-lmpi "
fi
if test "$MPI_LDFLAGS" == "" ; then
   MPI_LDFLAGS=" "
fi
if test "$EIGEN_CFLAGS" == "" ; then
   EIGEN_CFLAGS="-I$includedir/eigen3 "
fi

AC_SUBST(OFELI_LIBS)
AC_SUBST(OFELI_CFLAGS)
AC_SUBST(OFELI_LDFLAGS)
AC_SUBST(PETSC_CFLAGS)
AC_SUBST(PETSC_LIBS)
AC_SUBST(PETSC_LDFLAGS)
AC_SUBST(EIGEN_CFLAGS)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_LDFLAGS)

dnl Checks for programs.
AC_LANG_CPLUSPLUS
AC_PROG_MAKE_SET
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL

dnl Get the operating system and machine names
UNAME=`uname`
HOSTNAME=`hostname`


# AC_PROG_LN_S
AC_CHECK_TOOL(RANLIB, ranlib, :)
AC_CHECK_TOOL(AR, ar, :)
AC_PROG_RANLIB

dnl Check compiler features
AC_TYPE_SIZE_T
AC_C_CONST
AC_C_INLINE

AH_TEMPLATE(HAVE_INLINE,0)
AH_TEMPLATE(HIDE_INLINE_STATIC,0)
AH_TEMPLATE(HAVE_CLOCKS_PER_SEC,0)
AH_TEMPLATE(HAVE_EXIT_SUCCESS_AND_FAILURE,0)
#AH_TEMPLATE(HAVE_SPARCLINUX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_M68KLINUX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_PPCLINUX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_X86LINUX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_SUNOS4_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_SOLARIS_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_HPUX11_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_HPUX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_TRU64_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_IRIX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_AIX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_FREEBSD_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_OS2EMX_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_NETBSD_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_OPENBSD_IEEE_INTERFACE,0)
#AH_TEMPLATE(HAVE_DARWIN_IEEE_INTERFACE,0)
AH_TEMPLATE(HAVE_GUESSED_CLOCKS_PER_SEC,0)

dnl Check for "extern inline", using a modified version of the test
dnl for AC_C_INLINE from acspecific.mt
dnl
AC_CACHE_CHECK([for extern inline], ac_cv_c_extern_inline,
[ac_cv_c_extern_inline=no
AC_TRY_COMPILE([extern $ac_cv_c_inline double foo(double x);
extern $ac_cv_c_inline double foo(double x) { return x + 1.0 ; } ;
double foo (double x) { return x + 1.0 ; };],
[  foo(1.0)  ],
[ac_cv_c_extern_inline="yes"])
])

if test "$ac_cv_c_extern_inline" != no ; then
   AC_DEFINE(HAVE_INLINE,1)
   AC_SUBST(HAVE_INLINE)
fi

dnl Checks for header files.
dnl Checks for typedefs, structures, and compiler characteristics.
dnl Checks for library functions.

AC_CHECK_FUNCS(clock)
AC_CACHE_CHECK(for CLOCKS_PER_SEC,
ac_cv_decl_clocks_per_sec,
AC_EGREP_CPP(yes,
[
#include <time.h>
#ifdef CLOCKS_PER_SEC
yes
#endif
],
ac_cv_decl_clocks_per_sec=yes,
ac_cv_decl_clocks_per_sec=no)
)

if test "$ac_cv_decl_clocks_per_sec" = yes ; then
  AC_DEFINE(HAVE_CLOCKS_PER_SEC)
else
AC_MSG_CHECKING(whether CLOCKS_PER_SEC is known for $host_os)
  case "$host" in
   *-*-sunos4.*)
        AC_DEFINE(HAVE_GUESSED_CLOCKS_PER_SEC,1000000)
	AC_MSG_RESULT([yes (microseconds)])
        ;;
    *-*-nextstep*)
	dnl
	dnl Nextstep3 defines CLK_TCK 64 in /usr/include/ansi/time.h,
	dnl which seems to agree with observation of clock().
	dnl
   AC_DEFINE(HAVE_GUESSED_CLOCKS_PER_SEC,64)
	AC_MSG_RESULT([yes (15.625 milliseconds)])
        ;;
    *)
        AC_DEFINE(HAVE_GUESSED_CLOCKS_PER_SEC,1000000)
        AC_MSG_RESULT([no, assuming microseconds])
        ;;
  esac
fi

AC_CACHE_CHECK(for EXIT_SUCCESS and EXIT_FAILURE,
ac_cv_decl_exit_success_and_failure,
AC_EGREP_CPP(yes,
[
#include <stdlib.h>
#ifdef EXIT_SUCCESS
yes
#endif
],
ac_cv_decl_exit_success_and_failure=yes,
ac_cv_decl_exit_success_and_failure=no)
)

if test "$ac_cv_decl_exit_success_and_failure" = yes ; then
  AC_DEFINE(HAVE_EXIT_SUCCESS_AND_FAILURE)
fi ;

dnl Add copious amounts of debugging with gcc, egcs, etc
if test "${GXX}" = "yes" ; then
CPPFLAGS="${CPPFLAGS}"
fi
if test "${GCC}" = "yes" ; then
CFLAGS="${CFLAGS}"
fi

AC_ARG_VAR([OFELI_CFLAGS], [OFELI Compilation flags.
                            Default = "-I${includedir}/ofeli"])
AC_ARG_VAR([OFELI_LIBS], [OFELI Library option.
                          Default = "-lofeli"])
AC_ARG_VAR([OFELI_LDFLAGS], [OFELI Library location.
                             Default = "-L$libdir"])
AC_ARG_VAR([PETSC_CFLAGS], [Petsc Compilation flags.
                            Default = "-I${includedir}/petsc"])
AC_ARG_VAR([PETSC_LIBS], [Petsc Library option.
                          Default = "-lpetsc"])
AC_ARG_VAR([MPI_CFLAGS], [MPI Compilation flags.
                          Default = "-I${includedir}/mpi"])
AC_ARG_VAR([MPI_LIBS], [MPI library option.
                        Default = "-lmpi"])
AC_ARG_VAR([MPI_LDFLAGS], [MPI library location.
                           Default = "-L$libdir"])
AC_ARG_VAR([EIGEN_CFLAGS], [Eigen Compilation flags.
                            Default = "-I${includedir}/eigen3"])

dnl Debugging turned on
AC_MSG_CHECKING(for debugging)
AC_ARG_ENABLE(debug,[  --enable-debug          compile for debugging])
AC_MSG_RESULT([$enable_debug])
if test "x$enable_debug" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -g -DDEBUG"
CFLAGS="${CFLAGS} -g -DDEBUG"
fi

dnl Debugging and Tracing turned on
AC_MSG_CHECKING(for debugging and tracing)
AC_ARG_ENABLE(debug,[  --enable-trace-debug    compile for debugging and tracing])
AC_MSG_RESULT([$enable_trace_debug])
CPPFLAGS="${CPPFLAGS}"
CFLAGS="${CFLAGS}"
if test "x$enable_trace_debug" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -g -DDEBUG -D_OFELI_DEBUG"
CFLAGS="${CFLAGS} -g -DDEBUG -D_OFELI_DEBUG"
fi

dnl Building for release
AC_MSG_CHECKING(for release build)
AC_ARG_ENABLE(release,[  --enable-release        compile for release])
AC_MSG_RESULT([$enable_release])
if test "x$enable_release" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -O3 -DNDEBUG"
CFLAGS="${CFLAGS} -O3 -DNDEBUG"
if test "${GXX}" = "yes" ; then
CPPFLAGS="${CPPFLAGS} -ffast-math -fstrict-aliasing"
CFLAGS="${CFLAGS} -ffast-math -fstrict-aliasing"
fi
fi
 
dnl Building for fast
AC_MSG_CHECKING(for fast build)
AC_ARG_ENABLE(fast,[  --enable-fast           compile for fast configuration])
AC_MSG_RESULT([$enable_fast])
if test "x$enable_fast" = "xyes" ; then
CPPFLAGS="${CPPFLAGS} -Ofast -DNDEBUG"
CFLAGS="${CFLAGS} -Ofast -DNDEBUG"
if test "${GXX}" = "yes" ; then
CPPFLAGS="${CPPFLAGS} -ffast-math -fstrict-aliasing"
CFLAGS="${CFLAGS} -ffast-math -fstrict-aliasing"
fi
fi

dnl Building with PETSc
AC_MSG_CHECKING(for PETSc and dependent libraries)
AC_ARG_WITH(petsc,[  --with-petsc            use PETSc in conjunction with OFELI])
AC_MSG_RESULT([$with_petsc])
if test "x$with_petsc" = "xyes" ; then
   CPPFLAGS="${CPPFLAGS} -DUSE_PETSC ${PETSC_CFLAGS} ${MPI_CFLAGS}"
#   OFELI_LIB="-DUSE_PETSC ${PETSC_CFLAGS} ${MPI_CFLAGS}"
   LIBS="$LIBS -DUSE_PETSC $PETSC_LIBS $MPI_LIBS $PETSC_LDFLAGS $MPI_LDFLAGS"
#   AC_CHECK_LIB(lapack, sgeev, [PETSC_LIBS+="-llapack "], AC_MSG_ERROR(liblapack is required))
   AC_CHECK_LIB(petsc, send, [PETSC_LIBS], [])
fi

dnl Building with Eigen
AC_MSG_CHECKING(for Eigen and dependent libraries)
AC_ARG_WITH(eigen,[  --with-eigen            use eigen library in conjunction with OFELI])
AC_MSG_RESULT([$with_eigen])
if test "x$with_eigen" = "xyes" ; then
   CPPFLAGS="${CPPFLAGS} -DUSE_EIGEN ${EIGEN_CFLAGS}"
fi

# Restore CFLAGS
CFLAGS="$save_cflags "

# Variables that will be substituted in the Makefile.in's

AC_SUBST(PACKAGE)
OFELI_LIB="-lofeli"
inc_dir="$includedir/ofeli"
OFELI_FLAGS="-I\$(top_srcdir)/include"
AC_SUBST(OFELI_LIB)
AC_SUBST(OFELI_FLAGS)

#save_cflags="$CFLAGS"

# Now restore our (possibly augmented) CFLAGS, LIBS and CPPFLAGS.
CFLAGS="$save_cflags"
OFELI_CFLAGS="$OFELI_CFLAGS"
if test "x$with_petsc" = "xyes" ; then
OFELI_LIBS="$OFELI_LIBS $PETSC_LIBS $MPI_LIBS $OFELI_LDFLAGS $PETSC_LDFLAGS $MPI_LDFLAGS"
fi

# List of executables
util_exes="g2m cmesh cfield vmesh "
demo_exes="elliptic heat transport "
demo_exes+="laplace_demo1 laplace_demo2 laplace_demo3 laplace_demo4 "
demo_exes+="std3 stdc2 ttd2 beam truss contact lelas2d lelas3d tiff2 tins2 lh2d euler-2d Helmholtz ls_demo1 ls_demo2 "
demo_exes+="ode_demo1 ode_demo2 ode_demo3 ode_demo4 ts_demo1 ts_demo2 eigen_demo1 eigen_demo2 "
demo_exes+="opt_demo1 opt_demo2 nl_demo1 nl_demo2 nl_demo3 nl_demo4 "

AC_OUTPUT(
   Makefile
   include/Makefile
   include/equations/Makefile
   include/equations/cl/Makefile
   include/equations/electromagnetics/Makefile
   include/equations/fluid/Makefile
   include/equations/interface/Makefile
   include/equations/laplace/Makefile
   include/equations/solid/Makefile
   include/equations/therm/Makefile
   include/equations/porous/Makefile
   include/io/Makefile
   include/io/fparser/Makefile
   include/io/tclap/Makefile
   include/io/xmlsp/Makefile
   include/linear_algebra/Makefile
   include/linear_algebra/petsc/Makefile
   include/mesh/Makefile
   include/mesh/bamg/Makefile
   include/mesh/metis/Makefile
   include/post/Makefile
   include/shape_functions/Makefile
   include/solvers/Makefile
   include/util/Makefile
   src/Makefile
   src/equations/Makefile
   src/equations/cl/Makefile
   src/equations/electromagnetics/Makefile
   src/equations/fluid/Makefile
   src/equations/interface/Makefile
   src/equations/laplace/Makefile
   src/equations/porous/Makefile
   src/equations/solid/Makefile
   src/equations/therm/Makefile
   src/io/Makefile
   src/io/fparser/Makefile
   src/io/xmlsp/Makefile
   src/linear_algebra/Makefile
   src/mesh/Makefile
   src/mesh/bamg/Makefile
   src/mesh/metis/Makefile
   src/post/Makefile
   src/shape_functions/Makefile
   src/solvers/Makefile
   src/util/Makefile
   material/Makefile
   util/Makefile
   util/g2m/Makefile
   util/conv/Makefile
   util/conv/src/Makefile
   util/conv/examples/Makefile
   util/conv/matlab/Makefile
   util/vmesh/Makefile
   demos/Makefile
   demos/1D/Makefile
   demos/laplace/Makefile
   demos/thermal/Makefile
   demos/thermal/std3/Makefile
   demos/thermal/stdc2/Makefile
   demos/thermal/ttd2/Makefile
   demos/solid/Makefile
   demos/solid/beam/Makefile
   demos/solid/truss/Makefile
   demos/solid/contact/Makefile
   demos/solid/lelas2d/Makefile
   demos/solid/lelas3d/Makefile
   demos/fluid/Makefile
   demos/fluid/tiff2/Makefile
   demos/fluid/tins2/Makefile
   demos/fluid/lh2d/Makefile
   demos/fluid/euler-2d/Makefile
   demos/electromagnetics/Makefile
   demos/electromagnetics/Helmholtz/Makefile
   demos/solvers/Makefile
   demos/solvers/LS/Makefile
   demos/solvers/ODE/Makefile
   demos/solvers/EIGEN/Makefile
   demos/solvers/OPTIM/Makefile
   demos/solvers/NLAS/Makefile
   demos/adapt/Makefile
   doc/Makefile
   doc/html/Makefile
   doc/html/search/Makefile
   doc/demos/Makefile
   doc/im/Makefile
   doc/css/Makefile
   doc/presentations/Makefile
   tools/Makefile
   ofeli-config
   test_ofeli.sh
)

echo "********************************************************************"
echo "OFELI is configured for"
echo "  - OS                  : ${UNAME} on ${HOSTNAME}"
echo "  - C++ compiler        : ${CXX}"
echo "  - Exe linker          : ${CXX}"
echo "  - C++ flags           : ${CPPFLAGS}"
echo "  - Utility executables : ${util_exes}"
echo "  - Demo executables    : ${demo_exes}"
echo "********************************************************************"
echo ""
echo "Configure finished!"
echo "            Type: 'make' to compile OFELI, utilities and examples."
echo "            Then: 'make install' to install it."
echo "            Tests are not installed."
echo "            They can be run once installation is completed."
echo "            For this, type: make tests"
echo ""
