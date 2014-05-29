#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([kwikfit],[2014.4.5],[mkeith@pulsarastronomy.net])

if test -n "${CXXFLAGS}"; then
	user_set_cxxflags=yes
fi

AC_CONFIG_MACRO_DIR([autoconf])
AC_CONFIG_AUX_DIR([config])

# Enable large-file support
AC_SYS_LARGEFILE

# Determine the machine endian
AC_C_BIGENDIAN

# Checks for compilers
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_F77
AC_PROG_F77_C_O
AC_F77_LIBRARY_LDFLAGS
AC_PROG_FC

if test -n "${CXXFLAGS}"; then
	user_set_cxxflags=yes
fi
AC_PROG_CXX
if test X$user_set_cxxflags != Xyes; then
	CXXFLAGS="-O2"
fi


AC_PROG_INSTALL
AC_PROG_LIBTOOL


if test $can_build_shared == no ; then
	AC_WARN("Libtool cannot build shared libraries - ")
fi

# Check whether the 'cos' function can be found in library 'm'
# (standard C math library).
# Math library is required for the tempo2 predictor library
AC_CHECK_LIB(m, cos)

# look for third-party packages in various places (or not)
SWIN_LOCAL


SWIN_LIB_TEMPO2

if test x"$have_tempo2" != xyes; then
   AC_ERROR(The TEMPO2 library was not found. Cannot build.)
fi

SWIN_LIB_PGPLOT
SWIN_LIB_FFTW
SWIN_LIB_CFITSIO

# No longer need openmp.
AX_OPENMP([have_openmp=yes],[have_openmp=no])
if test X$have_openmp == Xyes; then
	CFLAGS="$CFLAGS $OPENMP_CFLAGS -DHAVE_OPENMP"
	CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS -DHAVE_OPENMP"
fi

#GSL basic stuff
AX_PATH_GSL(1.11.0, have_gsl=yes,have_gsl=no)
AM_CONDITIONAL(HAVE_GSL, [test x"$have_gsl" = xyes])

AC_F77_WRAPPERS

ACX_BLAS(have_blas=yes,have_blas=no)
AC_MSG_CHECKING([for BLAS])
AC_MSG_RESULT($have_blas)
AM_CONDITIONAL(HAVE_BLAS, [test x$have_blas = xyes])

# test for LAPACK acceleration
ACX_LAPACK(have_lapack=yes, have_lapack=no)
AC_MSG_CHECKING([for LAPACK])
AC_MSG_RESULT($have_lapack)
AM_CONDITIONAL(HAVE_LAPACK, [test x"$have_lapack" = xyes])


AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile src/pfits/Makefile ])

AC_OUTPUT

AC_MSG_NOTICE(Fortran compiler is $F77 $FFLAGS)
AC_MSG_NOTICE(C compiler is $CC $CFLAGS)
AC_MSG_NOTICE(C++ compiler is $CXX $CXXFLAGS)
AC_MSG_NOTICE(Global linking flags are $LDFLAGS)


# User warnings
if test X$have_pgplot != Xyes; then
	AC_MSG_NOTICE(Plugins: PGPLOT was not detected - graphical utilities will not be compiled)
else
	AC_MSG_NOTICE(PGPLOT compile flags are "$PGPLOT_CFLAGS")
	AC_MSG_NOTICE(PGPLOT linking flags are "$PGPLOT_LIBS")
fi

if test X$have_fftw3 != Xyes; then
	AC_MSG_NOTICE(Plugins: FFTW3 was not found - Thats ok)
fi

if test X$have_cfitsio != Xyes; then
   AC_MSG_NOTICE(Plugins: cfitsio was not found - Wont be able to read fits)
fi

if test X$have_lapack = Xyes; then
	AC_MSG_NOTICE(Using LAPACK library for Cholesky decomposition)
fi

if test X$have_blas = Xyes; then
	AC_MSG_NOTICE(Using BLAS library for matrix multiplication)
fi

echo ""
echo "****"
echo "Ok! Good to start compiling... try make && make install'"
echo "****"



