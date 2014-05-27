bin_PROGRAMS=kwikfit-ascii
lib_LTLIBRARIES=libkwikfit.la
include_HEADERS = kwikfit.h


headers = kwikfit.h
libkwikfit_la_SOURCES = $(headers) kwikfit.c

libkwikfit_la_CFLAGS = -fPIC @TEMPO2_CFLAGS@
libkwikfit_la_LIBADD =@TEMPO2_LIBS@

if HAVE_LAPACK
libkwikfit_la_LIBADD+=@LIBS_LAPACK@
endif

if HAVE_BLAS
libkwikfit_la_LIBADD+=@BLAS_LIBS@
endif



