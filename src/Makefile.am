bin_PROGRAMS=kwikfit-ascii
lib_LTLIBRARIES=libkwikfit.la
include_HEADERS = kwikfit.h

kwikfit_ascii_LDFLAGS=libkwikfit.la
kwikfit_ascii_SOURCES=kwikfit-ascii.cpp mjk_cmd.cpp

headers = kwikfit.h mjklog.h mjk_cmd.h
libkwikfit_la_SOURCES = $(headers) kwikfit.cpp vonMises.cpp

libkwikfit_la_CFLAGS = -fPIC @TEMPO2_CFLAGS@ 
libkwikfit_la_LIBADD =@TEMPO2_LIBS@ 

if HAVE_PGPLOT
libkwikfit_la_SOURCES += kwikfit_pgplot.cpp
libkwikfit_la_CFLAGS += @PGPLOT_CFLAGS@ 
libkwikfit_la_LIBADD +=@PGPLOT_LIBS@ 
endif

if HAVE_LAPACK
libkwikfit_la_LIBADD+=@LIBS_LAPACK@
endif

if HAVE_BLAS
libkwikfit_la_LIBADD+=@BLAS_LIBS@
endif



