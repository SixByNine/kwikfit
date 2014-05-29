if HAVE_CFITSIO
lib_LTLIBRARIES=libpfits.la
libpfits_la_CFLAGS=-fPIC @CFITSIO_CFLAGS@
libpfits_la_LIBADD=@CFITSIO_LIBS@
libpfits_la_SOURCES=pfits.c
endif

if HAVE_PGPLOT
bin_PROGRAMS=pfits_plotfold pfits_fv pfits_plot pfits_plotfold
pfits_fv_LDADD=libpfits.la @PGPLOT_LIBS@
pfits_plotfold_LDADD=libpfits.la @PGPLOT_LIBS@
pfits_plot_LDADD=libpfits.la @PGPLOT_LIBS@
pfits_plotfold_LDADD=libpfits.la @PGPLOT_LIBS@
pfits_bpass_LDADD=libpfits.la @PGPLOT_LIBS@
endif


