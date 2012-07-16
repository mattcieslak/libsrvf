# Paths for SRVF libs and headers in the build tree
SRVF_CPPFLAGS="-I$(top_srcdir)/../srvf/src"
SRVF_LDFLAGS="-L$(top_srcdir)/../srvf/src/.libs"
SRVF_LIBS=" -lsrvf"

AM_CPPFLAGS=$(SRVF_CPPFLAGS) $(FLTK_CFLAGS)
AM_LDFLAGS=$(SRVF_LDFLAGS) $(FLTK_LDFLAGS)
LIBS+=$(SRVF_LIBS)

lib_LTLIBRARIES=libsrvf-plot.la

#libsrvf_plot_la_LDFLAGS= -static

srvf_includedir=$(includedir)/srvf
srvf_include_HEADERS=plot.h plothelper.h plotwin.h render.h

libsrvf_plot_la_SOURCES=plot.cc plothelper.cc plotwin.cc render.cc 
