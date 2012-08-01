# Paths for SRVF libs and headers in the build tree
SRVF_CPPFLAGS="-I$(top_srcdir)/srvf/src"

AM_CPPFLAGS=$(SRVF_CPPFLAGS) $(FLTK_CFLAGS)
AM_LDFLAGS=$(FLTK_LDFLAGS)

lib_LTLIBRARIES=libsrvf-plot.la

#libsrvf_plot_la_LDFLAGS= -static

srvf_includedir=$(includedir)/srvf
srvf_include_HEADERS=plot.h plothelper.h plotwin.h render.h

libsrvf_plot_la_SOURCES=plot.cc plothelper.cc plotwin.cc render.cc 
libsrvf_plot_la_LIBADD=$(top_srcdir)/srvf/src/libsrvf.la
