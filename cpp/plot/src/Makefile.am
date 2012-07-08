AM_CPPFLAGS=-I$(top_srcdir)/../srvf/src $(FLTK_CFLAGS)
AM_LDFLAGS=-L$(top_srcdir)/../srvf/src $(FLTK_LDFLAGS)
LIBS+=-lsrvf

lib_LTLIBRARIES=libsrvf-plot.la

libsrvf_plot_la_LDFLAGS= -static

srvf_includedir=$(includedir)/srvf

srvf_include_HEADERS=plot.h plothelper.h plotwin.h render.h

libsrvf_plot_la_SOURCES=plot.cc plothelper.cc plotwin.cc render.cc 
