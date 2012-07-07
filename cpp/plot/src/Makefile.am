AM_CPPFLAGS=-I$(top_srcdir)/../src
AM_LDFLAGS=-L$(top_srcdir)/../src
LIBS+=-lsrvf

lib_LTLIBRARIES=libsrvf-plot.la

srvf_includedir=$(includedir)/srvf

srvf_include_HEADERS=plot.h plothelper.h plotwin.h render.h

libsrvf_plot_la_SOURCES=plot.cc plothelper.cc plotwin.cc render.cc
