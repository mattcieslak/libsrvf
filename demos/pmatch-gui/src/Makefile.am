bin_PROGRAMS=srvf_pmatch_gui

noinst_HEADERS=matchview.h ui.h

AM_CPPFLAGS=$(FLTK_CFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS=$(FLTK_LDFLAGS) $(OPENMP_CXXFLAGS)

srvf_pmatch_gui_SOURCES=main.cc ui.cc
srvf_pmatch_gui_LDADD=$(top_srcdir)/src/libsrvf.la $(top_srcdir)/plot/libsrvf-plot.la
srvf_pmatch_gui_LDFLAGS=-static
