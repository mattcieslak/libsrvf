bin_PROGRAMS=srvf_pmatch_gui

AM_CPPFLAGS=-I$(top_srcdir)/srvf/src -I$(top_srcdir)/plot/src $(FLTK_CFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS=$(FLTK_LDFLAGS) $(OPENMP_CXXFLAGS)

srvf_pmatch_gui_SOURCES=main.cc ui.cc
srvf_pmatch_gui_LDADD=$(top_srcdir)/srvf/src/libsrvf.la $(top_srcdir)/plot/src/libsrvf-plot.la
srvf_pmatch_gui_LDFLAGS=-static
