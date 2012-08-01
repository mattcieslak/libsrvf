check_PROGRAMS=runtests

AM_CPPFLAGS=-I$(top_srcdir)/srvf/src -I$(top_srcdir)/plot/src $(FLTK_CFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS=$(FLTK_LDFLAGS) $(OPENMP_CXXFLAGS)

LIBS+=-lboost_test_exec_monitor

runtests_SOURCES=fileio_tests.cc function_tests.cc interp_tests.cc kdtree_tests.cc main.cc matrix_tests.cc minmaxheap_tests.cc numeric_tests.cc opencurves_tests.cc plf_tests.cc plot_tests.cc plotwin_tests.cc point_tests.cc pointset_tests.cc qmap_tests.cc render_tests.cc reparam_tests.cc rotate_tests.cc srvf_tests.cc util_tests.cc
runtests_LDADD=$(top_srcdir)/srvf/src/libsrvf.la $(top_srcdir)/plot/src/libsrvf-plot.la

EXTRA_DIST=data
