check_PROGRAMS=runtests
AM_CPPFLAGS=-I$(top_srcdir)/src -I$(top_srcdir)/../plot/src $(FLTK_CFLAGS)
AM_LDFLAGS=-L$(top_srcdir)/src/.libs -L$(top_srcdir)/../plot/src/.libs $(FLTK_LDFLAGS)
LIBS+=-lsrvf -lsrvf-plot -lboost_test_exec_monitor -lgomp
runtests_SOURCES=fileio_tests.cc function_tests.cc interp_tests.cc kdtree_tests.cc main.cc matrix_tests.cc minmaxheap_tests.cc numeric_tests.cc opencurves_tests.cc plf_tests.cc plot_tests.cc plotwin_tests.cc point_tests.cc pointset_tests.cc qmap_tests.cc render_tests.cc reparam_tests.cc rotate_tests.cc srvf_tests.cc util_tests.cc
