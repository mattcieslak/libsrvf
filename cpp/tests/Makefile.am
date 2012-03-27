bin_PROGRAMS=test
AM_CPPFLAGS=-I $(top_srcdir)/src
AM_LDFLAGS=-L $(top_srcdir)/src
LIBS+=-lsrvf -lboost_test_exec_monitor
test_SOURCES=main.cc interp_tests.cc matrix_tests.cc util_tests.cc
test_SOURCES+=plf_tests.cc srvf_tests.cc qmap_tests.cc pointset_tests.cc
test_SOURCES+=fileio_tests.cc
