bin_PROGRAMS=test
AM_CPPFLAGS=-I$(top_srcdir)/../srvf/src -I../src
AM_LDFLAGS=-L$(top_srcdir)/../srvf/src/.libs -L../src/.libs
LIBS+=-lsrvf-fa -lsrvf -lboost_test_exec_monitor -lgomp
test_SOURCES=function_tests.cc
