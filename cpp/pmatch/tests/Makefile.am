check_PROGRAMS=test
AM_CPPFLAGS=-I$(top_srcdir)/../srvf/src -I$(top_srcdir)/src
AM_LDFLAGS=-L$(top_srcdir)/../srvf/src/.libs -L$(top_srcdir)/src/.libs
LIBS+=-lsrvf -lsrvf-pmatch -lboost_test_exec_monitor -lgomp
test_SOURCES=partialmatch_tests.cc main.cc
