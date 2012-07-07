AM_CPPFLAGS=-I$(top_srcdir)/../srvf/src
AM_LDFLAGS=-L$(top_srcdir)/../srvf/src/.libs

lib_LTLIBRARIES=libsrvf-pmatch.la

srvf_includedir=$(includedir)/srvf

LIBS+=-lsrvf -lgomp

libsrvf_pmatch_la_SOURCES=partialmatch.cc
