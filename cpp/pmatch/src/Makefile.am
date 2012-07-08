AM_CPPFLAGS=-I$(top_srcdir)/../srvf/src
AM_LDFLAGS=-L$(top_srcdir)/../srvf/src/.libs

lib_LTLIBRARIES=libsrvf-pmatch.la

libsrvf_pmatch_la_LDFLAGS= -static

srvf_includedir=$(includedir)/srvf
srvf_include_HEADERS=paretoset.h partialmatch.h

LIBS+=-lsrvf -lgomp

libsrvf_pmatch_la_SOURCES=partialmatch.cc
