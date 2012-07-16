# Paths for SRVF libs and headers in the build tree
SRVF_CPPFLAGS="-I$(top_srcdir)/../srvf/src"
SRVF_LDFLAGS="-L$(top_srcdir)/../srvf/src/.libs"
SRVF_LIBS=" -lsrvf"

AM_CPPFLAGS=$(SRVF_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS=$(SRVF_LDFLAGS)
LIBS+=$(SRVF_LIBS)

lib_LTLIBRARIES=libsrvf-pmatch.la

#libsrvf_pmatch_la_LDFLAGS= -static

srvf_includedir=$(includedir)/srvf
srvf_include_HEADERS=paretoset.h partialmatch.h pmatch_util.h rotselect.h

libsrvf_pmatch_la_SOURCES=partialmatch.cc pmatch_util.cc rotselect.cc
