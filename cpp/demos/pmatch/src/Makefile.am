bin_PROGRAMS=srvf_pmatch

AM_CPPFLAGS=-I$(top_srcdir)/srvf/src $(OPENMP_CXXFLAGS)
AM_LDFLAGS=$(SRVF_LDFLAGS) $(OPENMP_CXXFLAGS)

srvf_pmatch_SOURCES=main.cc
srvf_pmatch_LDADD=$(top_srcdir)/srvf/src/libsrvf.la
srvf_pmatch_LDFLAGS=-static
