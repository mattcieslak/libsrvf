bin_PROGRAMS=srvf_shapedist

AM_CPPFLAGS=$(OPENMP_CXXFLAGS)
AM_LDFLAGS=$(OPENMP_CXXFLAGS)

srvf_shapedist_SOURCES=shapedist.cc
srvf_shapedist_LDADD=$(top_srcdir)/src/libsrvf.la
srvf_shapedist_LDFLAGS=-static
