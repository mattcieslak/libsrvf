bin_PROGRAMS=srvf_pmatch_ui

AM_CPPFLAGS=$(SRVF_CPPFLAGS) $(FLTK_CFLAGS)
AM_LDFLAGS=$(SRVF_LDFLAGS) $(FLTK_LDFLAGS)

LIBS+=$(SRVF_LIBS) -lgomp

srvf_pmatch_ui_SOURCES=main.cc ui.cc

EXTRA_DIST=data
