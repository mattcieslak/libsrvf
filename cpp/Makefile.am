ACLOCAL_AMFLAGS='-Im4'

if COND_PLOT
MAYBE_PLOT=plot
endif

SUBDIRS=srvf $(MAYBE_PLOT) demos

EXTRA_DIST=testdata

bindistdir=$(PACKAGE)-$(VERSION)-bin
bindist: $(DISTFILES)
	mkdir -p $(bindistdir)
	mkdir -p $(bindistdir)/libs
	cp -d srvf/src/.libs/libsrvf*.so* $(bindistdir)/libs
	cp demos/pmatch/src/srvf_pmatch $(bindistdir)
if COND_PLOT
	cp -d plot/src/.libs/libsrvf*.so* $(bindistdir)/libs
	cp demos/pmatch-gui/src/srvf_pmatch_gui $(bindistdir)
endif
	tar -cvzf $(bindistdir).tar.gz $(bindistdir)
	rm -rf $(bindistdir)
