if COND_PLOT
MAYBE_PLOT=plot
endif

if COND_PMATCH
MAYBE_PMATCH=pmatch
endif

if COND_DEMOS
MAYBE_DEMOS=demos
endif

SUBDIRS=srvf $(MAYBE_PLOT) $(MAYBE_PMATCH) $(MAYBE_DEMOS)

EXTRA_DIST=testdata
bindistdir=$(PACKAGE)-$(VERSION)-bin

bindist: $(DISTFILES)
	mkdir -p $(bindistdir)
	mkdir -p $(bindistdir)/libs
	cp -d srvf/src/.libs/libsrvf*.so* $(bindistdir)/libs
if COND_PMATCH
	cp -d pmatch/src/.libs/libsrvf*.so* $(bindistdir)/libs
endif
if COND_PLOT
	cp -d plot/src/.libs/libsrvf*.so* $(bindistdir)/libs
endif
if COND_DEMOS
if COND_PMATCH
	cp demos/pmatch/src/srvf_pmatch $(bindistdir)
if COND_PLOT
	cp demos/pmatch-gui/src/srvf_pmatch_gui $(bindistdir)
endif	 # COND_PLOT
endif	 # COND_PMATCH
endif  # COND_DEMOS
	tar -cvzf $(bindistdir).tar.gz $(bindistdir)
	rm -rf $(bindistdir)
