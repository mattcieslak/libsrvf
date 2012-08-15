SRVF_LIBDIR=../src/.libs
SRVF_LIBS=-lsrvf
SRVF_INCDIR=..

# Change this to the command for compiling MEX files on your system.
MEXCMD=mkoctfile --mex

MEXFILES=srvf_fa_groupwise_reparam.mex srvf_fa_karcher_mean.mex srvf_fa_optimal_reparam.mex srvf_karcher_mean.mex srvf_optimal_reparam.mex srvf_pmatch_find_matches.mex

# Rule for building a .mex file from a .cc file.
# Links libsrvf.a in statically so that this will work even if the 
# library has been built but not installed.  Other libraries are 
# still linked dynamically.
%.mex: %.cc
	$(MEXCMD) -I$(SRVF_INCDIR) -L$(SRVF_LIBDIR) $< -Wl,-static $(SRVF_LIBS) -Wl,-Bdynamic -lgsl -lgslcblas

all: $(MEXFILES)

clean:
	rm -f $(MEXFILES)

