SRVF_LIBDIR=../build/src
SRVF_LIBS=-lsrvf
SRVF_INCDIR=../include

# Change this to the command for compiling MEX files on your system.
MEXCMD=mkoctfile --mex

MEXFILES=srvf_fa_optimal_reparam.mex

# Rule for building a .mex file from a .cc file.
# Links libsrvf.a in statically so that this will work even if the 
# library has been built but not installed.  Other libraries are 
# still linked dynamically.
%.mex: %.cc
	$(MEXCMD) -I$(SRVF_INCDIR) -L$(SRVF_LIBDIR) -Wl,-rpath=$(SRVF_LIBDIR) $< $(SRVF_LIBS) -lgsl -lgslcblas

all: $(MEXFILES)

clean:
	rm -f $(MEXFILES)
	rm -f *.o
