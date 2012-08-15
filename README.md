libsrvf
=======

A library for shape analysis of elastic curves, using the square root 
velocity framework


0. Required third-party libraries
------------------------------------------------------------------------------
libsrvf depends on the GNU Scientific Library (GSL) for its SVD routine.  
See www.gnu.org/software/gsl/ for instructions on getting this set up.

FLTK 1.1 and OpenGL 1.1 (or newer) are required if you want to enable the 
plotting features of libsrvf.  See www.fltk.org for instructions on 
installing FLTK.

If you want to run the unit tests, then you will need the Boost libraries.  
See www.boost.org for instructions on installing these.

Last but not least, you need a working C++ compiler and linker,  and a 
compatible 'make' program.  The configure script should check all of this 
for you.


I. Quick start guide
------------------------------------------------------------------------------

To build the project with plotting enabled and then install it under 
the directory /home/fred, do this:

> ./configure --prefix=/home/fred --enable-plot
> make
> make install

This will install the headers under /home/fred/include/srvf, the libraries 
will be installed under /home/fred/lib, and the demo programs will be 
installed in /home/fred/bin.  If you don't have FLTK and OpenGL, then omit 
the --enable-plot:

> ./configure --prefix=/home/fred
> make
> make install

If you have doxygen installed and want to generate the documentation, then do

> doxygen

The generated documentation will be in the doc/ directory.

If you want to run the unit tests, then do the following after you build 
the library:

> make check
> cd tests
> ./runtests

In order to run the unit tests, you need to have Boost, FLTK, and OpenGL 
installed, and you need to pass --enable-plot to the configure script so that 
the plotting features of the library get built.

To see the full list of supported configuration options, do

> ./configure --help


II. Demo programs
------------------------------------------------------------------------------

The demos/ directory contains several small programs using various 
features of the library.  These programs are built when the library 
is built, provided that the features they require (i.e. plotting) have 
been enabled.

1. srvf_shapedist
 - Description: computes the shape distance between two curves
 - Directory: demos/shapedist
 - Executable: srvf_shapedist
 - Required optional features: none

2. srvf_pmatch
 - Description: finds partial matches between two curves
 - Directory: demos/pmatch
 - Executable: srvf_pmatch
 - Required optional features: none

3. srvf_pmatch_gui
 - Description: a GUI version of srvf_pmatch
 - Directory: demos/pmatch-gui
 - Executable: srvf_pmatch_gui
 - Required optional features: --enable-plot


III. Matlab / Octave functions
------------------------------------------------------------------------------

The matlab/ directory contains Matlab / Octave implementations of most of 
the library features.  Some of the features are implemented as MEX files, 
and need to be compiled and linked against libsrvf.  Currently, this must 
be done as a separate step, after building the library.  To do this, first 
edit the Makefile in the matlab/ directory to reflect your system, then run 
make in the matlab/ directory.

