libsrvf
=======

A library for shape analysis of elastic curves, using the square root 
velocity framework


Required third-party libraries
------------------------------------------------------------------------------
If you want to use libsrvf for rotational alignment, define the symbol 
USE_GSL=1 during compilation.  libsrvf depends on the GNU Scientific Library (GSL)
for its SVD routine.  See www.gnu.org/software/gsl/ for instructions on getting this set up.

If you want to run the unit tests, then you will need the Boost libraries.  
See www.boost.org for instructions on installing these.

To build the library from source, you will need CMake, available at http://www.cmake.org/