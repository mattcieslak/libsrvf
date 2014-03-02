include_directories(${LIBSRVF_SOURCE_DIR}/include)

add_library(srvf dpnbhd.cc fileio.cc functions.cc interp.cc matchselect.cc matrix.cc opencurves.cc partialmatch.cc plf.cc qmap.cc reparam.cc rotate.cc rotselect.cc srvf.cc util.cc)

target_link_libraries(srvf gsl gslcblas)
