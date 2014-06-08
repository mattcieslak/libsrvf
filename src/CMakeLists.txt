include_directories(${LIBSRVF_SOURCE_DIR}/include)

add_library(srvf SHARED c_api.cc dpnbhd.cc fileio.cc functions.cc interp.cc matchselect.cc matrix.cc opencurves.cc partialmatch.cc plf.cc qmap.cc reparam.cc rotate.cc rotselect.cc srvf.cc util.cc)

if(GSL)
  target_link_libraries(srvf gsl gslcblas)
  set_target_properties(srvf PROPERTIES COMPILE_FLAGS "-DUSE_GSL")
endif()

install(TARGETS srvf LIBRARY DESTINATION lib)
