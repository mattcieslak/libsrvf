cmake_minimum_required(VERSION 2.6)
project(LIBSRVF)

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

add_subdirectory(src)
add_subdirectory(tests)
