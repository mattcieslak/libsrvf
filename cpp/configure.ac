AC_INIT([libsrvf], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Main srvf package is always enabled
AC_CONFIG_SUBDIRS([srvf])

# Optional plotting support
AC_ARG_ENABLE(plot,
  [AS_HELP_STRING([--enable-plot],
    [Enable plotting support (requires FLTK and OpenGL)])],
  [if test "$enableval" = yes; then AC_CONFIG_SUBDIRS([plot]) fi],
  [])

# Optional partial matching algorithms
AC_ARG_ENABLE(pmatch,
  [AS_HELP_STRING([--enable-pmatch],
    [Enable partial matching algorithms])],
  [if test "$enableval" = yes; then AC_CONFIG_SUBDIRS([pmatch]) fi],
  [])

# Enable unit tests?
enable_unit_tests=no
AC_ARG_ENABLE([unit-tests],
  [AS_HELP_STRING([--enable-unit-tests],
    [Enable unit tests])],
  [enable_unit_tests="yes"], [])
AM_CONDITIONAL([COND_UNIT_TESTS],[test "$enable_unit_tests" = "yes"])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
