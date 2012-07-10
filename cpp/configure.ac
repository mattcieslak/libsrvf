AC_INIT([libsrvf], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Main srvf package is always enabled
AC_CONFIG_SUBDIRS([srvf])

# Optional plotting support
plot_enabled="no"
AC_ARG_ENABLE(plot,
  [AS_HELP_STRING([--enable-plot],
    [Enable plotting support (requires >= FLTK 1.1 and OpenGL)])],
  [if test "$enableval" = yes; then 
    AC_CONFIG_SUBDIRS([plot]) 
    plot_enabled="yes"
  fi],
  [])
AM_CONDITIONAL([COND_PLOT],[test "$plot_enabled" = "yes"])


# Optional partial matching algorithms
pmatch_enabled="no"
AC_ARG_ENABLE(pmatch,
  [AS_HELP_STRING([--enable-pmatch],
    [Enable partial matching algorithms])],
  [if test "$enableval" = yes; then 
    AC_CONFIG_SUBDIRS([pmatch]) 
    pmatch_enabled="yes"
  fi],
  [])
AM_CONDITIONAL([COND_PMATCH],[test "$pmatch_enabled" = "yes"])


# Build the demo programs?
demos_enabled="no"
AC_ARG_ENABLE(demos,
  [AS_HELP_STRING([--enable-demos],
    [Build the demo programs])],
  [if test "$enableval" = yes; then 
    AC_CONFIG_SUBDIRS([demos]) 
    demos_enabled="yes"
  fi],
  [])
AM_CONDITIONAL([COND_DEMOS],[test "$demos_enabled" = "yes"])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
