AC_INIT([libsrvf], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_AUX_DIR([.])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Check whether plotting support is enabled in the library
libsrvf_plot_enabled="no"
AC_ARG_ENABLE(plot,
  [AS_HELP_STRING([--enable-plot],
    [Enable plotting support (requires >= FLTK 1.1 and OpenGL)])],
  [if test "$enableval" = yes; then 
    libsrvf_plot_enabled="yes"
  fi],
  [])

# Check whether partial matching is enabled in the library
libsrvf_pmatch_enabled="no"
AC_ARG_ENABLE(pmatch,
  [AS_HELP_STRING([--enable-pmatch],
    [Enable partial matching algorithms])],
  [if test "$enableval" = yes; then 
    libsrvf_pmatch_enabled="yes"
  fi],
  [])


# Conditionally enable each program if its dependencies are satisfied
demo_pmatch_enabled="no"
if test "$libsrvf_plot_enabled" = yes && 
   test "$libsrvf_pmatch_enabled" = yes; then :
  demo_pmatch_enabled="yes"
  AC_CONFIG_SUBDIRS([pmatch]) 
fi

AM_CONDITIONAL([COND_DEMO_PMATCH],[test "$demo_pmatch_enabled" = "yes"])


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
