#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ticcltools],[0.11],[lamasoftware@science.ru.nl]) #adapt version number in codemeta.json as well
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([.])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

if test x"${CXXFLAGS+set}" = xset; then
  # the user set CXXFLAGS; don't override it.
  cxx_flags_were_set=true
else
  cxx_flags_were_set=false
fi
if $cxx_flags_were_set; then
  CXXFLAGS=$CXXFLAGS
fi

# Checks for programs.
AC_PROG_CXX( [g++] )

# when running tests, use CXX
AC_LANG([C++])

# check OpenMP support
AC_OPENMP
if test "x$ac_cv_prog_cxx_openmp" != "x"; then
  if test "x$ac_cv_prog_cxx_openmp" != "xunsupported"; then
    CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
    AC_DEFINE([HAVE_OPENMP], [1] , [Define to 1 if you have OpenMP] )
  else
    AC_MSG_NOTICE([We don't have OpenMP. Multithreaded operation is disabled])
  fi
fi

# use libtool
LT_INIT

#checks for libraries.

# Checks for header files.
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_INT8_T
AC_C_INLINE

# Checks for library functions.

if test "x$prefix" = "xNONE"; then
   prefix="/usr/local"
fi

PKG_PROG_PKG_CONFIG
if test "x$PKG_CONFIG_PATH" = x; then
    export PKG_CONFIG_PATH="/usr/lib/pkgconfig:$prefix/lib/pkgconfig"
else
    export PKG_CONFIG_PATH="/usr/lib/pkgconfig:$PKG_CONFIG_PATH"
fi

AC_OSX_PKG( [icu4c] )

PKG_CHECK_MODULES([ICU], [icu-uc >= 56 icu-io] )
CXXFLAGS="$CXXFLAGS $ICU_CFLAGS"
LIBS="$ICU_LIBS $LIBS"

PKG_CHECK_MODULES([XML2], [libxml-2.0 >= 2.6.16] )
CXXFLAGS="$CXXFLAGS $XML2_CFLAGS"
LIBS="$XML2_LIBS $LIBS"

PKG_CHECK_MODULES([ticcutils], [ticcutils >= 0.32] )
CXXFLAGS="$ticcutils_CFLAGS $CXXFLAGS"
LIBS="$ticcutils_LIBS $LIBS"

AC_CONFIG_FILES([
  Makefile
  m4/Makefile
  src/Makefile
  include/Makefile
  include/ticcl/Makefile
  docs/Makefile
])
AC_OUTPUT
