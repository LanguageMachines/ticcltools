#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac 16569 2013-10-07 15:28:10Z sloot $
# $URL $

AC_PREREQ(2.61)
AC_INIT([ticcltools],[0.2],[ko.vandersloot@uvt.nl])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([.])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config.h])

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

AC_OPENMP
if test "x$ac_cv_prog_cxx_openmp" != "xunsupported"; then
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
   AC_DEFINE(HAVE_OPENMP, 1 , Define to 1 if you have OpenMP )
else
   AC_MSG_NOTICE(We don't have OpenMP. Multithreaded operation is disabled)
fi


# use libtool
AC_PROG_LIBTOOL
LT_INIT

#checks for libraries.

# Checks for header files.
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_INT8_T
AC_C_INLINE

# Checks for library functions.

# ugly hack when PKG_CONFIG_PATH isn't defined.
# couldn't get it to work otherwise
if test "x$PKG_CONFIG_PATH" = x; then
    export PKG_CONFIG_PATH="/usr/lib/pkgconfig"
fi

PKG_PROG_PKG_CONFIG
PKG_CHECK_MODULES([XML2], [libxml-2.0 >= 2.6.16] )
CXXFLAGS="$CXXFLAGS $XML2_CFLAGS"
LIBS="$XML2_LIBS $LIBS"

AX_ICU_CHECK( [4.2],
	      [CXXFLAGS="$CXXFLAGS $ICU_CPPSEARCHPATH"
	      LIBS="$ICU_LIBPATH $ICU_LIBS $ICU_IOLIBS $LIBS"],
	      [AC_MSG_FAILURE( "No ICU development environment found. Please check if libicu-dev
 or the like is installed" )] )

AC_ARG_WITH(ticcutils,
       [  --with-ticcutils=DIR       use ticcutils installed in <DIR>;
               note that you can install ticcutils in a non-default directory wi
th
               ./configure --prefix=<DIR> in the ticcutils installation director
y],
       [PKG_CONFIG_PATH="$withval/lib/pkgconfig:$PKG_CONFIG_PATH"],
       [PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$PKG_CONFIG_PATH"])
#AC_MSG_NOTICE( [pkg-config search path: $PKG_CONFIG_PATH] )
PKG_CHECK_MODULES([ticcutils], [ticcutils >= 0.6] )
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
