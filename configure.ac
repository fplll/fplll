#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
# don't forget to update library version in Makefile.am before a release
# https://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html#Updating-version-info
AC_INIT(fplll, 5.0.0)
AC_CONFIG_SRCDIR([fplll/fplll.cpp])

# cf http://comments.gmane.org/gmane.comp.sysutils.autoconf.general/15737

# First extract pieces from the version number string
FPLLL_MAJOR_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $1}'`
FPLLL_MINOR_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $2}'`
FPLLL_MICRO_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $3}'`
FPLLL_VERSION=$FPLLL_MAJOR_VERSION.$FPLLL_MINOR_VERSION.$FPLLL_MICRO_VERSION$FPLLL_MICRO_VERSION_SUFFIX
FPLLL_VERSION_NUMBER=`expr $FPLLL_MAJOR_VERSION \* 1000000 + $FPLLL_MINOR_VERSION \* 1000 + $FPLLL_MICRO_VERSION`

AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AM_SILENT_RULES([yes])

AC_CONFIG_HEADER(config.h)

DEFAULT_CXXFLAGS="-O3 -Wall"

CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

#AC_DISABLE_SHARED
AC_PROG_LIBTOOL

AC_PROG_INSTALL

# Checks for programs.
AC_LANG_CPLUSPLUS
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])

# Checks for libraries.
AC_ARG_WITH(gmp,
   AS_HELP_STRING([--with-gmp=@<:@=DIR@:>@], [GMP install directory]), [
      CPPFLAGS="$CPPFLAGS -I$withval/include"
      LDFLAGS="$LDFLAGS -L$withval/lib"
      gmp_lib_path="$withval/lib"
   ])

AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR(
[GNU MP not found, see http://gmplib.org])])

AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
[GMP version too old, need >= 4.2.0, see http://gmplib.org])])

AC_ARG_WITH(mpfr,
   AS_HELP_STRING([--with-mpfr=@<:@=DIR@:>@], [MPFR install directory]), [
      CPPFLAGS="$CPPFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib"
       mpfr_lib_path="$withval/lib"
   ])

AC_CHECK_LIB(mpfr, mpfr_add, , [AC_MSG_ERROR(
[MPFR not found, see http://www.mpfr.org])])

AC_CHECK_LIB(mpfr, mpfr_fms, , [AC_MSG_ERROR(
[MPFR version too old, need >= 2.3.0, see http://www.mpfr.org])])

AC_ARG_WITH(qd,
   AS_HELP_STRING([--with-qd=@<:@=DIR@:>@], [quaddouble install directory]), [
      CPPFLAGS="$CPPFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib"
       qd_lib_path="$withval/lib"
   ])

AC_CHECK_HEADER(qd/dd_real.h,
    have_libqd="yes",
    have_libqd="no")

if test "x${have_libqd}" = "xyes"; then
   LIBQD_LIBADD="-lqd"
   AC_DEFINE([FPLLL_WITH_QD], [1], [defined when libqd-dev is available])
fi

AC_SUBST(LIBQD_LIBADD)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h cstdio iostream string limits vector])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow rint sqrt strtol])

# enumeration dimension
max_enumeration_dimension=256

AC_ARG_WITH(max-enumeration-dimension,
        AS_HELP_STRING([--with-max-enumeration-dimension@<:@=DIM@:>@],
         [Maximum supported enumeration dimension (default: 256).]),
        [max_enumeration_dimension=$withval])

AC_DEFINE_UNQUOTED([FPLLL_MAX_ENUM_DIMENSION], $max_enumeration_dimension, [maximum supported enumeration dimension])

# Store version numbers in header

AC_DEFINE_UNQUOTED([FPLLL_MAJOR_VERSION],[$FPLLL_MAJOR_VERSION],[major version])
AC_DEFINE_UNQUOTED([FPLLL_MINOR_VERSION],[$FPLLL_MINOR_VERSION],[minor version])
AC_DEFINE_UNQUOTED([FPLLL_MICRO_VERSION],[$FPLLL_MICRO_VERSION],[micro version])
AC_DEFINE_UNQUOTED([FPLLL_VERSION],[$FPLLL_VERSION],[full version])
AC_DEFINE_UNQUOTED([FPLLL_VERSION_INFO],[$FPLLL_VERSION_INFO],[long version string])

AC_CONFIG_HEADERS([fplll/fplll_config.h])

AC_CONFIG_FILES([Makefile
                 fplll/Makefile
                 tests/Makefile
                 fplll.pc])
AC_OUTPUT

