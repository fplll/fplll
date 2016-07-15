#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(fplll,4.0.5)
AC_CONFIG_SRCDIR([fplll/fplll.cpp])

# cf http://comments.gmane.org/gmane.comp.sysutils.autoconf.general/15737

# First extract pieces from the version number string
FPLLL_MAJOR_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $1}'`
FPLLL_MINOR_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $2}'`
FPLLL_MICRO_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $3}'`
FPLLL_VERSION=$FPLLL_MAJOR_VERSION.$FPLLL_MINOR_VERSION.$FPLLL_MICRO_VERSION$FPLLL_MICRO_VERSION_SUFFIX
FPLLL_VERSION_NUMBER=`expr $FPLLL_MAJOR_VERSION \* 1000000 + $FPLLL_MINOR_VERSION \* 1000 + $FPLLL_MICRO_VERSION`
AGE=`expr $FPLLL_MAJOR_VERSION '*' 1000 + $FPLLL_MINOR_VERSION`
REVISION=$FPLLL_MICRO_VERSION
FPLLL_VERSION_INFO=$CURRENT:$REVISION:$AGE

AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_CONFIG_HEADER(config.h)

DEFAULT_CFLAGS="-O2 -Wall"

CFLAGS=${CFLAGS:-$DEFAULT_CFLAGS}
CXXFLAGS=${CXXFLAGS:-$DEFAULT_CFLAGS}

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


# Recursive enumeration
AC_ARG_ENABLE([recursive-enumation],
    AS_HELP_STRING([--disable-recursive-enumation], [Disable (faster) recursive enumeration implementation.]))

AS_IF([test "x$enable_recursive_enum" != "xno"],
      [have_recursive_enum=yes],
      [have_recursive_enum=no])

AS_IF([test "x$have_recursive_enum" = "xyes"],
      AC_DEFINE([FPLLL_WITH_RECURSIVE_ENUM], [1], [defined when recursive enumeration is used]))

# Store version numbers in header

AC_DEFINE_UNQUOTED([FPLLL_MAJOR_VERSION],[$FPLLL_MAJOR_VERSION],[major version])
AC_DEFINE_UNQUOTED([FPLLL_MINOR_VERSION],[$FPLLL_MINOR_VERSION],[minor version])
AC_DEFINE_UNQUOTED([FPLLL_MICRO_VERSION],[$FPLLL_MICRO_VERSION],[micro version])
AC_DEFINE_UNQUOTED([FPLLL_SONUM],[$FPLLL_SONUM],[sonum])
AC_DEFINE_UNQUOTED([FPLLL_VERSION],[$FPLLL_VERSION],[full version])
AC_DEFINE_UNQUOTED([FPLLL_VERSION_INFO],[$FPLLL_VERSION_INFO],[long version string])

AC_CONFIG_HEADERS([fplll/fplll_config.h])

AC_CONFIG_FILES([Makefile
                 fplll/Makefile
                 tests/Makefile
                 fplll.pc])
AC_OUTPUT

