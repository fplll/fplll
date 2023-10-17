# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([fplll],[5.4.5])
AC_CONFIG_SRCDIR([fplll/fplll.cpp])
AM_MAINTAINER_MODE([enable])
AC_CANONICAL_TARGET

# cf http://comments.gmane.org/gmane.comp.sysutils.autoconf.general/15737

# First extract pieces from the version number string
FPLLL_MAJOR_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $1}'`
FPLLL_MINOR_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $2}'`
FPLLL_MICRO_VERSION=`echo AC_PACKAGE_VERSION | awk -F. '{print $3}'`
FPLLL_VERSION=$FPLLL_MAJOR_VERSION.$FPLLL_MINOR_VERSION.$FPLLL_MICRO_VERSION$FPLLL_MICRO_VERSION_SUFFIX
FPLLL_VERSION_NUMBER=`expr $FPLLL_MAJOR_VERSION \* 1000000 + $FPLLL_MINOR_VERSION \* 1000 + $FPLLL_MICRO_VERSION`

AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])

DEFAULT_CXXFLAGS="-O3 -Wall"

CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

LT_INIT

# The version of the libtool library is of the form current:revision:age
#
# See http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
#
# When doing a release, they should be updated like this:
# 1. If no interfaces changed, only implementations: just increment
#    revision.
# 2. If interfaces were added, none removed: increment current, set
#    revision to zero and increment age.
# 3. If interfaces were removed (breaks backward compatibility): increment
#    current, and set both revision and age to zero.

FPLLL_LT_CURRENT=8
FPLLL_LT_REVISION=1
FPLLL_LT_AGE=0

AC_SUBST(FPLLL_LT_CURRENT)
AC_SUBST(FPLLL_LT_REVISION)
AC_SUBST(FPLLL_LT_AGE)

AC_PROG_INSTALL

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])

# Checks for libraries.
AC_ARG_ENABLE(mpir,
        AS_HELP_STRING([--enable-mpir],
         [Enable MPIR instead of GMP]))

AS_IF([test "x$enable_mpir" = "xyes"], [
	AC_ARG_WITH(mpir,
	   AS_HELP_STRING([--with-mpir=@<:@=DIR@:>@], [MPIR install directory]), [
	      CPPFLAGS="$CPPFLAGS -I$withval/include"
	      LDFLAGS="$LDFLAGS -L$withval/lib"
	      mpir_lib_path="$withval/lib"
	   ])
	AC_CHECK_LIB(mpir, __gmp_get_memory_functions, , [AC_MSG_ERROR(
	[MPIR version >= 1.0.0 needed, see http://mpir.org])])
])

AS_IF([test "x$enable_mpir" != "xyes"], [
	AC_ARG_WITH(gmp,
	   AS_HELP_STRING([--with-gmp=@<:@=DIR@:>@], [GMP install directory]), [
	      CPPFLAGS="$CPPFLAGS -I$withval/include"
	      LDFLAGS="$LDFLAGS -L$withval/lib"
	      mpir_lib_path="$withval/lib"
	   ])
	AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
	[GMP version >= 4.2.0 needed, see http://gmplib.org])])
])

AC_ARG_WITH(mpfr,
   AS_HELP_STRING([--with-mpfr=@<:@=DIR@:>@], [MPFR install directory]), [
      CPPFLAGS="$CPPFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib"
       mpfr_lib_path="$withval/lib"
   ])

AC_CHECK_LIB(mpfr, mpfr_fms, , [AC_MSG_ERROR(
[MPFR version >= 2.3.0 needed, see http://www.mpfr.org])])

AC_ARG_WITH(qd, AS_HELP_STRING([--with-qd=@<:@=DIR@:>@], [quaddouble install directory]),)

# Act as if --with-qd was passed, by default.
AS_IF([test -z "$with_qd"], [with_qd=yes])

# Passing --with-qd (without an argument) will use pkg-config to
# detect qd on the system, or fall back to manual detection of headers
# and libraries if that didn't work for some reason.
#
# We only show a warning if we cannot find libqd, it is optional but recommended
AS_IF([test "x$with_qd" = "xyes"], [
  PKG_CHECK_MODULES([LIBQD], [qd], [have_libqd="yes"], [
    # fall back to manual search.
    AC_SEARCH_LIBS(c_dd_add,
                   qd,
                   [
                     AC_CHECK_HEADER(qd/dd_real.h,
                       [have_libqd="yes"],
                       [AC_MSG_WARN([unable to find header qd/dd_real.h])])
                   ],
                   [ AC_MSG_WARN([unable to find c_dd_add() in libqd]) ]
    )]
  )],
  # else if
  [test "x$with_qd" != "xno"], [
  CPPFLAGS="$CPPFLAGS -I$with_qd/include"
  AC_CHECK_HEADER(qd/dd_real.h, [
    LDFLAGS="$LDFLAGS -L$with_qd/lib"
    qd_lib_path="$with_qd/lib"
    LIBQD_LIBS="-lqd"
    have_libqd="yes"],
    [ AC_MSG_WARN([unable to find header qd/dd_real.h]) ]
  )]
)

AS_IF([test "x${have_libqd}" = "xyes"], [
  AC_DEFINE([FPLLL_WITH_QD], [1], [defined when libqd is usable])
])

# This will be non-empty if we found qd using pkg-config, or if a path
# was passed to --with-qd=<path> and we set the variable ourselves. If
# we found libqd via a manual library search, on the other hand, then
# the linker flag will be found in LIBS instead because that's just how
# the AC_SEARCH_LIBS macro works.
AC_SUBST(LIBQD_LIBS)

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

AX_PTHREAD()

CHECK_ATOMIC()

# enumeration dimension
max_enum_dim=256

AC_ARG_WITH(max-enum-dim,
        AS_HELP_STRING([--with-max-enum-dim@<:@=DIM@:>@],
         [Maximum supported enumeration dimension (default: 256).]),
        [max_enum_dim=$withval])

AC_DEFINE_UNQUOTED([FPLLL_MAX_ENUM_DIM], $max_enum_dim, [maximum supported enumeration dimension])

AC_ARG_ENABLE(recursive-enum,
        AS_HELP_STRING([--disable-recursive-enum],
         [Disable recursive enumeration]))

AS_IF([test "x$enable_recursive_enum" != "xno"], [
        AC_DEFINE([FPLLL_WITH_RECURSIVE_ENUM], [1], [recursive enumeration enabled])])


max_parallel_enum_dim=120
AC_ARG_WITH(max-parallel-enum-dim,
        AS_HELP_STRING([--with-max-parallel-enum-dim@<:@=DIM@:>@],
         [Maximum supported parallel enumeration dimension (max=160, default=120, disable=0).]),
        [max_parallel_enum_dim=$withval])

AC_DEFINE_UNQUOTED([FPLLL_MAX_PARALLEL_ENUM_DIM], $max_parallel_enum_dim, [parallel enumeration enabled])
AM_CONDITIONAL(FPLLL_PARALLEL_ENUM, test "x$max_parallel_enum_dim" != "x0")

# Use external enumeration
AC_ARG_WITH(extenum-dir,
  AS_HELP_STRING([--with-extenum-dir@<:@=DIR@:>@], [Specify external enumeration library install directory]),
  [LDFLAGS="$LDFLAGS -L$withval"])
AC_ARG_WITH(extenum-lib,
  AS_HELP_STRING([--with-extenum-lib@<:@=FILE@:>@], [Specify external enumeration library name]),
  [extenum_lib=$withval])
AC_ARG_WITH(extenum-func,
  AS_HELP_STRING([--with-extenum-func@<:@=FUNC@:>@], [Specify external enumeration function name]),
  [extenum_func=$withval])
AS_IF([test "x$extenum_lib" != "x"],[
  AS_IF([test "x$extenum_func" != "x"],,[extenum_func=$extenum_lib])
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([
          #include "fplll/enum/enumerate_ext_api.h"
          extern extenum_fc_enumerate $extenum_func;
    ],[
          auto abc = & $extenum_func;
    ])]
    ,,AC_MSG_ERROR([Failed to find external enumeration]))
  EXTENUM_LIBS="-l$extenum_lib"
  AC_DEFINE_UNQUOTED([FPLLL_EXTENUM_FUNC], $extenum_func, [use external enumeration function])
])
AC_SUBST(EXTENUM_LIBS)

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

