include_fpllldir=$(includedir)/fplll
TOPSRCDIR = $(srcdir)/..
strategydir = $(pkgdatadir)/strategies
# see https://stackoverflow.com/questions/5867136/autoconf-how-to-get-installation-paths-into-config-h
AM_CPPFLAGS = -DFPLLL_DEFAULT_STRATEGY_PATH=\"$(strategydir)\" -DFPLLL_DEFAULT_STRATEGY=\"$(strategydir)/default.json\" -I$(TOPSRCDIR)
EXTRA_DIST = io/json.hpp ballvol.const factorial.const

CLEANFILES = *.gcov .libs/*.gcda .libs/*.gcno *.gcno *.gcda */*.gcno */*.gcda

nobase_include_fplll_HEADERS=defs.h fplll.h fplll_config.h \
	nr/dpe.h \
	nr/matrix.h \
	nr/matrix.cpp \
	nr/nr_FP_dd.inl \
	nr/nr_FP_d.inl \
	nr/nr_FP_dpe.inl \
	nr/nr_FP.inl \
	nr/nr_FP_ld.inl \
	nr/nr_FP_misc.inl \
	nr/nr_FP_mpfr.inl \
	nr/nr_FP_qd.inl \
	nr/nr.h \
	nr/nr_rand.inl \
	nr/nr_Z_d.inl \
	nr/nr_Z.inl \
	nr/nr_Z_l.inl \
	nr/nr_Z_misc.inl \
	nr/nr_Z_mpz.inl \
	nr/numvect.h \
	util.h \
	svpcvp.h bkz.h lll.h gso_interface.h gso_gram.h gso.h  \
	enum/evaluator.h \
	wrapper.h \
	bkz_param.h \
	enum/enumerate.h enum/enumerate_base.h enum/enumerate_ext.h enum/enumerate_ext_api.h \
	pruner/pruner.h pruner/pruner_simplex.h \
	householder.h hlll.h \
	threadpool.h io/thread_pool.hpp

bin_PROGRAMS=fplll latticegen
check_PROGRAMS=llldiff
EXTRA_PROGRAMS=fplll_dbg latticegen_dbg
lib_LTLIBRARIES=libfplll.la
EXTRA_LTLIBRARIES=libfplllv.la libfpllld.la

# fplll bin
fplll_SOURCES=main.cpp main.h
fplll_LDADD=libfplll.la
fplll_dbg_SOURCES=$(fplll_SOURCES)
fplll_dbg_CPPFLAGS=-DDEBUG $(AM_CPPFLAGS)
fplll_dbg_CXXFLAGS=-g $(AM_CXXFLAGS)
fplll_dbg_LDADD=libfpllld.la

# latticegen bin
latticegen_SOURCES=latticegen.cpp
latticegen_LDADD=libfplll.la
latticegen_dbg_SOURCES=$(latticegen_SOURCES)
latticegen_dbg_CPPFLAGS=-DDEBUG $(AM_CPPFLAGS)
latticegen_dbg_CXXFLAGS=-g $(AM_CXXFLAGS)
latticegen_dbg_LDADD=libfpllld.la

# llldiff bin
llldiff_SOURCES=llldiff.cpp
llldiff_LDADD=libfplll.la

# libfplll
libfplll_la_SOURCES=fplll.cpp fplll.h \
	util.cpp util.h \
	enum/topenum.cpp enum/topenum.h \
	enum/enumerate.cpp enum/enumerate.h \
	enum/enumerate_base.cpp enum/enumerate_base.h \
	enum/enumerate_ext.cpp enum/enumerate_ext.h enum/enumerate_ext_api.h \
	enum/evaluator.cpp enum/evaluator.h \
	lll.cpp lll.h \
	wrapper.cpp wrapper.h \
	bkz.cpp bkz.h \
	bkz_param.cpp bkz_param.h \
	gso_interface.cpp gso_interface.h gso_gram.cpp gso_gram.h gso.cpp gso.h \
	pruner/pruner.cpp \
	pruner/pruner.h \
	pruner/pruner_simplex.h \
	pruner/pruner_cost.cpp \
	pruner/pruner_optimize.cpp \
	pruner/pruner_optimize_tc.cpp \
	pruner/pruner_optimize_tp.cpp \
	pruner/pruner_prob.cpp \
	pruner/pruner_util.cpp \
	householder.cpp householder.h hlll.cpp hlll.h \
	io/json.hpp \
	threadpool.h threadpool.cpp io/thread_pool.hpp
libfplll_la_CXXFLAGS=$(PTHREAD_CFLAGS)

EXTRA_libfplll_la_SOURCES= svpcvp.cpp
libfplll_la_LIBADD=$(LIBQD_LIBS) $(PTHREAD_LIBS) $(EXTENUM_LIBS)
libfplll_la_LDFLAGS=-no-undefined -version-info @FPLLL_LT_CURRENT@:@FPLLL_LT_REVISION@:@FPLLL_LT_AGE@ $(PTHREAD_CFLAGS)

if FPLLL_PARALLEL_ENUM
libfplll_la_SOURCES+=\
	enum-parallel/enumeration.h enum-parallel/enumlib.h enum-parallel/fplll_types.h \
	enum-parallel/enumlib.cpp \
	enum-parallel/enumlib_dim.160.cpp enum-parallel/enumlib_dim.150.cpp enum-parallel/enumlib_dim.140.cpp enum-parallel/enumlib_dim.130.cpp \
	enum-parallel/enumlib_dim.120.cpp enum-parallel/enumlib_dim.110.cpp enum-parallel/enumlib_dim.100.cpp enum-parallel/enumlib_dim.90.cpp \
	enum-parallel/enumlib_dim.80.cpp enum-parallel/enumlib_dim.70.cpp enum-parallel/enumlib_dim.60.cpp enum-parallel/enumlib_dim.50.cpp \
	enum-parallel/enumlib_dim.40.cpp enum-parallel/enumlib_dim.30.cpp enum-parallel/enumlib_dim.20.cpp
EXTRA_libfplll_la_SOURCES+=enum-parallel/enumlib_dim.cpp
endif


# libfplll variants
libfplllv_la_SOURCES=$(libfplll_la_SOURCES)
EXTRA_libfplllv_la_SOURCES=$(EXTRA_libfplll_la_SOURCES)
libfplllv_la_CPPFLAGS=$(AM_CPPFLAGS)
libfplllv_la_LIBADD=$(libfplll_la_LIBADD)
libfplllv_la_LDFLAGS=$(libfplll_la_LDFLAGS)

libfpllld_la_SOURCES=$(libfplll_la_SOURCES)
EXTRA_libfpllld_la_SOURCES=$(EXTRA_libfplll_la_SOURCES)
libfpllld_la_CPPFLAGS=-DDEBUG $(AM_CPPFLAGS)
libfpllld_la_CXXFLAGS=-g $(AM_CXXFLAGS)
libfpllld_la_LIBADD=$(libfplll_la_LIBADD)
libfpllld_la_LDFLAGS=$(libfplll_la_LDFLAGS)

install-data-hook:
	echo "#include <fplll/fplll.h>" > "$(builddir)/fplll.h.root"
	$(INSTALL) -m 644 "$(builddir)/fplll.h.root" "${DESTDIR}$(includedir)/fplll.h"
	rm -f "$(builddir)/fplll.h.root"

uninstall-hook:
	if test -d "${DESTDIR}$(includedir)/fplll/enum"; then rmdir "${DESTDIR}$(includedir)/fplll/enum"; fi
	if test -d "${DESTDIR}$(includedir)/fplll/nr"; then rmdir "${DESTDIR}$(includedir)/fplll/nr"; fi
	if test -d "${DESTDIR}$(includedir)/fplll/pruner"; then rmdir "${DESTDIR}$(includedir)/fplll/pruner"; fi
	(test -d "${DESTDIR}$(includedir)/fplll" && rmdir "${DESTDIR}$(includedir)/fplll") || true
	rm -f "${DESTDIR}$(includedir)/fplll.h"

dbg: fplll_dbg latticegen_dbg

