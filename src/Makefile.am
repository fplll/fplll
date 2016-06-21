include_fpllldir=$(includedir)/fplll

nobase_include_fplll_HEADERS=defs.h fplll.h \
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
	svpcvp.h fplllv31.h bkz.h lll.h gso.h  \
	enum/evaluator.h \
	fplll_config.h wrapper.h \
  enum/enumerate.h

EXTRA_DIST=pruning.txt

bin_PROGRAMS=fplll latticegen
check_PROGRAMS=llldiff
EXTRA_PROGRAMS=fplll_dbg latticegen_dbg test 
lib_LTLIBRARIES=libfplll.la
EXTRA_LTLIBRARIES=libfplllv.la libfpllld.la

fplll_SOURCES=main.cpp main.h
fplll_LDADD=libfplll.la
fplll_dbg_SOURCES=$(fplll_SOURCES)
fplll_dbg_CPPFLAGS=-DDEBUG $(AM_CPPFLAGS)
fplll_dbg_CXXFLAGS=-g $(AM_CXXFLAGS)
fplll_dbg_LDADD=libfpllld.la

latticegen_SOURCES=latticegen.cpp
latticegen_LDADD=libfplll.la
latticegen_dbg_SOURCES=$(latticegen_SOURCES)
latticegen_dbg_CPPFLAGS=-DDEBUG $(AM_CPPFLAGS)
latticegen_dbg_CXXFLAGS=-g $(AM_CXXFLAGS)
latticegen_dbg_LDADD=libfpllld.la

llldiff_SOURCES=llldiff.cpp
llldiff_LDADD=libfplll.la

test_SOURCES=test.cpp
test_LDADD=libfplll.la

libfplll_la_SOURCES=fplll.cpp fplll.h util.cpp util.h enum/topenum.cpp enum/topenum.h \
	enum/evaluator.cpp enum/evaluator.h
EXTRA_libfplll_la_SOURCES=gso.cpp gso.h lll.cpp lll.h wrapper.cpp wrapper.h \
	bkz.cpp bkz.h enum/enumerate.cpp enum/enumerate.h svpcvp.cpp
libfplll_la_LIBADD=-lgmp -lmpfr
libfplll_la_LDFLAGS=-no-undefined -version-info 1:0:1

libfplllv_la_SOURCES=$(libfplll_la_SOURCES)
EXTRA_libfplllv_la_SOURCES=$(EXTRA_libfplll_la_SOURCES)
libfplllv_la_CPPFLAGS=$(AM_CPPFLAGS)
libfplllv_la_LIBADD=$(libfplll_la_LIBADD)

libfpllld_la_SOURCES=$(libfplll_la_SOURCES)
EXTRA_libfpllld_la_SOURCES=$(EXTRA_libfplll_la_SOURCES)
libfpllld_la_CPPFLAGS=-DDEBUG $(AM_CPPFLAGS)
libfpllld_la_CXXFLAGS=-g $(AM_CXXFLAGS)
libfpllld_la_LIBADD=$(libfplll_la_LIBADD)

install-data-hook:
	echo "#include <fplll/fplll.h>" > "$(srcdir)/fplll.h.root"
	echo "using namespace fplll;" >> "$(srcdir)/fplll.h.root"
	$(INSTALL) -m 644 "$(srcdir)/fplll.h.root" "${DESTDIR}$(includedir)/fplll.h"
	rm -f "$(srcdir)/fplll.h.root"

uninstall-hook:
	if test -d "${DESTDIR}$(includedir)/fplll/enum"; then rmdir "${DESTDIR}$(includedir)/fplll/enum"; fi
	if test -d "${DESTDIR}$(includedir)/fplll/nr"; then rmdir "${DESTDIR}$(includedir)/fplll/nr"; fi
	(test -d "${DESTDIR}$(includedir)/fplll" && rmdir "${DESTDIR}$(includedir)/fplll") || true
	rm -f "${DESTDIR}$(includedir)/fplll.h"

dbg: fplll_dbg latticegen_dbg

