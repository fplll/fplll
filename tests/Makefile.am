AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/../fplll
TOPBUILDDIR = $(builddir)/../fplll/
DEFINES =

include_fpllldir=$(includedir)/fplll

EXTRA_DIST=lattices/dim55_in lattices/example_in lattices/example_svp_in lattices/example_svp_out

# include TOPBUILDIR for fplll_config.h
AM_CPPFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR)

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lfplll
LDFLAGS = -no-install
LIBADD = $(LIBQD_LIBADD)

TESTS = test_lll test_svp test_bkz test_nr

test_lll_SOURCES = test_lll.cpp
test_svp_SOURCES = test_svp.cpp
test_bkz_SOURCES = test_bkz.cpp
test_nr_SOURCES = test_nr.cpp

check_PROGRAMS = $(TESTS)
