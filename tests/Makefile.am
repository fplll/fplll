AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/../src/
TOPBUILDDIR = $(builddir)/../src/
DEFINES =

include_fpllldir=$(includedir)/fplll

EXTRA_DIST=lattices/dim55_in lattices/example_in

# include TOPBUILDIR for fplll_config.h
AM_CPPFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR)

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lfplll
LDFLAGS = -no-install

TESTS = test_lll test_svp

test_lll_SOURCES = test_lll.cpp
test_svp_SOURCES = test_svp.cpp

check_PROGRAMS = $(TESTS)
