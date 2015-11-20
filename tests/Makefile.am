AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/../src/
TOPBUILDDIR = $(builddir)/../src/
DEFINES =

include_fpllldir=$(includedir)/fplll

EXTRA_DIST=lattices/dim55_in lattices/examples_in

# include TOPBUILDIR for fplll_config.h
AM_CPPFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR)

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lfplll
LDFLAGS = -no-install

TESTS = test_lll

test_lll_SOURCES = test_lll.cpp

check_PROGRAMS = $(TESTS)
