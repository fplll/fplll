AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/..
TOPBUILDDIR = $(builddir)/../fplll
DEFINES =

include_fpllldir=$(includedir)/fplll

EXTRA_DIST=lattices/dim55_in \
	lattices/example_in \
	lattices/example_svp_in \
	lattices/example_svp_out \
	lattices/example_dsvp_in \
	lattices/example_dsvp_out

# include TOPBUILDIR for fplll_config.h
AM_CPPFLAGS = -I$(TOPSRCDIR) -I$(TOPSRCDIR)/fplll -I$(TOPBUILDDIR) -DTESTDATADIR=\"$(TOPSRCDIR)/\"

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lfplll -no-install $(LIBQD_LIBADD)

TESTS = test_nr test_lll test_svp test_bkz test_pruner test_sieve

test_pruner_LDADD=-lgmp -lmpfr $(LIBQD_LIBADD)
test_sieve_LDADD=-lgmp -lmpfr $(LIBQD_LIBADD)

test_nr_SOURCES = test_nr.cpp
test_lll_SOURCES = test_lll.cpp
test_svp_SOURCES = test_svp.cpp
test_bkz_SOURCES = test_bkz.cpp
test_pruner_SOURCES = test_pruner.cpp
test_sieve_SOURCES = test_sieve.cpp

check_PROGRAMS = $(TESTS)
