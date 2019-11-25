SUBDIRS            = fplll tests
include_fpllldir   = $(includedir)/fplll
EXTRA_DIST         = README.md tools/reformat.pl tools/reformat_magma.pl tools/plot_gso_dump.py strategies/default.json
ACLOCAL_AMFLAGS    = -I m4

CLEANFILES = *.gcov .libs/*.gcda .libs/*.gcno *.gcno *.gcda

pkgconfigdir       = $(libdir)/pkgconfig
pkgconfig_DATA     = fplll.pc

strategydir        = $(pkgdatadir)/strategies
dist_strategy_DATA = strategies/default.json

CLANGFORMAT       ?= clang-format
.PHONY: check-style
check-style:
	-bash .check-m4.sh
	$(CLANGFORMAT) -i --style=file fplll/*.{cpp,h} fplll/*/*.{cpp,h,inl} tests/*.cpp
