PGM = r.hydro.hbv

LIBES = $(GISLIB) $(MATHLIB)
EXTRA_LIBS = $(OPENMP_LIBPATH) $(OPENMP_LIB)
DEPENDENCIES = $(GISDEP)
EXTRA_CFLAGS = $(OPENMP_CFLAGS)
EXTRA_INC = $(OPENMP_INCPATH)

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd

# Module.make's ETCFILES rule only mkdir's $(ETC)/$(PGM) itself, not
# nested subdirectories, so the bundled sample datasets are installed
# with a small extra rule instead of via ETCFILES.
cmd: install-data

install-data:
	$(MKDIR) $(ETC)/$(PGM)/data/original
	$(MKDIR) $(ETC)/$(PGM)/data/dicrim
	$(INSTALL_DATA) data/original/*.csv data/original/basin_ids.txt $(ETC)/$(PGM)/data/original/
	$(INSTALL_DATA) data/dicrim/*.csv data/dicrim/basin_ids.txt $(ETC)/$(PGM)/data/dicrim/

.PHONY: install-data
