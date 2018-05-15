PKG:=Modules

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE_F+= $(F_MODPATH)$(OD)

LIB_$(PKG):=$(LIBDIR)/libmod.a $(LIBDIR)/libminihealpix.a
MODOBJ:=planck_types.o ls_misc_utils.o linebufmod.o fitsmod2.o ls_paramfile_io.o focalplanemod.o dmc_io.o announce_mod.o ephemerides_f.o instrument_table.o
MHPOBJ:=minihealpix.o
MODOBJ:=$(MODOBJ:%=$(OD)/%)
MHPOBJ:=$(MHPOBJ:%=$(OD)/%)
MODULES:= planck_types ls_misc_utils linebufmod fitsmod2 ls_paramfile_io focalplanemod minihealpix dmc_io fts_io announce_mod ephemerides_f
MODULES:=$(MODULES:%=$(OD)/%.$(MOD))

$(MODOBJ) $(MHPOBJ): $(LIB_portability) | $(OD)_mkdir
$(LIBDIR)/libmod.a: $(MODOBJ)
$(LIBDIR)/libminihealpix.a: $(MHPOBJ)

$(OD)/fitsmod2.o: $(OD)/planck_types.o $(OD)/linebufmod.o $(OD)/ls_misc_utils.o
$(OD)/ls_paramfile_io.o: $(OD)/ls_misc_utils.o
$(OD)/focalplanemod.o: $(OD)/instrument_table.o $(OD)/ls_misc_utils.o $(OD)/dmc_io.o
$(OD)/instrument_table.o: $(OD)/ls_misc_utils.o $(OD)/dmc_io.o
$(OD)/ephemerides_f.o: $(OD)/ls_misc_utils.o $(OD)/dmc_io.o
$(OD)/announce_mod.o: $(OD)/ls_misc_utils.o $(OD)/dmc_io.o
$(OD)/dmc_io.o: $(OD)/planck_types.o $(OD)/ls_misc_utils.o $(OD)/ls_paramfile_io.o $(OD)/fitsmod2.o $(SD)/ddl.f90
$(OD)/minihealpix.o: $(OD)/ls_misc_utils.o

all_lib+=$(LIB_$(PKG))
all_mod+=$(MODULES)
