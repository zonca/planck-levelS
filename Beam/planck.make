PKG:=Beam

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE_F+= $(F_MODPATH)$(OD)

LIB_$(PKG):=$(LIBDIR)/libbeam.a
LIBOBJ:=beam_grid.o beam_cut.o beam_polar.o beam_square.o beam_bessel.o beam_alm.o beam_convert.o beam_interpolate.o beam_crosspol.o beam_transform.o beam.o beam2alm_main.o grasp2stokes_main.o crosspol_main.o gaussbeampol_main.o
BINOBJ:=grasp2stokes.o crosspol.o beam2alm.o gaussbeampol.o stokes_extract.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
BINOBJ:=$(BINOBJ:%=$(OD)/%)
ALLOBJ:=$(LIBOBJ) $(BINOBJ)
F90BIN:=grasp2stokes crosspol beam2alm gaussbeampol stokes_extract
MODULES:=beam2alm_main grasp2stokes_main crosspol_main gaussbeampol_main
MODULES:=$(MODULES:%=$(OD)/%.$(MOD))

ODEP:=$(LIB_Modules) $(LIB_libsharp_f) $(LIB_portability)
$(LIB_$(PKG)): $(LIBOBJ)
$(ALLOBJ): $(ODEP) | $(OD)_mkdir

BDEP:=$(LIB_$(PKG)) $(LIB_Modules) $(LIB_libsharp_f) $(LIB_portability) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)
$(BINOBJ): $(ODEP) $(LIB_$(PKG))

$(OD)/beam_alm.o: $(OD)/beam_bessel.o
$(OD)/beam_convert.o: $(OD)/beam_grid.o $(OD)/beam_cut.o $(OD)/beam_polar.o $(OD)/beam_square.o
$(OD)/beam_interpolate.o: $(OD)/beam_polar.o $(OD)/beam_square.o
$(OD)/beam_crosspol.o: $(OD)/beam_polar.o $(OD)/beam_square.o
$(OD)/beam_transform.o: $(OD)/beam_polar.o $(OD)/beam_alm.o
$(OD)/beam.o: $(OD)/beam_grid.o $(OD)/beam_cut.o $(OD)/beam_polar.o $(OD)/beam_square.o $(OD)/beam_bessel.o $(OD)/beam_alm.o $(OD)/beam_convert.o $(OD)/beam_interpolate.o $(OD)/beam_crosspol.o $(OD)/beam_transform.o
$(OD)/grasp2stokes_main.o $(OD)/beam2alm_main.o $(OD)/crosspol_main.o $(OD)/gaussbeampol_main.o: $(OD)/beam.o

F90BIN:=$(F90BIN:%=$(BINDIR)/%)
$(F90BIN): $(BINDIR)/%: $(OD)/%.o $(BDEP)

bin_$(PKG): $(F90BIN)

all_lib+=$(LIB_$(PKG))
all_f90bin+=$(F90BIN)
all_mod+=$(MODULES)

levels_bin+= $(F90BIN)
