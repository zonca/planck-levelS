PKG:=CAMB

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE_F+= $(F_MODPATH)$(OD)

LIB_$(PKG):=$(LIBDIR)/libcamb.a
LIBOBJ:=constants.o reionization.o subroutines.o inifile.o power_tilt.o recfast.o modules.o bessels.o equations.o lensing.o cmbmain.o camb.o writefits.o inidriver.o halofit.o utils.o reionization.o SeparableBispectrum.o
ALLOBJ:=$(LIBOBJ) camb_prog.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)
F90BIN:=camb
MODULES:=camb_main
MODULES:=$(MODULES:%=$(OD)/%.$(MOD))

ODEP:=$(LIB_Modules) $(LIB_portability)
$(LIB_$(PKG)): $(LIBOBJ)
$(ALLOBJ): $(ODEP) | $(OD)_mkdir

BDEP:=$(LIB_$(PKG)) $(LIB_Modules) $(LIB_portability) $(LIB_libcfitsio)

$(OD)/subroutines.o: $(OD)/utils.o
$(OD)/power_tilt.o: $(OD)/subroutines.o $(OD)/inifile.o
$(OD)/recfast.o: $(OD)/inifile.o $(OD)/constants.o $(OD)/subroutines.o
$(OD)/modules.o: $(OD)/constants.o $(OD)/subroutines.o $(OD)/recfast.o $(OD)/power_tilt.o $(OD)/inifile.o $(OD)/reionization.o
$(OD)/bessels.o: $(OD)/subroutines.o $(OD)/modules.o
$(OD)/SeparableBispectrum.o: $(OD)/bessels.o $(OD)/inifile.o $(OD)/modules.o $(OD)/lensing.o $(OD)/power_tilt.o
$(OD)/equations.o: $(OD)/subroutines.o $(OD)/modules.o
$(OD)/lensing.o: $(OD)/subroutines.o $(OD)/modules.o $(OD)/utils.o
$(OD)/cmbmain.o: $(OD)/subroutines.o $(OD)/bessels.o $(OD)/modules.o $(OD)/equations.o $(OD)/halofit.o $(OD)/SeparableBispectrum.o
$(OD)/camb.o: $(OD)/constants.o $(OD)/subroutines.o $(OD)/modules.o $(OD)/cmbmain.o $(OD)/lensing.o $(OD)/SeparableBispectrum.o
$(OD)/writefits.o: $(OD)/camb.o
$(OD)/inidriver.o: $(OD)/inifile.o $(OD)/camb.o $(OD)/modules.o $(OD)/SeparableBispectrum.o
$(OD)/halofit.o: $(OD)/modules.o $(OD)/equations.o
$(OD)/reionization.o: $(OD)/constants.o $(OD)/subroutines.o $(OD)/recfast.o $(OD)/inifile.o
$(OD)/camb_prog.o: $(OD)/inidriver.o $(OD)/SeparableBispectrum.o

F90BIN:=$(F90BIN:%=$(BINDIR)/%)
$(BINDIR)/camb: $(OD)/camb_prog.o $(BDEP)

all_lib+=$(LIB_$(PKG))
all_f90bin+=$(F90BIN)
all_mod+=$(MODULES)

levels_bin+=$(F90BIN)
