PKG:=cxxmod

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libcxxmod.a

OBJ:=detector_pointing.o interpolator.o interpolator2.o sat_info.o sampler.o dipole.o oofnoise.o satquat.o zle.o repointing.o wobble_correction.o ptcor.o
OBJ:=$(OBJ:%=$(OD)/%)

$(OBJ): $(HDR_$(PKG)) $(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libfftpack) $(HDR_c_utils)
$(OD)/zle.o: $(HDR_zlelib)
$(OBJ): | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
