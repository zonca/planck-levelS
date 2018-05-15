PKG:=pntsrcconvolver

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libPointSourceConvolver.a

OBJ:=BeamCharacteristic.o PointSourceConvolver.o Planet.o
OBJ:=$(OBJ:%=$(OD)/%)

$(OBJ): $(HDR_$(PKG)) $(HDR_cxxmod) $(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libfftpack) $(HDR_c_utils)
$(OBJ): | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
