PKG:=conviqt

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

LIB_$(PKG):=$(LIBDIR)/libconviqt_v3.a
CXXBIN:=conviqt_v3 conviqt_v4 ringset_compare
LIBOBJ:=conviqt_v3_module.o conviqt_v4_module.o
ALLOBJ:=$(LIBOBJ) conviqt_v3.o conviqt_v4.o ringset_compare.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)

ODEP:=$(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libsharp) $(HDR_libfftpack) $(HDR_c_utils)
BDEP:=$(LIB_$(PKG)) $(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)

$(LIB_$(PKG)): $(LIBOBJ)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)
$(CXXBIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_lib+=$(LIB_$(PKG))
all_cxxbin+=$(CXXBIN)

levels_bin+=$(CXXBIN)
