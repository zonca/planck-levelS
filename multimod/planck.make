PKG:=multimod

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libmultimod.a
CXXBIN:=multimod gapmaker
LIBOBJ:=writers.o multimod_module.o
ALLOBJ:=$(LIBOBJ) $(CXXBIN:%=%.o)
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_pntsrcconvolver) $(HDR_cxxmod) $(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libsharp) $(HDR_libfftpack) $(HDR_c_utils)
ifeq ($(MIXBIN_SUPPORTED),yes)
  BDEP:=$(LIB_$(PKG)) $(LIB_pntsrcconvolver) $(LIB_cxxmod) $(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_zlelib) $(LIB_portability) $(LIB_libcfitsio)
else
  BDEP:=$(LIB_$(PKG)) $(LIB_pntsrcconvolver) $(LIB_cxxmod) $(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)
endif
$(LIB_$(PKG)): $(LIBOBJ)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)
$(CXXBIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
ifeq ($(MIXBIN_SUPPORTED),yes)
  all_mixbin+=$(CXXBIN)
else
  all_cxxbin+=$(CXXBIN)
endif

levels_bin+=multimod gapmaker
