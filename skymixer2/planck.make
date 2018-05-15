PKG:=skymixer2

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libskymixer.a
CXXBIN:=almmixer skymixer3 pntsrcmixer
LIBOBJ:=skymixer3_module.o almmixer_module.o
ALLOBJ:=$(LIBOBJ) $(CXXBIN:%=%.o)
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libsharp) $(HDR_libfftpack) $(HDR_c_utils)
BDEP:=$(LIB_$(PKG)) $(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)

$(LIB_$(PKG)): $(LIBOBJ)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
$(OD)/skymixer3_module.o: $(SD)/spectra.cc $(SD)/detector_responses.cc $(SD)/emitters.cc
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)
$(CXXBIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
all_cxxbin+=$(CXXBIN)

levels_bin+=$(CXXBIN)
