PKG:=cxxtools

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

CXXBIN:=addMaps makeMap HPXconvert_cxx rotmap_cxx dmc_transfer beamsampler pntsrc2alm addTOI fpdbhelper focalplane ddlconvert sat2quat pointing_errors ringset2map pntsrc2map read_horizons object_inspector ephemeris_helper text2powspec rotalm_phi earl_importer ringset2mapset
ALLOBJ:=$(CXXBIN:%=$(OD)/%.o)

ODEP:=$(HDR_pntsrcconvolver) $(HDR_cxxmod) $(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libsharp) $(HDR_libfftpack) $(HDR_c_utils)
BDEP:=$(LIB_pntsrcconvolver) $(LIB_cxxmod) $(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)
$(CXXBIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_cxxbin+=$(CXXBIN)

levels_bin+=$(CXXBIN)
