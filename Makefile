ifndef LEVELS_TARGET
  LEVELS_TARGET:=$(error LEVELS_TARGET undefined. Please see README.compilation for help)UNDEFINED
endif

default: compile_levels
SRCROOT:=$(shell pwd)
include $(SRCROOT)/config/config.$(LEVELS_TARGET)
include $(SRCROOT)/config/rules.common

.SILENT:

all_hdr:=
all_mod:=
all_lib:=
all_cbin:=
all_cxxbin:=
all_f90bin:=
all_mixbin:=
levels_bin:=

FULL_INCLUDE:=
FULL_INCLUDE_F:=

include libcfitsio/planck.make
include c_utils/planck.make
include libfftpack/planck.make
include libsharp/planck.make
include portability/planck.make
include Modules/planck.make
include simmission3/planck.make
include libsharp_f/planck.make
include CAMB/planck.make
include Beam/planck.make
include cxxsupport/planck.make
include Healpix_cxx/planck.make
include skymixer2/planck.make
include zlelib/planck.make
include cxxmod/planck.make
include pntsrcconvolver/planck.make
include multimod/planck.make
include cxxtools/planck.make
include LFI-specific/planck.make
include conviqt/planck.make
include docsrc/planck.make

$(all_lib): %: | $(LIBDIR)_mkdir dump_config
	@echo "#  creating library $*"
	$(ARCREATE) $@ $^

$(all_cxxbin): %: | $(BINDIR)_mkdir dump_config
	@echo "#  linking C++ binary $*"
	$(CXXL) $(CXXLFLAGS) -o $@ $^ $(CXX_EXTRALIBS)

$(all_cbin): %: | $(BINDIR)_mkdir dump_config
	@echo "#  linking C binary $*"
	$(CXXL) $(CXXLFLAGS) -o $@ $^ $(CXX_EXTRALIBS)

$(all_f90bin): %: | $(BINDIR)_mkdir dump_config
	@echo "#  linking F90 binary $*"
	$(FL) $(FLFLAGS) -o $@ $^ $(F_EXTRALIBS)

$(all_mixbin): %: | $(BINDIR)_mkdir dump_config
	@echo "#  linking mixed C++/F90 binary $*"
	$(CXXL) $(CXXLFLAGS) -o $@ $^ $(CXX_EXTRALIBS) $(F_EXTRALIBS) $(MIX_EXTRALIBS)

compile_levels: $(levels_bin) hdrmodcopy

compile_all: $(all_f90bin) $(all_cbin) $(all_cxxbin) $(all_mixbin) hdrmodcopy

compile_f90: $(all_f90bin)

hdrclean:
	@if [ -d $(INCDIR) ]; then rm -rf $(INCDIR)/* ; fi

hdrmodcopy: $(all_f90bin) | hdrclean $(INCDIR)_mkdir
	@if [ "$(all_hdr)" ]; then cp -p $(all_hdr) $(INCDIR); fi
	@if [ "$(all_mod)" ]; then cp -p $(all_mod) $(INCDIR); fi

hdrcopy: | $(INCDIR)_mkdir
	@if [ "$(all_hdr)" ]; then cp -p $(all_hdr) $(INCDIR); fi

modcopy: $(all_f90bin) | $(INCDIR)_mkdir
	@if [ "$(all_mod)" ]; then cp -p $(all_mod) $(INCDIR); fi

$(notdir $(all_f90bin) $(all_cbin) $(all_cxxbin) $(all_mixbin)) : % : $(BINDIR)/%
$(notdir $(all_lib)) libcfitsio.a : % : $(LIBDIR)/%

ddl_rebuild: ddlconvert
	$(BINDIR)/ddlconvert cxxtools/ddl.txt cxxsupport/ddl.cc Modules/ddl.f90

cxxlibs: $(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio) hdrcopy
