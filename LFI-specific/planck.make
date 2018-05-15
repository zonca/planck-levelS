PKG:=LFI-specific

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

CXXBIN:=quantum ls2lfitoi ahf2satpt
ALLOBJ:=$(CXXBIN:%=$(OD)/%.o)

ODEP:=$(HDR_cxxsupport)
BDEP:=$(LIB_cxxsupport) $(LIB_libcfitsio)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)
$(CXXBIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_cxxbin+=$(CXXBIN)

levels_bin+=$(CXXBIN)
