PKG:=zlelib

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)
FULL_INCLUDE_F+= $(F_MODPATH)$(SD)

LIB_$(PKG):=$(LIBDIR)/libzlelib.a
HDR_$(PKG):=$(SD)/*.h

LIBOBJ:=zlelib_full.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)

$(LIBOBJ): $(HDR_$(PKG)) | $(OD)_mkdir
$(LIBDIR)/libzlelib.a: $(LIBOBJ)

$(OD)/zlelib_full.o: $(SD)/*.f90 $(LIB_portability)

all_lib+=$(LIB_$(PKG))
