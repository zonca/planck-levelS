PKG:=libsharp_f

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE_F+= $(F_MODPATH)$(OD)

LIB_$(PKG):=$(LIBDIR)/libsharp_f.a
OBJ:=sharp_f.o sharpf_mod.o
OBJ:=$(OBJ:%=$(OD)/%)

$(OBJ): $(SD)/sharp_f_inc.c

ODEP:=$(LIB_portability) $(HDR_c_utils) $(HDR_libsharp)

$(OBJ): $(ODEP) | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ)

all_lib+=$(LIB_$(PKG))
