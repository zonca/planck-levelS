PKG:=portability

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE_F+= $(F_MODPATH)$(OD)

LIB_$(PKG):=$(LIBDIR)/libportability.a
OBJ:=planck_config.o
OBJ:=$(OBJ:%=$(OD)/%)

$(OBJ): | $(OD)_mkdir
MODULES:=planck_config
MODULES:=$(MODULES:%=$(OD)/%.$(MOD))

$(LIB_$(PKG)): $(OBJ)

all_lib+=$(LIB_$(PKG))
all_mod+=$(MODULES)
