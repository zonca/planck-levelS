PKG:=libcfitsio

ifeq ($(EXTERNAL_CFITSIO),yes)

F_EXTRALIBS+=$(CFITSIO_EXT_LIB)
CXX_EXTRALIBS+=$(CFITSIO_EXT_LIB)
FULL_INCLUDE+= $(CFITSIO_EXT_INC)

else

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=fitsio.h longnam.h
HDR_$(PKG):=$(HDR_$(PKG):%=$(SD)/%)
LIB_$(PKG):=$(LIBDIR)/libcfitsio.a

$(LIB_$(PKG)): $(SD)/* | $(OD)_mkdir $(LIBDIR)_mkdir
	cd $(BLDROOT)/libcfitsio/ && \
	rm -rf * && \
	MAKE="$(MAKE)" FC="$(FC)" CC="$(CC)" CFLAGS="$(CCFLAGS_NO_C)" \
	$(SRCROOT)/libcfitsio/configure && \
	echo VPATH=$(SRCROOT)/libcfitsio > Makefile2 && \
	cat Makefile >> Makefile2 && \
	mv Makefile2 Makefile && \
	$(MAKE) && \
	cp -p libcfitsio.a $(LIBDIR) && \
	$(MAKE) clean

all_hdr+=$(HDR_$(PKG))

endif
