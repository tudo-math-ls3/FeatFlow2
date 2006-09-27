#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk
include ./Globals.mk

LIBNAME=hdf5
LIBVERSION=1.6.5

LIB=$(LIBNAME:%=$(LIBDIR)/lib%*.a)
PWD=$(shell pwd)

all: greet lib
	@echo "Done," $(LIBNAME) "is ready."
	@echo

greet:
	@echo "Compiling library" $(LIBNAME) "in"
	@pwd

lib: prepare configure $(LIB) 

prepare:
ifeq ($(shell test ! -d ./$(ID) && echo 1),1)
	@mkdir -p ./$(ID)
	@gunzip < $(LIBNAME)-$(LIBVERSION).tar.gz | tar -C ./$(ID) -xf -
endif

configure: prepare
ifeq ($(shell test ! -f ./$(ID)/$(LIBNAME)-$(LIBVERSION)/config.status && echo 1),1)
	@( cd ./$(ID)/$(LIBNAME)-$(LIBVERSION) &&\
        env FFLAGS="$(FFLAGS)" env CFLAGS="$(CFLAGS)" env CC=$(CC) \
	env F9X=$(FC) ./configure --enable-fortran --disable-cxx \
	--disable-shared --disable-fast-install --disable-debug --disable-parallel \
	--with-zlib=$(PWD)/../zlib/include,$(PWD)/$(LIBDIR) \
	--with-szlib=$(PWD)/../sz/include,$(PWD)/$(LIBDIR) )
endif

libhdf5: configure
ifeq ($(shell test ! -f ./$(ID)/$(LIBNAME)-$(LIBVERSION)/src/.libs/libhdf5.a && echo 1),1)
	@( cd ./$(ID)/$(LIBNAME)-$(LIBVERSION)/src &&\
	make lib )
endif

libhdf5_hl: configure
ifeq ($(shell test ! -f ./$(ID)/$(LIBNAME)-$(LIBVERSION)/hl/src/.libs/libhdf5_hl.a && echo 1),1)
	@( cd ./$(ID)/$(LIBNAME)-$(LIBVERSION)/hl/src &&\
	make lib )
endif

libhdf5_fortran: configure
ifeq ($(shell test ! -f ./$(ID)/$(LIBNAME)-$(LIBVERSION)/fortran/src/.libs/libhdf5_fortran.a && echo 1),1)
	@( cd ./$(ID)/$(LIBNAME)-$(LIBVERSION)/fortran/src &&\
	make lib )
endif

$(LIB): libhdf5 libhdf5_hl libhdf5_fortran
	@mkdir -p $(LIBDIR)
	@cp $(PWD)/$(ID)/$(LIBNAME)-$(LIBVERSION)/src/.libs/$(LIBNAME:%=lib%*.a) $(LIBDIR)
	@cp $(PWD)/$(ID)/$(LIBNAME)-$(LIBVERSION)/hl/src/.libs/$(LIBNAME:%=lib%*.a) $(LIBDIR)
	@cp $(PWD)/$(ID)/$(LIBNAME)-$(LIBVERSION)/fortran/src/.libs/$(LIBNAME:%=lib%*.a) $(LIBDIR)
	@cp $(PWD)/$(ID)/$(LIBNAME)-$(LIBVERSION)/fortran/src/*.mod $(LIBDIR)

clean:
	-rm -r -f ./$(ID)

purge: clean
	-rm -f $(LIB)
	-rm -f $(LIBDIR)/h5*.mod

purge_all: clean
	-rm -f $(FEATFLOW)/object/libraries/lib-*/$(LIBNAME:%=lib%*.a)

id: .id
help: .help
debug: all
.PHONY : prepare configure build

