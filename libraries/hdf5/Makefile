#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

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
	@echo "Please note that this library may cause some problems on different"
	@echo "architectures and therefore is not included in the list of standard"
	@echo "libraries. Due to the fact that this library strongly depends on the."
	@echo "underlying hardware, the original configure procedure has not been"
	@echo "replaced. In contrast, a pseudo installation is performed and the"
	@echo "required libraries are copied to their Feat-style directories."

lib: | prepare configure build $(LIB)

prepare:
ifeq ($(shell test ! -d ./$(ID) && echo 1),1)
	@mkdir -p ./$(ID)
	@gunzip < $(LIBNAME)-$(LIBVERSION).tar.gz | tar -C ./$(ID) -xf -
endif

configure: prepare
ifeq ($(shell test ! -f ./$(ID)/$(LIBNAME)-$(LIBVERSION)/config.status && echo 1),1)
	@( cd ./$(ID)/$(LIBNAME)-$(LIBVERSION) &&\
	env CC=$(CC) env F9X=$(FC) ./configure --enable-fortran --disable-cxx \
	--disable-shared --disable-fast-install --disable-debug --disable-parallel \
	--with-zlib=$(PWD)/../zlib/include,$(PWD)/$(LIBDIR) \
	--with-szlib=$(PWD)/../sz/include,$(PWD)/$(LIBDIR)  \
	--prefix=$(PWD)/$(ID)/tmpinstall )
endif

build: configure
ifeq ($(shell test ! -d ./$(ID)/tmpinstall && echo 1),1)
	@( cd ./$(ID)/$(LIBNAME)-$(LIBVERSION) && $(MAKE) lib && \
	cd ./src && $(MAKE) install && \
	cd ../fortran && $(MAKE) install )
endif

$(LIB): build
	@mkdir -p $(LIBDIR)
	@cp $(PWD)/$(ID)/tmpinstall/lib/$(LIBNAME:%=lib%*.a) $(LIBDIR)
	@cp $(PWD)/$(ID)/tmpinstall/lib/*.mod $(LIBDIR)

clean:
	-rm -r -f ./$(ID)

purge: clean
	-rm -f $(LIB)

purge_all: clean
	-rm -f $(FEATFLOW)/object/libraries/lib-*/$(LIBNAME:%=lib%*.a)

id: .id
help: .help
debug: all
.PHONY : prepare configure build

