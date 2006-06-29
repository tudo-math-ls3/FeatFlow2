#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

LIBNAME=hdf5
LIB=$(LIBNAME:%=$(LIBDIR)/lib%*.a)

PREFIX=$(shell pwd)/tmpinstall

all: greet lib clean
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

lib:
	@PWD=$(shell pwd)
	@( cd orig && env CC=$(CC) env F9X=$(FC) \
	  configure --enable-fortran --disable-cxx --disable-shared \
	  --disable-fast-install --disable-debug --disable-parallel \
	  --with-zlib=$(PWD)/../zlib/include,$(PWD)/$(LIBDIR) \
	  --with-szlib=$(PWD)/../sz/include,$(PWD)/$(LIBDIR) \
	  --prefix=$(PREFIX) && \
	  make && \
	  make install )
	mkdir -p $(LIBDIR)
	cp $(PREFIX)/lib/libhdf5.a $(LIBDIR)/libhdf5.a
	cp $(PREFIX)/lib/libhdf5_fortran.a $(LIBDIR)/libhdf5_fortran.a
	cp $(PREFIX)/lib/libhdf5_hl.a $(LIBDIR)/libhdf5_hl.a

clean:
	( cd orig && make clean )
	-rm -r -f $(PREFIX)

purge: clean
	-rm -f $(LIB)

purge_all: clean
	-rm -f $(FEATFLOW)/object/libraries/lib-*/$(LIBNAME:%=lib%*.a)

id: .id
help: .help
debug: all
