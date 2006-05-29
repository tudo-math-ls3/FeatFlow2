#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

LIBNAME=amd

INCDIR= -I./include -I./src

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
        amd_order amd_control amd_info amd_valid amd_preprocess

UMF4FLAGS ?= -DINT -DLP64  #fall back definition if not defined

ifdef ARC
AR=$(ARC)
endif

SRCLIST= $(addsuffix .c, $(AMD))
DEFS= $(UMF4FLAGS)

include $(FEATFLOW)/Rules_libs.mk

update:
	@echo "Based on UMFPACKv4.3."
	@echo "Check http://www.cise.ufl.edu/research/sparse/umfpack/ for updates."
	@(test -e ./UMFPACKv4.3.tar.gz || wget http://www.cise.ufl.edu/research/sparse/umfpack/v4.3/UMFPACKv4.3.tar.gz )
	@tar -xzf UMFPACKv4.3.tar.gz UMFPACKv4.3/AMD
	@rm -rf orig
	@mv ./UMFPACKv4.3/AMD ./orig
	@rm -r ./UMFPACKv4.3
	@cp ./orig/Source/* ./src/
	@cp ./orig/Include/* ./include/
