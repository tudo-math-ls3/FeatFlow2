#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

LIBNAME=umfpack4

INCDIR=-I./include -I./src -I../amd/include -I../amd/src

ifdef ARC
AR=$(ARC)
endif


UMF4= umf_2by2 umf_analyze umf_apply_order umf_assemble \
  umf_blas3_update umf_build_tuples umf_colamd \
  umf_create_element umf_dump umf_extend_front umf_free \
  umf_fsize umf_garbage_collection umf_get_memory \
  umf_grow_front umf_init_front umf_is_permutation umf_kernel \
  umf_kernel_init umf_kernel_wrapup umf_local_search \
  umf_lsolve umf_ltsolve umf_malloc umf_mem_alloc_element \
  umf_mem_alloc_head_block umf_mem_alloc_tail_block \
  umf_mem_free_tail_block umf_mem_init_memoryspace \
  umfpack_col_to_triplet umfpack_qsymbolic umfpack_tictoc \
  umfpack_defaults umf_store_lu umf_set_stats \
  umfpack_report_control umfpack_numeric umf_scale_column \
  umfpack_free_numeric umf_usolve umfpack_free_symbolic \
  umfpack_scale umfpack_get_lunz umfpack_report_numeric \
  umfpack_get_numeric umfpack_get_symbolic umfpack_report_info \
  umfpack_report_matrix umfpack_solve umfpack_report_symbolic \
  umfpack_report_perm umfpack_report_status \
  umfpack_report_triplet umfpack_report_vector \
  umf_symbolic_usage umfpack_symbolic umfpack_timer \
  umfpack_transpose umf_realloc umfpack_triplet_to_col \
  umf_report_perm umf_report_vector umf_row_search umf_scale \
  umf_singletons umf_start_front umfpack_save_numeric \
  umf_transpose umf_triplet umf_tuple_lengths umf_utsolve \
  umf_valid_numeric umf_valid_symbolic umf_solve \
  umfpack_save_symbolic umfpack_load_numeric umfpack_load_symbolic


# our additional files for easy compilation and fortran 77 wrapper
AUX= umf_multicompile_1 umf_multicompile_2 umf_multicompile_3 \
     umf_multicompile_4 umf_multicompile_5 umf_multicompile_6 \
     umf_multicompile_7 umf_multicompile_8 umf4_f77wrapper_port

# we make only the double/integer version
# no need for complex version

UMF4FLAGS ?= -DINT -DLP64  #fall back definition if not defined

SRCLIST= $(addsuffix .c, $(UMF4) $(AUX))
DEFS= $(UMF4FLAGS)

include $(FEATFLOW)/Rules_libs.mk

$(AUX:%=$(OBJDIR)/%.o): $(OBJDIR)/%.o: modsrc/%.c
$(AUX:%=src/%.c): src/%.c: modsrc/%.c
	cp $< $@

src/umf_config.h: modsrc/umf_config.h
	cp $< $@

update:
	@echo "Based on UMFPACKv4.3."
	@echo "Check http://www.cise.ufl.edu/research/sparse/umfpack/ for updates."
	@(/bin/test -e ./UMFPACKv4.3.tar.gz || wget http://www.cise.ufl.edu/research/sparse/umfpack/v4.3/UMFPACKv4.3.tar.gz )
	@tar -xvzf UMFPACKv4.3.tar.gz
	@rm -r ./orig
	@cp -r -L ./UMFPACKv4.3/UMFPACK ./orig
	@rm -r ./UMFPACKv4.3
	@cp ./orig/Source/* ./src/
	@cp ./orig/Include/* ./include/
	@cp ./modsrc/* ./src/

