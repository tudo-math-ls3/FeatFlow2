/* This header file is meant to load the shared library
 * libmeshadapt.so into Matlab and to call the functions and
 * subroutines from the Featflow2 kernel and the mesh adaptation
 * applications from Matlab.
 *
 * NOTE: This header file only provides access to a limited number of
 *       functions and subroutine and sometimes omits optional
 *       arguments or makes them mandatory.
 *
 *       This header file is adjusted to the calling conventions of
 *       the Intel compiler suite. That is, the full name of a
 *       subroutine/functions implemented in a module file reads:
 *
 *       modfile_mp_subroutine_name_()
 */

#include "types.h"

struct t_meshAdapt {
  struct t_boundary* rboundary;
  struct t_triangulation rtriangulation;
  struct t_hadapt rhadapt;
  int nrefmax;
  double dreftol;
  double dcrstol;
};

#ifdef __cplusplus
extern "C" {
#endif
  // fsystem.f90
  void fsystem_mp_sys_init_simple_();
  
  // genoutput.f90
  void genoutput_mp_output_init_simple_();
  
  // storage.f90
  void storage_mp_storage_done_(void*);
  void storage_mp_storage_init_(int*, int*, void*);

  // meshadaptbase.f90
  void meshadaptbase_mp_madapt_alloc_(struct t_meshAdapt**);
  void meshadaptbase_mp_madapt_dealloc_(struct t_meshAdapt**);
  void meshadaptbase_mp_madapt_init_(struct t_meshAdapt*, int*, char*, int);
  void meshadaptbase_mp_madapt_step_fromfile_(struct t_meshAdapt*, char*, int*, double*, double*, int*, int);
  void meshadaptbase_mp_madapt_step_sngl2_(struct t_meshAdapt*, int*, float*, int*, float*, float*, int*);
  void meshadaptbase_mp_madapt_step_dble2_(struct t_meshAdapt*, int*, double*, int*, double*, double*, int*);
  void meshadaptbase_mp_madapt_done_(struct t_meshAdapt*);
  int  meshadaptbase_mp_madapt_getnel_(struct t_meshAdapt*, int*);
  int  meshadaptbase_mp_madapt_getnvt_(struct t_meshAdapt*, int*);
  int  meshadaptbase_mp_madapt_getndim_(struct t_meshAdapt*, int*);
  int  meshadaptbase_mp_madapt_getnnve_(struct t_meshAdapt*, int*);
  void meshadaptbase_mp_madapt_getvertexcoords_(struct t_meshAdapt*, double*, int*);
  void meshadaptbase_mp_madapt_getverticesatelement_(struct t_meshAdapt*, int*, int*);
  void meshadaptbase_mp_madapt_getneighboursatelement_(struct t_meshAdapt*, int*, int*);
  void meshadaptbase_mp_madapt_getvertexage(struct t_meshAdapt*, int*);
  void meshadaptbase_mp_madapt_getelementage(struct t_meshAdapt*, int*, int*);

#ifdef __cplusplus
} /* end extern "C" */
#endif
