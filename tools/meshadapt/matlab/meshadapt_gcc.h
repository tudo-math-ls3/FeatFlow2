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
 *       the GNU compiler suite. That is, the full name of a
 *       subroutine/functions implemented in a module file reads:
 *
 *       __modfile_MOD_subroutine_name()
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
  void __fsystem_MOD_sys_init_simple();
  
  // genoutput.f90
  void __genoutput_MOD_output_init_simple();
  
  // storage.f90
  void __storage_MOD_storage_done(void*);
  void __storage_MOD_storage_init(int*, int*, void*);

  // meshadaptbase.f90
  void __meshadaptbase_MOD_madapt_alloc(struct t_meshAdapt**);
  void __meshadaptbase_MOD_madapt_dealloc(struct t_meshAdapt**);
  void __meshadaptbase_MOD_madapt_init(struct t_meshAdapt*, int*, char*, int);
  void __meshadaptbase_MOD_madapt_step_fromfile(struct t_meshAdapt*, char*, int*, double*, double*, int);
  void __meshadaptbase_MOD_madapt_step_dble2(struct t_meshAdapt*, int*, double*, int*, double*, double*);
  void __meshadaptbase_MOD_madapt_done(struct t_meshAdapt*);
  int  __meshadaptbase_MOD_madapt_getnel(struct t_meshAdapt*);
  int  __meshadaptbase_MOD_madapt_getnvt(struct t_meshAdapt*);
  int  __meshadaptbase_MOD_madapt_getndim(struct t_meshAdapt*);
  int  __meshadaptbase_MOD_madapt_getnnve(struct t_meshAdapt*);
  void __meshadaptbase_MOD_madapt_getvertexcoords(struct t_meshAdapt*, double*);
  void __meshadaptbase_MOD_madapt_getverticesatelement(struct t_meshAdapt*, int*);
  void __meshadaptbase_MOD_madapt_getneighboursatelement(struct t_meshAdapt*, int*);

#ifdef __cplusplus
} /* end extern "C" */
#endif
