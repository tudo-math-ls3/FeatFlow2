/* This header file is meant to load the shared library
 * libmeshadapt.so into Matlab and to call the functions and
 * subroutines from the Featflow2 kernel and the mesh adaptation
 * applications from Matlab.
 *
 * NOTE: This header file only provides access to a limited number of
 *       functions and subroutine and sometimes omits optional
 *       arguments or makes them mandatory.
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
  void __meshadaptbase_MOD_madapt_init(struct t_meshAdapt*, int*, char*);
  void __meshadaptbase_MOD_madapt_done(struct t_meshAdapt*);
  
#ifdef __cplusplus
} /* end extern "C" */
#endif
