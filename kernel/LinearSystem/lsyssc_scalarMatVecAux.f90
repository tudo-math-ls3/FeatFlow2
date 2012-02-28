module lsyssc_scalarMatVecAux

  use fsystem
  use linearalgebra
  use perfconfig

  implicit none

  private
  public :: lsyssc_lax79dbledble
  public :: lsyssc_lax79sngldble
  public :: lsyssc_lax79dblesngl
  public :: lsyssc_lax79snglsngl
  public :: lsyssc_lax9rowcdbledble
  public :: lsyssc_lax9rowcsngldble
  public :: lsyssc_lax9rowcdblesngl
  public :: lsyssc_lax9rowcsnglsngl
  public :: lsyssc_lax79intl1dbledble
  public :: lsyssc_lax79intl1sngldble
  public :: lsyssc_lax79intl1dblesngl
  public :: lsyssc_lax79intl1snglsngl
  public :: lsyssc_lax79intlddbledble
  public :: lsyssc_lax79intldsngldble
  public :: lsyssc_lax79intlddblesngl
  public :: lsyssc_lax79intldsnglsngl
  public :: lsyssc_latxddbledble
  public :: lsyssc_latxdsngldble
  public :: lsyssc_latxddblesngl
  public :: lsyssc_latxdsnglsngl
  public :: lsyssc_ltx79dbledble
  public :: lsyssc_ltx79sngldble
  public :: lsyssc_ltx79dblesngl
  public :: lsyssc_ltx79snglsngl
  
contains

  ! =-=-=-=-=- Double-valued matrix, double-valued vector -=-=-=-=-=

#define MatDT DOUBLE_PREC
#define VecDT DOUBLE_PREC
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

  ! =-=-=-=-=- Double-valued matrix, single-valued vector -=-=-=-=-=

#define MatDT DOUBLE_PREC
#define VecDT SINGLE_PREC
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

  ! =-=-=-=-=- Single-valued matrix, double-valued vector -=-=-=-=-=

#define MatDT SINGLE_PREC
#define VecDT DOUBLE_PREC
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

  ! =-=-=-=-=- Single-valued matrix, single-valued vector -=-=-=-=-=

#define MatDT SINGLE_PREC
#define VecDT SINGLE_PREC
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

end module lsyssc_scalarMatVecAux
