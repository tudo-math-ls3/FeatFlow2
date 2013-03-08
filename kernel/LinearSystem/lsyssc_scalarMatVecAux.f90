module lsyssc_scalarMatVecAux

  use fsystem
  use linearalgebra
  use perfconfig

  implicit none

  private
  public :: lsyssc_lax79DPDP
  public :: lsyssc_lax79SPDP
  public :: lsyssc_lax79DPSP
  public :: lsyssc_lax79SPSP
  public :: lsyssc_lax9rowcDPDP
  public :: lsyssc_lax9rowcSPDP
  public :: lsyssc_lax9rowcDPSP
  public :: lsyssc_lax9rowcSPSP
  public :: lsyssc_lax79intl1DPDP
  public :: lsyssc_lax79intl1SPDP
  public :: lsyssc_lax79intl1DPSP
  public :: lsyssc_lax79intl1SPSP
  public :: lsyssc_lax79intldDPDP
  public :: lsyssc_lax79intldSPDP
  public :: lsyssc_lax79intldDPSP
  public :: lsyssc_lax79intldSPSP
  public :: lsyssc_latxdDPDP
  public :: lsyssc_latxdSPDP
  public :: lsyssc_latxdDPSP
  public :: lsyssc_latxdSPSP
  public :: lsyssc_ltx79DPDP
  public :: lsyssc_ltx79SPDP
  public :: lsyssc_ltx79DPSP
  public :: lsyssc_ltx79SPSP

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
