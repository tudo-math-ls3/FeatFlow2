module lsyssc_scalarMatVecAux

  use fsystem
  use linearalgebra
  use perfconfig

  implicit none

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
