module lsyssc_scalarMatVecAux

  use fsystem
  use linearalgebra
  use perfconfig

  implicit none

contains

  ! =-=-=-=-=- Double-valued matrix, double-valued vector -=-=-=-=-=

#define MatDT double
#define VecDT double
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

  ! =-=-=-=-=- Double-valued matrix, single-valued vector -=-=-=-=-=

#define MatDT double
#define VecDT single
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

  ! =-=-=-=-=- Single-valued matrix, double-valued vector -=-=-=-=-=

#define MatDT single
#define VecDT double
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

  ! =-=-=-=-=- Single-valued matrix, single-valued vector -=-=-=-=-=

#define MatDT single
#define VecDT single
#include "lsyssc_scalarMatVecAux.h"
#undef MatDT
#undef VecDT

end module lsyssc_scalarMatVecAux
