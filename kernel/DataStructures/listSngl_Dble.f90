module listSngl_Dble

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          Sngl
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#define D          Dble
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(DP)

#include "list.h"

end module
