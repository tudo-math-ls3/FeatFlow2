module listSP_SP

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#define D          SP
#define D_STORAGE  ST_SINGLE
#define D_TYPE     real(SP)

#include "kernel/DataStructures/list.h"

end module
