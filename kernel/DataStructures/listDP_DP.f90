module listDP_DP

!$ use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#define D          DP
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(DP)

#include "kernel/DataStructures/list.h"

end module
