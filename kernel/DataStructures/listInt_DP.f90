module listInt_DP

!$ use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer

#define D          DP
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(DP)

#include "kernel/DataStructures/list.h"

end module
