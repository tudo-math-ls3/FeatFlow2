module stackDP

!$ use omp_lib
  use fsystem
  use genoutput
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#include "kernel/DataStructures/stack.h"

end module
