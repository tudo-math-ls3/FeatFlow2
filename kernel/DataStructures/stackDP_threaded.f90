module stackDP_threaded

!$ use omp_lib
  use fsystem
  use genoutput
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)
#define T_THREAD_SAFE

#include "kernel/DataStructures/stack.h"

end module
