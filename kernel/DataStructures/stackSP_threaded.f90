module stackSP_threaded

!$ use omp_lib
  use fsystem
  use genoutput
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)
#define T_THREAD_SAFE

#include "kernel/DataStructures/stack.h"

end module
