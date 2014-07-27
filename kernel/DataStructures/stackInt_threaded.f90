module stackInt_threaded

!$ use omp_lib
  use fsystem
  use genoutput
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer
#define T_THREAD_SAFE

#include "kernel/DataStructures/stack.h"

end module
