module stackInt

!$use omp_lib
  use fsystem
  use genoutput
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer

#include "stack.h"

end module
