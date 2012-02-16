module stackDble

!$use omp_lib
  use fsystem
  use genoutput
  use storage

#define T          Dble
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#include "stack.h"

end module
