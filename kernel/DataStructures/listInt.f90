module listInt

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer
#undef  T_MODULE

#include "list.h"

end module
