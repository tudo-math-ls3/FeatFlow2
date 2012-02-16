module listDble_Int

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          Dble
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#define D          Int
#define D_STORAGE  ST_INT
#define D_TYPE     integer

#include "list.h"

end module
