module listSP

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#include "list.h"

end module
