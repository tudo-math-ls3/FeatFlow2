module listSngl

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          Sngl
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)
#undef  T_MODULE

#include "list.h"

end module
