module listSngl_Int

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          Sngl
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#define D          Int
#define D_STORAGE  ST_INT
#define D_TYPE     integer

#include "list.h"

end module
