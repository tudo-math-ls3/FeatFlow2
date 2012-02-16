module listInt_Sngl

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer

#define D          Sngl
#define D_STORAGE  ST_SINGLE
#define D_TYPE     real(SP)

#include "list.h"

end module
