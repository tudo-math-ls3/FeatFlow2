module arraylistSP_Int

!$ use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#define D          Int
#define D_STORAGE  ST_INT
#define D_TYPE     integer

#include "kernel/DataStructures/arraylist.h"

end module
