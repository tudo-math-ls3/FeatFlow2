module arraylistInt_SP

!$ use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer

#define D          SP
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(SP)

#include "kernel/DataStructures/arraylist.h"

end module
