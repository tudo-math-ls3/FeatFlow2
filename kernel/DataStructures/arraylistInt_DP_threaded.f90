module arraylistInt_DP

!$ use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer
#define T_THREAD_SAFE

#define D          DP
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(DP)

#include "kernel/DataStructures/arraylist.h"

end module
