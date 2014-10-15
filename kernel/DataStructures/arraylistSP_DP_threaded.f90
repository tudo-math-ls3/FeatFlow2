module arraylistSP_DP

!$ use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)
#define T_THREAD_SAFE

#define D          DP
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(DP)

#include "kernel/DataStructures/arraylist.h"

end module
