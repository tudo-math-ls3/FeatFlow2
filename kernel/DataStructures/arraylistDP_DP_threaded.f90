module arraylistDP_DP

!$ use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)
#define T_THREAD_SAFE

#define D          DP
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(DP)

#include "kernel/DataStructures/arraylist.h"

end module
