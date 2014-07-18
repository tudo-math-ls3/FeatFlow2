module mapDP_SP

!$ use omp_lib
  use mapbase
  use fsystem
  use genoutput
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#define D          SP
#define D_STORAGE  ST_SINGLE
#define D_TYPE     real(SP)

#include "kernel/DataStructures/map.h"

end module
