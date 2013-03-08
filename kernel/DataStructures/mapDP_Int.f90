module mapDP_Int

!$use omp_lib
  use mapbase
  use fsystem
  use genoutput
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#define D          Int
#define D_STORAGE  ST_INT
#define D_TYPE     integer

#include "map.h"

end module
