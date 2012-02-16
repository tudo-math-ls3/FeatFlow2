module mapDble_Dble

!$use omp_lib
  use mapbase
  use fsystem
  use genoutput
  use storage

#define T          Dble
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#define D          Dble
#define D_STORAGE  ST_DOUBLE
#define D_TYPE     real(DP)

#include "map.h"

end module
