module mapSngl_Sngl

!$use omp_lib
  use mapbase
  use fsystem
  use genoutput
  use storage

#define T          Sngl
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#define D          Sngl
#define D_STORAGE  ST_SINGLE
#define D_TYPE     real(SP)

#include "map.h"

end module
