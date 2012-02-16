module mapInt_Int

!$use omp_lib
  use mapbase
  use fsystem
  use genoutput
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer

#define D          Int
#define D_STORAGE  ST_INT
#define D_TYPE     integer

#include "map.h"

end module
