module octreeSngl

!$use omp_lib
  use fsystem
  use genoutput
  use octreebase
  use storage

#define T          Sngl
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#include "octree.h"

end module
