module octreeSP

!$ use omp_lib
  use fsystem
  use genoutput
  use octreebase
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#include "kernel/DataStructures/octree.h"

end module
