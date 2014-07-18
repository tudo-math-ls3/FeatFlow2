module octreeDP

!$ use omp_lib
  use fsystem
  use genoutput
  use octreebase
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#include "kernel/DataStructures/octree.h"

end module
