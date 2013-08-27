module quadtreeDP

!$use omp_lib
  use fsystem
  use genoutput
  use quadtreebase
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#include "kernel/DataStructures/quadtree.h"

end module
