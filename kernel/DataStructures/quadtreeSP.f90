module quadtreeSP

!$use omp_lib
  use fsystem
  use genoutput
  use quadtreebase
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#include "kernel/DataStructures/quadtree.h"

end module
