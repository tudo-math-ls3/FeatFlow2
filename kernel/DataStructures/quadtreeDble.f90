module quadtreeDble

!$use omp_lib
  use fsystem
  use genoutput
  use quadtreebase
  use storage

#define T          Dble
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#include "quadtree.h"

end module
