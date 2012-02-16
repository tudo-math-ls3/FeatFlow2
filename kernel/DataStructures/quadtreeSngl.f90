module quadtreeSngl

!$use omp_lib
  use fsystem
  use genoutput
  use quadtreebase
  use storage

#define T          Sngl
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#include "quadtree.h"

end module
