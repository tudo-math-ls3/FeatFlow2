module arraylistSP

!$use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          SP
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#include "arraylist.h"

end module
