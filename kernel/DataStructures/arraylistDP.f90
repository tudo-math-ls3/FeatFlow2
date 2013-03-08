module arraylistDP

!$use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          DP
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#include "arraylist.h"

end module
