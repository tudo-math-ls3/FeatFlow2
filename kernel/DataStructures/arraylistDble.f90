module arraylistDble

!$use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          Dble
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#include "arraylist.h"

end module
