module arraylistDble_Sngl

!$use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          Dble
#define T_STORAGE  ST_DOUBLE
#define T_TYPE     real(DP)

#define D          Sngl
#define D_STORAGE  ST_SINGLE
#define D_TYPE     real(SP)

#include "arraylist.h"

end module
