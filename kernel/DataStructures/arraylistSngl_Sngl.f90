module arraylistSngl_Sngl

!$use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          Sngl
#define T_STORAGE  ST_SINGLE
#define T_TYPE     real(SP)

#define D          Sngl
#define D_STORAGE  ST_SINGLE
#define D_TYPE     real(SP)

#include "arraylist.h"

end module
