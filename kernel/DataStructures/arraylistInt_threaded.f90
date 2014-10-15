module arraylistInt

!$ use omp_lib
  use arraylistbase
  use fsystem
  use genoutput
  use storage

#define T          Int
#define T_STORAGE  ST_INT
#define T_TYPE     integer
#define T_THREAD_SAFE

#include "kernel/DataStructures/arraylist.h"

end module
