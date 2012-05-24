/*#############################################################################
 ******************************************************************************
 * <name> coproc_transpose </name>
 ******************************************************************************
 *
 * <purpose>
 * This file provides routines for transposing multi-dimensional arrays.
 * </purpose>
 *
 *#############################################################################
 */

#include <stdint.h>
#include "coproc_core.h"
#include "coproc_transpose.h"

using namespace std;

/*******************************************************************************
 * Transpose two-dimensional array
 ******************************************************************************/
template<typename T>
int coproc_transpose2d(const void *ptrSrc, void *ptrDest, int n1, int n2)
{
  T *src = (T*) ptrSrc;
  T *dest = (T*) ptrDest;
  
  for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      dest[n2*i1+i2] = src[n1*i2+i1];
}
/******************************************************************************/
template<typename T,int num_elmt>
int coproc_transpose2d(const void *ptrSrc, void *ptrDest, int n1, int n2)
{
  T *src = (T*) ptrSrc;
  T *dest = (T*) ptrDest;
  
  for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
#pragma unroll
      for (int ielmt=0; ielmt<num_elmt; ++ielmt)
	dest[sizeof(T)*(n2*i1+i2)+ielmt] = src[sizeof(T)*(n1*i2+i1)+ielmt];
}

/******************************************************************************/
int coproc_transposeOnHost2d(void *ptrSrc, void *ptrDest,
			     int bytes_per_elmt, int n1, int n2)
{
  switch (bytes_per_elmt) {
  
  case sizeof(uint8_t): // 1 byte
    return coproc_transpose2d<uint8_t>(ptrSrc, ptrDest, n1, n2);

  case sizeof(uint16_t): // 2 byte
    return coproc_transpose2d<uint16_t>(ptrSrc, ptrDest, n1, n2);

  case (3):              // 3 bytes
    return coproc_transpose2d<uint8_t,3>(ptrSrc, ptrDest, n1, n2);
    
  case sizeof(uint32_t): // 4 byte
    return coproc_transpose2d<uint32_t>(ptrSrc, ptrDest, n1, n2);
    
  case (5):              // 5 bytes
    return coproc_transpose2d<uint8_t,5>(ptrSrc, ptrDest, n1, n2);
    
  case (6):              // 6 bytes
    return coproc_transpose2d<uint16_t,3>(ptrSrc, ptrDest, n1, n2);
    
  case (7):              // 7 bytes
    return coproc_transpose2d<uint8_t,7>(ptrSrc, ptrDest, n1, n2);

  case sizeof(uint64_t): // 8 byte
    return coproc_transpose2d<uint64_t>(ptrSrc, ptrDest, n1, n2);

  case (9):              // 9 bytes
    return coproc_transpose2d<uint8_t,9>(ptrSrc, ptrDest, n1, n2);

  case (10):             // 10 bytes
    return coproc_transpose2d<uint16_t,5>(ptrSrc, ptrDest, n1, n2);

  case (11):             // 11 bytes
    return coproc_transpose2d<uint8_t,11>(ptrSrc, ptrDest, n1, n2);

  case (12):              // 12 bytes
    return coproc_transpose2d<uint32_t,3>(ptrSrc, ptrDest, n1, n2);

  case (13):             // 13 bytes
    return coproc_transpose2d<uint8_t,13>(ptrSrc, ptrDest, n1, n2);

  case (14):             // 14 bytes
    return coproc_transpose2d<uint16_t,7>(ptrSrc, ptrDest, n1, n2);
    
  case (15):             // 15 bytes
    return coproc_transpose2d<uint8_t,15>(ptrSrc, ptrDest, n1, n2);

  case (16):             // 16 bytes
    return coproc_transpose2d<uint64_t,2>(ptrSrc, ptrDest, n1, n2);

  case (24):             // 24 bytes
    return coproc_transpose2d<uint64_t,3>(ptrSrc, ptrDest, n1, n2);

  case (32):             // 32 bytes
    return coproc_transpose2d<uint64_t,4>(ptrSrc, ptrDest, n1, n2);

  case (64):             // 64 bytes
    return coproc_transpose2d<uint64_t,8>(ptrSrc, ptrDest, n1, n2);

  case (128):            // 128 bytes
    return coproc_transpose2d<uint64_t,16>(ptrSrc, ptrDest, n1, n2);

  case (256):            // 256 bytes
    return coproc_transpose2d<uint64_t,32>(ptrSrc, ptrDest, n1, n2);

  case (512):            // 512 bytes
    return coproc_transpose2d<uint64_t,64>(ptrSrc, ptrDest, n1, n2);

  default:
    char *src = (char*) ptrSrc;
    char *dest = (char*) ptrDest;
    
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
	for (int ibyte=0; ibyte<bytes_per_elmt; ++ibyte)
	  dest[bytes_per_elmt*(n2*i1+i2)+ibyte] = src[bytes_per_elmt*(n1*i2+i1)+ibyte];
    return 0;
    break;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_transposeonhost2d)(void **ptrSrc, void **ptrDest,
					__INT *bytes_per_elmt, __INT *n1, __INT *n2)
  {
    return (__INT) coproc_transposeOnHost2d(*ptrSrc, *ptrDest, *bytes_per_elmt, *n1, *n2);
  }
};
    
/*******************************************************************************
 * Transpose three-dimensional array
 ******************************************************************************/
template<int bytes_per_elmt>
int coproc_transposeOnHost3d(void **ptrSrc, void **ptrDest, int n1, int n2, int n3)
{
}

/******************************************************************************/
int coproc_transposeOnHost3d(void **ptrSrc, void **ptrDest,
			     int bytes_per_elmt, int n1, int n2, int n3)
{
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_transposeonhost3d)(void **ptrSrc, void **ptrDest,
					__INT bytes_per_elmt, __INT n1, __INT n2, __INT n3)
  {
  }
};
