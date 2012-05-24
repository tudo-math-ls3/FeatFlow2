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
#include <iostream>
#include "coproc_core.h"
#include "coproc_transpose.h"

using namespace std;

/*******************************************************************************
 * Transpose two-dimensional array
 ******************************************************************************/
template<typename T>
int coproc_transpose2d(const void *__restrict__ ptrSrc, 
		       void *__restrict__ ptrDest,
		       int nelmt1, int nelmt2)
{
  T *src = (T*)(ptrSrc);
  T *dest = (T*)(ptrDest);
  
  for (int i2=0; i2<nelmt2; ++i2)
    for (int i1=0; i1<nelmt1; ++i1)
      dest[nelmt2*i1+i2] = src[nelmt1*i2+i1];

  return 0;
}
/******************************************************************************/
template<typename T>
int coproc_transpose2d(int num_elmt, const void *__restrict__ ptrSrc, 
		       void *__restrict__ ptrDest, int nelmt1, int nelmt2)
{
  T *src = (T*)(ptrSrc);
  T *dest = (T*)(ptrDest);
  
  for (int i2=0; i2<nelmt2; ++i2)
    for (int i1=0; i1<nelmt1; ++i1)
      for (int ielmt=0; ielmt<num_elmt; ++ielmt)
	dest[num_elmt*(nelmt2*i1+i2)+ielmt] = src[num_elmt*(nelmt1*i2+i1)+ielmt];
  
  return 0;
}
/******************************************************************************/
template<typename T,int num_elmt>
int coproc_transpose2d(const void *__restrict__ ptrSrc,
		       void *__restrict__ ptrDest, int nelmt1, int nelmt2)
{
  T *src = (T*)(ptrSrc);
  T *dest = (T*)(ptrDest);
  
  for (int i2=0; i2<nelmt2; ++i2)
    for (int i1=0; i1<nelmt1; ++i1)
#pragma unroll
      for (int ielmt=0; ielmt<num_elmt; ++ielmt)
	dest[num_elmt*(nelmt2*i1+i2)+ielmt] = src[num_elmt*(nelmt1*i2+i1)+ielmt];
  return 0;

}
/******************************************************************************/
int coproc_transposeOnHost2d(const void *__restrict__ ptrSrc,
			     void *__restrict__ ptrDest,
			     int bytes_per_elmt, int nelmt1, int nelmt2)
{
  switch (bytes_per_elmt) {
  
    //
    // Check for 1,2,4,8 bytes
    //

  case sizeof(uint8_t): // 1 byte
    return coproc_transpose2d<uint8_t>(ptrSrc, ptrDest, nelmt1, nelmt2);

  case sizeof(uint16_t): // 2 byte
    return coproc_transpose2d<uint16_t>(ptrSrc, ptrDest, nelmt1, nelmt2);
    
  case sizeof(uint32_t): // 4 byte
    return coproc_transpose2d<uint32_t>(ptrSrc, ptrDest, nelmt1, nelmt2);
    
  case sizeof(uint64_t): // 8 byte
    return coproc_transpose2d<uint64_t>(ptrSrc, ptrDest, nelmt1, nelmt2);

    //
    // Check for special cases bytes=2^n
    //
    
  case 16:
    return coproc_transpose2d<uint64_t,2>(ptrSrc, ptrDest, nelmt1, nelmt2);
  case 32:
    return coproc_transpose2d<uint64_t,4>(ptrSrc, ptrDest, nelmt1, nelmt2);
  case 64:
    return coproc_transpose2d<uint64_t,8>(ptrSrc, ptrDest, nelmt1, nelmt2);    
  case 128:
    return coproc_transpose2d<uint64_t,16>(ptrSrc, ptrDest, nelmt1, nelmt2);
  case 256:
    return coproc_transpose2d<uint64_t,32>(ptrSrc, ptrDest, nelmt1, nelmt2);
  case 512:
    return coproc_transpose2d<uint64_t,64>(ptrSrc, ptrDest, nelmt1, nelmt2);
  case 1024:
    return coproc_transpose2d<uint64_t,128>(ptrSrc, ptrDest, nelmt1, nelmt2);

  default:
    //
    // Default: check if we have a multiple of 8,4,2 bytes
    //
    if ((bytes_per_elmt&7) == 0) {
      return coproc_transpose2d<uint64_t>((bytes_per_elmt>>3),
					  ptrSrc, ptrDest, nelmt1, nelmt2);
    }
    else if ((bytes_per_elmt&3) == 0) {
      return coproc_transpose2d<uint32_t>((bytes_per_elmt>>2),
					  ptrSrc, ptrDest, nelmt1, nelmt2);
    }
    else if ((bytes_per_elmt&1) == 0) {
      return coproc_transpose2d<uint16_t>((bytes_per_elmt>>1),
					  ptrSrc, ptrDest, nelmt1, nelmt2);
    }
    else {
      return coproc_transpose2d<uint8_t>(bytes_per_elmt,
					 ptrSrc, ptrDest, nelmt1, nelmt2);
    }
  }
}
/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_transposeonhost2d)(const void *__restrict__ ptrSrc,
					void *__restrict__ ptrDest,
					__INT *bytes_per_elmt,
					__INT *nelmt1, __INT *nelmt2)
  {
    return (__INT) coproc_transposeOnHost2d(ptrSrc, ptrDest, *bytes_per_elmt,
					    *nelmt1, *nelmt2);
  }
};
    
/*******************************************************************************
 * Transpose three-dimensional array
 ******************************************************************************/
template<typename T>
int coproc_transpose3d(const void *__restrict__ ptrSrc,
		       void *__restrict__ ptrDest,
		       int nelmt1, int nelmt2, int nelmt3)
{
  T *src = (T*)(ptrSrc);
  T *dest = (T*)(ptrDest);
  
  for (int i3=0; i3<nelmt3; ++i3)
    for (int i2=0; i2<nelmt2; ++i2)
      for (int i1=0; i1<nelmt1; ++i1)
	dest[nelmt2*nelmt3*i1+nelmt3*i2+i3] = src[nelmt1*nelmt2*i3+nelmt1*i2+i1];
  
  return 0;
}
/******************************************************************************/
template<typename T>
int coproc_transpose3d(int num_elmt, const void *__restrict__ ptrSrc, 
		       void *__restrict__ ptrDest, int nelmt1, int nelmt2, int nelmt3)
{
  T *src = (T*)(ptrSrc);
  T *dest = (T*)(ptrDest);
  
  for (int i3=0; i3<nelmt3; ++i3)
    for (int i2=0; i2<nelmt2; ++i2)
      for (int i1=0; i1<nelmt1; ++i1)
	for (int ielmt=0; ielmt<num_elmt; ++ielmt)
	  dest[num_elmt*(nelmt2*nelmt3*i1+nelmt3*i2+i3)+ielmt] =
	    src[num_elmt*(nelmt1*nelmt2*i3+nelmt1*i2+i1)+ielmt];

  return 0;
}
/******************************************************************************/
template<typename T,int num_elmt>
int coproc_transpose3d(const void *__restrict__ ptrSrc,
		       void *__restrict__ ptrDest, int nelmt1, int nelmt2, int nelmt3)
{
  T *src = (T*)(ptrSrc);
  T *dest = (T*)(ptrDest);
  
  for (int i3=0; i3<nelmt3; ++i3)
    for (int i2=0; i2<nelmt2; ++i2)
      for (int i1=0; i1<nelmt1; ++i1)
#pragma unroll
	for (int ielmt=0; ielmt<num_elmt; ++ielmt)
	  dest[num_elmt*(nelmt2*nelmt3*i1+nelmt3*i2+i3)+ielmt] =
	    src[num_elmt*(nelmt1*nelmt2*i3+nelmt1*i2+i1)+ielmt];

  return 0;
}
/******************************************************************************/
int coproc_transposeOnHost3d(const void *__restrict__ ptrSrc,
			     void *__restrict__ ptrDest,
			     int bytes_per_elmt, int nelmt1, int nelmt2, int nelmt3)
{
  switch (bytes_per_elmt) {
  
    //
    // Check for 1,2,4,8 bytes
    //

  case sizeof(uint8_t): // 1 byte
    return coproc_transpose3d<uint8_t>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);

  case sizeof(uint16_t): // 2 byte
    return coproc_transpose3d<uint16_t>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
    
  case sizeof(uint32_t): // 4 byte
    return coproc_transpose3d<uint32_t>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
    
  case sizeof(uint64_t): // 8 byte
    return coproc_transpose3d<uint64_t>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);

    //
    // Check for special cases bytes=2^n
    //
    
  case 16:
    return coproc_transpose3d<uint64_t,2>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
  case 32:
    return coproc_transpose3d<uint64_t,4>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
  case 64:
    return coproc_transpose3d<uint64_t,8>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
  case 128:
    return coproc_transpose3d<uint64_t,16>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
  case 256:
    return coproc_transpose3d<uint64_t,32>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
  case 512:
    return coproc_transpose3d<uint64_t,64>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
  case 1024:
    return coproc_transpose3d<uint64_t,128>(ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);

  default:
    //
    // Default: check if we have a multiple of 8,4,2 bytes
    //
    if ((bytes_per_elmt&7) == 0) {
      return coproc_transpose3d<uint64_t>((bytes_per_elmt>>3),
					  ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
    }
    else if ((bytes_per_elmt&3) == 0) {
      return coproc_transpose3d<uint32_t>((bytes_per_elmt>>2),
					  ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
    }
    else if ((bytes_per_elmt&1) == 0) {
      return coproc_transpose3d<uint16_t>((bytes_per_elmt>>1),
					  ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
    }
    else {
      return coproc_transpose3d<uint8_t>(bytes_per_elmt,
					 ptrSrc, ptrDest, nelmt1, nelmt2, nelmt3);
    }
  }
}
/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_transposeonhost3d)(const void *__restrict__ ptrSrc,
					void *__restrict__ ptrDest,
					__INT *bytes_per_elmt,
					__INT *nelmt1, __INT *nelmt2, __INT *nelmt3)
  {
    return (__INT) coproc_transposeOnHost3d(ptrSrc, ptrDest, *bytes_per_elmt,
					    *nelmt1, *nelmt2, *nelmt3);
  }
};
