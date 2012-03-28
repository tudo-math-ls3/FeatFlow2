/*#############################################################################
 ******************************************************************************
 * <name> coproc_core_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This file provides the core routines for handling coprocessor
 * support in CUDA.
 * </purpose>
 *
 *#############################################################################
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "coproc_core.h"

/*******************************************************************************
 * Check for CUDA errors
 *******************************************************************************
 */
void coproc_checkErrors(__CHAR *label)
{
  #ifdef ENABLE_PARAMETER_CHECK
  cudaError_t err;
  
  err = cudaThreadSynchronize();
  if (err != cudaSuccess)
    {
      char *e = (__CHAR*) cudaGetErrorString(err);
      printf("CUDA Error: %s (at %s)\n", e, label);
    }
  
  err = cudaGetLastError();
  if (err != cudaSuccess)
    {
      char *e = (__CHAR*) cudaGetErrorString(err);
      printf("CUDA Error: %s (at %s)\n", e, label);
    }
  #endif
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_checkerrors)(__CHAR *label)
  {
    coproc_checkErrors(label);
  }
}

/*******************************************************************************
 * Initialisation of the CUDA subsystem
 *******************************************************************************
 */
int coproc_init(int deviceNumber=0)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0)
    {
      __coproc__error__("No CUDA device found");
      return 1;
    }

  int i = 0;
  printf("Available devices:\n");
  for (i=0 ; i < deviceCount ; ++i)
    {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      printf("%d: %s\n", i, prop.name);
      printf("   Total global memory:     %ld bytes\n", (long)prop.totalGlobalMem);
      printf("   Total constant memory:   %ld bytes\n", (long)prop.totalConstMem);
      printf("   Shared memory per block: %ld bytes\n", (long)prop.sharedMemPerBlock);
      printf("   Clock frequency:         %d khz\n",    prop.clockRate);
      printf("   Compute capability;      %d.%d\n",     prop.major,prop.minor);
      printf("   Concurrent memory copy:  %s\n",        prop.deviceOverlap == 1 ?
	                                                "supported" : "not supported");
      printf("   Mapping of host memory:  %s\n",        prop.canMapHostMemory == 1 ?
	                                                "supported" : "not supported");
    }
  
  deviceNumber = max(deviceNumber,0);
  if (deviceNumber >= deviceCount)
    {
      fprintf(stderr, "Choose device ID between 0 and %d!\n", deviceCount-1);
      return 1;
    }
  else
    cudaSetDevice(deviceNumber);
  
  printf("Selected device: %d\n", deviceNumber);
  coproc_checkErrors("coproc_init");
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_init)(__INT *deviceNumber)
  {
    return (__INT) coproc_init((__INT)*deviceNumber);
  }
}

/*******************************************************************************
 * Check for Fortran <-> C interoperability
 *******************************************************************************
 */
size_t coproc_getSizeOf(int cdatatype)
{
  // Check if __SIZET and size_t coincide
  if (sizeof(size_t) != sizeof(__SIZET)) {
    __coproc__error__("coproc_getSizeOf");
    return (size_t) 0;
  }
  
  switch (cdatatype)
    {
    case 1: // ST_SINGLE
      return sizeof(__SP);
    case 2: // ST_DOUBLE
      return sizeof(__DP);
    case 3: // ST_QUAD
      return sizeof(__QP);
    case 4: // ST_INT
      return sizeof(__INT);
    case 5: // ST_INT8
      return sizeof(__I8);
    case 6: // ST_INT16
      return sizeof(__I16);
    case 7: // ST_INT32
      return sizeof(__I32);
    case 8: // ST_INT64
      return sizeof(__I64);
    case 9: // ST_LOGICAL
      return sizeof(__LOGICAL);
    case 10: // ST_CHAR
      return sizeof(__CHAR);
    default:
      __coproc__error__("coproc_getSizeOf");
      return (size_t) 0;
    }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_getsizeof)(__INT *cdatatype, __SIZET *isize)
  {
    *isize = (__SIZET) coproc_getSizeOf((int)*cdatatype);
    return (__INT) 0;
  }
}

/*******************************************************************************
 * Create asynchroneous stream
 *******************************************************************************
 */
int coproc_createStream(cudaStream_t *stream)
{
  if (cudaStreamCreate(stream) != cudaSuccess) {
    __coproc__error__("coproc_createStream");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_createstream)(__I64 *stream)
  {
    return coproc_createStream((cudaStream_t*)stream);
  }
}

/*******************************************************************************
 * Destroy asynchroneous stream
 *******************************************************************************
 */
int coproc_destroyStream(cudaStream_t stream)
{
  if (cudaStreamDestroy(stream) != cudaSuccess) {
    __coproc__error__("coproc_destroyStream");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_destroystream)(__I64 *stream)
  {
    return coproc_destroyStream((cudaStream_t)*stream);
  }
}

/*******************************************************************************
 * Synchronize asynchroneous stream
 *******************************************************************************
 */
int coproc_synchronizeStream(cudaStream_t stream)
{
  if (cudaStreamSynchronize(stream) != cudaSuccess) {
    __coproc__error__("coproc_synchronizeStream");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_synchronizestream)(__I64 *stream)
  {
    return coproc_synchronizeStream((cudaStream_t)*stream);
  }
}
