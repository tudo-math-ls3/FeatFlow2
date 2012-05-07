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

#include <cmath>
#include <iostream>
#include "coproc_core.h"

using namespace std;

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
      cout << "CUDA Error: " << e << "(at " << label << ")" << endl;
    }
  
  err = cudaGetLastError();
  if (err != cudaSuccess)
    {
      char *e = (__CHAR*) cudaGetErrorString(err);
      cout << "CUDA Error: " << e << "(at " << label << ")" << endl;
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

  cout << "Available devices:" << endl;
  for (int i=0 ; i < deviceCount ; ++i)
    {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      cout << i << ": " << prop.name << endl;
      cout << "   Compute capability;      " << prop.major << "." << prop.minor << endl;
      cout << "   Total global memory:     " << prop.totalGlobalMem << " bytes" << endl;
      cout << "   Total constant memory:   " << prop.totalConstMem << " bytes" << endl;
      cout << "   Shared memory per block: " << prop.sharedMemPerBlock << " bytes" << endl;
      cout << "   Multiprocessors:         " << prop.multiProcessorCount << endl;
      cout << "   Blocks per grid:         "
	   << prop.maxGridSize[0]*prop.maxGridSize[1]*prop.maxGridSize[2]
	   << " (" << prop.maxGridSize[0]
	   << "," << prop.maxGridSize[1]
	   << "," << prop.maxGridSize[2] << ")" << endl;
      cout << "   Registers per block:     " << prop.regsPerBlock << endl;
      cout << "   Threads per block:       "
	   << prop.maxThreadsPerBlock
	   << " (" << prop.maxThreadsDim[0]
	   << "," << prop.maxThreadsDim[1]
	   << "," << prop.maxThreadsDim[2] << ")" << endl;
      cout << "   Threads per warp:        " << prop.warpSize << endl;
      cout << "   Clock frequency:         " << prop.clockRate << " kHz" << endl;
      cout << "   Concurrent memory copy:  "
	   << (prop.deviceOverlap == 1 ? "supported" : "not supported") << endl;
      cout << "   Mapping of host memory:  "
	<< (prop.canMapHostMemory == 1 ? "supported" : "not supported") << endl;
      cout << "   Concurrent kernels:      "
	   << (prop.concurrentKernels == 1 ? "supported" : "not supported") << endl;
    }
  
  // Autoselect device?
  if (deviceNumber < 0) {
    cudaDeviceProp prop;
    memset(&prop, 0, sizeof(cudaDeviceProp));
    
    // Set major compute capability
#if defined(HAS_CUDA10) || defined(HAS_CUDA11) || defined(HAS_CUDA12) || defined(HAS_CUDA13)
    prop.major = 1;
#endif

#if defined(HAS_CUDA20) || defined(HAS_CUDA21)
    prop.major = 2;
#endif

#if defined(HAS_CUDA30)
    prop.major = 3;
#endif

    // Set minor compute capability
#if defined(HAS_CUDA10) || defined(HAS_CUDA20) || defined(HAS_CUDA30)
    prop.minor = 0;
#endif

#if defined(HAS_CUDA11) || defined(HAS_CUDA21)
    prop.minor = 1;
#endif

#if defined(HAS_CUDA12)
    prop.minor = 2;
#endif

#if defined(HAS_CUDA13)
    prop.minor = 3;
#endif

    // Choose device according to compute capability
    int dev;
    cudaChooseDevice(&dev, &prop);
    cudaSetDevice(dev);
    cout << "Selected device: " << dev << endl;
  } 
  else if (deviceNumber >= deviceCount) {
    cerr << "Choose device ID between 0 and " << deviceCount-1 << "!" << endl;
    return 1;
  }
  else {
    cudaSetDevice(deviceNumber);
    cout << "Selected device: " << deviceNumber << endl;
  }
  
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
