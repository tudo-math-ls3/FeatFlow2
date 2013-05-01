/*#############################################################################
 ******************************************************************************
 * <name> coproc_core </name>
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

static cudaDeviceProp *DeviceProp;

/*******************************************************************************
 * Check for CUDA error
 ******************************************************************************/
void coproc_checkError(const char *label)
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
  void FNAME(coproc_checkerror)(__CHAR *label)
  {
    coproc_checkError(label);
  }
}

/*******************************************************************************
 * Initialisation of the CUDA subsystem
 ******************************************************************************/
void coproc_init(int deviceNumber=0)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0)
    {
      coproc_throwError("No CUDA device found");
    }
  
  // Initialise device property array
  DeviceProp = new cudaDeviceProp[deviceCount];

  cout << "Available devices:" << endl;
  for (int i=0 ; i < deviceCount ; ++i)
    {
      cudaGetDeviceProperties(&DeviceProp[i], i);
      cout << i << ": " << DeviceProp[i].name << endl;
      cout << "   Compute capability;      " << DeviceProp[i].major << "." << DeviceProp[i].minor << endl;
      cout << "   Total global memory:     " << DeviceProp[i].totalGlobalMem << " bytes" << endl;
      cout << "   Total constant memory:   " << DeviceProp[i].totalConstMem << " bytes" << endl;
      cout << "   Shared memory per block: " << DeviceProp[i].sharedMemPerBlock << " bytes" << endl;
      cout << "   Multiprocessors:         " << DeviceProp[i].multiProcessorCount << endl;
      cout << "   Blocks per grid:         "
		   << "(" << DeviceProp[i].maxGridSize[0]
		   << "," << DeviceProp[i].maxGridSize[1]
		   << "," << DeviceProp[i].maxGridSize[2] << ")" << endl;
      cout << "   Registers per block:     " << DeviceProp[i].regsPerBlock << endl;
      cout << "   Threads per block:       "
		   << DeviceProp[i].maxThreadsPerBlock
		   << " (" << DeviceProp[i].maxThreadsDim[0]
		   << "," << DeviceProp[i].maxThreadsDim[1]
		   << "," << DeviceProp[i].maxThreadsDim[2] << ")" << endl;
      cout << "   Threads per warp:        " << DeviceProp[i].warpSize << endl;
      cout << "   Clock frequency:         " << DeviceProp[i].clockRate << " kHz" << endl;
      cout << "   Concurrent memory copy:  "
		   << (DeviceProp[i].deviceOverlap == 1 ? "supported" : "not supported") << endl;
      cout << "   Mapping of host memory:  "
		   << (DeviceProp[i].canMapHostMemory == 1 ? "supported" : "not supported") << endl;
      cout << "   Concurrent kernels:      "
		   << (DeviceProp[i].concurrentKernels == 1 ? "supported" : "not supported") << endl;
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
    exit(-1);
  }
  else {
    cudaSetDevice(deviceNumber);
    cout << "Selected device: " << deviceNumber << endl;
  }
  
  coproc_checkError("coproc_init");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_init)(__INT *deviceNumber)
  {
    coproc_init((__INT)*deviceNumber);
  }
}

/*******************************************************************************
 * Finalisation of the CUDA subsystem
 ******************************************************************************/
void coproc_done()
{
  delete[] DeviceProp;
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_done)()
  {
    coproc_done();
  }
}

/*******************************************************************************
 * Get property list of device deviceNumber
 ******************************************************************************/
extern "C" {
  const cudaDeviceProp* coproc_getDeviceProp(int deviceNumber)
  {
    return &DeviceProp[deviceNumber];
  }
}

/*******************************************************************************
 * Get property list of current device
 ******************************************************************************/
extern "C" {
  const cudaDeviceProp* coproc_getCurrentDeviceProp()
  {
    int deviceNumber;
    cudaGetDevice(&deviceNumber);
    return &DeviceProp[deviceNumber];
  }
}

/*******************************************************************************
 * Check for Fortran <-> C interoperability
 ******************************************************************************/
size_t coproc_getSizeOf(int cdatatype)
{
  // Check if __SIZET and size_t coincide
  if (sizeof(size_t) != sizeof(__SIZET)) {
    coproc_throwError("coproc_getSizeOf");
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
      coproc_throwError("coproc_getSizeOf");
      exit(-1);
    }
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_getsizeof)(__INT *cdatatype, __SIZET *isize)
  {
    *isize = (__SIZET) coproc_getSizeOf((int)*cdatatype);
  }
}

/*******************************************************************************
 * Create asynchroneous stream
 ******************************************************************************/
void coproc_createStream(cudaStream_t *stream)
{
  cudaStreamCreate(stream);
  coproc_checkError("coproc_createStream");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_createstream)(__I64 *stream)
  {
    coproc_createStream((cudaStream_t*)stream);
  }
}

/*******************************************************************************
 * Destroy asynchroneous stream
 ******************************************************************************/
void coproc_destroyStream(cudaStream_t stream)
{
  cudaStreamDestroy(stream);
  coproc_checkError("coproc_destroyStream");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_destroystream)(__I64 *stream)
  {
    coproc_destroyStream((cudaStream_t)*stream);
  }
}

/*******************************************************************************
 * Synchronize asynchroneous stream
 ******************************************************************************/
void coproc_synchronizeStream(cudaStream_t stream)
{
  cudaStreamSynchronize(stream);
  coproc_checkError("coproc_synchronizeStream");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_synchronizestream)(__I64 *stream)
  {
    coproc_synchronizeStream((cudaStream_t)*stream);
  }
}
