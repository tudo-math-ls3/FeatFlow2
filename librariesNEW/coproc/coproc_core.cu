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
 * Device to use in case there is more than one
 *******************************************************************************
 */
static int selectedDevice = 0;

/*******************************************************************************
 * Initialisation of the CUDA subsystem
 *******************************************************************************
 */
int coproc_init()
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0)
    {
      fprintf(stderr, "No CUDA device found!");
      return 1;
    }
  if (selectedDevice >= deviceCount)
    {
      fprintf(stderr, "Choose device ID between 0 and %d!\n", deviceCount-1);
      return 1;
    }
  cudaSetDevice(selectedDevice);
  
  int i = 0;
  printf("Available devices:\n");
  for (i=0 ; i < deviceCount ; ++i)
    {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      printf("%d: %s\n", i, prop.name);
      printf("    memory: %d bytes\n", prop.totalGlobalMem);
      printf("    clockrate: %d khz\n", prop.clockRate);
    }
  printf("Selected device: %d\n", selectedDevice);
  coproc_checkErrors("coproc_init");  
  return 0;
}

int FNAME(coproc_init)()
{
  return coproc_init();
}

/*******************************************************************************
 * Check for CUDA errors
 *******************************************************************************
 */
void coproc_checkErrors(char *label)
{
#ifdef ENABLE_PARAMETER_CHECK
  cudaError_t err;
  
  err = cudaThreadSynchronize();
  if (err != cudaSuccess)
    {
      char *e = (char*) cudaGetErrorString(err);
      printf("CUDA Error: %s (at %s)\n", e, label);
    }
  
  err = cudaGetLastError();
  if (err != cudaSuccess)
    {
      char *e = (char*) cudaGetErrorString(err);
      printf("CUDA Error: %s (at %s)\n", e, label);
    }
#endif
}
