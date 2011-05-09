/*#############################################################################
 ******************************************************************************
 * <name> coproc_storage_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This file provides the basic routines for handling coprocessor
 * support in CUDA. Routines for creating and deleting storage on the
 * device and transfering data between host and device memory are
 * provided.
 * </purpose>
 *
 *#############################################################################
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "coproc_core.h"
#include "coproc_storage_cuda.h"


/*******************************************************************************
 * Create new double data on device
 *******************************************************************************
 */
int coproc_newDoubleOnDevice(unsigned long * h_Ddata, int * n)
{
  double * d_Ddata = 0;

  cudaMalloc((void**)&d_Ddata, (*n)*sizeof(double));
  *h_Ddata = (unsigned long) d_Ddata;
  coproc_checkErrors("coproc_newDoubleOnDevice");
  return 0;
}

int FNAME(coproc_newdoubleondevice)(unsigned long * h_Ddata, int * n)
{
  return coproc_newDoubleOnDevice(h_Ddata, n);
}

/*******************************************************************************
 * Create new single data on device
 *******************************************************************************
 */
int coproc_newSingleOnDevice(unsigned long * h_Fdata, int * n)
{
  float * d_Fdata = 0;
  
  cudaMalloc((void**)&d_Fdata, (*n)*sizeof(float));
  *h_Fdata = (unsigned long) d_Fdata;
  coproc_checkErrors("coproc_newSingleOnDevice");
  return 0;
}

int FNAME(coproc_newsingleondevice)(unsigned long * h_Fdata, int * n)
{
  return coproc_newSingleOnDevice(h_Fdata, n);
}

/*******************************************************************************
 * Create new integer data on device
 *******************************************************************************
 */
int coproc_newIntOnDevice(unsigned long * h_Idata, int * n)
{
  int * d_Idata = 0;

  cudaMalloc((void**)&d_Idata, (*n)*sizeof(int));
  *h_Idata = (unsigned long) d_Idata;
  coproc_checkErrors("coproc_newIntOnDevice");
  return 0;
}

int FNAME(coproc_newintondevice)(unsigned long * h_Idata, int * n)
{
  return coproc_newIntOnDevice(h_Idata, n);
}

/*******************************************************************************
 * Free existing double data on device
 *******************************************************************************
 */
int coproc_freeDoubleOnDevice(unsigned long * h_Ddata)
{
  double * d_Ddata = (double*)(*h_Ddata);

  cudaFree(d_Ddata);
  *h_Ddata = 0;
  coproc_checkErrors("coproc_freeDoubleOnDevice");
  return 0;
}

int FNAME(coproc_freedoubleondevice)(unsigned long * h_Ddata)
{
  return coproc_freeDoubleOnDevice(h_Ddata);
}

/*******************************************************************************
 * Free existing single data on device
 *******************************************************************************
 */
int coproc_freeSingleOnDevice(unsigned long * h_Fdata)
{
  float * d_Fdata = (float*)(*h_Fdata);

  cudaFree(d_Fdata);
  *h_Fdata = 0;
  coproc_checkErrors("coproc_freeSingleOnDevice");
  return 0;
}

int FNAME(coproc_freesingleondevice)(unsigned long * h_Fdata)
{
  return coproc_freeSingleOnDevice(h_Fdata);
}

/*******************************************************************************
 * Free existing integer data on device
 *******************************************************************************
 */
int coproc_freeIntOnDevice(unsigned long * h_Idata)
{
  int * d_Idata = (int*)(*h_Idata);

  cudaFree(d_Idata);
  *h_Idata = 0;
  coproc_checkErrors("coproc_freeIntOnDevice");
  return 0;
}

int FNAME(coproc_freeintondevice)(unsigned long * h_Idata)
{
  return coproc_freeIntOnDevice(h_Idata);
}

/*******************************************************************************
 * Clear double data on device
 *******************************************************************************
 */
int coproc_clearDoubleOnDevice(unsigned long * h_Ddata, int * n)
{
  cudaMemset((void*)(*h_Ddata), 0, (*n)*sizeof(double));
  coproc_checkErrors("coproc_clearDoubleOnDevice");
  return 0;
}

int FNAME(coproc_cleardoubleondevice)(unsigned long * h_Ddata, int * n)
{
  return coproc_clearDoubleOnDevice(h_Ddata, n);
}

/*******************************************************************************
 * Clear single data on device
 *******************************************************************************
 */
int coproc_clearSingleOnDevice(unsigned long * h_Fdata, int * n)
{
  cudaMemset((void*)(*h_Fdata), 0, (*n)*sizeof(float));
  coproc_checkErrors("coproc_clearSingleOnDevice");
  return 0;
}

int FNAME(coproc_clearsingleondevice)(unsigned long * h_Fdata, int * n)
{
  return coproc_clearSingleOnDevice(h_Fdata, n);
}

/*******************************************************************************
 * Clear integer data on device
 *******************************************************************************
 */
int coproc_clearIntOnDevice(unsigned long * h_Idata, int * n)
{
  cudaMemset((void*)(*h_Idata), 0, (*n)*sizeof(int));
  coproc_checkErrors("coproc_clearIntOnDevice");
  return 0;
}

int FNAME(coproc_clearintondevice)(unsigned long * h_Idata, int * n)
{
  return coproc_clearIntOnDevice(h_Idata, n);
}

/*******************************************************************************
 * Copy host double data to device
 *******************************************************************************
 */
int coproc_copyDoubleToDevice(double * Ddata,  unsigned long * h_Ddata, int * n)
{
  double * d_Ddata = (double*)(*h_Ddata);

  cudaMemcpy(d_Ddata, Ddata, (*n)*sizeof(double), cudaMemcpyHostToDevice);
  coproc_checkErrors("coproc_copyDoubleToDevice");
  return 0;
}

int FNAME(coproc_copydoubletodevice)(double * Ddata, unsigned long * h_Ddata, int * n)
{
  return coproc_copyDoubleToDevice(Ddata, h_Ddata, n);
}

/********************************************************************************
 * Copy host single data to device
 ********************************************************************************
 */
int coproc_copySingleToDevice(float * Fdata, unsigned long * h_Fdata, int * n)
{
  float * d_Fdata = (float*)(*h_Fdata);

  cudaMemcpy(d_Fdata, Fdata, (*n)*sizeof(float), cudaMemcpyHostToDevice);
  coproc_checkErrors("coproc_copySingleToDevice");
  return 0;
}

int FNAME(coproc_copysingletodevice)(float * Fdata, unsigned long * h_Fdata, int * n)
{
  return coproc_copySingleToDevice(Fdata, h_Fdata, n);
}

/********************************************************************************
 * Copy host integer data to device
 ********************************************************************************
 */
int coproc_copyIntToDevice(int * Idata, unsigned long * h_Idata, int * n)
{
  int * d_Idata = (int*)(*h_Idata);

  cudaMemcpy(d_Idata, Idata, (*n)*sizeof(int), cudaMemcpyHostToDevice);
  coproc_checkErrors("coproc_copyIntToDevice");
  return 0;
}

int FNAME(coproc_copyinttodevice)(int * Idata, unsigned long * h_Idata, int * n)
{
  return coproc_copyIntToDevice(Idata, h_Idata, n);
}

/*******************************************************************************
 * Copy device double data to host
 *******************************************************************************
 */
int coproc_copyDoubleToHost(unsigned long * h_Ddata, double * Ddata, int * n)
{
  double * d_Ddata = (double*)(*h_Ddata);

  cudaMemcpy(Ddata, d_Ddata, (*n)*sizeof(double), cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  coproc_checkErrors("coproc_copyDoubleToHost");
  return 0;
}

int FNAME(coproc_copydoubletohost)(unsigned long * h_Ddata, double * Ddata, int * n)
{
  return coproc_copyDoubleToHost(h_Ddata, Ddata, n);
}

/*******************************************************************************
 * Copy device single data to host
 *******************************************************************************
 */
int coproc_copySingleToHost(unsigned long * h_Fdata, float * Fdata, int * n)
{
  float * d_Fdata = (float*)(*h_Fdata);

  cudaMemcpy(Fdata, d_Fdata, (*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  coproc_checkErrors("coproc_copySingleToHost");
  return 0;
}

int FNAME(coproc_copysingletohost)(unsigned long * h_Fdata, float * Fdata, int * n)
{
  return coproc_copySingleToHost(h_Fdata, Fdata, n);
}

/*******************************************************************************
 * Copy device integer data to host
 *******************************************************************************
 */
int coproc_copyIntToHost(unsigned long * h_Idata, int * Idata, int * n)
{
  int * d_Idata = (int*)(*h_Idata);

  cudaMemcpy(Idata, d_Idata, (*n)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  coproc_checkErrors("coproc_copyIntToHost");
  return 0;
}

int FNAME(coproc_copyinttohost)(unsigned long * h_Idata, int * Idata, int * n)
{
  return coproc_copyIntToHost(h_Idata, Idata, n);
}

