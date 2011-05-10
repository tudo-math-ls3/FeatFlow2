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
 * Create new memory on device
 *******************************************************************************
 */
int coproc_newMemoryOnDevice(unsigned long * p_MemoryBlock,
			     unsigned long * imemBytes)
{
  void * d_MemoryBlock = 0;

  cudaMalloc((void**)&d_MemoryBlock, *imemBytes);
  *p_MemoryBlock = (unsigned long) d_MemoryBlock;
  coproc_checkErrors("coproc_newMemoryOnDevice");
  return 0;
}

int FNAME(coproc_newmemoryondevice)(unsigned long * p_MemoryBlock,
				    unsigned long * imemBytes)
{
  return coproc_newMemoryOnDevice(p_MemoryBlock, imemBytes);
}


/*******************************************************************************
 * Free existing memory on device
 *******************************************************************************
 */
int coproc_freeMemoryOnDevice(unsigned long * p_MemoryBlock)
{
  void * d_MemoryBlock = (void*)(*p_MemoryBlock);

  cudaFree(d_MemoryBlock);
  *p_MemoryBlock = 0;
  coproc_checkErrors("coproc_freeMemoryOnDevice");
  return 0;
}

int FNAME(coproc_freememoryondevice)(unsigned long * p_MemoryBlock)
{
  return coproc_freeMemoryOnDevice(p_MemoryBlock);
}


/*******************************************************************************
 * Clear memory on device
 *******************************************************************************
 */
int coproc_clearMemoryOnDevice(unsigned long * p_MemoryBlock,
			       unsigned long * imemBytes)
{
  cudaMemset((void*)(*p_MemoryBlock), 0, *imemBytes);
  coproc_checkErrors("coproc_clearMemoryOnDevice");
  return 0;
}

int FNAME(coproc_clearmemoryondevice)(unsigned long * p_MemoryBlock,
				      unsigned long * imemBytes)
{
  return coproc_clearMemoryOnDevice(p_MemoryBlock, imemBytes);
}


/*******************************************************************************
 * Copy host memory to device memory
 *******************************************************************************
 */
int coproc_copyMemoryHostToDevice(unsigned long * p_MemoryBlockOnHost, 
				  unsigned long * p_MemoryBlockOnDevice,
				  unsigned long * imemBytes)
{
  void * d_MemoryBlockOnDevice = (void*)(*p_MemoryBlockOnDevice);

  cudaMemcpy(d_MemoryBlockOnDevice, p_MemoryBlockOnHost,
	     *imemBytes, cudaMemcpyHostToDevice);
  coproc_checkErrors("coproc_copyMemoryHostToDev");
  return 0;
}

int FNAME(coproc_copymemoryhosttodevice)(unsigned long * p_MemoryBlockOnHost,
					 unsigned long * p_MemoryBlockOnDevice,
					 unsigned long * imemBytes)
{
  return coproc_copyMemoryHostToDevice(p_MemoryBlockOnHost,
				       p_MemoryBlockOnDevice, imemBytes);
}


/*******************************************************************************
 * Copy device memory to host memory
 *******************************************************************************
 */
int coproc_copyMemoryDeviceToHost(unsigned long * p_MemoryBlockOnDevice,
				  unsigned long * p_MemoryBlockOnHost,
				  unsigned long * imemBytes)
{
  void * d_MemoryBlockOnDevice = (void*)(*p_MemoryBlockOnDevice);

  cudaMemcpy(p_MemoryBlockOnHost, d_MemoryBlockOnDevice,
	     *imemBytes, cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  coproc_checkErrors("coproc_copyMemoryDeviceToHost");
  return 0;
}

int FNAME(coproc_copymemorydevicetohost)(unsigned long * p_MemoryBlockOnDevice,
					 unsigned long * p_MemoryBlockOnHost,
					 unsigned long * imemBytes)
{
  return coproc_copyMemoryDeviceToHost(p_MemoryBlockOnDevice,
				       p_MemoryBlockOnHost, imemBytes);
}


/*******************************************************************************
 * Copy device memory data to device memory
 *******************************************************************************
 */
int coproc_copyMemoryDeviceToDevice(unsigned long * p_MemoryBlockSrc,
				    unsigned long * p_MemoryBlockDest,
				    unsigned long * imemBytes)
{
  void * d_MemoryBlockSrc  = (void*)(*p_MemoryBlockSrc);
  void * d_MemoryBlockDest = (void*)(*p_MemoryBlockDest);
  
  if (d_MemoryBlockDest != d_MemoryBlockSrc) {
    cudaMemcpy(d_MemoryBlockDest, d_MemoryBlockSrc,
	       *imemBytes, cudaMemcpyDeviceToDevice);
    cudaThreadSynchronize();
    coproc_checkErrors("coproc_copyMemoryDeviceToDevice");
  }
  return 0;
}

int FNAME(coproc_copymemorydevicetodevice)(unsigned long * p_MemoryBlockSrc,
					   unsigned long * p_MemoryBlockDest,
					   unsigned long * imemBytes)
{
  return coproc_copyMemoryDeviceToDevice(p_MemoryBlockSrc,
					 p_MemoryBlockDest, imemBytes);
}

/*******************************************************************************
 * Add two memory blocks in device memory
 *******************************************************************************
 */
int coproc_addMemoryOnDevice(unsigned long * p_MemoryBlock1,
			     unsigned long * p_MemoryBlock2,
			     unsigned long * p_MemoryBlockDest,
			     unsigned long * imemBytes)
{
  void * d_MemoryBlock1    = (void*)(*p_MemoryBlock1);
  void * d_MemoryBlock2    = (void*)(*p_MemoryBlock2);
  void * d_MemoryBlockDest = (void*)(*p_MemoryBlockDest);

  // Here, we need to write a CUDA kernel.

  return 0;
}

int FNAME(coproc_addmemoryondevice)(unsigned long * p_MemoryBlock1,
				    unsigned long * p_MemoryBlock2,
				    unsigned long * p_MemoryBlockDest,
				    unsigned long * imemBytes)
{
  return coproc_addMemoryOnDevice(p_MemoryBlock1, p_MemoryBlock2,
				  p_MemoryBlockDest, imemBytes);
}
