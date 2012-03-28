/*#############################################################################
 ******************************************************************************
 * <name> coproc_storage_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This file provides the basic routines for creating and deleting
 * storage on the device and transfering data between host and device
 * memory are provided.
 * </purpose>
 *
 *#############################################################################
 */

#include <stdio.h>
#include <math.h>
#include "coproc_core.h"
#include "coproc_storage_cuda.h"

/*******************************************************************************
 * Wrappers for malloc/free using cudaMallocHost and cudaHostFree
 *******************************************************************************
 */

int coproc_malloc(void **ptr, size_t size)
{
  if (cudaMallocHost(ptr, size) != cudaSuccess) {
    __coproc__error__("coproc_malloc");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/

int coproc_free(void **ptr)
{
  if (cudaFreeHost(*ptr) != cudaSuccess) {
    __coproc__error__("coproc_free");
    return 1;
  } else {
    *ptr = NULL;
    return 0;
  }
}

/*******************************************************************************
 ***
 ***  HOST CODE
 ***
 *******************************************************************************
 */

/*******************************************************************************
 * Allocate new memory on host
 *******************************************************************************
 */
int coproc_newMemoryOnHost(void **h_ptr,
			   size_t size,
			   unsigned int flags)
{
  if (cudaHostAlloc(h_ptr, size, flags) != cudaSuccess) {
    __coproc__error__("coproc_newMemoryOnHost");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_newmemoryonhost)(void **h_ptr,
				      __SIZET *size,
				      __INT *flags=cudaHostAllocDefault)
  {
    return (__INT) coproc_newMemoryOnHost(h_ptr, *size, *flags);
  }
}

/*******************************************************************************
 * Free memory on host
 *******************************************************************************
 */
int coproc_freeMemoryOnHost(void *h_ptr)
{
  if (cudaFreeHost(h_ptr) != cudaSuccess) {
    __coproc__error__("coproc_freeMemoryOnHost");
    return 1;
  } else {
    h_ptr = NULL;
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_freememoryonhost)(void **h_ptr)
  {
    return (__INT) coproc_freeMemoryOnHost(*h_ptr);
  }
}

/*******************************************************************************
 * Clear memory on host
 *******************************************************************************
 */
int coproc_clearMemoryOnHost(void *h_ptr,
			     size_t size)
{
  memset(h_ptr, 0, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_clearmemoryonhost)(void **h_ptr,
					__SIZET *size)
  {
    return (__INT) coproc_clearMemoryOnHost(*h_ptr, *size);
  }
}

/*******************************************************************************
 ***
 ***  DEVICE CODE
 ***
 *******************************************************************************
 */

/*******************************************************************************
 * Allocate new memory on device
 *******************************************************************************
 */
int coproc_newMemoryOnDevice(void **d_ptr,
			     size_t size)
{
  if (cudaMalloc(d_ptr, size) != cudaSuccess) {
    __coproc__error__("coproc_newMemoryOnDevice");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_newmemoryondevice)(void **d_ptr,
					__SIZET *size)
  {
    return (__INT) coproc_newMemoryOnDevice(d_ptr, *size);
  }
}

/*******************************************************************************
 * Free existing memory on device
 *******************************************************************************
 */
int coproc_freeMemoryOnDevice(void *d_ptr)
{
  if (cudaFree(d_ptr) != cudaSuccess) {
    __coproc__error__("coproc_freeMemoryOnDevice");
    return 1;
  } else {
    d_ptr = NULL;
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_freememoryondevice)(void **d_ptr)
  {
    return (__INT) coproc_freeMemoryOnDevice(*d_ptr);
  }
}

/*******************************************************************************
 * Clear memory on device
 *******************************************************************************
 */
int coproc_clearMemoryOnDevice(void *d_ptr,
			       size_t size)
{
  if (cudaMemset(d_ptr, 0, size) != cudaSuccess) {
    __coproc__error__("coproc_clearMemoryOnDevice");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_clearmemoryondevice)(void **d_ptr,
					  __SIZET *size)
  {
    return (__INT) coproc_clearMemoryOnDevice(*d_ptr, *size);
  }
}

/*******************************************************************************
 * Copy host memory to host memory synchroneously
 *******************************************************************************
 */
int coproc_memcpyHostToHost(void *h_ptrSrc, 
			    void *h_ptrDest,
			    size_t size)
{
  if (cudaMemcpy(h_ptrDest, h_ptrSrc,
		 size, cudaMemcpyHostToHost) != cudaSuccess) {
    __coproc__error__("coproc_memcpyHostToHost");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpyhosttohost)(void **h_ptrSrc,
				       void **h_ptrDest,
				       __SIZET *size)
  {
    return (__INT) coproc_memcpyHostToHost(*h_ptrSrc, *h_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy host memory to host memory asynchroneously
 *******************************************************************************
 */
int coproc_memcpyHostToHostAsync(void *h_ptrSrc, 
				 void *h_ptrDest,
				 size_t size,
				 cudaStream_t stream=0)
{
  cudaMemcpyAsync(h_ptrDest, h_ptrSrc,
		  size, cudaMemcpyHostToHost, stream);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpyhosttohosteasync)(void **h_ptrSrc,
					     void **h_ptrDest,
					     __SIZET *size,
					     __I64 *stream)
  {
    return (__INT) coproc_memcpyHostToHostAsync(*h_ptrSrc, *h_ptrDest, *size,
						(cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Copy host memory to device memory synchroneously
 *******************************************************************************
 */
int coproc_memcpyHostToDevice(void *h_ptrSrc, 
			      void *d_ptrDest,
			      size_t size)
{
  if (cudaMemcpy(d_ptrDest, h_ptrSrc,
		 size, cudaMemcpyHostToDevice) != cudaSuccess) {
    __coproc__error__("coproc_memcpyHostToDevice");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpyhosttodevice)(void *h_ptrSrc,
					 void **d_ptrDest,
					 __SIZET *size)
  {
    return (__INT) coproc_memcpyHostToDevice(h_ptrSrc, *d_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy host memory to device memory asynchroneously
 *******************************************************************************
 */
int coproc_memcpyHostToDeviceAsync(void *h_ptrSrc, 
				   void *d_ptrDest,
				   size_t size,
				   cudaStream_t stream=0)
{
  cudaMemcpyAsync(d_ptrDest, h_ptrSrc,
		  size, cudaMemcpyHostToDevice, stream);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpyhosttodeviceasync)(void **h_ptrSrc,
					      void **d_ptrDest,
					      __SIZET *size,
					      __I64 *stream)
  {
    return (__INT) coproc_memcpyHostToDeviceAsync(*h_ptrSrc, *d_ptrDest, *size,
						  (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Copy device memory to host memory synchroneously
 *******************************************************************************
 */
int coproc_memcpyDeviceToHost(void *d_ptrSrc,
			      void *h_ptrDest,
			      size_t size)
{
  if (cudaMemcpy(h_ptrDest, d_ptrSrc,
		 size, cudaMemcpyDeviceToHost) != cudaSuccess) {
    __coproc__error__("coproc_memcpyDeviceToHost");
    return 1;
  } else {
    return 0;
  }
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpydevicetohost)(void **d_ptrSrc,
					 void *h_ptrDest,
					 __SIZET *size)
  {
    return (__INT) coproc_memcpyDeviceToHost(*d_ptrSrc, h_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy device memory to host memory asynchroneously
 *******************************************************************************
 */
int coproc_memcpyDeviceToHostAsync(void *d_ptrSrc,
				   void *h_ptrDest,
				   size_t size,
				   cudaStream_t stream=0)
{
  cudaMemcpyAsync(h_ptrDest, d_ptrSrc,
		  size, cudaMemcpyDeviceToHost,stream);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpydevicetohostasync)(void **d_ptrSrc,
					      void **h_ptrDest,
					      __SIZET *size,
					      __I64 *stream)
  {
    return (__INT) coproc_memcpyDeviceToHostAsync(*d_ptrSrc, *h_ptrDest, *size,
						  (cudaStream_t)*stream);
  }
}

/*******************************************************************************
 * Copy device memory data to device memory synchroneously
 *******************************************************************************
 */
int coproc_memcpyDeviceToDevice(void *d_ptrSrc,
				void *d_ptrDest,
				size_t size)
{
  if (d_ptrDest != d_ptrSrc) {
    if (cudaMemcpy(d_ptrDest, d_ptrSrc,
		   size, cudaMemcpyDeviceToDevice) != cudaSuccess) {
      __coproc__error__("coproc_memcpyDeviceToDevice");
      return 1;
    } else {
      return 0;
    }
  }
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpydevicetodevice)(void **d_ptrSrc,
					   void **d_ptrDest,
					   __SIZET *size)
  {
    return (__INT) coproc_memcpyDeviceToDevice(*d_ptrSrc, *d_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy device memory data to device memory asynchroneously
 *******************************************************************************
 */
int coproc_memcpyDeviceToDeviceAsync(void *d_ptrSrc,
				     void *d_ptrDest,
				     size_t size,
				     cudaStream_t stream=0)
{
  if (d_ptrDest != d_ptrSrc) {
    cudaMemcpyAsync(d_ptrDest, d_ptrSrc,
		    size, cudaMemcpyDeviceToDevice, stream);
    return 0;
  }
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_memcpydevicetodeviceasync)(void **d_ptrSrc,
						void **d_ptrDest,
						__SIZET *size,
						__I64 *stream)
  {
    return (__INT) coproc_memcpyDeviceToDeviceAsync(*d_ptrSrc, *d_ptrDest, *size,
						    (cudaStream_t)*stream);
  }
}

/*******************************************************************************
 * Combine two generic memory blocks in device memory
 *******************************************************************************
 */
template<typename T>
__global__ void combineOnDevice_knl(T *ptr1,
				    T *ptrSrc2,
				    T *ptrDest,
				    size_t size)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx<size)
    {
      ptrDest[idx] = ptr1[idx] + ptrSrc2[idx];
    }
}

template<>
__global__ void combineOnDevice_knl(bool *ptr1,
				    bool *ptrSrc2,
				    bool *ptrDest,
				    size_t size)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx<size)
    {
      ptrDest[idx] = ptr1[idx] || ptrSrc2[idx];
    }
}

/*******************************************************************************
 * Combine two single memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineSingleOnDevice(void *d_ptrSrc1,
				 void *d_ptrSrc2,
				 void *d_ptrDest,
				 size_t size,
				 cudaStream_t stream=0)
{
  __SP *ptrSrc1 = (__SP*)(d_ptrSrc1);
  __SP *ptrSrc2 = (__SP*)(d_ptrSrc2);
  __SP *ptrDest = (__SP*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combinesingleondevice)(void **d_ptrSrc1,
					    void **d_ptrSrc2,
					    void **d_ptrDest,
					    __SIZET *size,
					    __I64 *stream)
  {
    return (__INT) coproc_combineSingleOnDevice(*d_ptrSrc1, *d_ptrSrc2,
						*d_ptrDest, *size,
						(cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two double memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineDoubleOnDevice(void *d_ptrSrc1,
				 void *d_ptrSrc2,
				 void *d_ptrDest,
				 size_t size,
				 cudaStream_t stream=0)
{
  __DP *ptrSrc1 = (__DP*)(d_ptrSrc1);
  __DP *ptrSrc2 = (__DP*)(d_ptrSrc2);
  __DP *ptrDest = (__DP*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size); 
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combinedoubleondevice)(void **d_ptrSrc1,
					    void **d_ptrSrc2,
					    void **d_ptrDest,
					    __SIZET *size,
					    __I64 *stream)
  {
    return (__INT) coproc_combineDoubleOnDevice(*d_ptrSrc1, *d_ptrSrc2,
						*d_ptrDest, *size,
						(cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two quadrupel memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineQuadOnDevice(void *d_ptrSrc1,
			       void *d_ptrSrc2,
			       void *d_ptrDest,
			       size_t size,
			       cudaStream_t stream=0)
{
  __QP *ptrSrc1 = (__QP*)(d_ptrSrc1);
  __QP *ptrSrc2 = (__QP*)(d_ptrSrc2);
  __QP *ptrDest = (__QP*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size); 
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combinequadondevice)(void **d_ptrSrc1,
					  void **d_ptrSrc2,
					  void **d_ptrDest,
					  __SIZET *size,
					  __I64 *stream)
  {
    return (__INT) coproc_combineQuadOnDevice(*d_ptrSrc1, *d_ptrSrc2,
					      *d_ptrDest, *size,
					      (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two integer memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineIntegerOnDevice(void *d_ptrSrc1,
				  void *d_ptrSrc2,
				  void *d_ptrDest,
				  size_t size,
				  cudaStream_t stream=0)
{
  __INT *ptrSrc1 = (__INT*)(d_ptrSrc1);
  __INT *ptrSrc2 = (__INT*)(d_ptrSrc2);
  __INT *ptrDest = (__INT*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combineintegerondevice)(void **d_ptrSrc1,
					     void **d_ptrSrc2,
					     void **d_ptrDest,
					     __SIZET *size,
					     __I64 *stream)
  {
    return (__INT) coproc_combineIntegerOnDevice(*d_ptrSrc1, *d_ptrSrc2,
						 *d_ptrDest, *size,
						 (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int8 memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineInt8OnDevice(void *d_ptrSrc1,
			       void *d_ptrSrc2,
			       void *d_ptrDest,
			       size_t size,
			       cudaStream_t stream=0)
{
  __I8 *ptrSrc1 = (__I8*)(d_ptrSrc1);
  __I8 *ptrSrc2 = (__I8*)(d_ptrSrc2);
  __I8 *ptrDest = (__I8*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combineint8ondevice)(void **d_ptrSrc1,
					  void **d_ptrSrc2,
					  void **d_ptrDest,
					  __SIZET *size,
					  __I64 *stream)
  {
    return (__INT) coproc_combineInt8OnDevice(*d_ptrSrc1, *d_ptrSrc2,
					      *d_ptrDest, *size,
					      (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int16 memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineInt16OnDevice(void *d_ptrSrc1,
				void *d_ptrSrc2,
				void *d_ptrDest,
				size_t size,
				cudaStream_t stream=0)
{
  __I16 *ptrSrc1 = (__I16*)(d_ptrSrc1);
  __I16 *ptrSrc2 = (__I16*)(d_ptrSrc2);
  __I16 *ptrDest = (__I16*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combineint16ondevice)(void **d_ptrSrc1,
					   void **d_ptrSrc2,
					   void **d_ptrDest,
					   __SIZET *size,
					   __I64 *stream)
  {
    return (__INT) coproc_combineInt16OnDevice(*d_ptrSrc1, *d_ptrSrc2,
					       *d_ptrDest, *size,
					       (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int32 memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineInt32OnDevice(void *d_ptrSrc1,
				void *d_ptrSrc2,
				void *d_ptrDest,
				size_t size,
				cudaStream_t stream=0)
{
  __I32 *ptrSrc1 = (__I32*)(d_ptrSrc1);
  __I32 *ptrSrc2 = (__I32*)(d_ptrSrc2);
  __I32 *ptrDest = (__I32*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combineint32ondevice)(void **d_ptrSrc1,
					   void **d_ptrSrc2,
					   void **d_ptrDest,
					   __SIZET *size,
					   __I64 *stream)
  {
    return (__INT) coproc_combineInt32OnDevice(*d_ptrSrc1, *d_ptrSrc2,
					       *d_ptrDest, *size,
					       (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int64 memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineInt64OnDevice(void *d_ptrSrc1,
				void *d_ptrSrc2,
				void *d_ptrDest,
				size_t size,
				cudaStream_t stream=0)
{
  __I64 *ptrSrc1 = (__I64*)(d_ptrSrc1);
  __I64 *ptrSrc2 = (__I64*)(d_ptrSrc2);
  __I64 *ptrDest = (__I64*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combineint64ondevice)(void **d_ptrSrc1,
					   void **d_ptrSrc2,
					   void **d_ptrDest,
					   __SIZET *size,
					   __I64 *stream)
  {
    return (__INT) coproc_combineInt64OnDevice(*d_ptrSrc1, *d_ptrSrc2,
					       *d_ptrDest, *size,
					       (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two logical memory blocks in device memory
 *******************************************************************************
 */
int coproc_combineLogicalOnDevice(void *d_ptrSrc1,
				  void *d_ptrSrc2,
				  void *d_ptrDest,
				  size_t size,
				  cudaStream_t stream=0)
{
  __LOGICAL *ptrSrc1 = (__LOGICAL*)(d_ptrSrc1);
  __LOGICAL *ptrSrc2 = (__LOGICAL*)(d_ptrSrc2);
  __LOGICAL *ptrDest = (__LOGICAL*)(d_ptrDest);
  
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil(size/(double)(block.x));
  combineOnDevice_knl<<<grid, block, 0, stream>>>(ptrSrc1, ptrSrc2,
						  ptrDest, size);
  return 0;
}

/******************************************************************************/
extern "C" {
  __INT FNAME(coproc_combinelogicalondevice)(void **d_ptrSrc1,
					     void **d_ptrSrc2,
					     void **d_ptrDest,
					     __SIZET *size,
					     __I64 *stream)
  {
    return (__INT) coproc_combineLogicalOnDevice(*d_ptrSrc1, *d_ptrSrc2,
						 *d_ptrDest, *size,
						 (cudaStream_t)(*stream));
  }
}
