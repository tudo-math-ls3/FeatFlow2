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

#include <cmath>
#include <iostream>
#include "coproc_core.h"
#include "coproc_storage_cuda.h"
#include "coproc_transpose.h"

/*******************************************************************************
 * Wrappers for malloc/free using cudaMallocHost and cudaHostFree
 ******************************************************************************/

void coproc_malloc(void **ptr, size_t size)
{
  cudaMallocHost(ptr, size);
  coproc_checkError("coproc_malloc");
}

/******************************************************************************/

void coproc_free(void **ptr)
{
  cudaFreeHost(*ptr);
  coproc_checkError("coproc_free");
}

/*******************************************************************************
 ***
 ***  HOST CODE
 ***
 ******************************************************************************/

/*******************************************************************************
 * Allocate new memory on host
 ******************************************************************************/
void coproc_newMemoryOnHost(void **h_ptr,
			    size_t size,
			    unsigned int flags)
{
  cudaHostAlloc(h_ptr, size, flags);
  coproc_checkError("coproc_newMemoryOnHost");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_newmemoryonhost)(void **h_ptr,
				     __SIZET *size,
				     __INT *flags=cudaHostAllocDefault)
  {
    coproc_newMemoryOnHost(h_ptr, *size, *flags);
  }
}

/*******************************************************************************
 * Free memory on host
 ******************************************************************************/
void coproc_freeMemoryOnHost(void *h_ptr)
{
  cudaFreeHost(h_ptr);
  coproc_checkError("coproc_freeMemoryOnHost");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_freememoryonhost)(void **h_ptr)
  {
    coproc_freeMemoryOnHost(*h_ptr);
  }
}

/*******************************************************************************
 * Clear memory on host
 ******************************************************************************/
void coproc_clearMemoryOnHost(void * __restrict__ h_ptr,
			      size_t size)
{
  memset(h_ptr, 0, size);
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_clearmemoryonhost)(void **h_ptr,
				       __SIZET *size)
  {
    coproc_clearMemoryOnHost(*h_ptr, *size);
  }
}

/*******************************************************************************
 ***
 ***  DEVICE CODE
 ***
 ******************************************************************************/

/*******************************************************************************
 * Allocate new memory on device
 ******************************************************************************/
void coproc_newMemoryOnDevice(void **d_ptr,
			      size_t size)
{
  cudaMalloc(d_ptr, size);
  coproc_checkError("coproc_newMemoryOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_newmemoryondevice)(void **d_ptr,
				       __SIZET *size)
  {
    coproc_newMemoryOnDevice(d_ptr, *size);
  }
}

/*******************************************************************************
 * Free existing memory on device
 ******************************************************************************/
void coproc_freeMemoryOnDevice(void *d_ptr)
{
  cudaFree(d_ptr);
  coproc_checkError("coproc_freeMemoryOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_freememoryondevice)(void **d_ptr)
  {
    coproc_freeMemoryOnDevice(*d_ptr);
  }
}

/*******************************************************************************
 * Clear memory on device
 ******************************************************************************/
void coproc_clearMemoryOnDevice(void * __restrict__ d_ptr,
				size_t size)
{
  cudaMemset(d_ptr, 0, size);
  coproc_checkError("coproc_clearMemoryOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_clearmemoryondevice)(void **d_ptr,
					 __SIZET *size)
  {
    coproc_clearMemoryOnDevice(*d_ptr, *size);
  }
}

/*******************************************************************************
 * Copy host memory to host memory synchroneously
 ******************************************************************************/
void coproc_memcpyHostToHost(const void * __restrict__ h_ptrSrc, 
			     void * __restrict__ h_ptrDest,
			     size_t size)
{
  cudaMemcpy(h_ptrDest, h_ptrSrc, size, cudaMemcpyHostToHost);
  coproc_checkError("coproc_memcpyHostToHost");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpyhosttohost)(const void **h_ptrSrc,
				      void **h_ptrDest,
				      __SIZET *size)
  {
    coproc_memcpyHostToHost(*h_ptrSrc, *h_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy host memory to host memory asynchroneously
 ******************************************************************************/
void coproc_memcpyHostToHostAsync(const void * __restrict__ h_ptrSrc, 
				  void * __restrict__ h_ptrDest,
				  size_t size,
				  cudaStream_t stream)
{
  cudaMemcpyAsync(h_ptrDest, h_ptrSrc,
		  size, cudaMemcpyHostToHost, stream);
  coproc_checkError("coproc_memcpyHostToHostAsync");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpyhosttohosteasync)(const void **h_ptrSrc,
					    void **h_ptrDest,
					    __SIZET *size,
					    __I64 *stream)
  {
    coproc_memcpyHostToHostAsync(*h_ptrSrc, *h_ptrDest, *size,
				 (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Copy host memory to device memory synchroneously
 ******************************************************************************/
void coproc_memcpyHostToDevice(const void * __restrict__ h_ptrSrc, 
			       void * __restrict__ d_ptrDest,
			       size_t size)
{
  cudaMemcpy(d_ptrDest, h_ptrSrc, size, cudaMemcpyHostToDevice);
  coproc_checkError("coproc_memcpyHostToDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpyhosttodevice)(const void *h_ptrSrc,
					void **d_ptrDest,
					__SIZET *size)
  {
    coproc_memcpyHostToDevice(h_ptrSrc, *d_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy host memory to device memory asynchroneously
 ******************************************************************************/
void coproc_memcpyHostToDeviceAsync(const void * __restrict__ h_ptrSrc, 
				    void * __restrict__ d_ptrDest,
				    size_t size,
				    cudaStream_t stream)
{
  cudaMemcpyAsync(d_ptrDest, h_ptrSrc,
		  size, cudaMemcpyHostToDevice, stream);
  coproc_checkError("coproc_memcpyHostToDeviceAsync");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpyhosttodeviceasync)(const void **h_ptrSrc,
					     void **d_ptrDest,
					     __SIZET *size,
					     __I64 *stream)
  {
    coproc_memcpyHostToDeviceAsync(*h_ptrSrc, *d_ptrDest, *size,
				   (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Copy device memory to host memory synchroneously
 ******************************************************************************/
void coproc_memcpyDeviceToHost(const void * __restrict__  d_ptrSrc,
			       void * __restrict__ h_ptrDest,
			       size_t size)
{
  cudaMemcpy(h_ptrDest, d_ptrSrc, size, cudaMemcpyDeviceToHost);
  coproc_checkError("coproc_memcpyDeviceToHost");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpydevicetohost)(const void **d_ptrSrc,
					void *h_ptrDest,
					__SIZET *size)
  {
    coproc_memcpyDeviceToHost(*d_ptrSrc, h_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy device memory to host memory asynchroneously
 ******************************************************************************/
void coproc_memcpyDeviceToHostAsync(const void * __restrict__ d_ptrSrc,
				    void * __restrict__ h_ptrDest,
				    size_t size,
				    cudaStream_t stream)
{
  cudaMemcpyAsync(h_ptrDest, d_ptrSrc,
		  size, cudaMemcpyDeviceToHost,stream);
  coproc_checkError("coproc_memcpyDeviceToHostAsync");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpydevicetohostasync)(const void **d_ptrSrc,
					     void **h_ptrDest,
					     __SIZET *size,
					     __I64 *stream)
  {
    coproc_memcpyDeviceToHostAsync(*d_ptrSrc, *h_ptrDest, *size,
				   (cudaStream_t)*stream);
  }
}

/*******************************************************************************
 * Copy device memory data to device memory synchroneously
 ******************************************************************************/
void coproc_memcpyDeviceToDevice(const void * __restrict__ d_ptrSrc,
				 void * __restrict__ d_ptrDest,
				 size_t size)
{
  if (d_ptrDest != d_ptrSrc) {
    cudaMemcpy(d_ptrDest, d_ptrSrc, size, cudaMemcpyDeviceToDevice);
    coproc_checkError("coproc_memcpyDeviceToDevice");
  }
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpydevicetodevice)(const void **d_ptrSrc,
					  void **d_ptrDest,
					  __SIZET *size)
  {
    coproc_memcpyDeviceToDevice(*d_ptrSrc, *d_ptrDest, *size);
  }
}

/*******************************************************************************
 * Copy device memory data to device memory asynchroneously
 ******************************************************************************/
void coproc_memcpyDeviceToDeviceAsync(const void * __restrict__ d_ptrSrc,
				      void * __restrict__ d_ptrDest,
				      size_t size,
				      cudaStream_t stream)
{
  if (d_ptrDest != d_ptrSrc) {
    cudaMemcpyAsync(d_ptrDest, d_ptrSrc,
		    size, cudaMemcpyDeviceToDevice, stream);
    coproc_checkError("coproc_memcpyDeviceToDeviceAsync");
  }
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_memcpydevicetodeviceasync)(const void **d_ptrSrc,
					       void **d_ptrDest,
					       __SIZET *size,
					       __I64 *stream)
  {
    coproc_memcpyDeviceToDeviceAsync(*d_ptrSrc, *d_ptrDest, *size,
				     (cudaStream_t)*stream);
  }
}

/*******************************************************************************
 * Copy host memory to host memory synchroneously (using transposition)
 ******************************************************************************/

//
// 2d version
//
void coproc_tmemcpyHostToHost2d(const void * __restrict__ h_ptrSrc, 
				void * __restrict__ h_ptrDest,
				size_t size, int nelmt1, int nelmt2,
				int packsize)
{
  const int n = nelmt1/packsize;
  const int bytes_per_elmt = (int)(size/nelmt2/n);
  coproc_transposeOnHost2d(h_ptrSrc, h_ptrDest,
			   bytes_per_elmt, n, nelmt2);
}

//
// 3d version
//
void coproc_tmemcpyHostToHost3d(const void * __restrict__ h_ptrSrc, 
				void * __restrict__ h_ptrDest,
				size_t size, int nelmt1, int nelmt2,
				int nelmt3, int packsize)
{
  const int n = nelmt1/packsize;
  const int bytes_per_elmt = (int)(size/nelmt2/nelmt3/n);
  coproc_transposeOnHost3d(h_ptrSrc, h_ptrDest,
			   bytes_per_elmt, n, nelmt2, nelmt3);
}

/******************************************************************************/
extern "C" {

  //
  // 2d version
  //
  void FNAME(coproc_tmemcpyhosttohost2d)(const void **h_ptrSrc,
					 void **h_ptrDest,
					 __SIZET *size,
					 __INT *nelmt1, __INT *nelmt2,
					 __INT *packsize)
  {
    coproc_tmemcpyHostToHost2d(*h_ptrSrc, *h_ptrDest, *size,
			       *nelmt1, *nelmt2, *packsize);
  }

  //
  // 3d version
  //
  void FNAME(coproc_tmemcpyhosttohost3d)(const void **h_ptrSrc,
					 void **h_ptrDest,
					 __SIZET *size,
					 __INT *nelmt1, __INT *nelmt2,
					 __INT *nelmt3, __INT *packsize)
  {
    coproc_tmemcpyHostToHost3d(*h_ptrSrc, *h_ptrDest, *size,
			       *nelmt1, *nelmt2, *nelmt3, *packsize);
  }
}

/*******************************************************************************
 * Copy host memory to device memory synchroneously (using transposition)
 ******************************************************************************/

//
// 2d version
//
void coproc_tmemcpyHostToDevice2d(const void * __restrict__ h_ptrSrc, 
				  void * __restrict__ d_ptrDest,
				  size_t size, int nelmt1, int nelmt2,
				  int packsize)
{
  void *h_ptr;
  coproc_newMemoryOnHost(&h_ptr, size, cudaHostAllocDefault);
  
  const int n = nelmt1/packsize;
  const int bytes_per_elmt = (int)(size/nelmt2/n);
  coproc_transposeOnHost2d(h_ptrSrc, h_ptr, bytes_per_elmt, n, nelmt2);
  
  coproc_memcpyHostToDevice(h_ptr, d_ptrDest, size);
  coproc_freeMemoryOnHost(h_ptr);
}

//
// 3d version
//
void coproc_tmemcpyHostToDevice3d(const void * __restrict__ h_ptrSrc, 
				  void * __restrict__ d_ptrDest,
				  size_t size, int nelmt1, int nelmt2,
				  int nelmt3, int packsize)
{
  void *h_ptr;
  coproc_newMemoryOnHost(&h_ptr, size, cudaHostAllocDefault);
  
  const int n = nelmt1/packsize;
  const int bytes_per_elmt = (int)(size/nelmt2/nelmt3/n);
  coproc_transposeOnHost3d(h_ptrSrc, h_ptr, bytes_per_elmt, n, nelmt2, nelmt3);
  
  coproc_memcpyHostToDevice(h_ptr, d_ptrDest, size);
  coproc_freeMemoryOnHost(h_ptr);
}

/******************************************************************************/
extern "C" {

  //
  // 2d version
  //
  void FNAME(coproc_tmemcpyhosttodevice2d)(const void *h_ptrSrc,
					   void **d_ptrDest,
					   __SIZET *size,
					   __INT *nelmt1, __INT *nelmt2,
					   __INT *packsize)
  {
    coproc_tmemcpyHostToDevice2d(h_ptrSrc, *d_ptrDest, *size,
				 *nelmt1, *nelmt2, *packsize);
  }

  //
  // 3d version
  //
  void FNAME(coproc_tmemcpyhosttodevice3d)(const void *h_ptrSrc,
					   void **d_ptrDest,
					   __SIZET *size,
					   __INT *nelmt1, __INT *nelmt2,
					   __INT *nelmt3, __INT *packsize)
  {
    coproc_tmemcpyHostToDevice3d(h_ptrSrc, *d_ptrDest, *size,
				 *nelmt1, *nelmt2, *nelmt3, *packsize);
  }
}

/*******************************************************************************
 * Copy device memory to host memory synchroneously (using transposition)
 ******************************************************************************/

//
// 2d version
//
void coproc_tmemcpyDeviceToHost2d(const void * __restrict__ d_ptrSrc, 
				  void * __restrict__ h_ptrDest,
				  size_t size, int nelmt1, int nelmt2,
				  int packsize)
{
  void *d_ptr;
  coproc_newMemoryOnDevice(&d_ptr, size);
  coproc_clearMemoryOnDevice(d_ptr, size);
  
  const int n = nelmt1/packsize;
  const int bytes_per_elmt = (int)(size/nelmt2/n);
  coproc_transposeOnDevice2d(d_ptrSrc, d_ptr, bytes_per_elmt, n, nelmt2);
  
  coproc_memcpyDeviceToHost(d_ptr, h_ptrDest, size);
  coproc_freeMemoryOnDevice(d_ptr);
}

//
// 3d version
//
void coproc_tmemcpyDeviceToHost3d(const void * __restrict__ d_ptrSrc, 
				  void * __restrict__ h_ptrDest,
				  size_t size, int nelmt1, int nelmt2,
				  int nelmt3, int packsize)
{
  /*
  void *h_ptr;
  coproc_newMemoryOnHost(&h_ptr, size, cudaHostAllocDefault);
  
  const int n = nelmt1/packsize;
  const int bytes_per_elmt = (int)(size/nelmt2/nelmt3/n);
  coproc_transposeOnHost3d(h_ptrSrc, h_ptr, bytes_per_elmt, n, nelmt2, nelmt3);
  
  coproc_memcpyHostToDevice(h_ptr, d_ptrDest, size);
  coproc_freeMemoryOnHost(h_ptr);
  */
}

/******************************************************************************/
extern "C" {

  //
  // 2d version
  //
  void FNAME(coproc_tmemcpydevicetohost2d)(const void **d_ptrSrc,
					   void *h_ptrDest,
					   __SIZET *size,
					   __INT *nelmt1, __INT *nelmt2,
					   __INT *packsize)
  {
    coproc_tmemcpyDeviceToHost2d(*d_ptrSrc, h_ptrDest, *size,
				 *nelmt1, *nelmt2, *packsize);
  }

  //
  // 3d version
  //
  void FNAME(coproc_tmemcpydevicetohost3d)(const void **d_ptrSrc,
					   void *h_ptrDest,
					   __SIZET *size,
					   __INT *nelmt1, __INT *nelmt2,
					   __INT *nelmt3, __INT *packsize)
  {
    coproc_tmemcpyDeviceToHost3d(*d_ptrSrc, h_ptrDest, *size,
				 *nelmt1, *nelmt2, *nelmt3, *packsize);
  }
}

/*******************************************************************************
 * Copy device memory to device memory synchroneously (using transposition)
 ******************************************************************************/

//
// 2d version
//
void coproc_tmemcpyDeviceToDevice2d(const void * __restrict__ d_ptrSrc, 
				    void * __restrict__ d_ptrDest,
				    size_t size, int nelmt1, int nelmt2,
				    int packsize)
{
  
}

//
// 3d version
//
void coproc_tmemcpyDeviceToDevice3d(const void * __restrict__ d_ptrSrc, 
				    void * __restrict__ d_ptrDest,
				    size_t size, int nelmt1, int nelmt2,
				    int nelmt3, int packsize)
{
  
}

/******************************************************************************/
extern "C" {

  //
  // 2d version
  //
  void FNAME(coproc_tmemcpydevicetodevice2d)(const void **d_ptrSrc,
					     void **d_ptrDest,
					     __SIZET *size,
					     __INT *nelmt1, __INT *nelmt2,
					     __INT *packsize)
  {
    coproc_tmemcpyDeviceToDevice2d(*d_ptrSrc, *d_ptrDest, *size,
				   *nelmt1, *nelmt2, *packsize);
  }
  
  //
  // 3d version
  //
  void FNAME(coproc_tmemcpydevicetodevice3d)(const void **d_ptrSrc,
					     void **d_ptrDest,
					     __SIZET *size,
					     __INT *nelmt1, __INT *nelmt2,
					     __INT *nelmt3, __INT *packsize)
  {
    coproc_tmemcpyDeviceToDevice3d(*d_ptrSrc, *d_ptrDest, *size,
				   *nelmt1, *nelmt2, *nelmt3, *packsize);
  }
}

/*******************************************************************************
 * Combine two generic memory blocks in device memory
 ******************************************************************************/
template<typename T>
__global__ void combineOnDevice_knl(const T *ptrSrc1,
				    const T *ptrSrc2,
				    T *ptrDest,
				    size_t size)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx<size)
    {
      ptrDest[idx] = ptrSrc1[idx] + ptrSrc2[idx];
    }
}

template<>
__global__ void combineOnDevice_knl(const bool *ptrSrc1,
				    const bool *ptrSrc2,
				    bool *ptrDest,
				    size_t size)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx<size)
    {
      ptrDest[idx] = ptrSrc1[idx] || ptrSrc2[idx];
    }
}

/*******************************************************************************
 * Combine two single memory blocks in device memory
 ******************************************************************************/
void coproc_combineSingleOnDevice(const void *d_ptrSrc1,
				  const void *d_ptrSrc2,
				  void *d_ptrDest,
				  size_t size,
				  cudaStream_t stream)
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
  coproc_checkError("coproc_combineSingleOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combinesingleondevice)(const void **d_ptrSrc1,
					   const void **d_ptrSrc2,
					   void **d_ptrDest,
					   __SIZET *size,
					   __I64 *stream)
  {
    coproc_combineSingleOnDevice(*d_ptrSrc1, *d_ptrSrc2,
				 *d_ptrDest, *size,
				 (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two double memory blocks in device memory
 ******************************************************************************/
void coproc_combineDoubleOnDevice(const void *d_ptrSrc1,
				  const void *d_ptrSrc2,
				  void *d_ptrDest,
				  size_t size,
				  cudaStream_t stream)
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
  coproc_checkError("coproc_combineDoubleOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combinedoubleondevice)(const void **d_ptrSrc1,
					   const void **d_ptrSrc2,
					   void **d_ptrDest,
					   __SIZET *size,
					   __I64 *stream)
  {
    coproc_combineDoubleOnDevice(*d_ptrSrc1, *d_ptrSrc2,
				 *d_ptrDest, *size,
				 (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two quadrupel memory blocks in device memory
 ******************************************************************************/
void coproc_combineQuadOnDevice(const void *d_ptrSrc1,
				const void *d_ptrSrc2,
				void *d_ptrDest,
				size_t size,
				cudaStream_t stream)
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
  coproc_checkError("coproc_combineQuadOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combinequadondevice)(const void **d_ptrSrc1,
					 const void **d_ptrSrc2,
					 void **d_ptrDest,
					 __SIZET *size,
					 __I64 *stream)
  {
    coproc_combineQuadOnDevice(*d_ptrSrc1, *d_ptrSrc2,
			       *d_ptrDest, *size,
			       (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two integer memory blocks in device memory
 ******************************************************************************/
void coproc_combineIntegerOnDevice(const void *d_ptrSrc1,
				   const void *d_ptrSrc2,
				   void *d_ptrDest,
				   size_t size,
				   cudaStream_t stream)
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
  coproc_checkError("coproc_combineIntegerOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combineintegerondevice)(const void **d_ptrSrc1,
					    const void **d_ptrSrc2,
					    void **d_ptrDest,
					    __SIZET *size,
					    __I64 *stream)
  {
    coproc_combineIntegerOnDevice(*d_ptrSrc1, *d_ptrSrc2,
				  *d_ptrDest, *size,
				  (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int8 memory blocks in device memory
 ******************************************************************************/
void coproc_combineInt8OnDevice(const void *d_ptrSrc1,
				const void *d_ptrSrc2,
				void *d_ptrDest,
				size_t size,
				cudaStream_t stream)
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
  coproc_checkError("coproc_combineInt8OnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combineint8ondevice)(const void **d_ptrSrc1,
					 const void **d_ptrSrc2,
					 void **d_ptrDest,
					 __SIZET *size,
					 __I64 *stream)
  {
    coproc_combineInt8OnDevice(*d_ptrSrc1, *d_ptrSrc2,
			       *d_ptrDest, *size,
			       (cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int16 memory blocks in device memory
 ******************************************************************************/
void coproc_combineInt16OnDevice(const void *d_ptrSrc1,
				 const void *d_ptrSrc2,
				 void *d_ptrDest,
				 size_t size,
				 cudaStream_t stream)
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
  coproc_checkError("coproc_combineInt16OnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combineint16ondevice)(const void **d_ptrSrc1,
					  const void **d_ptrSrc2,
					  void **d_ptrDest,
					  __SIZET *size,
					  __I64 *stream)
  {
    coproc_combineInt16OnDevice(*d_ptrSrc1, *d_ptrSrc2,
				*d_ptrDest, *size,
				(cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int32 memory blocks in device memory
 ******************************************************************************/
void coproc_combineInt32OnDevice(const void *d_ptrSrc1,
				 const void *d_ptrSrc2,
				 void *d_ptrDest,
				 size_t size,
				 cudaStream_t stream)
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
  coproc_checkError("coproc_combineInt32OnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combineint32ondevice)(const void **d_ptrSrc1,
					  const void **d_ptrSrc2,
					  void **d_ptrDest,
					  __SIZET *size,
					  __I64 *stream)
  {
    coproc_combineInt32OnDevice(*d_ptrSrc1, *d_ptrSrc2,
				*d_ptrDest, *size,
				(cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two int64 memory blocks in device memory
 ******************************************************************************/
void coproc_combineInt64OnDevice(const void *d_ptrSrc1,
				 const void *d_ptrSrc2,
				 void *d_ptrDest,
				 size_t size,
				 cudaStream_t stream)
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
  coproc_checkError("coproc_combineInt64OnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combineint64ondevice)(const void **d_ptrSrc1,
					  const void **d_ptrSrc2,
					  void **d_ptrDest,
					  __SIZET *size,
					  __I64 *stream)
  {
    coproc_combineInt64OnDevice(*d_ptrSrc1, *d_ptrSrc2,
				*d_ptrDest, *size,
				(cudaStream_t)(*stream));
  }
}

/*******************************************************************************
 * Combine two logical memory blocks in device memory
 ******************************************************************************/
void coproc_combineLogicalOnDevice(const void *d_ptrSrc1,
				   const void *d_ptrSrc2,
				   void *d_ptrDest,
				   size_t size,
				   cudaStream_t stream)
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
  coproc_checkError("coproc_combineLogicalOnDevice");
}

/******************************************************************************/
extern "C" {
  void FNAME(coproc_combinelogicalondevice)(const void **d_ptrSrc1,
					    const void **d_ptrSrc2,
					    void **d_ptrDest,
					    __SIZET *size,
					    __I64 *stream)
  {
    coproc_combineLogicalOnDevice(*d_ptrSrc1, *d_ptrSrc2,
				  *d_ptrDest, *size,
				  (cudaStream_t)(*stream));
  }
}
