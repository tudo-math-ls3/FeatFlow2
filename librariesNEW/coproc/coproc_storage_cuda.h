#ifndef _COPROC_STORAGE_CUDA_H_
#define _COPROC_STORAGE_CUDA_H_

/*#############################################################################
 ******************************************************************************
 * <name> coproc_storage_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This header file provides declarations for the file coproc_storage_cuda.cu
 * </purpose>
 *
 *#############################################################################
 */

extern "C"
{
  void coproc_malloc(void**, size_t);
  void coproc_newMemoryOnHost(void**, size_t);
  void coproc_newMemoryOnDevice(void**, size_t); 

  void coproc_free(void**);
  void coproc_freeMemoryOnHost(void*);
  void coproc_freeMemoryOnDevice(void*);

  void coproc_clearMemoryOnHost(void*, size_t);
  void coproc_clearMemoryOnDevice(void*, size_t);

  void coproc_memcpyHostToHost(const void*, void*, size_t);
  void coproc_memcpyHostToHostAsync(const void*, void*, size_t, cudaStream_t=0);
  void coproc_memcpyHostToDevice(const void*, void*, size_t);
  void coproc_memcpyHostToDeviceAsync(const void*, void*, size_t, cudaStream_t=0);
  void coproc_memcpyDeviceToHost(const void*, void*, size_t);
  void coproc_memcpyDeviceToHostAsync(const void*, void*, size_t, cudaStream_t=0);
  void coproc_memcpyDeviceToDevice(const void*, void*, size_t);
  void coproc_memcpyDeviceToDeviceAsync(const void*, void*, size_t, cudaStream_t=0);

  void coproc_tmemcpyHostToHost2d(const void*, void*, size_t, int, int, int);
  void coproc_tmemcpyHostToDevice2d(const void*, void*, size_t, int, int, int);
  void coproc_tmemcpyDeviceToHost2d(const void*, void*, size_t, int, int, int);
  void coproc_tmemcpyDeviceToDevice2d(const void*, void*, size_t, int, int, int);
  void coproc_tmemcpyHostToHost3d(const void*, void*, size_t, int, int, int, int);
  void coproc_tmemcpyHostToDevice3d(const void*, void*, size_t, int, int, int, int);
  void coproc_tmemcpyDeviceToHost3d(const void*, void*, size_t, int, int, int, int);
  void coproc_tmemcpyDeviceToDevice3d(const void*, void*, size_t, int, int, int, int);
  
  void coproc_combinesingleOnDevice(const void*, const void*, void*, size_t);
  void coproc_combineDoubleOnDevice(const void*, const void*, void*, size_t);
  void coproc_combineQuadOnDevice(const void*, const void*, void*, size_t);
  void coproc_combineIntegerOnDevice(const void*, const void*, void*, size_t);
  void coproc_combineInt8OnDevice(const void*, const void*, void*, size_t);
  void coproc_combineInt16OnDevice(const void*, const void*, void*, size_t);
  void coproc_combineInt32OnDevice(const void*, const void*, void*, size_t);
  void coproc_combineInt64OnDevice(const void*, const void*, void*, size_t);
  void coproc_combineLogicalOnDevice(const void*, const void*, void*, size_t);
}

#endif
