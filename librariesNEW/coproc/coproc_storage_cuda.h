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
  int coproc_malloc(void **ptr, size_t size);

  int coproc_free(void *ptr);  

  int coproc_newMemoryOnHost(void **h_ptr,
			     size_t size);

  int coproc_freeMemoryOnHost(void *h_ptr);

  int coproc_clearMemoryOnHost(void *h_ptr,
			       size_t size);
  
  int coproc_newMemoryOnDevice(void **d_ptr,
			       size_t size);
  
  int coproc_freeMemoryOnDevice(void *d_ptr);
  
  int coproc_clearMemoryOnDevice(void *d_ptr,
				 size_t size);

  int coproc_memcpyHostToDevice(void *h_ptrSrc, 
				void *d_ptrDest,
				size_t size);

  int coproc_memcpyHostToDeviceAsync(void *h_ptrSrc, 
				     void *d_ptrDest,
				     size_t size,
				     cudaStream_t stream=0);
  
  int coproc_memcpyDeviceToHost(void *d_ptrSrc,
				void *h_ptrDest,
				size_t size);
  
  int coproc_memcpyDeviceToHostAsync(void *d_ptrSrc,
				     void *h_ptrDest,
				     size_t size,
				     cudaStream_t stream=0);
  
  int coproc_memcpyDeviceToDevice(void *d_ptrSrc,
				  void *d_ptrDest,
				  size_t size);

  int coproc_memcpyDeviceToDeviceAsync(void *d_ptrSrc,
				       void *d_ptrDest,
				       size_t size,
				       cudaStream_t stream=0);
  
  int coproc_combinesingleOnDevice(void *d_ptrSrc1,
				   void *d_ptrSrc2,
				   void *d_ptrDest,
				   size_t size);
  
  int coproc_combineDoubleOnDevice(void *d_ptrSrc1,
				   void *d_ptrSrc2,
				   void *d_ptrDest,
				   size_t size);
  
  int coproc_combineQuadOnDevice(void *d_ptrSrc1,
				 void *d_ptrSrc2,
				 void *d_ptrDest,
				 size_t size);

  int coproc_combineIntegerOnDevice(void *d_ptrSrc1,
				    void *d_ptrSrc2,
				    void *d_ptrDest,
				    size_t size);
  
  int coproc_combineInt8OnDevice(void *d_ptrSrc1,
				 void *d_ptrSrc2,
				 void *d_ptrDest,
				 size_t size);

  int coproc_combineInt16OnDevice(void *d_ptrSrc1,
				  void *d_ptrSrc2,
				  void *d_ptrDest,
				  size_t size);
  
  int coproc_combineInt32OnDevice(void *d_ptrSrc1,
				  void *d_ptrSrc2,
				  void *d_ptrDest,
				  size_t size);
  
  int coproc_combineInt64OnDevice(void *d_ptrSrc1,
				  void *d_ptrSrc2,
				  void *d_ptrDest,
				  size_t size);
  
  int coproc_combineLogicalOnDevice(void *d_ptrSrc1,
				    void *d_ptrSrc2,
				    void *d_ptrDest,
				    size_t size);
}

#endif
