#ifndef _COPROC_TRANSPOSE_H_
#define _COPROC_TRANSPOSE_H_

/*#############################################################################
 ******************************************************************************
 * <name> coproc_storage_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This header file provides declarations for the file coproc_transpose_cuda.cu
 * </purpose>
 *
 *#############################################################################
 */

extern "C" {
  void coproc_transposeOnHost2d(const void *, void *, int, int, int);
  void coproc_transposeOnHost3d(const void *, void *, int, int, int, int);

  void coproc_transposeOnDevice2d(const void *, void *, int, int, int, cudaStream_t=0);
  void coproc_transposeOnDevice3d(const void *, void *, int, int, int, int, cudaStream_t=0);
}

#endif
