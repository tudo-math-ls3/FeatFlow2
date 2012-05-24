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
  int coproc_transposeOnHost2d(void *, void *, int, int, int);
  int coproc_transposeOnHost3d(void *, void *, int, int, int, int);
}

#endif
