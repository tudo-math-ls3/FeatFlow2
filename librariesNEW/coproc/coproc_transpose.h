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
  int coproc_transpose2d(void **ptrSrc, void **ptrDest,
			 int bytes_per_elmt, int n1, int n2);
  int coproc_transpose3d(void **ptrSrc, void **ptrDest,
			 int bytes_per_elmt, int n1, int n2, int n3);
}

#endif
