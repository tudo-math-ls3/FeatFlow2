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
  int coproc_newDoubleOnDevice(unsigned long * h_Ddata, int * n);
  int FNAME(coproc_newdoubleondevice)(unsigned long * h_Ddata, int * n);
  int coproc_newSingleOnDevice(unsigned long * h_Fdata, int * n);
  int FNAME(coproc_newsingleondevice)(unsigned long * h_Fdata, int * n);
  int coproc_newIntOnDevice(unsigned long * h_Idata, int * n);
  int FNAME(coproc_newintondevice)(unsigned long * h_Idata, int * n);

  int coproc_freeDoubleOnDevice(unsigned long * h_Ddata);
  int FNAME(coproc_freedoubleondevice)(unsigned long * h_Ddata);
  int coproc_freeSingleOnDevice(unsigned long * h_Fdata);
  int FNAME(coproc_freesingleondevice)(unsigned long * h_Fdata);
  int coproc_freeIntOnDevice(unsigned long * h_Idata);
  int FNAME(coproc_freeintondevice)(unsigned long * h_Idata);

  int coproc_clearDoubleOnDevice(unsigned long * h_Ddata, int * n);
  int FNAME(coproc_cleardoubleondevice)(unsigned long * h_Ddata, int * n);
  int coproc_clearSingleOnDevice(unsigned long * h_Fdata, int * n);
  int FNAME(coproc_clearsingleondevice)(unsigned long * h_Fdata, int * n);
  int coproc_clearIntOnDevice(unsigned long * h_Idata, int * n);
  int FNAME(coproc_clearintondevice)(unsigned long * h_Idata, int * n);

  int coproc_copyDoubleToDevice(double * Ddata, unsigned long * h_Ddata, int * n);
  int FNAME(coproc_copydoubletodevice)(double * Ddata, unsigned long * h_Ddata, int * n);
  int coproc_copySingleToDevice(float * Fdata, unsigned long * h_Fdata, int * n);
  int FNAME(coproc_copysingletodevice)(float * Fdata, unsigned long * h_Fdata, int * n);
  int coproc_copyIntToDevice(int * Idata, unsigned long * h_Idata, int * n);
  int FNAME(coproc_copyinttodevice)(int * Idata, unsigned long * h_Idata, int * n);

  int coproc_copyDoubleToHost(unsigned long * h_Ddata, double * Ddata, int * n);
  int FNAME(coproc_copydoubletohost)(unsigned long * h_Ddata, double * Ddata, int * n);
  int coproc_copySingleToHost(unsigned long * h_Fdata, float * Fdata, int * n);
  int FNAME(coproc_copysingletohost)(unsigned long * h_Fdata, float * Fdata, int * n);
  int coproc_copyIntToHost(unsigned long * h_Idata, int * Idata, int * n);
  int FNAME(coproc_copyinttohost)(unsigned long * h_Idata, int * Idata, int * n);
}

#endif
