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
  int coproc_newMemoryOnHost(unsigned long * p_MemoryBlock,
			     unsigned long * imemBytes);
  int FNAME(coproc_newmemoryonhost)(unsigned long * p_MemoryBlock,
				    unsigned long * imemBytes);

  int coproc_freeMemoryOnHost(unsigned long * p_MemoryBlock);
  int FNAME(coproc_freememoryonhost)(unsigned long * p_MemoryBlock);

  int coproc_newMemoryOnDevice(unsigned long * p_MemoryBlock,
			       unsigned long * imemBytes);
  int FNAME(coproc_newmemoryondevice)(unsigned long * p_MemoryBlock,
				      unsigned long * imemBytes);

  int coproc_freeMemoryOnDevice(unsigned long * p_MemoryBlock);
  int FNAME(coproc_freememoryondevice)(unsigned long * p_MemoryBlock);

  int coproc_clearMemoryOnDevice(unsigned long * p_MemoryBlock,
				 unsigned long * imemBytes);
  int FNAME(coproc_clearmemoryondevice)(unsigned long * p_MemoryBlock,
					unsigned long * imemBytes);

  int coproc_copyMemoryHostToDevice(unsigned long * p_MemoryBlockOnHost, 
				    unsigned long * p_MemoryBlockOnDevice,
				    unsigned long * imemBytes);
  int FNAME(coproc_copymemoryhosttodevice)(unsigned long * p_MemoryBlockOnHost, 
					   unsigned long * p_MemoryBlockOnDevice,
					   unsigned long * imemBytes);

  int coproc_copyMemoryDeviceToHost(unsigned long * p_MemoryBlockOnDevice,
				    unsigned long * p_MemoryBlockOnHost,
				    unsigned long * imemBytes);
  int FNAME(coproc_copymemorydevicetohost)(unsigned long * p_MemoryBlockOnDevice,
					   unsigned long * p_MemoryBlockOnHost,
					   unsigned long * imemBytes);

  int coproc_copyMemoryDeviceToDevice(unsigned long * p_MemoryBlockSrc,
				      unsigned long * p_MemoryBlockDest,
				      unsigned long * imemBytes);
  int FNAME(coproc_copymemorydevicetodevice)(unsigned long * p_MemoryBlockSrc,
					     unsigned long * p_MemoryBlockDest,
					     unsigned long * imemBytes);

  int coproc_addsingleOnDevice(unsigned long * p_MemoryBlock1,
			       unsigned long * p_MemoryBlock2,
			       unsigned long * p_MemoryBlockDest,
			       unsigned long * imemBytes);
  int FNAME(coproc_addsingleondevice)(unsigned long * p_MemoryBlock1,
				      unsigned long * p_MemoryBlock2,
				      unsigned long * p_MemoryBlockDest,
				      unsigned long * imemBytes);

  int coproc_addDoubleOnDevice(unsigned long * p_MemoryBlock1,
			       unsigned long * p_MemoryBlock2,
			       unsigned long * p_MemoryBlockDest,
			       unsigned long * imemBytes);
  int FNAME(coproc_adddoubleondevice)(unsigned long * p_MemoryBlock1,
				      unsigned long * p_MemoryBlock2,
				      unsigned long * p_MemoryBlockDest,
				      unsigned long * imemBytes);

  int coproc_addIntegerOnDevice(unsigned long * p_MemoryBlock1,
				unsigned long * p_MemoryBlock2,
				unsigned long * p_MemoryBlockDest,
				unsigned long * imemBytes);
  int FNAME(coproc_addintegerondevice)(unsigned long * p_MemoryBlock1,
				       unsigned long * p_MemoryBlock2,
				       unsigned long * p_MemoryBlockDest,
				       unsigned long * imemBytes);

  int coproc_addLogicalOnDevice(unsigned long * p_MemoryBlock1,
				unsigned long * p_MemoryBlock2,
				unsigned long * p_MemoryBlockDest,
				unsigned long * imemBytes);
  int FNAME(coproc_addlogicalondevice)(unsigned long * p_MemoryBlock1,
				       unsigned long * p_MemoryBlock2,
				       unsigned long * p_MemoryBlockDest,
				       unsigned long * imemBytes);
}

#endif
