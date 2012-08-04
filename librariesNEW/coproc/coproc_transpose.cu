/*#############################################################################
******************************************************************************
* <name> coproc_transpose </name>
******************************************************************************
*
* <purpose>
* This file provides routines for transposing multi-dimensional arrays.
* </purpose>
*
*#############################################################################
*/

#include <stdint.h>
#include <iostream>
#include <algorithm>
#include "coproc_core.h"
#include "coproc_transpose.h"

using namespace std;


/*#############################################################################
 * Transpose two-dimensional arrays in host memory
 *#############################################################################
 */

/******************************************************************************
 * 2D-transpose in host memory:
 * Each element is assumed to have a single data item
 ******************************************************************************/
template<typename T>
void coproc_transposeOnHost2d(const void *__restrict__ h_ptrSrc, 
			      void *__restrict__ h_ptrDest,
			      int nelmt1, int nelmt2)
{
  T *src = (T*)(h_ptrSrc);
  T *dest = (T*)(h_ptrDest);
  
  for (int i2=0; i2<nelmt2; ++i2) {
    for (int i1=0; i1<nelmt1; ++i1) {
      dest[nelmt2*i1+i2] = src[nelmt1*i2+i1];
    }
  }
}

/******************************************************************************
 * 2D-transpose in host memory:
 * Each element is assumed to have num_elmt data items given at run-time
 ******************************************************************************/
template<typename T>
void coproc_transposeOnHost2d(int num_elmt,
			      const void *__restrict__ h_ptrSrc, 
			      void *__restrict__ h_ptrDest,
			      int nelmt1, int nelmt2)
{
  T *src = (T*)(h_ptrSrc);
  T *dest = (T*)(h_ptrDest);
  
  for (int i2=0; i2<nelmt2; ++i2) {
    for (int i1=0; i1<nelmt1; ++i1) {
      for (int ielmt=0; ielmt<num_elmt; ++ielmt) {
	dest[num_elmt*(nelmt2*i1+i2)+ielmt] = src[num_elmt*(nelmt1*i2+i1)+ielmt];
      }
    }
  }
}

/******************************************************************************
 * 2D-transpose in host memory:
 * Each element is assumed to have num_elmt data items known at compile-time
 ******************************************************************************/
template<typename T, int num_elmt>
void coproc_transposeOnHost2d(const void *__restrict__ h_ptrSrc,
			      void *__restrict__ h_ptrDest,
			      int nelmt1, int nelmt2)
{
  T *src = (T*)(h_ptrSrc);
  T *dest = (T*)(h_ptrDest);
  
  for (int i2=0; i2<nelmt2; ++i2) {
    for (int i1=0; i1<nelmt1; ++i1) {
#pragma unroll
      for (int ielmt=0; ielmt<num_elmt; ++ielmt) {
	dest[num_elmt*(nelmt2*i1+i2)+ielmt] = src[num_elmt*(nelmt1*i2+i1)+ielmt];
      }
    }
  }
}

/******************************************************************************
 * Wrapper routine in C++
 ******************************************************************************/
void coproc_transposeOnHost2d(const void *__restrict__ h_ptrSrc,
			      void *__restrict__ h_ptrDest,
			      int bytes_per_elmt, int nelmt1, int nelmt2)
{
  switch (bytes_per_elmt) {
  
    //
    // Check for 1,2,4,8 bytes
    //

  case sizeof(uint8_t): // 1 byte
    coproc_transposeOnHost2d<uint8_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;

  case sizeof(uint16_t): // 2 byte
    coproc_transposeOnHost2d<uint16_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;
    
  case sizeof(uint32_t): // 4 byte
    coproc_transposeOnHost2d<uint32_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;
    
  case sizeof(uint64_t): // 8 byte
    coproc_transposeOnHost2d<uint64_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;

    //
    // Check for special cases bytes=2^n
    //
    
  case 16:
    coproc_transposeOnHost2d<uint64_t,2>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;
    
  case 32:
    coproc_transposeOnHost2d<uint64_t,4>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;
    
  case 64:
    coproc_transposeOnHost2d<uint64_t,8>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);    
    break;
    
  case 128:
    coproc_transposeOnHost2d<uint64_t,16>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;

  case 256:
    coproc_transposeOnHost2d<uint64_t,32>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;

  case 512:
    coproc_transposeOnHost2d<uint64_t,64>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;

  case 1024:
    coproc_transposeOnHost2d<uint64_t,128>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    break;

  default:
    //
    // Default: check if we have a multiple of 8,4,2 bytes
    //
    if ((bytes_per_elmt&7) == 0) {
      coproc_transposeOnHost2d<uint64_t>((bytes_per_elmt>>3),
					 h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    }
    else if ((bytes_per_elmt&3) == 0) {
      coproc_transposeOnHost2d<uint32_t>((bytes_per_elmt>>2),
					 h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    }
    else if ((bytes_per_elmt&1) == 0) {
      coproc_transposeOnHost2d<uint16_t>((bytes_per_elmt>>1),
					 h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    }
    else {
      coproc_transposeOnHost2d<uint8_t>(bytes_per_elmt,
					h_ptrSrc, h_ptrDest, nelmt1, nelmt2);
    }
    break;
  }
  coproc_checkError("coproc_transposeOnHost2d");
}

/******************************************************************************
 * Wrapper routines to be called from Fortran
 ******************************************************************************/
extern "C" {
  void FNAME(coproc_transposeonhost2d)(const void *__restrict__ h_ptrSrc,
				       void *__restrict__ h_ptrDest,
				       __INT *bytes_per_elmt,
				       __INT *nelmt1, __INT *nelmt2)
  {
    coproc_transposeOnHost2d(h_ptrSrc, h_ptrDest, *bytes_per_elmt,
			     *nelmt1, *nelmt2);
  }
};



/*#############################################################################
 * Transpose three-dimensional arrays in host memory
 *#############################################################################
 */



/******************************************************************************
 * Each element is assumed to have a single data item
 ******************************************************************************/
template<typename T>
void coproc_transposeOnHost3d(const void *__restrict__ h_ptrSrc,
			      void *__restrict__ h_ptrDest,
			      int nelmt1, int nelmt2, int nelmt3)
{
  T *src = (T*)(h_ptrSrc);
  T *dest = (T*)(h_ptrDest);
  
  for (int i3=0; i3<nelmt3; ++i3) {
    for (int i2=0; i2<nelmt2; ++i2) {
      for (int i1=0; i1<nelmt1; ++i1) {
	dest[nelmt2*nelmt3*i1+nelmt3*i2+i3] = src[nelmt1*nelmt2*i3+nelmt1*i2+i1];
      }
    }
  }
}

/******************************************************************************
 * Each element is assumed to have num_elmt data items given at run-time
 ******************************************************************************/
template<typename T>
void coproc_transposeOnHost3d(int num_elmt,
			      const void *__restrict__ h_ptrSrc, 
			      void *__restrict__ h_ptrDest,
			      int nelmt1, int nelmt2, int nelmt3)
{
  T *src = (T*)(h_ptrSrc);
  T *dest = (T*)(h_ptrDest);
  
  for (int i3=0; i3<nelmt3; ++i3) {
    for (int i2=0; i2<nelmt2; ++i2) {
      for (int i1=0; i1<nelmt1; ++i1) {
	for (int ielmt=0; ielmt<num_elmt; ++ielmt) {
	  dest[num_elmt*(nelmt2*nelmt3*i1+nelmt3*i2+i3)+ielmt] =
	    src[num_elmt*(nelmt1*nelmt2*i3+nelmt1*i2+i1)+ielmt];
	}
      }
    }
  }
}

/******************************************************************************
 * Each element is assumed to have num_elmt data items known at compile-time
 ******************************************************************************/
template<typename T, int num_elmt>
void coproc_transposeOnHost3d(const void *__restrict__ h_ptrSrc,
			      void *__restrict__ h_ptrDest,
			      int nelmt1, int nelmt2, int nelmt3)
{
  T *src = (T*)(h_ptrSrc);
  T *dest = (T*)(h_ptrDest);
  
  for (int i3=0; i3<nelmt3; ++i3) {
    for (int i2=0; i2<nelmt2; ++i2) {
      for (int i1=0; i1<nelmt1; ++i1) {
#pragma unroll
	for (int ielmt=0; ielmt<num_elmt; ++ielmt) {
	  dest[num_elmt*(nelmt2*nelmt3*i1+nelmt3*i2+i3)+ielmt] =
	    src[num_elmt*(nelmt1*nelmt2*i3+nelmt1*i2+i1)+ielmt];
	}
      }
    }
  }
}

/******************************************************************************
 * Wrapper routine in C++
 ******************************************************************************/
void coproc_transposeOnHost3d(const void *__restrict__ h_ptrSrc,
			      void *__restrict__ h_ptrDest,
			      int bytes_per_elmt, int nelmt1, int nelmt2, int nelmt3)
{
  switch (bytes_per_elmt) {
  
    //
    // Check for 1,2,4,8 bytes
    //

  case sizeof(uint8_t): // 1 byte
    coproc_transposeOnHost3d<uint8_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

  case sizeof(uint16_t): // 2 byte
    coproc_transposeOnHost3d<uint16_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;
    
  case sizeof(uint32_t): // 4 byte
    coproc_transposeOnHost3d<uint32_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;
    
  case sizeof(uint64_t): // 8 byte
    coproc_transposeOnHost3d<uint64_t>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

    //
    // Check for special cases bytes=2^n
    //
    
  case 16:
    coproc_transposeOnHost3d<uint64_t,2>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

  case 32:
    coproc_transposeOnHost3d<uint64_t,4>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

  case 64:
    coproc_transposeOnHost3d<uint64_t,8>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

  case 128:
    coproc_transposeOnHost3d<uint64_t,16>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;
    
  case 256:
    coproc_transposeOnHost3d<uint64_t,32>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

  case 512:
    coproc_transposeOnHost3d<uint64_t,64>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

  case 1024:
    coproc_transposeOnHost3d<uint64_t,128>(h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    break;

  default:
    //
    // Default: check if we have a multiple of 8,4,2 bytes
    //
    if ((bytes_per_elmt&7) == 0) {
      coproc_transposeOnHost3d<uint64_t>((bytes_per_elmt>>3),
					 h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    }
    else if ((bytes_per_elmt&3) == 0) {
      coproc_transposeOnHost3d<uint32_t>((bytes_per_elmt>>2),
					 h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    }
    else if ((bytes_per_elmt&1) == 0) {
      coproc_transposeOnHost3d<uint16_t>((bytes_per_elmt>>1),
					 h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    }
    else {
      coproc_transposeOnHost3d<uint8_t>(bytes_per_elmt,
					h_ptrSrc, h_ptrDest, nelmt1, nelmt2, nelmt3);
    }
    break;
  }
  coproc_checkError("coproc_transposeOnHost3d");
}

/******************************************************************************
 * Wrapper routine to be called from Fortran
 ******************************************************************************/
extern "C" {
  void FNAME(coproc_transposeonhost3d)(const void *__restrict__ h_ptrSrc,
				       void *__restrict__ h_ptrDest,
				       __INT *bytes_per_elmt,
				       __INT *nelmt1, __INT *nelmt2, __INT *nelmt3)
  {
    coproc_transposeOnHost3d(h_ptrSrc, h_ptrDest, *bytes_per_elmt,
			     *nelmt1, *nelmt2, *nelmt3);
  }
};



/*#############################################################################
 * Transpose two-dimensional arrays in device memory
 *#############################################################################
 */

/******************************************************************************
 * Naive 2D-transpose CUDA kernel:
 * Each element is assumed to have a single data item.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_naive_knl(T *Src, T *Dest,
				      int nelmt1, int nelmt2) {

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int index_in  = xIndex + nelmt1 * yIndex;
  int index_out = yIndex + nelmt2 * xIndex;

  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
    if (index_out+i < nelmt1*nelmt2)
      Dest[index_out+i] = Src[index_in+i*nelmt1];
  }
}

/******************************************************************************
 * Naive 2D-transpose CUDA kernel:
 * Each element is assumed to have a single data item.
 * Remark: same as before but this kernel can handle multiple items per thread.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_naive_mult_knl(T *Src, T *Dest,
					   int nelmt1, int nelmt2,
					   int nitem1, int nitem2) {

  for (int k1=0; k1<nitem1; k1++) {
    for (int k2=0; k2<nitem2; k2++) {
      
      int xIndex = (k1*gridDim.x + blockIdx.x) * TILE_DIM + threadIdx.x;
      int yIndex = (k2*gridDim.y + blockIdx.y) * TILE_DIM + threadIdx.y;
      
      int index_in  = xIndex + nelmt1 * yIndex;
      int index_out = yIndex + nelmt2 * xIndex;
      
      for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
	if (index_out+i < nelmt1*nelmt2)
	  Dest[index_out+i] = Src[index_in+i*nelmt1];
      }
    }
  }
}

/******************************************************************************
 * Coalesced 2D-transpose CUDA kernel:
 * Each element is assumed to have a single data item.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_knl(T *Src, T *Dest,
				int nelmt1, int nelmt2) {
  
  __shared__ T tile[TILE_DIM][TILE_DIM+1]; 
  
  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index  = xIndex + yIndex * nelmt1;
  
  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) 
    tile[threadIdx.y+i][threadIdx.x] = Src[index+i*nelmt1]; 
  
  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  index  = xIndex + yIndex * nelmt2;
  
  __syncthreads();
  
  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) 
    Dest[index+i*nelmt2] = tile[threadIdx.x][threadIdx.y+i];
}

/******************************************************************************
 * Coalesced 2D-transpose CUDA kernel:
 * Each element is assumed to have a single data item.
 * Remark: same as before but this kernel can handle multiple items per thread.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_mult_knl(T *Src, T *Dest,
				     int nelmt1, int nelmt2,
				     int nitem1, int nitem2) {
  
  __shared__ T tile[TILE_DIM][TILE_DIM+1]; 
  
  for (int k1=0; k1<nitem1; k1++) {
    for (int k2=0; k2<nitem2; k2++) {
      
      int xIndex = (k1*gridDim.x + blockIdx.x) * TILE_DIM + threadIdx.x;
      int yIndex = (k2*gridDim.y + blockIdx.y) * TILE_DIM + threadIdx.y;
      int index  = xIndex + yIndex * nelmt1;
      
      for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) 
	tile[threadIdx.y+i][threadIdx.x] = Src[index+i*nelmt1]; 
      
      xIndex = (k2*gridDim.y + blockIdx.y) * TILE_DIM + threadIdx.x;
      yIndex = (k1*gridDim.x + blockIdx.x) * TILE_DIM + threadIdx.y;
      index  = xIndex + yIndex * nelmt2;
      
      __syncthreads();
      
      for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) 
	Dest[index+i*nelmt2] = tile[threadIdx.x][threadIdx.y+i];

      __syncthreads();
    }
  }
}

/******************************************************************************
 * Coalesced 2D-transpose CUDA kernel with offset and bound checking:
 * Each element is assumed to have a single data item.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_offset_knl(T *Src, T *Dest,
				       int nelmt1, int nelmt2,
				       int offset1, int offset2) {
  
  __shared__ T tile[TILE_DIM][TILE_DIM+1]; 
  
  int xIndex = offset1 + blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = offset2 + blockIdx.y * TILE_DIM + threadIdx.y;
  int index  = xIndex + yIndex * nelmt1;
  
  if (xIndex < nelmt1) {
    for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
      if (index+i*nelmt1 < nelmt1*nelmt2) {
	tile[threadIdx.y+i][threadIdx.x] = Src[index+i*nelmt1];
      }
    }
  }
  
  xIndex = offset2 + blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = offset1 + blockIdx.x * TILE_DIM + threadIdx.y;
  index  = xIndex + yIndex * nelmt2;
  
  __syncthreads();
  
  if (xIndex < nelmt2) {
    for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
      if (index+i*nelmt2 < nelmt1*nelmt2) {
	Dest[index+i*nelmt2] = tile[threadIdx.x][threadIdx.y+i];
      }
    }
  }
}

/******************************************************************************
 * Coalesced 2D-transpose CUDA kernel with offset and bound checking:
 * Each element is assumed to have a single data item.
 * Remark: same as before but this kernel can handle multiple items per thread.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_mult_offset_knl(T *Src, T *Dest,
					    int nelmt1, int nelmt2,
					    int nitem1, int nitem2,
					    int offset1, int offset2) {
  
  __shared__ T tile[TILE_DIM][TILE_DIM+1]; 
  
  for (int k1=0; k1<nitem1; k1++) {
    for (int k2=0; k2<nitem2; k2++) {
      
      int xIndex = offset1 + (k1*gridDim.x + blockIdx.x) * TILE_DIM + threadIdx.x;
      int yIndex = offset2 + (k2*gridDim.y + blockIdx.y) * TILE_DIM + threadIdx.y;
      int index  = xIndex + yIndex * nelmt1;
      
      if (xIndex < nelmt1) {
	for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
	  if (index+i*nelmt1 < nelmt1*nelmt2) {
	    tile[threadIdx.y+i][threadIdx.x] = Src[index+i*nelmt1];
	  }
	}
      }
      
      xIndex = offset2 + (k2*gridDim.y + blockIdx.y) * TILE_DIM + threadIdx.x;
      yIndex = offset1 + (k1*gridDim.x + blockIdx.x) * TILE_DIM + threadIdx.y;
      index  = xIndex + yIndex * nelmt2;
      
      __syncthreads();
      
      if (xIndex < nelmt2) {
	for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
	  if (index+i*nelmt2 < nelmt1*nelmt2) {
	    Dest[index+i*nelmt2] = tile[threadIdx.x][threadIdx.y+i];
	  }
	}
      }

      __syncthreads();
    }
  }
}

/******************************************************************************
 * Coalesced 2D-transpose CUDA kernel:
 * Each element is assumed to have a single data item
 ******************************************************************************/

template<int TILE_DIM, int FIXED_DIM, class T>
__global__ void transpose2d_rect_knl(T *Src, T *Dest,
				     int nelmt1, int nelmt2) {
  
  if (FIXED_DIM == 1) {
    // First dimension is fixed, i.e. not larger than TILE_DIM
    __shared__ T tile[TILE_DIM*(TILE_DIM-1)];
    
    // All threads read nelmt1*TILE_DIM data items from the source
    // array from contiguous memory positions
    int index = blockIdx.x * TILE_DIM * nelmt1;

    for (int i=0; i<nelmt1; i++) {
      int tid = threadIdx.x + i * TILE_DIM;
      tile[tid] = Src[index+tid];
    }
    
    __syncthreads();
    
    // All threads write TILE_DIM contiguous items from the same
    // column of the shared memory block and store them in the same
    // row of the destination array; Repeat until all nelmt1<TILE_DIM
    // rows have been processed.
    index = blockIdx.x * TILE_DIM + threadIdx.x;
    
    for (int i=0; i<nelmt1; i++)
      Dest[index + i*nelmt2] = tile[threadIdx.x*nelmt2+i];
  }
  else if (FIXED_DIM == 2) {
    // Second dimension is fixed, i.e. not larger than TILE_DIM
    __shared__ T tile[TILE_DIM*(TILE_DIM-1)];
    
    // All threads read TILE_DIM contiguous items from the same row
    // and store them in the same column of the shared memory block;
    // Repeat until all nelmt2<TILE_DIM rows have been processed.
    int index = blockIdx.x * TILE_DIM + threadIdx.x;
    
    for (int i=0; i<nelmt2; i++)
      tile[threadIdx.x*nelmt2+i] = Src[index + i*nelmt1];
    
    __syncthreads();
    
    // All threads write nelmt2*TILE_DIM data items from the shared
    // memory block into the destination array at contiguous positions
    index = blockIdx.x * TILE_DIM * nelmt2;
    
    for (int i=0; i<nelmt2; i++) {
      int tid = threadIdx.x + i * TILE_DIM;
      Dest[index+tid] = tile[tid];
    }
  }
}

/******************************************************************************
 * 2D-transpose in device memory:
 * Each element is assumed to have a single data item
 ******************************************************************************/
template<typename T>
void coproc_transposeOnDevice2d(const void *__restrict__ d_ptrSrc, 
				void *__restrict__ d_ptrDest,
				int nelmt1, int nelmt2, cudaStream_t stream)
{
  T *ptrSrc = (T*)(d_ptrSrc);
  T *ptrDest = (T*)(d_ptrDest);

  const cudaDeviceProp *devProp = coproc_getCurrentDeviceProp();
  const int TILE_DIM   = 32;
  const int BLOCK_ROWS = 2;

  if (nelmt1 >= TILE_DIM && nelmt2 >= TILE_DIM) {
    //
    // Transpose matrix using TILE_DIM x TILE_DIM tiles
    //

    // Compute number of tiles and grid size
    const int n1 = nelmt1/TILE_DIM;
    const int n2 = nelmt2/TILE_DIM;
    const int m1 = (n1+devProp->maxGridSize[0]-1)/(devProp->maxGridSize[0]);
    const int m2 = (n2+devProp->maxGridSize[1]-1)/(devProp->maxGridSize[1]);
  
    // Coalesced transpose using quadratic tiles
    dim3 grid(n1/m1, n2/m2);
    dim3 threads(TILE_DIM,BLOCK_ROWS);
    
    if (m1*m2 == 1) {
      transpose2d_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
	(ptrSrc, ptrDest, nelmt1, nelmt2);
    }
    else {
      transpose2d_mult_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
	(ptrSrc, ptrDest, nelmt1, nelmt2, m1, m2);
    }

    const int l1 = m1*grid.x*TILE_DIM;
    const int l2 = m2*grid.y*TILE_DIM;
    
    if (nelmt1 > l1) {
      const int nn1 = (nelmt1-l1+TILE_DIM-1)/TILE_DIM;
      const int nn2 = (nelmt2+TILE_DIM-1)/TILE_DIM;
      const int mm1 = (nn1+devProp->maxGridSize[0]-1)/(devProp->maxGridSize[0]);
      const int mm2 = (nn2+devProp->maxGridSize[1]-1)/(devProp->maxGridSize[1]);
      
      // Coalesced transpose of last few (less than TILE_DIM) rows
      grid = dim3((nn1+mm1-1)/mm1, (nn2+mm2-1)/mm2);

      if (mm1*mm2 == 1) {
	transpose2d_offset_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
	  (ptrSrc, ptrDest, nelmt1, nelmt2, l1, 0);
      }
      else {
	transpose2d_mult_offset_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
	  (ptrSrc, ptrDest, nelmt1, nelmt2, mm1, mm2, l1, 0);
      }
    }
    
    if (nelmt2 > l2) {
      const int nn1 = (nelmt1+TILE_DIM-1)/TILE_DIM;
      const int nn2 = (nelmt2-l2+TILE_DIM-1)/TILE_DIM;
      const int mm1 = (nn1+devProp->maxGridSize[0]-1)/(devProp->maxGridSize[0]);
      const int mm2 = (nn2+devProp->maxGridSize[1]-1)/(devProp->maxGridSize[1]);

      // Coalesced transpose of last few (less than TILE_DIM) columns
      grid = dim3((nn1+mm1-1)/mm1, (nn2+mm2-1)/mm2);
      
      if (mm1*mm2 == 1) {
	transpose2d_offset_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
	  (ptrSrc, ptrDest, nelmt1, nelmt2, 0, l2);
      }
      else {
	transpose2d_mult_offset_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
	  (ptrSrc, ptrDest, nelmt1, nelmt2, mm1, mm2, 0, l2);
      }
    }
  }
  else if (nelmt1 > TILE_DIM) {
    //
    // Transpose array with first dimension larger than TILE_DIM
    //
    
    // Compute number of tiles and grid size
    const int n1 = nelmt1/TILE_DIM;
    const int m1 = (n1+devProp->maxGridSize[0]-1)/(devProp->maxGridSize[0]);
    
    dim3 grid(n1/m1,1);
    dim3 threads(TILE_DIM,1);
    transpose2d_rect_knl<TILE_DIM,2><<<grid, threads, 0, stream>>>
      (ptrSrc, ptrDest, nelmt1, nelmt2);
    
    const int l1 = m1*grid.x*TILE_DIM;
    if (nelmt1 > l1) {
      const int nn1 = (nelmt1-l1+TILE_DIM-1)/TILE_DIM;
      const int mm1 = (nn1+devProp->maxGridSize[0]-1)/(devProp->maxGridSize[0]);
      
      // Coalesced transpose of last few (less than TILE_DIM) rows
      grid = dim3((nn1+mm1-1)/mm1,1);
      dim3 threads(TILE_DIM,BLOCK_ROWS);
      transpose2d_offset_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
	(ptrSrc, ptrDest, nelmt1, nelmt2, l1, 0);
    }
  }
  else if (nelmt2 > TILE_DIM) {
    //
    // Transpose array with second dimension larger than TILE_DIM
    //
    // Compute number of tiles and grid size
    const int n2 = nelmt2/TILE_DIM;
    const int m2 = (n2+devProp->maxGridSize[1]-1)/(devProp->maxGridSize[1]);
    
    dim3 grid(1,n2/m2);
    dim3 threads(TILE_DIM,1);
    transpose2d_rect_knl<TILE_DIM,1><<<grid, threads, 0, stream>>>
      (ptrSrc, ptrDest, nelmt1, nelmt2);
    
    const int l2 = m2*grid.y*TILE_DIM;
    if (nelmt2 > l2) {
      const int nn2 = (nelmt2-l2+TILE_DIM-1)/TILE_DIM;
      const int mm2 = (nn2+devProp->maxGridSize[1]-1)/(devProp->maxGridSize[1]);
      
      // Coalesced transpose of last few (less than TILE_DIM) rows
      grid = dim3(1,(nn2+mm2-1)/mm2);
      dim3 threads(TILE_DIM,BLOCK_ROWS);
      transpose2d_offset_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
    	(ptrSrc, ptrDest, nelmt1, nelmt2, 0, l2);
    }
  }
  else {
    //
    // Transpose array with both dimensions smaller than TILE_DIM
    //
    dim3 grid(1,1);
    dim3 threads(TILE_DIM,BLOCK_ROWS);
    transpose2d_offset_knl<TILE_DIM,BLOCK_ROWS><<<grid, threads, 0, stream>>>
      (ptrSrc, ptrDest, nelmt1, nelmt2, 0, 0);
  }
}

/******************************************************************************
 * Naive 2D-transpose CUDA kernel:
 * Each element is assumed to have num_elmt data items given at run-time.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_naive_knl(int num_elmt,
				      T *Src, T *Dest,
				      int nelmt1, int nelmt2) {

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int index_in  = xIndex + nelmt1 * yIndex;
  //int index_out = yIndex + nelmt2 * (xIndex);

  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
    if (index_in+i*nelmt1 < nelmt1*nelmt2) {
      int xIdx = index_in/nelmt1;
      int yIdx = index_in-xIdx*nelmt1;
      
      xIndex = xIdx/num_elmt;
      yIndex = yIndex*num_elmt + xIdx%num_elmt;
      int index_out = yIndex + nelmt2 * xIndex;
      
      Dest[index_out] = Src[index_in+i*nelmt1];
    }
  }
}

/******************************************************************************
 * Naive 2D-transpose CUDA kernel:
 * Each element is assumed to have num_elmt data items given at run-time.
 * Remark: same as before but this kernel can handle multiple items per thread.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, class T>
__global__ void transpose2d_naive_mult_knl(int num_elmt,
					   T *Src, T *Dest,
					   int nelmt1, int nelmt2,
					   int nitem1, int nitem2) {

  for (int k1=0; k1<nitem1; k1++) {
    for (int k2=0; k2<nitem2; k2++) {
      
      int xIndex = (k1*gridDim.x + blockIdx.x) * TILE_DIM + threadIdx.x;
      int yIndex = (k2*gridDim.y + blockIdx.y) * TILE_DIM + threadIdx.y;
      
      int index_in  = xIndex + nelmt1 * yIndex;
      int index_out = yIndex + nelmt2 * xIndex;
      
      for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
	if (index_out+i < nelmt1*nelmt2)
	  Dest[index_out+i] = Src[index_in+i*nelmt1];
      }
    }
  }
}

/******************************************************************************
 * 2D-transpose in device memory:
 * Each element is assumed to have num_elmt data items given at run-time
 ******************************************************************************/
template<typename T>
void coproc_transposeOnDevice2d(int num_elmt,
				const void *__restrict__ d_ptrSrc, 
				void *__restrict__ d_ptrDest,
				int nelmt1, int nelmt2, cudaStream_t stream)
{
  T *ptrSrc = (T*)(d_ptrSrc);
  T *ptrDest = (T*)(d_ptrDest);
  
  const cudaDeviceProp *devProp = coproc_getCurrentDeviceProp();
  const int TILE_DIM   = 32;
  const int BLOCK_ROWS = 2;

  // Compute number of tiles and grid size
  const int n1 = (nelmt1+TILE_DIM-1)/TILE_DIM;
  const int n2 = (nelmt2+TILE_DIM-1)/TILE_DIM;
  const int m1 = (n1+devProp->maxGridSize[0]-1)/(devProp->maxGridSize[0]);
  const int m2 = (n2+devProp->maxGridSize[1]-1)/(devProp->maxGridSize[1]);
  
  dim3 grid((n1+m1-1)/m1, (n2+m2-1)/m2);
  dim3 threads(TILE_DIM,BLOCK_ROWS);
  
  // Naive transpose of rectangular arrays
  if (m1*m2 == 1) {
    transpose2d_naive_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
      (num_elmt, ptrSrc, ptrDest, nelmt1, nelmt2);
  }
  else {
    transpose2d_naive_mult_knl<TILE_DIM,BLOCK_ROWS,T><<<grid, threads, 0, stream>>>
      (num_elmt, ptrSrc, ptrDest, nelmt1, nelmt2, m1, m2);
  }
}

/******************************************************************************
 * Naive 2D-transpose CUDA kernel:
 * Each element is assumed to have num_elmt data items given at run-time.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, int num_elmt, class T>
__global__ void transpose2d_naive_knl(T *Src, T *Dest,
				      int nelmt1, int nelmt2) {

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int index_in  = xIndex + nelmt1 * yIndex;
  //int index_out = yIndex + nelmt2 * (xIndex);

  for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
    int index = index_in+i*nelmt1;
    if (index < nelmt1*nelmt2) {
      int yIdx = index/nelmt1;
      int xIdx = index-yIdx*nelmt1;
      
      xIndex = xIdx/num_elmt;
      yIndex = yIndex*num_elmt + xIdx%num_elmt;
      int index_out = yIdx + nelmt2 * xIdx;
      
      Dest[index_out] = Src[index];
    }
  }
}

/******************************************************************************
 * Naive 2D-transpose CUDA kernel:
 * Each element is assumed to have num_elmt data items given at run-time.
 * Remark: same as before but this kernel can handle multiple items per thread.
 ******************************************************************************/

template<int TILE_DIM, int BLOCK_ROWS, int num_elmt, class T>
__global__ void transpose2d_naive_mult_knl(T *Src, T *Dest,
					   int nelmt1, int nelmt2,
					   int nitem1, int nitem2) {

  for (int k1=0; k1<nitem1; k1++) {
    for (int k2=0; k2<nitem2; k2++) {
      
      int xIndex = (k1*gridDim.x + blockIdx.x) * TILE_DIM + threadIdx.x;
      int yIndex = (k2*gridDim.y + blockIdx.y) * TILE_DIM + threadIdx.y;
      
      int index_in  = xIndex + nelmt1 * yIndex;
      int index_out = yIndex + nelmt2 * xIndex;
      
      for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
	if ((index_out+i < nelmt1*nelmt2) &&
	    (index_in+i*nelmt1 < nelmt1*nelmt2))
	  Dest[index_out+i] = Src[index_in+i*nelmt1];
      }
    }
  }
}

/******************************************************************************
 * 2D-transpose in device memory:
 * Each element is assumed to have num_elmt data items known at compile-time
 ******************************************************************************/
template<typename T, int num_elmt>
void coproc_transposeOnDevice2d(const void *__restrict__ d_ptrSrc,
				void *__restrict__ d_ptrDest,
				int nelmt1, int nelmt2, cudaStream_t stream)
{
  T *ptrSrc = (T*)(d_ptrSrc);
  T *ptrDest = (T*)(d_ptrDest);
  
  const cudaDeviceProp *devProp = coproc_getCurrentDeviceProp();
  const int TILE_DIM   = 32;
  const int BLOCK_ROWS = 2;

  // Compute number of tiles and grid size
  const int n1 = (nelmt1*num_elmt+TILE_DIM-1)/TILE_DIM;
  const int n2 = (nelmt2+TILE_DIM-1)/TILE_DIM;
  const int m1 = (n1+devProp->maxGridSize[0]-1)/(devProp->maxGridSize[0]);
  const int m2 = (n2+devProp->maxGridSize[1]-1)/(devProp->maxGridSize[1]);
  
  dim3 grid((n1+m1-1)/m1, (n2+m2-1)/m2);
  dim3 threads(TILE_DIM,BLOCK_ROWS);
  
  cout << nelmt1 << "," << num_elmt << ":" << nelmt2 << endl;

  // Naive transpose of rectangular arrays
  if (m1*m2 == 1) {
    transpose2d_naive_knl<TILE_DIM,BLOCK_ROWS,num_elmt,T><<<grid, threads, 0, stream>>>
      (ptrSrc, ptrDest, nelmt1*num_elmt, nelmt2);
  }
  else {
    transpose2d_naive_mult_knl<TILE_DIM,BLOCK_ROWS,num_elmt,T><<<grid, threads, 0, stream>>>
      (ptrSrc, ptrDest, nelmt1*num_elmt, nelmt2, m1, m2);
  }
}

/******************************************************************************
 * Wrapper routine in C++
 ******************************************************************************/
void coproc_transposeOnDevice2d(const void *__restrict__ d_ptrSrc,
				void *__restrict__ d_ptrDest,
				int bytes_per_elmt, int nelmt1, int nelmt2,
				cudaStream_t stream)
{
  cudaEvent_t start,stop;
  float inTime;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,stream);

  switch (bytes_per_elmt) {
    
    //
    // Check for 1,2,4,8 bytes
    //
    
  case sizeof(uchar1): // 1 byte
    coproc_transposeOnDevice2d<uchar1>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case sizeof(uchar2): // 2 byte
    coproc_transposeOnDevice2d<uchar2>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case sizeof(uchar3): // 3 byte
    coproc_transposeOnDevice2d<uchar3>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case sizeof(float1): // 4 byte
    coproc_transposeOnDevice2d<float1>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case sizeof(float2): // 8 byte
    coproc_transposeOnDevice2d<float2>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case sizeof(float3): // 12 byte
    coproc_transposeOnDevice2d<float3>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case sizeof(float4): // 16 byte
    coproc_transposeOnDevice2d<float4>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
    //
    // Multiples of float3 : n*12 byte
    //
    
  case 24:
    cout << "THIS IS IT" << endl;
    coproc_transposeOnDevice2d<float3,2>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 60:
    coproc_transposeOnDevice2d<float3,5>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case 72:
    coproc_transposeOnDevice2d<float3,6>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 84:
    coproc_transposeOnDevice2d<float3,7>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 108:
    coproc_transposeOnDevice2d<float3,9>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 120:
    coproc_transposeOnDevice2d<float3,10>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 132:
    coproc_transposeOnDevice2d<float3,11>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 156:
    coproc_transposeOnDevice2d<float3,13>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 168:
    coproc_transposeOnDevice2d<float3,14>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 180:
    coproc_transposeOnDevice2d<float3,15>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 204:
    coproc_transposeOnDevice2d<float3,17>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 216:
    coproc_transposeOnDevice2d<float3,18>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case 228:
    coproc_transposeOnDevice2d<float3,19>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 252:
    coproc_transposeOnDevice2d<float3,21>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
    //
    // Multiples of float4 : n * 16 byte
    //
    
  case 32:
    coproc_transposeOnDevice2d<float4,2>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 48:
    coproc_transposeOnDevice2d<float4,3>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 64:
    coproc_transposeOnDevice2d<float4,4>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 80:
    coproc_transposeOnDevice2d<float4,5>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 96:
    coproc_transposeOnDevice2d<float4,6>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 112:
    coproc_transposeOnDevice2d<float4,7>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case 128:
    coproc_transposeOnDevice2d<float4,8>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 144:
    coproc_transposeOnDevice2d<float4,9>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 160:
    coproc_transposeOnDevice2d<float4,10>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case 176:
    coproc_transposeOnDevice2d<float4,11>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 192:
    coproc_transposeOnDevice2d<float4,12>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 208:
    coproc_transposeOnDevice2d<float4,13>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case 224:
    coproc_transposeOnDevice2d<float4,14>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 240:
    coproc_transposeOnDevice2d<float4,15>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;

  case 256:
    coproc_transposeOnDevice2d<float4,16>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case 512:
    coproc_transposeOnDevice2d<float4,32>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
  case 1024:
    coproc_transposeOnDevice2d<float4,64>(d_ptrSrc, d_ptrDest, nelmt1, nelmt2, stream);
    break;
    
    
  default:
    //
    // Default: check if we have a multiple of 16,8,4,2 bytes
    //
    if ((bytes_per_elmt&15) == 0) {
      coproc_transposeOnDevice2d<float4>((bytes_per_elmt>>4), d_ptrSrc, d_ptrDest,
					 nelmt1, nelmt2, stream);
    }
    if ((bytes_per_elmt&7) == 0) {
      coproc_transposeOnDevice2d<float2>((bytes_per_elmt>>3), d_ptrSrc, d_ptrDest,
					 nelmt1, nelmt2, stream);
    }
    else if ((bytes_per_elmt&3) == 0) {
      coproc_transposeOnDevice2d<float1>((bytes_per_elmt>>2), d_ptrSrc, d_ptrDest,
					 nelmt1, nelmt2, stream);
    }
    else if ((bytes_per_elmt&1) == 0) {
      coproc_transposeOnDevice2d<char2>((bytes_per_elmt>>1), d_ptrSrc, d_ptrDest,
					nelmt1, nelmt2, stream);
    }
    else {
      coproc_transposeOnDevice2d<char1>(bytes_per_elmt, d_ptrSrc, d_ptrDest,
					nelmt1, nelmt2, stream);
    }
    break;
  }
  coproc_checkError("coproc_transposeOnDevice2d");

  cudaEventRecord(stop,stream);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&inTime, start, stop);
  
  cout << "Bandwidth: " << 2.0f*1000.0f*nelmt1*nelmt2*bytes_per_elmt/(1024*1024*1024)/(inTime) << "GB/s" << endl;
}

/******************************************************************************
 * Wrapper routines to be called from Fortran
 ******************************************************************************/
extern "C" {
  void FNAME(coproc_transposeondevice2d)(const void *__restrict__ d_ptrSrc,
					 void *__restrict__ d_ptrDest,
					 __INT *bytes_per_elmt,
					 __INT *nelmt1, __INT *nelmt2,
					 __I64 *stream)
  {
    coproc_transposeOnDevice2d(d_ptrSrc, d_ptrDest, *bytes_per_elmt,
			       *nelmt1, *nelmt2, (cudaStream_t)(*stream));
  }
};



/*#############################################################################
 * Transpose three-dimensional arrays in device memory
 *#############################################################################
 */

/******************************************************************************
 * Naive 3D-transpose CUDA kernel:
 * Each element is assumed to have a single data item
 ******************************************************************************/

template<class T>
__global__ void transpose3d_knl(T *Src, T *Dest,
				int nelmt1, int nelmt2, int nelmt3,
				int offset1=0, int offset2=0, int offset3=0) {
  
  int xIndex = offset1 + blockIdx.x + threadIdx.x;
  int yIndex = offset2 + blockIdx.y + threadIdx.y;
  int zIndex = offset3 + blockIdx.z + threadIdx.z;
  
  int index_in  = xIndex + nelmt1 * yIndex + nelmt1 * nelmt2 * zIndex;
  int index_out = zIndex + nelmt3 * yIndex + nelmt2 * nelmt3 * xIndex;
  
  Dest[index_out] = Src[index_in];
}

/******************************************************************************
 * Naive 3D-transpose CUDA kernel:
 * Each element is assumed to have num_elmt data items given at run-time
 ******************************************************************************/

template<class T>
__global__ void transpose3d_knl(int num_elmt,
				T *Src, T *Dest,
				int nelmt1, int nelmt2, int nelmt3,
				int offset1=0, int offset2=0, int offset3=0) {
  
  int xIndex = offset1 + blockIdx.x + threadIdx.x;
  int yIndex = offset2 + blockIdx.y + threadIdx.y;
  int zIndex = offset3 + blockIdx.z + threadIdx.z;
  
  int index_in  = (xIndex + nelmt1 * yIndex + nelmt1 * nelmt2 * zIndex) * num_elmt;
  int index_out = (zIndex + nelmt3 * yIndex + nelmt2 * nelmt3 * xIndex) * num_elmt;
  
  for (int ielmt=0; ielmt<num_elmt; ++ielmt)
    Dest[index_out+ielmt] = Src[index_in+ielmt];
}

/******************************************************************************
 * Naive 3D-transpose CUDA kernel:
 * Each element is assumed to have num_elmt data items known at compile-time
 ******************************************************************************/

template<class T, int num_elmt>
__global__ void transpose3d_knl(T *Src, T *Dest,
				int nelmt1, int nelmt2, int nelmt3,
				int offset1=0, int offset2=0, int offset3=0) {
  
  int xIndex = offset1 + blockIdx.x + threadIdx.x;
  int yIndex = offset2 + blockIdx.y + threadIdx.y;
  int zIndex = offset3 + blockIdx.z + threadIdx.z;
  
  int index_in  = (xIndex + nelmt1 * yIndex + nelmt1 * nelmt2 * zIndex) * num_elmt;
  int index_out = (zIndex + nelmt3 * yIndex + nelmt2 * nelmt3 * xIndex) * num_elmt;
  
  for (int ielmt=0; ielmt<num_elmt; ++ielmt)
    Dest[index_out+ielmt] = Src[index_in+ielmt];
}
