/*#############################################################################
 ******************************************************************************
 * <name> hydro_calcOperator2d_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This file provides CUDA kernels to compute the operator for the low-order
 * scheme in 2D using different types if artificial viscosities.
 * </purpose>
 *
 *#############################################################################/
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "coproc_core.h"
#include "coproc_storage_cuda.h"

#define LANGUAGE LANGUAGE_C
#include "../../../../../kernel/System/idxmanager.h"

#define HYDRO_NDIM 2
#include "hydro.h"

/*******************************************************************************
 * External C functions which can be called from the Fortran code
 *******************************************************************************/

extern "C"
{
  int FNAME(hydro_calcmatdiag2d_cuda)(unsigned long * h_CoeffsAtDiag,
				      unsigned long * h_IdiagList,
				      unsigned long * h_Dx,
				      unsigned long * h_Da[NVAR2D][NVAR2D],
				      double * scale,
				      bool * bclear,
				      int * nblocks,
				      int * na,
				      int * neq,
				      int * nvar,
				      int * ncoeff);

  int FNAME(hydro_calcmatscdiss2d_cuda)(unsigned long * h_CoeffsAtEdge,
					unsigned long * h_IedgeList,
					unsigned long * h_Dx,
					unsigned long * h_Da[NVAR2D][NVAR2D],
					double * scale,
					bool * bclear,
					int * nblocks,
					int * na,
					int * neq,
					int * nvar,
					int * nedge,
					int * ncoeff,
					int * nedges,
					int * iedgeset);
}


/*******************************************************************************
 * CUDA kernels
 *******************************************************************************/

template <int isystemformat>
struct gather_DataAtNode
{ 
};

/*******************************************************************************/

template <>
struct gather_DataAtNode<SYSTEM_SEGREGATED>
{
  template <typename T>
  __host__ __device__ inline
  static void eval (T * DataAtNode, 
		    T * Dx,
		    int ieq, 
		    int neq)
  {
    // Solution vector is stored in interleaved format
    IDX2(DataAtNode,1,1,NVAR2D,2) = IDX2_REVERSE(Dx,1,ieq,NVAR2D,neq);
    IDX2(DataAtNode,2,1,NVAR2D,2) = IDX2_REVERSE(Dx,2,ieq,NVAR2D,neq);
    IDX2(DataAtNode,3,1,NVAR2D,2) = IDX2_REVERSE(Dx,3,ieq,NVAR2D,neq);
    IDX2(DataAtNode,4,1,NVAR2D,2) = IDX2_REVERSE(Dx,4,ieq,NVAR2D,neq);
  }
};

/*******************************************************************************/

template <>
struct gather_DataAtNode<SYSTEM_ALLCOUPLED>
{
  template <typename T>
  __host__ __device__ inline
  static void eval (T * DataAtNode, 
		    T * Dx,
		    int ieq, 
		    int neq)
  {
    // Solution vector is stored in block format
    IDX2(DataAtNode,1,1,NVAR2D,2) = IDX2_FORWARD(Dx,1,ieq,NVAR2D,neq);
    IDX2(DataAtNode,2,1,NVAR2D,2) = IDX2_FORWARD(Dx,2,ieq,NVAR2D,neq);
    IDX2(DataAtNode,3,1,NVAR2D,2) = IDX2_FORWARD(Dx,3,ieq,NVAR2D,neq);
    IDX2(DataAtNode,4,1,NVAR2D,2) = IDX2_FORWARD(Dx,4,ieq,NVAR2D,neq);
  }
};

/*******************************************************************************/

template <int isystemformat>
struct gather_DataAtEdge
{ 
};

/*******************************************************************************/

template <>
struct gather_DataAtEdge<SYSTEM_SEGREGATED>
{
  template <typename T>
  __host__ __device__ inline
  static void eval (T * DataAtEdge, 
		    T * Dx,
		    int i, 
		    int j, 
		    int neq)
  {
    // Solution vector is stored in interleaved format
    IDX2(DataAtEdge,1,1,NVAR2D,2) = IDX2_REVERSE(Dx,1,i,NVAR2D,neq);
    IDX2(DataAtEdge,2,1,NVAR2D,2) = IDX2_REVERSE(Dx,2,i,NVAR2D,neq);
    IDX2(DataAtEdge,3,1,NVAR2D,2) = IDX2_REVERSE(Dx,3,i,NVAR2D,neq);
    IDX2(DataAtEdge,4,1,NVAR2D,2) = IDX2_REVERSE(Dx,4,i,NVAR2D,neq);
    
    IDX2(DataAtEdge,1,2,NVAR2D,2) = IDX2_REVERSE(Dx,1,j,NVAR2D,neq);
    IDX2(DataAtEdge,2,2,NVAR2D,2) = IDX2_REVERSE(Dx,2,j,NVAR2D,neq);
    IDX2(DataAtEdge,3,2,NVAR2D,2) = IDX2_REVERSE(Dx,3,j,NVAR2D,neq);
    IDX2(DataAtEdge,4,2,NVAR2D,2) = IDX2_REVERSE(Dx,4,j,NVAR2D,neq);
  }
};

/*******************************************************************************/

template <>
struct gather_DataAtEdge<SYSTEM_ALLCOUPLED>
{
  template <typename T>
  __host__ __device__ inline
  static void eval (T * DataAtEdge,	
		    T * Dx,
		    int i, 
		    int j, 
		    int neq)
  {
    // Solution vector is stored in block format
    IDX2(DataAtEdge,1,1,NVAR2D,2) = IDX2_FORWARD(Dx,1,i,NVAR2D,neq);
    IDX2(DataAtEdge,2,1,NVAR2D,2) = IDX2_FORWARD(Dx,2,i,NVAR2D,neq);
    IDX2(DataAtEdge,3,1,NVAR2D,2) = IDX2_FORWARD(Dx,3,i,NVAR2D,neq);
    IDX2(DataAtEdge,4,1,NVAR2D,2) = IDX2_FORWARD(Dx,4,i,NVAR2D,neq);
    
    IDX2(DataAtEdge,1,2,NVAR2D,2) = IDX2_FORWARD(Dx,1,j,NVAR2D,neq);
    IDX2(DataAtEdge,2,2,NVAR2D,2) = IDX2_FORWARD(Dx,2,j,NVAR2D,neq);
    IDX2(DataAtEdge,3,2,NVAR2D,2) = IDX2_FORWARD(Dx,3,j,NVAR2D,neq);
    IDX2(DataAtEdge,4,2,NVAR2D,2) = IDX2_FORWARD(Dx,4,j,NVAR2D,neq);
  }
};

/*******************************************************************************/

template <int isystemformat,bool bclear>
struct scatter_CoefficientsAtNode
{
};

/*******************************************************************************/

template <int isystemformat,bool bclear>
struct scatter_CoefficientsAtEdge
{
};

/*******************************************************************************/

template <typename T, int isystemformat, int idissipationtype, bool bclear>
__global__ void hydro_calcMatDiag2d_knl(T * CoeffsAtDiag,
					int * IdiagList,
					T * Dx,
					T (* Da)[NVAR2D][NVAR2D],
					T scale,
					int na,
					int neq,
					int nvar,
					int ncoeff)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx<neq)
    {
    }
}

/*******************************************************************************
 * External C functions which can be called from the Fortran code
 *******************************************************************************/

template <typename T>
inline
int hydro_calcMatDiag2d_cuda(unsigned long * h_CoeffsAtDiag,
			     unsigned long * h_IdiagList,
			     unsigned long * h_Dx,
			     unsigned long * h_Da[NVAR2D][NVAR2D],
			     T * scale,
			     bool * bclear,
			     int * nblocks,
			     int * na,
			     int * neq,
			     int * nvar,
			     int * ncoeff)
{
  T * d_Dx = (T*)(*h_Dx);
  T * d_CoeffsAtDiag = (T*)(*h_CoeffsAtDiag);
  int * d_IdiagList = (int*)(*h_IdiagList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((*neq)/(double)(block.x));
  
  if (*nblocks == 1) {
    //    printf("%lu\n",d_Da);
    printf("Hello 1\n");
  } else {
    printf("Hello 2\n");
  }
  coproc_checkErrors("hydro_calcMatDiag2d_cuda");
  return 0;
};

/*******************************************************************************/

#ifdef HAS_CUDADOUBLEPREC
int FNAME(hydro_calcmatdiag2d_cuda)(unsigned long * h_CoeffsAtDiag,
				    unsigned long * h_IdiagList,
				    unsigned long * h_Dx,
				    unsigned long * h_Da[NVAR2D][NVAR2D],
				    double * scale,
				    bool * bclear,
				    int * nblocks,
				    int * na,
				    int * neq,
				    int * nvar,
				    int * ncoeff)
{
  return hydro_calcMatDiag2d_cuda(h_CoeffsAtDiag, h_IdiagList, h_Dx, h_Da, scale,
				  bclear, nblocks, na, neq, nvar, ncoeff);
}
#else
int FNAME(hydro_calcmatdiag2d_cuda)(unsigned long * h_CoeffsAtDiag,
				    unsigned long * h_IdiagList,
				    unsigned long * h_Dx,
				    unsigned long * h_Da[NVAR2D][NVAR2D],
				    float * scale,
				    bool * bclear,
				    int * nblocks,
				    int * na,
				    int * neq,
				    int * nvar,
				    int * ncoeff)
{
  return hydro_calcMatDiag2d_cuda(h_CoeffsAtDiag, h_IdiagList, h_Dx, h_Da, scale,
				  bclear, nblocks, na, neq, nvar, ncoeff);
}
#endif

/*******************************************************************************/

template <typename T>
inline
int hydro_calcMatScDiss2d_cuda(unsigned long * h_CoeffsAtEdge,
			       unsigned long * h_IedgeList,
			       unsigned long * h_Dx,
			       unsigned long * h_Da[NVAR2D][NVAR2D],
			       T * scale, bool * bclear,
			       int * nblocks, int * na,
			       int * neq, int * nvar,
			       int * nedge, int * ncoeff,
			       int * nedges, int * iedgeset)
{
  coproc_checkErrors("hydro_calcMatScDiss2d_cuda");
  return 0;
}

/*******************************************************************************/

#ifdef HAS_CUDADOUBLEPREC
int FNAME(hydro_calcmatscdiss2d_cuda)(unsigned long * h_CoeffsAtEdge,
				      unsigned long * h_IedgeList,
				      unsigned long * h_Dx,
				      unsigned long * h_Da[NVAR2D][NVAR2D],
				      double * scale, bool * bclear,
				      int * nblocks, int * na,
				      int * neq, int * nvar,
				      int * nedge, int * ncoeff,
				      int * nedges, int * iedgeset)
{
  return hydro_calcMatScDiss2d_cuda(h_CoeffsAtEdge, h_IedgeList, h_Dx, h_Da, scale,
				    bclear, nblocks, na, neq, nvar, nedge, ncoeff,
				    nedges, iedgeset);
}
#else
int FNAME(hydro_calcmatscdiss2d_cuda)(unsigned long * h_CoeffsAtEdge,
				      unsigned long * h_IedgeList,
				      unsigned long * h_Dx,
				      unsigned long * h_Da[NVAR2D][NVAR2D],
				      float * scale, bool * bclear,
				      int * nblocks, int * na,
				      int * neq, int * nvar,
				      int * nedge, int * ncoeff,
				      int * nedges, int * iedgeset)
{
  return hydro_calcMatScDiss2d_cuda(h_CoeffsAtEdge, h_IedgeList, h_Dx, h_Da, scale,
				    bclear, nblocks, na, neq, nvar, nedge, ncoeff,
				    nedges, iedgeset);
}
#endif
