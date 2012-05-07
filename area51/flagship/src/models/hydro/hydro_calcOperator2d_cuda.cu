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
#include <cmath>
#include <cfloat>
#include <iostream>
#include "coproc_core.h"
#include "coproc_storage_cuda.h"

#define LANGUAGE LANGUAGE_C
#include "../../../../../kernel/System/idxmanager.h"

#define HYDRO_NDIM 2
#include "hydro.h"

#ifdef ENABLE_COPROC_SHMEM
#define SHMEM_IDX idx
#define SHMEM_BLOCKSIZE blockDim.x
#else
#define SHMEM_IDX 1
#define SHMEM_BLOCKSIZE 1
#endif

// Memory pool in constant device memory
__device__ __constant__ __SIZET constMemPool[NVAR2D*NVAR2D];

namespace hydro2d_cuda
{

  /*****************************************************************************
   * CUDA kernels for hydrydynamic model in 2D
   ****************************************************************************/

  /*****************************************************************************
   * Thus CUDA kernel collects the nodal solution data from the global
   * solution vector.
   ****************************************************************************/

  template <int isystemformat>
  struct gather_DataAtNode
  { 
  };

  /*****************************************************************************
   * Input:  solution vector Dx stored in interleaved format
   * Output: DataAtNode vector
   ****************************************************************************/

  template <>
  struct gather_DataAtNode<SYSTEM_SCALAR>
  {
    template <typename TdSrc,
	      typename TdDest,
	      typename Ti>
    __device__ inline
    static void eval (TdDest *DataAtNode, 
		      TdSrc *Dx,
		      Ti ieq, 
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in interleaved format
      IDX2(DataAtNode,1,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,1,ieq,NVAR2D,neq);
      IDX2(DataAtNode,2,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,2,ieq,NVAR2D,neq);
      IDX2(DataAtNode,3,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,3,ieq,NVAR2D,neq);
      IDX2(DataAtNode,4,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,4,ieq,NVAR2D,neq);
    }
  };

  /*****************************************************************************
   * Input:  solution vector Dx stored in block format
   * Output: DataAtNode vector
   ****************************************************************************/

  template <>
  struct gather_DataAtNode<SYSTEM_BLOCK>
  {
    template <typename TdSrc,
	      typename TdDest,
	      typename Ti>
    __device__ inline
    static void eval (TdDest *DataAtNode,
		      TdSrc *Dx,
		      Ti ieq, 
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in block format
      IDX2(DataAtNode,1,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,1,ieq,NVAR2D,neq);
      IDX2(DataAtNode,2,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,2,ieq,NVAR2D,neq);
      IDX2(DataAtNode,3,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,3,ieq,NVAR2D,neq);
      IDX2(DataAtNode,4,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,4,ieq,NVAR2D,neq);
    }
  };

  /*****************************************************************************
   * This CUDA kernel collects the nodal solution data at the two
   * endpoints of the given edge from the global solution vector.
   ****************************************************************************/

  template <int isystemformat>
  struct gather_DataAtEdge
  { 
  };

  /*****************************************************************************
   * Input:  solution vector Dx stored in interleaved format
   * Output: DataAtEdge vector
   ****************************************************************************/

  template <>
  struct gather_DataAtEdge<SYSTEM_SCALAR>
  {
    template <typename TdSrc,
	      typename TdDest,
	      typename Ti>
    __device__ inline
    static void eval (TdDest *DataAtEdge, 
		      TdSrc *Dx,
		      Ti i, 
		      Ti j, 
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in interleaved format
      
      // Gather solution data at first end point i
      IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,1,i,NVAR2D,neq);
      IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,2,i,NVAR2D,neq);
      IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,3,i,NVAR2D,neq);
      IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,4,i,NVAR2D,neq);

      // Gather solution data at second end point j
      IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,1,j,NVAR2D,neq);
      IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,2,j,NVAR2D,neq);
      IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,3,j,NVAR2D,neq);
      IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,4,j,NVAR2D,neq);
    }
  };

  /*****************************************************************************
   * Input:  solution vector Dx stored in block format
   * Output: DataAtEdge vector
   ****************************************************************************/

  template <>
  struct gather_DataAtEdge<SYSTEM_BLOCK>
  {
    template <typename TdSrc,
	      typename TdDest,
	      typename Ti>
    __device__ inline
    static void eval (TdDest *DataAtEdge,
		      TdSrc *Dx,
		      Ti i, 
		      Ti j, 
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in block format

      // Gather solution data at first end point i
      IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,1,i,NVAR2D,neq);
      IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,2,i,NVAR2D,neq);
      IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,3,i,NVAR2D,neq);
      IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,4,i,NVAR2D,neq);
    
      // Gather solution data at second end point j
      IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,1,j,NVAR2D,neq);
      IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,2,j,NVAR2D,neq);
      IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,3,j,NVAR2D,neq);
      IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,4,j,NVAR2D,neq);
    }
  };

  /*****************************************************************************
   * This CUDA kernel scatters the nodal matrix coefficients into the
   * global operator.
   ****************************************************************************/

  template <int isystemformat, int isystemcoupling>
  struct scatter_CoefficientsAtNode
  {
  };

  /*****************************************************************************
   * Input:  MatrixAtNode
   * Output: block diagonal global operator Da stored in interleaved format
   ****************************************************************************/
  
  template <>
  struct scatter_CoefficientsAtNode<SYSTEM_SCALAR,SYSTEM_SEGREGATED>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtNode,
		     Td *Da,
		     Ti ia,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in interleaved format
      for (int i=1; i <= NVAR2D; i++)
	IDX2_REVERSE(Da,i,ia,NVAR2D,na) +=
	  IDX2(MatrixAtNode,i,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE);
    }
  };
  
  /*****************************************************************************
   * Input:  MatrixAtNode
   * Output: block diagonal global operator Da stored in block format
   ****************************************************************************/
  
  template <>
  struct scatter_CoefficientsAtNode<SYSTEM_BLOCK,SYSTEM_SEGREGATED>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtNode,
		     Td *Da,
		     Ti ia,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory constMemPool
      for (int i=1; i <= NVAR2D; i++)
	IDX1(((Td*)constMemPool[i-1]),ia) +=
	  IDX2(MatrixAtNode,i,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE);
    }
  };
  
  /*****************************************************************************
   * Input:  MatrixAtNode
   * Output: full global operator Da stored in interleaved format
   ****************************************************************************/
  
  template <>
  struct scatter_CoefficientsAtNode<SYSTEM_SCALAR,SYSTEM_ALLCOUPLED>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtNode,
		     Td *Da,
		     Ti ia,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in interleaved format
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
      	IDX2_REVERSE(Da,i,ia,NVAR2D*NVAR2D,na) +=
	  IDX2(MatrixAtNode,i,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE);
    }
  };
  
  /*****************************************************************************
   * Input:  MatrixAtNode
   * Output: full global operator Da stored in block format
   ****************************************************************************/
  
  template <>
  struct scatter_CoefficientsAtNode<SYSTEM_BLOCK,SYSTEM_ALLCOUPLED>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtNode,
		     Td *Da,
		     Ti ia,
		     Ti na,
		     Ti idx)
    {      
      // Operator is full matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	IDX1(((Td*)constMemPool[i-1]),ia) +=
	  IDX2(MatrixAtNode,i,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE);
    }
  };
  
  /*****************************************************************************
   * This CUDA kernel scatters the matrix coefficients for a given
   * edge into the global operator.
   ****************************************************************************/

  template <int isystemformat, int isystemcoupling, 
	    bool bstabilise, bool blumped>
  struct scatter_CoefficientsAtEdge
  {
  };

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in interleaved format
   * Flags:  no stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_SEGREGATED,false,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in interleaved format;
      // no stabilisation and no lumping is applied; 
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ij,NVAR2D,na) = 
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,ji,NVAR2D,na) = 
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 
  
  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in block format
   * Flags:  no stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_SEGREGATED,false,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // no stabilisation and no lumping is applied; 
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ij) =
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),ji) =
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in interleaved format
   * Flags:  no stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_SEGREGATED,false,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in interleaved format;
      // no stabilisation is applied but matrix is lumped;
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ii,NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,jj,NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in block format
   * Flags:  no stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_SEGREGATED,false,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // no stabilisation is applied but matrix is lumped;
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ii) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),jj) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in interleaved format
   * Flags:  stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_SEGREGATED,true,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in interleaved format;
      // no lumping is applied but stabilisation is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ii,NVAR2D,na) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,jj,NVAR2D,na) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,ij,NVAR2D,na) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,ji,NVAR2D,na) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 
  
  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in block format
   * Flags:  stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_SEGREGATED,true,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // no lumping is applied but stabilisation is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ii) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),jj) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),ij) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),ji) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in interleaved format
   * Flags:  stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_SEGREGATED,true,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in interleaved format;
      // stabilisation and lumping is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      // but it it cancelled due to lumping
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ii,NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,jj,NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: block-diagonal global operator Da stored in block format
   * Flags:  stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_SEGREGATED,true,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is block-diagonal matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // stabilisation and lumping is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      // but it it cancelled due to lumping
      for (int i=1; i<=NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ii) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),jj) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in interleaved format
   * Flags:  no stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_ALLCOUPLED,false,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in interleaved format;
      // no stabilisation and no lumping is applied;
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ij,NVAR2D*NVAR2D,na) =
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,ji,NVAR2D*NVAR2D,na) =
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in block format
   * Flags:  no stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_ALLCOUPLED,false,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {


      // Operator is full matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // no stabilisation and no lumping is applied;
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ij) =
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),ji) =
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 
  
  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in interleaved format
   * Flags:  no stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_ALLCOUPLED,false,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in interleaved format;
      // no stabilisation is applied but matrix is lumped;
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ii,NVAR2D*NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,jj,NVAR2D*NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in block format
   * Flags:  no stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_ALLCOUPLED,false,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // no stabilisation is applied but matrix is lumped;
      // Galerkin coefficients are stored in positions 2 and 3
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ii) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),jj) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in interleaved format
   * Flags:  stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_ALLCOUPLED,true,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in interleaved format;
      // no lumping is applied but stabilisation is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ii,NVAR2D*NVAR2D,na) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSize);
	  IDX2_REVERSE(Da,i,jj,NVAR2D*NVAR2D,na) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,ij,NVAR2D*NVAR2D,na) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,ji,NVAR2D*NVAR2D,na) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in block format
   * Flags:  stabilisation; no lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_ALLCOUPLED,true,false>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // no lumping is applied but stabilisation is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ii) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),jj) -=
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),ij) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),ji) =
	    IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE)+
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 
  
  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in interleaved format
   * Flags:  stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_SCALAR,SYSTEM_ALLCOUPLED,true,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in interleaved format;
      // stabilisation and lumping is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      // but it it cancelled due to lumping
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX2_REVERSE(Da,i,ii,NVAR2D*NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX2_REVERSE(Da,i,jj,NVAR2D*NVAR2D,na) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 
  
  /*****************************************************************************
   * Input:  MatrixAtEdge
   * Output: full global operator Da stored in block format
   * Flags:  stabilisation; lumping
   ****************************************************************************/

  template <>
  struct scatter_CoefficientsAtEdge<SYSTEM_BLOCK,SYSTEM_ALLCOUPLED,true,true>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Td *Da,
		     Ti ii,
		     Ti ij,
		     Ti ji,
		     Ti jj,
		     Ti na,
		     Ti idx)
    {
      // Operator is full matrix stored in block format;
      // thus we need to fetch the starting address of each memory
      // block from constant device memory
      // stabilisation and lumping is applied;
      // Galerkin coefficients are stored in positions 2 and 3;
      // artificial diffusion coefficient is stored in position 1
      // but it it cancelled due to lumping
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	{
	  IDX1(((Td*)constMemPool[i-1]),ii) +=
	    IDX3(MatrixAtEdge,i,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	  IDX1(((Td*)constMemPool[i-1]),jj) +=
	    IDX3(MatrixAtEdge,i,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE);
	}
    }
  }; 

  /*****************************************************************************
   * This CUDA kernel calculates the flux Jacobian matrix at the given node.
   ****************************************************************************/

  template <int isystemformat>
  struct calc_FluxJacobiMatrixAtNode
  { 
  };

  /*****************************************************************************
   * Flags: block-diagonal flux Jacobian matrix
   ****************************************************************************/

  template <>
  struct calc_FluxJacobiMatrixAtNode<SYSTEM_SEGREGATED>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtNode,
		     Tc *CoeffsAtNode,
		     Td scale,
		     Td ui,
		     Td vi, 
		     Ti ieq,
		     Ti neq,
		     Ti ncoeff,
		     Ti idx)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX2(MatrixAtNode,1,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,2,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,3,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,4,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
#else
      // Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX2(MatrixAtNode,1,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,2,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,3,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,4,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
#endif
    }
  };

  /*****************************************************************************
   * Flags: full flux Jacobian matrix
   ****************************************************************************/

  template <>
  struct calc_FluxJacobiMatrixAtNode<SYSTEM_ALLCOUPLED>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtNode,
		     Tc *CoeffsAtNode,
		     Td scale,
		     Td ui,
		     Td vi,
		     Td Ei, 
		     Ti ieq,
		     Ti neq,
		     Ti ncoeff,
		     Ti idx)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtNode,1,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,2,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX21(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,3,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX31(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,4,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX41(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,5,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX12(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,6,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,7,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX32(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,8,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX42(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,9,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX13(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,10,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX23(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,11,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,12,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX43(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,13,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX14(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,14,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX24(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,15,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX34(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,16,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
#else
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtNode,1,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,2,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX21(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,3,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX31(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,4,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX41(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,5,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX12(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,6,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,7,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX32(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,8,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX42(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,9,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX13(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,10,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX23(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,11,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,12,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX43(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,13,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX14(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,14,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX24(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,15,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX34(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,16,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE) = -
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
#endif
    }
  };
  
  /*****************************************************************************
   * This CUDA kernel calculates the flux Jacobian matrix at the given edge.
   ****************************************************************************/

  template <int isystemformat>
  struct calc_FluxJacobiMatrixAtEdge
  { 
  };

  /*****************************************************************************
   * Flags: block-diagonal flux Jacobian matrix
   ****************************************************************************/
  
  template <>
  struct calc_FluxJacobiMatrixAtEdge<SYSTEM_SEGREGATED>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Tc *CoeffsAtEdge,
		     Td scale,
		     Td ui,
		     Td uj,
		     Td vi,
		     Td vj,
		     Ti iedge,
		     Ti nedge,
		     Ti ncoeff,
		     Ti idx)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(MatrixAtEdge,1,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,2,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,3,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,4,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);

      // Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(MatrixAtEdge,1,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,2,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,3,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,4,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
#else
      // Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(MatrixAtEdge,1,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,2,2,SHMEM_IDX,NVAR2D,2,1,) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,3,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,4,2,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      
      // Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(MatrixAtEdge,1,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,2,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,3,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,4,3,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
#endif
    }
  };

  /*****************************************************************************
   * Flags: full flux Jacobian matrix
   ****************************************************************************/
  
  template <>
  struct calc_FluxJacobiMatrixAtEdge<SYSTEM_ALLCOUPLED>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Tc *CoeffsAtEdge,
		     Td scale,
		     Td ui,
		     Td uj,
		     Td vi,
		     Td vj,
		     Td Ei,
		     Td Ej,
		     Ti iedge, 
		     Ti nedge,
		     Ti ncoeff,
		     Ti idx)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(MatrixAtEdge,1,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,2,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX21(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,3,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX31(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,4,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX41(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,5,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX12(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,6,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,7,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX32(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,8,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX42(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,9,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX13(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,10,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX23(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,11,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,12,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX43(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,13,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX14(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,14,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX24(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,15,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX34(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,16,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);

      // Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(MatrixAtEdge,1,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,2,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX21(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,3,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX31(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,4,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX41(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,5,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX12(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,6,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,7,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX32(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,8,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX42(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,9,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX13(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,10,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX23(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,11,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,12,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX43(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,13,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX14(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,14,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX24(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,15,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX34(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,16,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
#else
      // Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(MatrixAtEdge,1,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,2,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX21(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,3,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX31(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,4,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX41(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,5,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX12(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,6,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,7,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX32(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,8,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX42(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,9,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX13(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,10,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX23(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,11,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,12,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX43(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,13,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX14(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,14,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX24(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,15,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX34(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,16,2,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      
      // Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(MatrixAtEdge,1,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,2,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX21(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,3,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX31(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,4,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX41(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,5,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX12(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,6,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,7,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX32(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,8,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX42(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,9,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX13(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,10,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX23(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,11,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,12,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX43(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,13,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX14(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,14,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX24(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,15,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX34(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,16,3,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
#endif
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the artificial dissipation at the given edge.
   ****************************************************************************/

  template <int isystemcoupling,int idissipationtype>
  struct calc_DissipationAtEdge
  { 
  };
  
  /*****************************************************************************
   * Zero artificial dissipation, aka standard Galerkin approach
   * Flags: block-diagonal flux Jacobian matrix
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<SYSTEM_SEGREGATED,DISSIPATION_ZERO>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Tc *CoeffsAtEdge,
		     Td *DataAtEdge,
		     Td scale,
		     Td ui,
		     Td uj,
		     Td vi,
		     Td vj,
		     Ti iedge, 
		     Ti nedge,
		     Ti ncoeff,
		     Ti idx)
    {
      for (int i=1; i<=NVAR2D; i++)
	IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D,3,SHMEM_BLOCKSIZE) = 0.0;
    }
  };

  /*****************************************************************************
   * Zero artificial dissipation, aka standard Galerkin approach
   * Flags: full flux Jacobian matrix
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<SYSTEM_ALLCOUPLED,DISSIPATION_ZERO>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *MatrixAtEdge,
		     Tc *CoeffsAtEdge,
		     Td *DataAtEdge,
		     Td scale,
		     Td ui,
		     Td uj,
		     Td vi,
		     Td vj,
		     Td Ei,
		     Td Ej,
		     Ti iedge, 
		     Ti nedge,
		     Ti ncoeff,
		     Ti idx)
    {
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	IDX3(MatrixAtEdge,i,1,SHMEM_IDX,NVAR2D*NVAR2D,3,SHMEM_BLOCKSIZE) = 0.0;
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the diagonal entries of the
   * block-diagonal global operator.
   ****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat>
  __global__ void hydro_calcMatDiagMatD2d_knl(Tc *CoeffsAtDiag,
					      Ti *IdiagList,
					      TdSrc *Dx,
					      TdDest *Da,
					      TdDest scale,
					      Ti neq,
					      Ti na,
					      Ti ncoeff)
  {
#ifdef ENABLE_COPROC_SHMEM
    // Use shared memory
    extern __shared__ TdDest shmemData[];
#endif

    Ti idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx<neq)
      {
	// Get actual equation number
	Ti ieq = IDX2(IdiagList,1,idx+1,2,neq);
		
#ifdef ENABLE_COPROC_SHMEM
	// Local solution data at node from shared memory
	TdDest *DataAtNode = shmemData;
#else
	// Local solution data at node from local memory
	TdDest DataAtNode[NVAR2D];
#endif
	
	// Get solution values at node
	gather_DataAtNode<isystemformat>::
	  eval(DataAtNode,Dx,ieq,neq,idx);

	// Compute velocities
	TdDest ui = XVELOCITY2(DataAtNode,IDX2,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE);
	TdDest vi = YVELOCITY2(DataAtNode,IDX2,SHMEM_IDX,NVAR2D,SHMEM_BLOCKSIZE);
		
	// Compute Galerkin coefficient $K_ii$
	calc_FluxJacobiMatrixAtNode<SYSTEM_SEGREGATED>::
	  eval(DataAtNode,CoeffsAtDiag,scale,ui,vi,ieq,neq,ncoeff,idx);

	// Get diagonal position in the global matrix
	Ti ia  = IDX2(IdiagList,2,idx+1,2,neq);

	// Build coefficients into global operator
	scatter_CoefficientsAtNode<isystemformat,SYSTEM_SEGREGATED>::
	  eval(DataAtNode,Da,ia,na,idx);
      }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the diagonal entries of the
   * full global operator.
   ****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat>
  __global__ void hydro_calcMatDiag2d_knl(Tc *CoeffsAtDiag,
					  Ti *IdiagList,
					  TdSrc *Dx,
					  TdDest *Da,
					  TdDest scale,
					  Ti neq,
					  Ti na,
					  Ti ncoeff)
  {
#ifdef ENABLE_COPROC_SHMEM
    // Use shared memory
    extern __shared__ TdDest shmemData[];
#endif

    Ti idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx<neq)
      {
	// Get actual equation
	Ti ieq = IDX2(IdiagList,1,idx+1,2,neq);

#ifdef ENABLE_COPROC_SHMEM
	// Local solution data at node from shared memory
	TdDest *DataAtNode = shmemData;
#else
	// Local solution data at node from local memory
	TdDest DataAtNode[NVAR2D*NVAR2D];
#endif
      
	// Get solution values at node
	gather_DataAtNode<isystemformat>::
	  eval(DataAtNode,Dx,ieq,neq,idx);

	// Compute velocities and energy
	TdDest ui = XVELOCITY2(DataAtNode,IDX2,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE);
	TdDest vi = YVELOCITY2(DataAtNode,IDX2,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE);
	TdDest Ei = SPECIFICTOTALENERGY2(DataAtNode,IDX2,SHMEM_IDX,NVAR2D*NVAR2D,SHMEM_BLOCKSIZE);

	// Compute Galerkin coefficient $K_ii$
	calc_FluxJacobiMatrixAtNode<SYSTEM_ALLCOUPLED>::
	  eval(DataAtNode,CoeffsAtDiag,scale,ui,vi,Ei,ieq,neq,ncoeff,idx);

	// Get diagonal position in the global matrix
	Ti ia  = IDX2(IdiagList,2,idx+1,2,neq);

	// Build coefficients into global operator
	scatter_CoefficientsAtNode<isystemformat,SYSTEM_ALLCOUPLED>::
	  eval(DataAtNode,Da,ia,na,idx);
      }
  }
  
  /*****************************************************************************
   * This CUDA kernel calculates the off-diagonal entries of the
   * block-diagonal global operator and assembles the artificial
   * dissipation tensor if required.
   ****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipation,
	    bool bstabilised,
	    bool blumped>
  __global__ void hydro_calcMatrixMatD2d_knl(Tc *CoeffsAtEdge,
					     Ti *IedgeList,
					     TdSrc *Dx,
					     TdDest *Da,
					     TdDest scale,
					     Ti neq,
					     Ti na,
					     Ti nedge,
					     Ti ncoeff,
					     Ti nedges,
					     Ti iedgeset)
  {
#ifdef ENABLE_COPROC_SHMEM
    // Use shared memory
    extern __shared__ TdDest shmemData[];
#endif

    Ti idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx<nedges)
      {
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(IedgeList,1,iedgeset+idx,6,nedge);
	Ti j = IDX2(IedgeList,2,iedgeset+idx,6,nedge);

#ifdef ENABLE_COPROC_SHMEM
	// Local solution data at edge from shared memory
	TdDest *DataAtEdge = shmemData;
#else
	// Local solution data at edge from local memory
	TdDest DataAtEdge[2*NVAR2D];
#endif
	
	// Get solution values at edge endpoints
	gather_DataAtEdge<isystemformat>::
	  eval(DataAtEdge,Dx,i,j,neq,idx);

	// Compute velocities
	TdDest ui = XVELOCITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
	TdDest vi = YVELOCITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
      
	TdDest uj = XVELOCITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
	TdDest vj = YVELOCITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
      
#ifdef ENABLE_COPROC_SHMEM
	// Local matrix data at edge from shared memory
	TdDest *MatrixAtEdge = &shmemData[2*NVAR2D*SHMEM_BLOCKSIZE];
#else
	// Local matrix data at edge from local memory
	TdDest MatrixAtEdge[3*NVAR2D];
#endif

	// Compute Galerkin coefficient $K_ij$ and $K_ji$
	calc_FluxJacobiMatrixAtEdge<SYSTEM_SEGREGATED>::
	  eval(MatrixAtEdge,CoeffsAtEdge,scale,
	       ui,uj,vi,vj,iedgeset+idx,nedge,ncoeff,idx);

	// Compute contribution of artificial diffusion
	calc_DissipationAtEdge<SYSTEM_SEGREGATED,idissipation>::
	  eval(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
	       scale,ui,uj,vi,vj,iedgeset+idx,nedge,ncoeff,idx);

	// Get positions of matrix positions (idx starts at zero)
	Ti ij = IDX2(IedgeList,3,iedgeset+idx,6,nedge);
	Ti ji = IDX2(IedgeList,4,iedgeset+idx,6,nedge);
	Ti ii = IDX2(IedgeList,5,iedgeset+idx,6,nedge);
	Ti jj = IDX2(IedgeList,6,iedgeset+idx,6,nedge);

	// Build matrix coefficients into global operator
	scatter_CoefficientsAtEdge
	  <isystemformat,SYSTEM_SEGREGATED,bstabilised,blumped>::
	  eval(MatrixAtEdge,Da,ii,ij,ji,jj,na,idx);
      }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the off-diagonal entries of the full
   * global operator and assembles the artificial dissipation tensor
   * if required.
   ****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipation,
	    bool bstabilised,
	    bool blumped>
  __global__ void hydro_calcMatrix2d_knl(Tc *CoeffsAtEdge,
					 Ti *IedgeList,
					 TdSrc *Dx,
					 TdDest *Da,
					 TdDest scale,
					 Ti neq,
					 Ti na,
					 Ti nedge,
					 Ti ncoeff,
					 Ti nedges,
					 Ti iedgeset)
  {
#ifdef ENABLE_COPROC_SHMEM
    // Use shared memory
    extern __shared__ TdDest shmemData[];
#endif

    Ti idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx<nedges)
      {
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(IedgeList,1,iedgeset+idx,6,nedge);
	Ti j = IDX2(IedgeList,2,iedgeset+idx,6,nedge);
	
#ifdef ENABLE_COPROC_SHMEM
	// Local solution data at edge from shared memory
	TdDest *DataAtEdge = shmemData;
#else
	// Local solution data at edge from local memory
	TdDest DataAtEdge[2*NVAR2D];
#endif
	
	// Get solution values at edge endpoints
	gather_DataAtEdge<isystemformat>::
	  eval(DataAtEdge,Dx,i,j,neq,idx);
	
	// Compute velocities
	TdDest ui = XVELOCITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
	TdDest vi = YVELOCITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
	
	TdDest uj = XVELOCITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
	TdDest vj = YVELOCITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
	
	// Compute specific energies
	TdDest Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);
	TdDest Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR2D,2,SHMEM_BLOCKSIZE);

#ifdef ENABLE_COPROC_SHMEM
	// Local matrix data at edge from shared memory
	TdDest *MatrixAtEdge = &shmemData[2*NVAR2D*SHMEM_BLOCKSIZE];
#else
	// Local matrix data at edge from local memory
	TdDest MatrixAtEdge[3*NVAR2D*NVAR2D];
#endif
	
	// Compute Galerkin coefficient $K_ij$ and $K_ji$
	calc_FluxJacobiMatrixAtEdge<SYSTEM_ALLCOUPLED>::
	  eval(MatrixAtEdge,CoeffsAtEdge,scale,
	       ui,uj,vi,vj,Ei,Ej,iedgeset+idx,nedge,ncoeff,idx);
	
	// Compute contribution of artificial diffusion
	calc_DissipationAtEdge<SYSTEM_ALLCOUPLED,idissipation>::
	  eval(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
	       scale,ui,uj,vi,vj,Ei,Ej,iedgeset+idx,nedge,ncoeff,idx);

	// Get positions of matrix positions (idx starts at zero)
	Ti ij = IDX2(IedgeList,3,iedgeset+idx,6,nedge);
	Ti ji = IDX2(IedgeList,4,iedgeset+idx,6,nedge);
	Ti ii = IDX2(IedgeList,5,iedgeset+idx,6,nedge);
	Ti jj = IDX2(IedgeList,6,iedgeset+idx,6,nedge);
	
	// Build matrix coefficients into global operator
	scatter_CoefficientsAtEdge
	  <isystemformat,SYSTEM_ALLCOUPLED,bstabilised,blumped>::
	  eval(MatrixAtEdge,Da,ii,ij,ji,jj,na,idx);
      }
  }
  
  /*****************************************************************************
   * Internal C++ functions which invoke the CUDA kernels
   *****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcMatDiagMatD2d_cuda(__SIZET *d_CoeffsAtDiag,
				   __SIZET *d_IdiagList,
				   __SIZET *d_Dx,
				   __SIZET *d_Da,
				   TdDest scale,
				   Ti nblocks,
				   Ti neq,
				   Ti na,
				   Ti ncoeff,
				   cudaStream_t stream=0)
  {
    TdSrc *Dx = (TdSrc*)(*d_Dx);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((neq)/(double)(block.x));
  
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      TdDest *Da = (TdDest*)(*d_Da);
      
      hydro_calcMatDiagMatD2d_knl
       	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR>
      	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
       				     IdiagList,
				     Dx, Da, scale,
				     neq, na, ncoeff);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      __SIZET Da[NVAR2D];
      for (int i=0; i<NVAR2D; i++)
	Da[i] = d_Da[i*(NVAR2D+1)];
      
      cudaMemcpyToSymbolAsync("constMemPool", Da,
			      sizeof(__SIZET)*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      hydro_calcMatDiagMatD2d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK>
	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
				     IdiagList,
				     Dx, NULL, scale,
				     neq, na, ncoeff);
    }
    coproc_checkErrors("hydro_calcMatDiagMatD2d_cuda");
    return 0;
  };

  /*****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcMatDiag2d_cuda(__SIZET *d_CoeffsAtDiag,
			       __SIZET *d_IdiagList,
			       __SIZET *d_Dx,
			       __SIZET *d_Da,
			       TdDest scale,
			       Ti nblocks,
			       Ti neq,
			       Ti na,
			       Ti ncoeff,
			       cudaStream_t stream=0)
  {
    TdSrc *Dx = (TdSrc*)(*d_Dx);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((neq)/(double)(block.x));
  
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      TdDest *Da = (TdDest*)(*d_Da);
      
      hydro_calcMatDiag2d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR>
	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
				     IdiagList,
				     Dx, Da, scale,
				     neq, na, ncoeff);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      cudaMemcpyToSymbolAsync("constMemPool", d_Da,
			      sizeof(__SIZET)*NVAR2D*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      hydro_calcMatDiag2d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK>
	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
				     IdiagList,
				     Dx, NULL, scale,
				     neq, na, ncoeff);
    }
    coproc_checkErrors("hydro_calcMatDiag2d_cuda");
    return 0;
  };

  /*****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    bool blumped>
  inline
  int hydro_calcMatGalMatD2d_cuda(__SIZET *d_CoeffsAtEdge,
				  __SIZET *d_IedgeList,
				  __SIZET *d_Dx,
				  __SIZET *d_Da,
				  TdDest scale,
				  Ti nblocks,
				  Ti neq,
				  Ti na,
				  Ti nedge,
				  Ti ncoeff,
				  Ti nedges,
				  Ti iedgeset,
				  cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
    
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      TdDest *Da = (TdDest*)(*d_Da);

      hydro_calcMatrixMatD2d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_ZERO,false,blumped>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Da, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      __SIZET Da[NVAR2D];
      for (int i=0; i<NVAR2D; i++)
	Da[i] = d_Da[i*(NVAR2D+1)];
      
      cudaMemcpyToSymbolAsync("constMemPool", Da,
			      sizeof(__SIZET)*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);

      hydro_calcMatrixMatD2d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_ZERO,false,blumped>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, NULL, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcMatGalMatD2d_cuda");
    return 0;
  };

  /*****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    bool blumped>
  inline
  int hydro_calcMatGalerkin2d_cuda(__SIZET *d_CoeffsAtEdge,
				   __SIZET *d_IedgeList,
				   __SIZET *d_Dx,
				   __SIZET *d_Da,
				   TdDest scale,
				   Ti nblocks,
				   Ti neq,
				   Ti na,
				   Ti nedge,
				   Ti ncoeff,
				   Ti nedges,
				   Ti iedgeset,
				   cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
    
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));

    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      TdDest *Da = (TdDest*)(*d_Da);

      hydro_calcMatrix2d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_ZERO,false,blumped>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Da, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      cudaMemcpyToSymbolAsync("constMemPool", d_Da,
			      sizeof(__SIZET)*NVAR2D*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      hydro_calcMatrix2d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_ZERO,false,blumped>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, NULL, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcMatGalerkin2d_cuda");
    return 0;
  };
  
  /*****************************************************************************
   * External C functions which can be called from the Fortran code
   *****************************************************************************/

  extern "C" {
    __INT FNAME(hydro_calcmatdiagmatd2d_cuda)(__SIZET *d_CoeffsAtDiag,
					      __SIZET *d_IdiagList,
					      __SIZET *d_Dx,
					      __SIZET *d_Da,
					      __DP *scale,
					      __INT *nblocks,
					      __INT *neq,
					      __INT *na,
					      __INT *ncoeff,
					      __I64 *stream)
    {
      return (__INT) hydro_calcMatDiagMatD2d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtDiag, d_IdiagList, d_Dx, d_Da,
			       *scale, *nblocks, *neq, *na, *ncoeff,
			       (cudaStream_t)(*stream));
    }

    /***************************************************************************/
    
    __INT FNAME(hydro_calcmatdiag2d_cuda)(__SIZET *d_CoeffsAtDiag,
					  __SIZET *d_IdiagList,
					  __SIZET *d_Dx,
					  __SIZET *d_Da,
					  __DP *scale,
					  __INT *nblocks,
					  __INT *neq,
					  __INT *na,
					  __INT *ncoeff,
					  __I64 *stream)
    {
      return (__INT) hydro_calcMatDiag2d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtDiag, d_IdiagList, d_Dx, d_Da,
			       *scale, *nblocks, *neq, *na, *ncoeff,
			       (cudaStream_t)(*stream));
    }

    /***************************************************************************/
    
    __INT FNAME(hydro_calcmatgalmatd2d_cuda)(__SIZET *d_CoeffsAtEdge,
					     __SIZET *d_IedgeList,
					     __SIZET *d_Dx,
					     __SIZET *d_Da,
					     __DP *scale,
					     __INT *nblocks,
					     __INT *neq,
					     __INT *na,
					     __INT *nedge,
					     __INT *ncoeff,
					     __INT *nedges,
					     __INT *iedgeset,
					     __INT *cconstrType,
					     __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatGalMatD2d_cuda
	  <__DP,__DP,__DP,__INT,false>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Da,
				       *scale, *nblocks, *neq, *na, *nedge,
				       *ncoeff, *nedges, *iedgeset,
				       (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatGalMatD2d_cuda
	  <__DP,__DP,__DP,__INT,true>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Da,
				      *scale, *nblocks, *neq, *na, *nedge,
				      *ncoeff, *nedges, *iedgeset,
				      (cudaStream_t)(*stream));
    }
    
    /***************************************************************************/
    
    __INT FNAME(hydro_calcmatgalerkin2d_cuda)(__SIZET *d_CoeffsAtEdge,
					      __SIZET *d_IedgeList,
					      __SIZET *d_Dx,
					      __SIZET *d_Da,
					      __DP *scale,
					      __INT *nblocks,
					      __INT *neq,
					      __INT *na,
					      __INT *nedge,
					      __INT *ncoeff,
					      __INT *nedges,
					      __INT *iedgeset,
					      __INT *cconstrType,
					      __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatGalerkin2d_cuda
	  <__DP,__DP,__DP,__INT,false>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Da,
				       *scale, *nblocks, *neq, *na, *nedge,
				       *ncoeff, *nedges, *iedgeset,
				       (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatGalerkin2d_cuda
	  <__DP,__DP,__DP,__INT,true>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Da,
				      *scale, *nblocks, *neq, *na, *nedge,
				      *ncoeff, *nedges, *iedgeset,
				      (cudaStream_t)(*stream));
    }
  };
    
}
