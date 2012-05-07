/*#############################################################################
 ******************************************************************************
 * <name> cudaGather </name>
 ******************************************************************************
 *
 * <purpose>
 * This header file contains CUDA kernels for gathering/scattering data.
 * </purpose>
 *
 *############################################################################
 */

#include "coproc_core.h"
#include "coproc_storage_cuda.h"
#include "cudaDMA.h"

#define LANGUAGE LANGUAGE_C
#include "../../../kernel/System/idxmanager.h"

/*******************************************************************************
 * VectorBase - one-dimensional array (basic facilities and specialisations)
 ******************************************************************************/
template <int nelem, int isystemformat>
  struct VectorBase
  {
  };

/*******************************************************************************
 * VectorBase: Vector is stored in interleaved format.
 ******************************************************************************/
template <int nelem>
  struct VectorBase<nelem,SYSTEM_SCALAR>
{ 
  /*****************************************************************************
   * Gather nodal data for neqsim nodes with nelem elements per node
   ****************************************************************************/
  template <int neqsim,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti ieq, 
				Ti neq,
				Ti tid)
	    {
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2(DataAtNode,i,tid,nelem,neqsim) = IDX2_REVERSE(vec,i,ieq,nelem,neq);
	    }
  
  /*****************************************************************************
   * Scatter nodal data for neqsim nodes with nelem elements per node
   ****************************************************************************/
  template <int neqsim,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti ieq, 
				 Ti neq,
				 Ti tid)
	    {
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2_REVERSE(vec,i,ieq,nelem,neq) += IDX2(DataAtNode,i,tid,nelem,neqsim);
	    }

  /*****************************************************************************
   * Gather edge-wise data for nedgesim edges with nelem elements per node
   ****************************************************************************/
  template <int nedgesim,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge, 
				Tv *vec,
				Ti ieq, 
				Ti jeq, 
				Ti neq,
				Ti tid)
	    {
	      // Gather data at first end point ieq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim) = IDX2_REVERSE(vec,i,ieq,nelem,neq);
	      
	      // Gather data at second end point jeq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim) = IDX2_REVERSE(vec,i,jeq,nelem,neq);
	    }

  /*****************************************************************************
   * Scatter edge-wise data for nedgesim edges with nelem elements per node
   ****************************************************************************/
  template <int nedgesim,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti ieq,
				 Ti jeq,
				 Ti neq,
				 Ti tid)
	    {
	      // Scatter data to first end point ieq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2_REVERSE(vec,i,ieq,nelem,neq) += IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
	      
	      // Scatter data to second end point jeq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2_REVERSE(vec,i,jeq,nelem,neq) += IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
	    }  
};

/*****************************************************************************
 * VectorBase: Vector is stored in block format.
 ****************************************************************************/
template <int nelem>
  struct VectorBase<nelem,SYSTEM_BLOCK>
{
  /*****************************************************************************
   * Gather nodal data for neqsim nodes with nelem elements per node
   ****************************************************************************/
  template <int neqsim,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti ieq, 
				Ti neq,
				Ti tid)
	    {
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2(DataAtNode,i,tid,nelem,neqsim) = IDX2_FORWARD(vec,i,ieq,nelem,neq);
	    }

  /*****************************************************************************
   * Scatter nodal data for neqsim nodes with nelem elements per node
   ****************************************************************************/
  template <int neqsim,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti ieq, 
				 Ti neq,
				 Ti tid)
	    {
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2_FORWARD(vec,i,ieq,nelem,neq) += IDX2(DataAtNode,i,tid,nelem,neqsim);
	    }
  
  /*****************************************************************************
   * Gather edge-wise data for nedgesim edges with nelem elements per node
   ****************************************************************************/
  template <int nedgesim,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge,
				Tv *vec,
				Ti ieq, 
				Ti jeq, 
				Ti neq,
				Ti tid)
	    {
	      // Gather data at first end point ieq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim) = IDX2_FORWARD(vec,i,ieq,nelem,neq);
	      
	      // Gather data at second end point jeq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim) = IDX2_FORWARD(vec,i,jeq,nelem,neq);
	    }

  /*****************************************************************************
   * Scatter edge-wise data for nedgesim edges with nelem elements per node
   ****************************************************************************/
  template <int nedgesim,
            typename Td,
            typename Tv,
            typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti ieq, 
				 Ti jeq, 
				 Ti neq,
				 Ti tid)
	    {
	      // Scatter data to first end point ieq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2_FORWARD(vec,i,ieq,nelem,neq) += IDX3(DataAtEdge,1,i,tid,nelem,2,nedgesim);
	      
	      // Scatter data to second end point jeq
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2_FORWARD(vec,i,jeq,nelem,neq) += IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
	    }
};


/*******************************************************************************
 * Vector - one-dimensional array
 ******************************************************************************/
template <int nelem, int isystemformat>
  struct Vector : public VectorBase<nelem,isystemformat>
{
  // Enable use of inherited functions
  using VectorBase<nelem,isystemformat>::gatherNodeData;
  using VectorBase<nelem,isystemformat>::scatterNodeData;
  using VectorBase<nelem,isystemformat>::gatherEdgeData;
  using VectorBase<nelem,isystemformat>::scatterEdgeData;

  /*****************************************************************************
   * Gather nodal data for a single nodes with nelem elements per node
   ****************************************************************************/
  template <typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti ieq,
				Ti neq)
	    {
	      VectorBase<nelem,isystemformat>::gatherNodeData<1>
		(DataAtNode,vec,ieq,neq,1);
	    }
  
  /*****************************************************************************
   * Scatter nodal data for a single nodes with nelem elements per node
   ****************************************************************************/
  template <typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti ieq,
				 Ti neq)
	    {
	      VectorBase<nelem,isystemformat>::scatterNodeData<1>
		(vec,DataAtNode,ieq,neq,1);
	    }
  
  /*****************************************************************************
   * Gather nodal data for a single edge with nelem elements per edge
   ****************************************************************************/
  template <typename Tv, typename Td, typename Ti>
     __device__ __forceinline__
     static void gatherEdgeData (Td *DataAtEdge,
				 Tv *vec,
				 Ti ieq,
				 Ti jeq,
				 Ti neq)
	    {
	      VectorBase<nelem,isystemformat>::gatherEdgeData<1>
		(DataAtEdge,vec,ieq,jeq,neq,1);
	    }

  /*****************************************************************************
   * Scatter nodal data for a single edge with nelem elements per edge
   ****************************************************************************/
  template <typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti ieq,
				 Ti jeq,
				 Ti neq)
	    {
	      VectorBase<nelem,isystemformat>::scatterEdgeData<1>
		(vec,DataAtEdge,ieq,jeq,neq,1);
	    }
};

/*******************************************************************************
 * MatrixBase - two dimensional array (basic facilities and specialisations)
 ******************************************************************************/

template <int nelem, int isystemformat>
  struct MatrixBase
  { 
  };

/*******************************************************************************
 * MatrixBase: Matrix is stored in interleaved format.
 ******************************************************************************/

template <int nelem>
  struct MatrixBase<nelem,SYSTEM_SCALAR>
{
  /*****************************************************************************
   * Scatter data for neqsim nodes with nelem elements per node
   ****************************************************************************/
  template <int neqsim,
            typename Tm,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tm *mat,
				 Td *DataAtNode,
				 Ti ia,
				 Ti na,
				 Ti tid)
	    {
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX2_REVERSE(mat,i,ia,nelem,na) += IDX2(DataAtNode,i,tid,nelem,neqsim);
	    }
  
  /*****************************************************************************
   * Scatter data for nedgesim edges with nelem elements per node;
   * off-diagonal matrix positions are given explicitly
   ****************************************************************************/
  template <int nedgesim,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti ii,
				 Ti ij,
				 Ti ji,
				 Ti jj,
				 Ti na,
				 Ti tid)
	    {
	      if (bstabilise)
		{
		  if (blumping)
		    {
		      // perform lumping and apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim);
			}
		      
		    } else
		    {
		      // perform no lumping but apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,ij,nelem,na)  = IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim)+
			                                     IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,ji,nelem,na)  = IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim)+
			                                     IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			}
		    }
		} else
		{
		  if (blumping)
		    {
		      // perform lumping but do not apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    } else
		    {
		      // perform no lumping and do not apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) = IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX2_REVERSE(mat,i,ji,nelem,na) = IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    }
		}
	    }
  
  /*****************************************************************************
   * Scatter data for nedgesim edges with nelem elements per node;
   * matrix positions are retrieved from edge list
   ****************************************************************************/
  template <int nedgesim,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti iedge,
				 Ti na,
				 Ti nedge,
				 Ti tid)
	    {
	      if (bstabilise)
		{
		  if (blumping)
		    {
		      // perform lumping and apply stabilisation
		      Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		      Ti jj = IDX2(IedgeList,6,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim);
			}
		      
		    } else
		    {
		      // perform no lumping but apply stabilisation
		      Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		      Ti ji = IDX2(IedgeList,4,iedge,6,nedge);
		      Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		      Ti jj = IDX2(IedgeList,6,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,ij,nelem,na)  = IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim)+
			                                     IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX2_REVERSE(mat,i,ji,nelem,na)  = IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim)+
			                                     IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			}
		    }
		} else
		{
		  if (blumping)
		    {
		      // perform lumping but do not apply stabilisation
		      Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		      Ti jj = IDX2(IedgeList,6,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    } else
		    {
		      // perform no lumping and do not apply stabilisation
		      Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		      Ti ji = IDX2(IedgeList,4,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) = IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX2_REVERSE(mat,i,ji,nelem,na) = IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    }
		}
	    }
};

/*******************************************************************************
 * MatrixBase: Matrix is stored in block format.
 ******************************************************************************/

template <int nelem>
  struct MatrixBase<nelem,SYSTEM_BLOCK>
{
  /*****************************************************************************
   * Scatter data for neqsim node with nelem elements per node
   ****************************************************************************/
  template <int neqsim,
            typename Tm,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tm *mat,
				 Td *DataAtNode,
				 Ti ia,
				 Ti na,
				 Ti tid)
	    {
#pragma unroll
	      for (int i=1; i<=nelem; i++)
		IDX1((((Tm**)mat)[i-1]),ia) += IDX2(DataAtNode,i,tid,nelem,neqsim);
	    }

  /*****************************************************************************
   * Scatter data for nedgesim edges with nelem elements per node;
   * off-diagonal matrix positions are given explicitly
   ****************************************************************************/
  template <int nedgesim,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti ii,
				 Ti ij,
				 Ti ji,
				 Ti jj,
				 Ti na,
				 Ti tid)
	    {
	      if (bstabilise)
		{
		  if (blumping)
		    {
		      // perform lumping and apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim);
			}
		      
		    } else
		    {
		      // perform no lumping but apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),ij)  = IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim)+
			                                 IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),ji)  = IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim)+
			                                 IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			}
		    }
		} else
		{
		  if (blumping)
		    {
		      // perform lumping but do not apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    } else
		    {
		      // perform no lumping and do not apply stabilisation
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) = IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),ji) = IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    }
		}
	    }
  
  /*****************************************************************************
   * Scatter data for nedgesim edges with nelem elements per node;
   * matrix positions are retrieved from edge list
   ****************************************************************************/
  template <int nedgesim,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti iedge,
				 Ti na,
				 Ti nedge,
				 Ti tid)
	    {
	      if (bstabilise)
		{
		  if (blumping)
		    {
		      // perform lumping and apply stabilisation
		      Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		      Ti jj = IDX2(IedgeList,6,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim);
			}
		      
		    } else
		    {
		      // perform no lumping but apply stabilisation
		      Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		      Ti ji = IDX2(IedgeList,4,iedge,6,nedge);
		      Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		      Ti jj = IDX2(IedgeList,6,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),ij)  = IDX3(DataAtEdge,i,1,tid,nelem,3,nedgesim)+
			                                 IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),ji)  = IDX3(DataAtEdge,i,2,tid,nelem,3,nedgesim)+
			                                 IDX3(DataAtEdge,i,3,tid,nelem,3,nedgesim);
			}
		    }
		} else
		{
		  if (blumping)
		    {
		      // perform lumping but do not apply stabilisation
		      Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		      Ti jj = IDX2(IedgeList,6,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    } else
		    {
		      // perform no lumping and do not apply stabilisation
		      Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		      Ti ji = IDX2(IedgeList,4,iedge,6,nedge);
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) = IDX3(DataAtEdge,i,1,tid,nelem,2,nedgesim);
			  IDX1((((Tm**)mat)[i-1]),ji) = IDX3(DataAtEdge,i,2,tid,nelem,2,nedgesim);
			}
		    }
		}
	    }
};

/*******************************************************************************
 * Matrix - two-dimensional array
 ******************************************************************************/
template <int nelem, int isystemformat>
  struct Matrix : public MatrixBase<nelem,isystemformat>
{
  // Enable use of inherited functions
  using MatrixBase<nelem,isystemformat>::scatterNodeData;
  using MatrixBase<nelem,isystemformat>::scatterEdgeData;
  
  /*****************************************************************************
   * Scatter nodal data for a single nodes with nelem elements per node
   ****************************************************************************/
  template <typename Tm, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tm *mat,
				 Td *DataAtNode,
				 Ti ieq,
				 Ti neq)
  {
    MatrixBase<nelem,isystemformat>::scatterNodeData<1>
      (mat,DataAtNode,ieq,neq,1);
  }

  /*****************************************************************************
   * Scatter data for a single edges with nelem elements per node;
   * off-diagonal matrix positions are given explicitly
   ****************************************************************************/
  template <bool bstabilise, bool blumping, typename Tm, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti ii,
				 Ti ij,
				 Ti ji,
				 Ti jj,
				 Ti na)
  {
    MatrixBase<nelem,isystemformat>::scatterEdgeData<1,bstabilise,blumping>
      (mat,DataAtEdge,ii,ij,ji,jj,na,1);
  }

  /*****************************************************************************
   * Scatter data for a single edges with nelem elements per node;
   * matrix positions are retrieved from edge list
   ****************************************************************************/
  template <bool bstabilise, bool blumping, typename Tm, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti iedge,
				 Ti na,
				 Ti nedge)
  {
    MatrixBase<nelem,isystemformat>::scatterEdgeData<1,bstabilise,blumping>
      (mat,DataAtEdge,IedgeList,iedge,na,nedge,1);
  }
};
