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

#include <coproc_core.h>
#include <coproc_storage_cuda.h>
#include "cudaDMA.h"

#define LANGUAGE LANGUAGE_C
#include "flagship.h"

/*******************************************************************************
 * VectorBase - one-dimensional array (basic facilities and specialisations)
 *
 * Template parameter nelem denotes the number of elements per vector
 * entry. By default, a vector is stored as an array-of-structures:
 *
 * <1,2,...,nelem>_1, <1,2,...,nelem>_2, ..., <1,2,...,nelem>_neq
 *
 * If the template parameter soa is true, then the vector is
 * assumed to be stored as structure-of-arrays:
 *
 * <1,2,...,neq>_1, <1,2,...,neq>_2, ..., <1,2,...neq>_nelem
 ******************************************************************************/
template <int nelem, bool soa>
  struct VectorBase
  {
  };

/*******************************************************************************
 * VectorBase: Vector is stored as array-of-structures
 ******************************************************************************/
template <int nelem>
  struct VectorBase<nelem,false>
{ 
  /*****************************************************************************
   * Gather nodal data for node ieq and store it at position ipos in
   * local data array DataAtNode
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti ipos,
				Ti ieq, 
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2T(DataAtNode,i,ipos,nelem,datalen) = IDX2_REVERSE(vec,i,ieq,nelem,neq);
		} 
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2(DataAtNode,i,ipos,nelem,datalen) = IDX2_REVERSE(vec,i,ieq,nelem,neq);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2T(DataAtNode,i,ipos,nelem,datalen) += IDX2_REVERSE(vec,i,ieq,nelem,neq);
		} 
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2(DataAtNode,i,ipos,nelem,datalen) += IDX2_REVERSE(vec,i,ieq,nelem,neq);
		}
	      }
	    }
  
  /*****************************************************************************
   * Scatter nodal data stored at position ipos in local data array
   * DataAtNode to node ieq into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti ipos,
				 Ti ieq, 
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) = IDX2(DataAtNode,i,ipos,nelem,datalen);
		} 
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) = IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) += IDX2(DataAtNode,i,ipos,nelem,datalen);
		} 
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) += IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	    }
  
  /*****************************************************************************
   * Gather nodal data for all nodes in IeqList and store it at
   * positions 1...datalen in local data array DataAtNode.
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti *IeqList,
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2T(DataAtNode,i,ipos+1,nelem,datalen) =
			  IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		} 
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2(DataAtNode,i,ipos+1,nelem,datalen) =
			  IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2T(DataAtNode,i,ipos+1,nelem,datalen) +=
			  IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		} 
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2(DataAtNode,i,ipos+1,nelem,datalen) +=
			  IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		}
	      }
	    }
  
  /*****************************************************************************
   * Scatter all nodal data stored at positions 1...datalen in local
   * data array DataAtNode to all nodes in IeqList into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti *IeqList, 
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq) =
			  IDX2T(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		} 
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq) =
			  IDX2(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq) +=
			  IDX2T(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		} 
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IeqList[ipos*incr],nelem,neq) +=
			  IDX2(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		}
	      }
	    }
  
  /*****************************************************************************
   * Gather edge-wise data for edge (ieq,jeq) and store it at position
   * ipos in local data array DataAtEdge
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge, 
				Tv *vec,
				Ti ipos,
				Ti ieq, 
				Ti jeq, 
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen) = IDX2_REVERSE(vec,i,ieq,nelem,neq);
		  
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen) = IDX2_REVERSE(vec,i,jeq,nelem,neq);
		} 
		else {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen) = IDX2_REVERSE(vec,i,ieq,nelem,neq);
		  
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen) = IDX2_REVERSE(vec,i,jeq,nelem,neq);
		}
	      }
	      else {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen) += IDX2_REVERSE(vec,i,ieq,nelem,neq);
		  
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen) += IDX2_REVERSE(vec,i,jeq,nelem,neq);
		} 
		else {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen) += IDX2_REVERSE(vec,i,ieq,nelem,neq);
		  
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen) += IDX2_REVERSE(vec,i,jeq,nelem,neq);
		}
	      }
	    }

  /*****************************************************************************
   * Scatter edge-wise data stored at position ipos in local data
   * array DataAtEdge to edge (ieq,jeq) into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti ipos,
				 Ti ieq,
				 Ti jeq,
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) = IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
		  
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,jeq,nelem,neq) = IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
		else {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) = IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
		  
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,jeq,nelem,neq) = IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
		  
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,jeq,nelem,neq) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
		else {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,ieq,nelem,neq) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
		  
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(vec,i,jeq,nelem,neq) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
	      }
	    }  

  /*****************************************************************************
   * Gather edge-wise data for all edges in IedgeList and store it at
   * positions 1...datalen in local data array DataAtEdge
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge, 
				Tv *vec,
				Ti *IedgeList, 
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen) =
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq);
		      
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen) =
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen) =
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq);
		      
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen) =
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen) +=
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq);
		      
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen) +=
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen) +=
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq);
		      
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen) +=
			  IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
	      }
	    }

  /*****************************************************************************
   * Scatter all edge-wise data stored at positions 1...datalen in local data
   * array DataAtEdge to all edges in IedgeList into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq) =
			  IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		      
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq) =
			  IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
	      }
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq) =
			  IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		      
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq) =
			  IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq) +=
			  IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		      
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq) +=
			  IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
	      }
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr],nelem,neq) +=
			  IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		      
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_REVERSE(vec,i,IedgeList[ipos*incr+1],nelem,neq) +=
			  IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
		}
	      }
	    }
};
  
/*****************************************************************************
 * VectorBase: Vector is stored as structure-of-arrays
 ****************************************************************************/
template <int nelem>
  struct VectorBase<nelem,true>
{
  /*****************************************************************************
   * Gather nodal data for node ieq and store it at position ipos in
   * local data array DataAtNode
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti ipos,
				Ti ieq, 
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2T(DataAtNode,i,ipos,nelem,datalen) = IDX2_FORWARD(vec,i,ieq,nelem,neq);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2(DataAtNode,i,ipos,nelem,datalen) = IDX2_FORWARD(vec,i,ieq,nelem,neq);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2T(DataAtNode,i,ipos,nelem,datalen) += IDX2_FORWARD(vec,i,ieq,nelem,neq);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2(DataAtNode,i,ipos,nelem,datalen) += IDX2_FORWARD(vec,i,ieq,nelem,neq);
		}
	      }
	    }

  /*****************************************************************************
   * Scatter nodal data stored at position ipos in local data array
   * DataAtNode to node ieq into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti ipos,
				 Ti ieq, 
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) = IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) = IDX2(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) += IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) += IDX2(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	    }
  
  /*****************************************************************************
   * Gather nodal data for all nodes in IeqList and store it at
   * positions 1...datalen in local data array DataAtNode
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti *IeqList, 
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2T(DataAtNode,i,ipos+1,nelem,datalen) =
			  IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2(DataAtNode,i,ipos+1,nelem,datalen) =
			  IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2T(DataAtNode,i,ipos+1,nelem,datalen) +=
			  IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2(DataAtNode,i,ipos+1,nelem,datalen) +=
			  IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq);
		    }
		}
	      }
	    }

  /*****************************************************************************
   * Scatter nodal data stored at positions 1...datalen in local data
   * array DataAtNode to all nodes in IeqList into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti *IeqList, 
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq) =
			  IDX2T(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq) =
			  IDX2(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq) +=
			  IDX2T(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IeqList[ipos*incr],nelem,neq) +=
			  IDX2(DataAtNode,i,ipos+1,nelem,datalen);
		    }
		}
	      }
	    }
  
  /*****************************************************************************
   * Gather edge-wise data for edge (ieq,jeq) and store it at position
   * ipos in local data array DataAtEdge
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge,
				Tv *vec,
				Ti ipos,
				Ti ieq, 
				Ti jeq, 
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen) = IDX2_FORWARD(vec,i,ieq,nelem,neq);
		
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen) = IDX2_FORWARD(vec,i,jeq,nelem,neq);
		}
		else {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen) = IDX2_FORWARD(vec,i,ieq,nelem,neq);
		
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen) = IDX2_FORWARD(vec,i,jeq,nelem,neq);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen) += IDX2_FORWARD(vec,i,ieq,nelem,neq);
		
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen) += IDX2_FORWARD(vec,i,jeq,nelem,neq);
		}
		else {
		  // Gather data at first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen) += IDX2_FORWARD(vec,i,ieq,nelem,neq);
		
		  // Gather data at second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen) += IDX2_FORWARD(vec,i,jeq,nelem,neq);
		}
	      }
	    }

  /*****************************************************************************
   * Scatter edge-wise data stored at position ipos in local data
   * array DataAtEdge to edge (ieq,jeq) into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Td,
            typename Tv,
            typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti ipos,
				 Ti ieq, 
				 Ti jeq, 
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) = IDX3T(DataAtEdge,1,i,ipos,nelem,2,datalen);
		
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,jeq,nelem,neq) = IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
		else {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) = IDX3(DataAtEdge,1,i,ipos,nelem,2,datalen);
		
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,jeq,nelem,neq) = IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) += IDX3T(DataAtEdge,1,i,ipos,nelem,2,datalen);
		
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,jeq,nelem,neq) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
		else {
		  // Scatter data to first end point ieq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,ieq,nelem,neq) += IDX3(DataAtEdge,1,i,ipos,nelem,2,datalen);
		
		  // Scatter data to second end point jeq
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_FORWARD(vec,i,jeq,nelem,neq) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		}
	      }
	    }

  /*****************************************************************************
   * Gather edge-wise data for all edges in IedgeList and store it at
   * positions 1...datalen in local data array DataAtEdge
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge, 
				Tv *vec,
				Ti *IedgeList, 
				Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen) =
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq);
		    
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen) =
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen) =
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq);
		    
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen) =
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen) +=
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq);
		    
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen) +=
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Gather data at first end point
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen) +=
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq);
		    
		      // Gather data at second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen) +=
			  IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq);
		    }
		}
	      }
	    }

  /*****************************************************************************
   * Scatter all edge-wise data stored at positions 1...datalen in local data
   * array DataAtEdge to all edges in IedgeList into global vector
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            int incr,
            typename Tv,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti neq)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq) =
			  IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		    
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq) =
			  IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq) =
			  IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		    
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq) =
			  IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq) +=
			  IDX3T(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		    
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq) +=
			  IDX3T(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
		}
		else {
		  for (int ipos=0; ipos<datalen; ipos++)
		    {
		      // Scatter data to first end point ieq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr],nelem,neq) +=
			  IDX3(DataAtEdge,i,1,ipos+1,nelem,2,datalen);
		    
		      // Scatter data to second end point jeq
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			IDX2_FORWARD(vec,i,IedgeList[ipos*incr+1],nelem,neq) +=
			  IDX3(DataAtEdge,i,2,ipos+1,nelem,2,datalen);
		    }  
		}
	      }
	    }
};


/*******************************************************************************
 * Vector - one-dimensional array
 ******************************************************************************/
template <int nelem, bool soa>
  struct Vector : public VectorBase<nelem,soa>
{
  // Enable use of inherited functions
  using VectorBase<nelem,soa>::gatherNodeData;
  using VectorBase<nelem,soa>::scatterNodeData;
  using VectorBase<nelem,soa>::gatherEdgeData;
  using VectorBase<nelem,soa>::scatterEdgeData;
  
  /*****************************************************************************
   * Gather nodal data for node ieq and store it at position ipos in
   * local data array DataAtNode
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti ieq,
				Ti neq)
  {
    VectorBase<nelem,soa>::gatherNodeData<1,false,boverwrite>
      (DataAtNode,vec,1,ieq,neq);
  }
  
  /*****************************************************************************
   * Scatter nodal data stored at position ipos in local data array
   * DataAtNode to node ieq into global vector
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti ieq,
				 Ti neq)
  {
    VectorBase<nelem,soa>::scatterNodeData<1,false,boverwrite>
      (vec,DataAtNode,1,ieq,neq);
  }
  
  /*****************************************************************************
   * Gather nodal data for first node in IeqList and store it at
   * position ipos in local data array DataAtNode
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void gatherNodeData (Td *DataAtNode,
				Tv *vec,
				Ti *IeqList,
				Ti neq)
  {
    VectorBase<nelem,soa>::gatherNodeData<1,false,boverwrite,1>
      (DataAtNode,vec,IeqList,neq);
  }

  /*****************************************************************************
   * Scatter nodal data stored at position ipos in local data array
   * DataAtNode to first node in IeqList into global vector
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tv *vec,
				 Td *DataAtNode,
				 Ti *IeqList,
				 Ti neq)
  {
    VectorBase<nelem,soa>::scatterNodeData<1,false,boverwrite,1>
      (vec,DataAtNode,1,IeqList,neq);
  }
  
  /*****************************************************************************
   * Gather edge-wise data for edge (ieq,jeq) and store it at position
   * ipos in local data array DataAtEdge
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge,
				Tv *vec,
				Ti ieq,
				Ti jeq,
				Ti neq)
  {
    VectorBase<nelem,soa>::gatherEdgeData<1,false,boverwrite>
      (DataAtEdge,vec,1,ieq,jeq,neq);
  }
 
  /*****************************************************************************
   * Scatter edge-wise data stored at position ipos in local data
   * array DataAtEdge to edge (ieq,jeq) into global vector
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti ieq,
				 Ti jeq,
				 Ti neq)
  {
    VectorBase<nelem,soa>::scatterEdgeData<1,false,boverwrite>
      (vec,DataAtEdge,1,ieq,jeq,neq);
  }

  /*****************************************************************************
   * Gather edge-wise data for first edge in IedgeList and store it at
   * position ipos in local data array DataAtEdge
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void gatherEdgeData (Td *DataAtEdge,
				Tv *vec,
				Ti *IedgeList,
				Ti neq)
  {
    VectorBase<nelem,soa>::gatherEdgeData<1,false,boverwrite,1>
      (DataAtEdge,vec,1,IedgeList,neq);
  }

  /*****************************************************************************
   * Scatter edge-wise data stored at position ipos in local data
   * array DataAtEdge to first edge in IedgeList into global vector
   ****************************************************************************/
  template <bool boverwrite, typename Tv, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tv *vec,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti neq)
  {
    VectorBase<nelem,soa>::scatterEdgeData<1,false,boverwrite>
      (vec,DataAtEdge,1,IedgeList,neq);
  }
};

/*******************************************************************************
 * MatrixBase - two dimensional array (basic facilities and specialisations)
 *
 * Template parameter nelem denotes the number of elements per matrix
 * entry. By default, a matrix is stored as an array-of-structures:
 *
 * <1,2,...,nelem>_1, <1,2,...,nelem>_2, ..., <1,2,...,nelem>_na
 *
 * If the template parameter soa is true, then the matrix is
 * assumed to be stored as structure-of-arrays:
 *
 * <1,2,...,na>_1, <1,2,...,na>_2, ..., <1,2,...na>_nelem

 ******************************************************************************/
template <int nelem, bool soa>
  struct MatrixBase
  { 
  };

/*******************************************************************************
 * MatrixBase: Matrix is stored as array-of-structures
 ******************************************************************************/
template <int nelem>
  struct MatrixBase<nelem,false>
{
  /*****************************************************************************
   * Scatter data for datalen nodes with nelem elements per node
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tm,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tm *mat,
				 Td *DataAtNode,
				 Ti ipos,
				 Ti ia,
				 Ti na)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(mat,i,ia,nelem,na) = IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(mat,i,ia,nelem,na) = IDX2(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(mat,i,ia,nelem,na) += IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX2_REVERSE(mat,i,ia,nelem,na) += IDX2(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	    }
  
  /*****************************************************************************
   * Scatter data for datalen edges with nelem elements per node;
   * off-diagonal matrix positions are given explicitly
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti ipos,
				 Ti ii,
				 Ti ij,
				 Ti ji,
				 Ti jj,
				 Ti na)
	    {
	      if (bstabilise) {
		if (blumping) {
		  // perform lumping and apply stabilisation
		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		    
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		  }
		}
		else {
		  // perform no lumping but apply stabilisation
		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na)  = IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na)  = IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na)  = IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
  			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na)  = IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
  			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		}
	      }
	      else {
		if (blumping) {
		  // perform lumping but do not apply stabilisation
		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  }
		}
		else {
		  // perform no lumping and do not apply stabilisation
		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) = IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) = IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) = IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) = IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    
		    // perform no lumping and do not apply stabilisation
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		}
	      }
	    }
  
  /*****************************************************************************
   * Scatter data for datalen edges with nelem elements per node;
   * matrix positions are retrieved from edge list
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti ipos,
				 Ti iedge,
				 Ti nedge,
				 Ti na)
	    {
	      if (bstabilise) {
		if (blumping) {
		  // perform lumping and apply stabilisation
		  Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		  Ti jj = IDX2(IedgeList,6,iedge,6,nedge);
		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		  }
		}
		else {
		  // perform no lumping but apply stabilisation
		  Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		  Ti ji = IDX2(IedgeList,4,iedge,6,nedge);
		  Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		  Ti jj = IDX2(IedgeList,6,iedge,6,nedge);

		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na)  = IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na)  = IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na)  = IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na)  = IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ii,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,jj,nelem,na) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                     IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		}
	      }
	      else {
		if (blumping) {
		  // perform lumping but do not apply stabilisation
		  Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		  Ti jj = IDX2(IedgeList,6,iedge,6,nedge);

		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX2_REVERSE(mat,i,ii,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX2_REVERSE(mat,i,jj,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  }
		}
		else {
		  // perform no lumping and do not apply stabilisation
		  Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		  Ti ji = IDX2(IedgeList,4,iedge,6,nedge);

		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) = IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) = IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) = IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) = IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX2_REVERSE(mat,i,ij,nelem,na) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX2_REVERSE(mat,i,ji,nelem,na) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		}
	      }
	    }
};

/*******************************************************************************
 * MatrixBase: Matrix is stored structure-of-arrays
 ******************************************************************************/

template <int nelem>
  struct MatrixBase<nelem,true>
{
  /*****************************************************************************
   * Scatter data for datalen node with nelem elements per node
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            typename Tm,
            typename Td,
            typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tm *mat,
				 Td *DataAtNode,
				 Ti ipos,
				 Ti ia,
				 Ti na)
	    {
	      if (boverwrite) {
		//
		// Overwrite existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX1((((Tm**)mat)[i-1]),ia) = IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX1((((Tm**)mat)[i-1]),ia) = IDX2(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	      else {
		//
		// Keep existing data
		//
		if (btranspose) {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX1((((Tm**)mat)[i-1]),ia) += IDX2T(DataAtNode,i,ipos,nelem,datalen);
		}
		else {
#pragma unroll
		  for (int i=1; i<=nelem; i++)
		    IDX1((((Tm**)mat)[i-1]),ia) += IDX2(DataAtNode,i,ipos,nelem,datalen);
		}
	      }
	    }
  
  /*****************************************************************************
   * Scatter data for datalen edges with nelem elements per node;
   * off-diagonal matrix positions are given explicitly
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti ipos,
				 Ti ii,
				 Ti ij,
				 Ti ji,
				 Ti jj,
				 Ti na)
	    {
	      if (bstabilise) {
		if (blumping) {
		  // perform lumping and apply stabilisation
		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		    
		  }
		}
		else {
		  // perform no lumping but apply stabilisation
		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij)  = IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji)  = IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij)  = IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji)  = IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		}
	      } else {
		if (blumping) {
		  // perform lumping but do not apply stabilisation
		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  }
		}
		else {
		  // perform no lumping and do not apply stabilisation
		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) = IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) = IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) = IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) = IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		}
	      }
	    }
  
  /*****************************************************************************
   * Scatter data for datalen edges with nelem elements per node;
   * matrix positions are retrieved from edge list
   ****************************************************************************/
  template <int datalen,
            bool btranspose,
            bool boverwrite,
            bool bstabilise,
            bool blumping,
            typename Tm,
            typename Td,
	    typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti ipos,
				 Ti iedge,
				 Ti nedge,
				 Ti na)
	    {
	      if (bstabilise) {
		if (blumping) {
		  // perform lumping and apply stabilisation
		  Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		  Ti jj = IDX2(IedgeList,6,iedge,6,nedge);

		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen);
		      }
		    
		  }
		} 
		else {
		  // perform no lumping but apply stabilisation
		  Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		  Ti ji = IDX2(IedgeList,4,iedge,6,nedge);
		  Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		  Ti jj = IDX2(IedgeList,6,iedge,6,nedge);

		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij)  = IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji)  = IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij)  = IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji)  = IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3T(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3T(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3T(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ii) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),jj) -= IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3(DataAtEdge,i,1,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3(DataAtEdge,i,2,ipos,nelem,3,datalen)+
			                                 IDX3(DataAtEdge,i,3,ipos,nelem,3,datalen);
			}
		    }
		  }
		}
	      }
	      else {
		if (blumping) {
		  // perform lumping but do not apply stabilisation
		  Ti ii = IDX2(IedgeList,5,iedge,6,nedge);
		  Ti jj = IDX2(IedgeList,6,iedge,6,nedge);

		  if (btranspose) {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  }
		  else {
#pragma unroll
		    for (int i=1; i<=nelem; i++)
		      {
			IDX1((((Tm**)mat)[i-1]),ii) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			IDX1((((Tm**)mat)[i-1]),jj) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
		      }
		  } 
		}
		else {
		  // perform no lumping and do not apply stabilisation
		  Ti ij = IDX2(IedgeList,3,iedge,6,nedge);
		  Ti ji = IDX2(IedgeList,4,iedge,6,nedge);
		  
		  if (boverwrite) {
		    //
		    // Overwrite existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) = IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) = IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) = IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) = IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		  else {
		    //
		    // Keep existing data
		    //
		    if (btranspose) {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3T(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3T(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		    else {
#pragma unroll
		      for (int i=1; i<=nelem; i++)
			{
			  IDX1((((Tm**)mat)[i-1]),ij) += IDX3(DataAtEdge,i,1,ipos,nelem,2,datalen);
			  IDX1((((Tm**)mat)[i-1]),ji) += IDX3(DataAtEdge,i,2,ipos,nelem,2,datalen);
			}
		    }
		  }
		}
	      }
	    }
};

/*******************************************************************************
 * Matrix - two-dimensional array
 ******************************************************************************/
template <int nelem, bool soa>
  struct Matrix : public MatrixBase<nelem,soa>
{
  // Enable use of inherited functions
  using MatrixBase<nelem,soa>::scatterNodeData;
  using MatrixBase<nelem,soa>::scatterEdgeData;
  
  /*****************************************************************************
   * Scatter nodal data for a single nodes with nelem elements per node
   ****************************************************************************/
  template <bool boverwrite, typename Tm, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterNodeData (Tm *mat,
				 Td *DataAtNode,
				 Ti ia,
				 Ti na)
  {
    MatrixBase<nelem,soa>::scatterNodeData<1,false,boverwrite>
      (mat,DataAtNode,1,ia,na);
  }

  /*****************************************************************************
   * Scatter data for a single edges with nelem elements per node;
   * off-diagonal matrix positions are given explicitly
   ****************************************************************************/
  template <bool boverwrite, bool bstabilise, bool blumping, typename Tm, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti ii,
				 Ti ij,
				 Ti ji,
				 Ti jj,
				 Ti na)
  {
    MatrixBase<nelem,soa>::scatterEdgeData<1,false,boverwrite,bstabilise,blumping>
      (mat,DataAtEdge,1,ii,ij,ji,jj,na);
  }

  /*****************************************************************************
   * Scatter data for a single edges with nelem elements per node;
   * matrix positions are retrieved from edge list
   ****************************************************************************/
  template <bool boverwrite, bool bstabilise, bool blumping, typename Tm, typename Td, typename Ti>
    __device__ __forceinline__
    static void scatterEdgeData (Tm *mat,
				 Td *DataAtEdge,
				 Ti *IedgeList,
				 Ti iedge,
				 Ti nedge,
				 Ti na)
  {
    MatrixBase<nelem,soa>::scatterEdgeData<1,false,boverwrite,bstabilise,blumping>
      (mat,DataAtEdge,IedgeList,1,iedge,na,nedge);
  }
};
