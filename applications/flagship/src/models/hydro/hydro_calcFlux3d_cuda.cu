/*#############################################################################
******************************************************************************
* <name> hydro_calcFlux3d_cuda </name>
******************************************************************************
*
* <purpose>
* This file provides CUDA kernels to compute the fluxes for the low-order
* scheme in 3D using different types if artificial viscosities.
* </purpose>
*
*#############################################################################
*/

#include <stdio.h>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <coproc_core.h>
#include <coproc_storage_cuda.h>
#include "../../cudaGatherScatter.h"
#ifdef HAS_INLINE_PTX
#include "../../cudaDMA.h"
#endif

#define LANGUAGE LANGUAGE_C
#include "../../flagship.h"
#include "../../cudaMacros.h"

#include "hydro.h"

// Define CUDA kernel which does not make use of the CUDADMA library
// and is applied to the remaining edges which are not processed in groups
// #define BASELINE_KERNEL hydro_calcFlux3d_shmem
#define BASELINE_KERNEL hydro_calcFlux3d_baseline

// Defines for baseline implementation
#define BASELINE_THREADS_PER_CTA  32*2
#define BASELINE_NEDGE_PER_THREAD 1

// Defines for shared memory implementation
#define SHMEM_DATA_TRANSPOSE   true
#define SHMEM_DATA_IDX3        IDX3T
#define SHMEM_NEDGE_PER_THREAD BASELINE_NEDGE_PER_THREAD

#ifdef HAS_INLINE_PTX
// Define CUDA kernel which makes use of the CUDADMA library to achive
// higher throughput between global and shared memory on the device
// #define CUDADMA_PREFETCH_SINGLE
#endif

// Defines for cudaDMA implementation without warp specialisation
#ifdef CUDADMA_NOSPEC
#define CUDADMA_KERNEL                  hydro_calcFlux3d_cudaDMA_nospec
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*2
#define CUDADMA_THREADS_PER_LD          0
#define CUDADMA_NEDGE_PER_THREAD        1
#define CUDADMA_DMA_LDS_IND             0
#define CUDADMA_DMA_LDS_SRC             0
#define CUDADMA_DMA_LDS_DEST            0
#define CUDADMA_DMA_LDS_COEFF           0
#define CUDADMA_DMA_LDS                 0
#endif

// Defines for cudaDMA single buffer implementation with prefetching of indices
#ifdef CUDADMA_PREFETCH_SINGLE
#define CUDADMA_KERNEL                  hydro_calcFlux3d_cudaDMA_prefetch_single
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*4
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        1*1
#define CUDADMA_DMA_LDS_IND             0
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           1
#define CUDADMA_DMA_LDS                 (CUDADMA_DMA_LDS_IND +	\
                                         CUDADMA_DMA_LDS_SRC +	\
                                         CUDADMA_DMA_LDS_DEST +	\
										 CUDADMA_DMA_LDS_COEFF)
#endif

// Defines for cudaDMA double buffer implementation with prefetching of indices
#ifdef CUDADMA_PREFETCH_DOUBLE
#define CUDADMA_KERNEL                  hydro_calcFlux3d_cudaDMA_prefetch_double
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*4
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        4*1
#define CUDADMA_DMA_LDS_IND             0
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           1
#define CUDADMA_DMA_LDS                 (CUDADMA_DMA_LDS_IND +	\
                                         CUDADMA_DMA_LDS_SRC +	\
                                         CUDADMA_DMA_LDS_DEST +	\
										 CUDADMA_DMA_LDS_COEFF)
#endif

// Defines for cudaDMA double buffer implementation
#ifdef CUDADMA_DOUBLE
#define CUDADMA_KERNEL                  hydro_calcFlux3d_cudaDMA_double
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*2
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        6*1
#define CUDADMA_DMA_LDS_IND             1
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           0
#define CUDADMA_DMA_LDS                 (3*CUDADMA_DMA_LDS_IND +	\
                                         2*CUDADMA_DMA_LDS_SRC +	\
                                         2*CUDADMA_DMA_LDS_DEST)
#endif

// Defines for cudaDMA manual buffer implementation
#ifdef CUDADMA_MANUAL
#define CUDADMA_KERNEL                  hydro_calcFlux3d_cudaDMA_manual
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*2
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        6*1
#define CUDADMA_DMA_LDS_IND             1
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           1
#define CUDADMA_DMA_LDS                 (CUDADMA_DMA_LDS_IND +	\
                                         CUDADMA_DMA_LDS_SRC +	\
                                         CUDADMA_DMA_LDS_DEST +	\
										 CUDADMA_DMA_LDS_COEFF)
#endif

// Defines for empty cudaDMA implementation
#ifndef CUDADMA_KERNEL
#define CUDADMA_COMPUTE_THREADS_PER_CTA 0
#define CUDADMA_THREADS_PER_LD          0
#define CUDADMA_NEDGE_PER_THREAD        0
#define CUDADMA_DMA_LDS_IND             0
#define CUDADMA_DMA_LDS_SRC             0
#define CUDADMA_DMA_LDS_DEST            0
#define CUDADMA_DMA_LDS_COEFF           0
#define CUDADMA_DMA_LDS                 0
#endif

using namespace std;

namespace hydro3d_cuda
{

  /*****************************************************************************
   * FluxBase: Compute inviscid fluxes for nedgesim edges
   ****************************************************************************/

  struct InviscidFluxBase
  {
    /*
     * Calculate the inviscid flux for x-direction (not skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcFluxXdir(Td *Fxi,
							 Td *Fxj,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc)
    {
      if (boverwrite) {
		// Overwrite destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fxi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2T(Fxj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fxi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2T(Fxj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fxi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2(Fxj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fxi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2(Fxj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
      }
      else {
		// Keep content of destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fxi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2T(Fxj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fxi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2T(Fxi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2T(Fxj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fxj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fxi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2(Fxj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fxi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
			IDX2(Fxi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi);
	    
			IDX2(Fxj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fxj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
      }
    }

    /*
     * Calculate the inviscid flux for y-direction (not skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcFluxYdir(Td *Fyi,
							 Td *Fyj,
							 Td *DataAtEdge,
							 Td vi,
							 Td vj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc)
    {
      if (boverwrite) {
		// Overwrite destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fyi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2T(Fyj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fyi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2T(Fyj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fyi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2(Fyj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fyi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2(Fyj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
      }
      else {
		// Keep content of destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fyi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2T(Fyj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fyi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2T(Fyi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2T(Fyj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fyj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fyi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2(Fyj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fyi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
			IDX2(Fyi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi);
	
			IDX2(Fyj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fyj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
      }
    }

    /*
     * Calculate the inviscid flux for z-direction (not skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcFluxZdir(Td *Fzi,
							 Td *Fzj,
							 Td *DataAtEdge,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc)
    {
      if (boverwrite) {
		// Overwrite destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fzi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2T(Fzj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fzi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2T(Fzj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fzi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2(Fzj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fzi,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2(Fzj,1,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,2,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,3,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,4,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,5,iposDest,NVAR3D,nedgesimDest) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
      }
      else {
		// Keep content of destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fzi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2T(Fzj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fzi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2T(Fzi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2T(Fzj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fzj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fzi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2(Fzj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fzi,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
			IDX2(Fzi,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi);
	
			IDX2(Fzj,1,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,2,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,3,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,4,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fzj,5,iposDest,NVAR3D,nedgesimDest) += INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
      }
    }

    /*
     * Calculate the inviscid flux for x-direction (skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcFluxXdir(Td *Fx_ij,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc)
    {
      if (boverwrite) {
		// Overwrite destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3t,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
      }
      else {
		// Keep content of destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3t,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2T(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fx_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
			IDX2(Fx_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
			  INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj);
		  }
		}
      }
    }
    
    /*
     * Calculate the inviscid flux for y-direction (skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcFluxYdir(Td *Fy_ij,
							 Td *DataAtEdge,
							 Td vi,
							 Td vj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc)
    {
      if (boverwrite) {
		// Overwrite destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
      }
      else {
		// Keep content of destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2T(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fy_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
			IDX2(Fy_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
			  INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj);
		  }
		}
      }
    }

    /*
     * Calculate the inviscid flux for z-direction (skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcFluxZdir(Td *Fz_ij,
							 Td *DataAtEdge,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc)
    {
      if (boverwrite) {
		// Overwrite destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) =
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
      }
      else {
		// Keep content of destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX2T(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Destination vector is transposed
			IDX2T(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2T(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX2(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		  else {
			// Both vectors are not transposed
			IDX2(Fz_ij,1,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,2,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,3,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,4,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
			IDX2(Fz_ij,5,iposDest,NVAR3D,nedgesimDest) +=
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
			  INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj);
		  }
		}
      }
    }

    /*
     * Calculate the inviscid fluxes in all directions (not skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *Fx_ij,
							 Td *Fy_ij,
							 Td *Fz_ij,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc)
    {
      // Compute Galerkin flux difference for x-direction
      InviscidFluxBase::
		calcFluxXdir<nedgesimDest,nedgesimSrc,btransposeDest,btransposeSrc,boverwrite>
		(Fx_ij,DataAtEdge,ui,uj,pi,pj,iposDest,iposSrc);
      
      // Compute Galerkin flux difference for y-direction
      InviscidFluxBase::
		calcFluxYdir<nedgesimDest,nedgesimSrc,btransposeDest,btransposeSrc,boverwrite>
		(Fy_ij,DataAtEdge,vi,vj,pi,pj,iposDest,iposSrc);

      // Compute Galerkin flux difference for z-direction
      InviscidFluxBase::
		calcFluxZdir<nedgesimDest,nedgesimSrc,btransposeDest,btransposeSrc,boverwrite>
		(Fz_ij,DataAtEdge,wi,wj,pi,pj,iposDest,iposSrc);
    }
   
    /*
     * Calculate the inviscid fluxes in all directions (skew-symmetric)
     * and multiply them by the precomputed finite element coefficients
     */
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  bool boverwrite, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *FluxesAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Td scale,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge, 
							 Ti nedge,
							 Ti ncoeff)
    {
      Td aux;
      
#ifdef HYDRO_USE_IBP
      // Calculate skew-symmetric inviscid fluxes

      // Flux component 1
      if (btransposeSrc) {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }
      else {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }

      if (boverwrite) {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3T(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
		else {
		  IDX3(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
      }
      else {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3T(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
		else {
		  IDX3(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
      }

      // Flux component 2
      if (btransposeSrc) {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }
      else {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }
      
      if (boverwrite) {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3T(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
		else {
		  IDX3(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
      }
      else {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3T(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
		else {
		  IDX3(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
      }

      // Flux component 3
      if (btransposeSrc) {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }
      else {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }
      
      if (boverwrite) {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3T(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
		else {
		  IDX3(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
      }
      else {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3T(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
		else {
		  IDX3(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
      }

      // Flux component 4
      if (btransposeSrc) {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }
      else {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,3,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }

      if (boverwrite) {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3T(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
		else {
		  IDX3(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
      }
      else {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3T(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
		else {
		  IDX3(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
      }

      // Flux component 5
      if (btransposeSrc) {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }
      else {
		aux = scale *
		  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj)
		   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)
		   -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
		   INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,3,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi));
      }

      if (boverwrite) {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3T(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
		else {
		  IDX3(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) =  aux;
		  IDX3(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) = -aux;
		}
      }
      else {
		if (btransposeDest) {
		  IDX3T(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3T(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
		else {
		  IDX3(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) += aux;
		  IDX3(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) -= aux;
		}
      }
#else
      // Calculate inviscid fluxes (not skew-symmetric)

      if (boverwrite) {
		// Overwrite destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX3T(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));


			IDX3T(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		  else {
			// Destination vector is transposed
			IDX3T(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));


			IDX3T(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX3(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));


			IDX3(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		  else {
			// Both vectors are not transposed
			IDX3(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));


			IDX3(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		}
      }
      else {
		// Keep content of destination vector
		if (btransposeDest) {
		  if (btransposeSrc) {
			// Both source and destination vector are transposed
			IDX3T(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));


			IDX3T(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		  else {
			// Destination vector is transposed
			IDX3T(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3T(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));


			IDX3T(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3T(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		}
		else {
		  if (btransposeSrc) {
			// Source vector is transposed
			IDX3(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));


			IDX3(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		  else {
			// Both vectors are not transposed
			IDX3(FluxesAtEdge,1,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,1,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,2,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,2,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,3,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,3,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
	    
			IDX3(FluxesAtEdge,4,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,4,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));

	    
			IDX3(FluxesAtEdge,5,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
	    
			IDX3(FluxesAtEdge,5,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,ui,pi)-
				INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,uj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,vi,pi)-
				INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,vj,pj))
			   +IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*
			   (INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc,wi,pi)-
				INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc,wj,pj)));
		  }
		}
      }
#endif
    }
  };

  /*****************************************************************************
   * InviscidFlux
   ****************************************************************************/

  struct InviscidFlux : public InviscidFluxBase
  {
    // Enable use of inherited functions
    using InviscidFluxBase::calcFluxXdir;
    using InviscidFluxBase::calcFluxYdir;
    using InviscidFluxBase::calcFluxZdir;
    using InviscidFluxBase::calcEdgeData;

    /**************************************************************************
     * Wrapper routine for processing a single edge
     *************************************************************************/

    /*
     * Calculate the inviscid flux for x-direction (not skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcFluxXdir(Td *Fxi,
							 Td *Fxj,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcFluxXdir<1,1,false,false,boverwrite>
		(Fxi,Fxj,DataAtEdge,ui,uj,pi,pj,1,1);
    }
    
    /*
     * Calculate the inviscid flux for y-direction (not skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcFluxYdir(Td *Fyi,
							 Td *Fyj,
							 Td *DataAtEdge,
							 Td vi,
							 Td vj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcFluxYdir<1,1,false,false,boverwrite>
		(Fyi,Fyj,DataAtEdge,vi,vj,pi,pj,1,1);
    }

    /*
     * Calculate the inviscid flux for z-direction (not skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcFluxZdir(Td *Fzi,
							 Td *Fzj,
							 Td *DataAtEdge,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcFluxZdir<1,1,false,false,boverwrite>
		(Fzi,Fzj,DataAtEdge,wi,wj,pi,pj,1,1);
    }

    /*
     * Calculate the inviscid flux for x-direction (skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcFluxXdir(Td *Fx_ij,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcFluxXdir<1,1,false,false,boverwrite>
		(Fx_ij,DataAtEdge,ui,uj,pi,pj,1,1);
    }
    
    /*
     * Calculate the inviscid flux for y-direction (skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcFluxYdir(Td *Fy_ij,
							 Td *DataAtEdge,
							 Td vi,
							 Td vj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcFluxYdir<1,1,false,false,boverwrite>
		(Fy_ij,DataAtEdge,vi,vj,pi,pj,1,1);
    }

    /*
     * Calculate the inviscid flux for z-direction (skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcFluxZdir(Td *Fz_ij,
							 Td *DataAtEdge,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcFluxZdir<1,1,false,false,boverwrite>
		(Fz_ij,DataAtEdge,wi,wj,pi,pj,1,1);
    }

    /*
     * Calculate the inviscid fluxes in all directions (not skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcEdgeData(Td *Fxi,
							 Td *Fxj,
							 Td *Fyi,
							 Td *Fyj,
							 Td *Fzi,
							 Td *Fzj,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcEdgeData<1,1,false,false,boverwrite>
		(Fxi,Fxj,Fyi,Fyj,Fzi,Fzj,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,1,1);
    }

    /*
     * Calculate the inviscid fluxes in all directions (skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcEdgeData(Td *Fx_ij,
							 Td *Fy_ij,
							 Td *Fz_ij,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj)
    {
      InviscidFluxBase::calcEdgeData<1,1,false,false,boverwrite>
		(Fx_ij,Fy_ij,Fz_ij,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,1,1);
    }

    /*
     * Calculate the inviscid fluxes in all directions (not skew-symmetric)
     * and multiply them by the precomputed finite element coefficients
     */
    template <bool boverwrite, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *FluxesAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Td scale,
							 Ti iedge, 
							 Ti nedge,
							 Ti ncoeff)
    {
      InviscidFluxBase::calcEdgeData<1,1,false,false,boverwrite>
		(FluxesAtEdge,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,
		 scale,1,1,iedge,nedge,ncoeff);
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase (basic functionality individual specialisations)
   ****************************************************************************/

  template <int idissipationtype>
  struct InviscidFluxDissipationBase
  {
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing zero
   * artificial dissipation, aka standard Galerkin
   ****************************************************************************/

  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_ZERO>
  {
    template <int nedgesimDest, int nedgesimSrc,
			  bool btransposeDest, bool btransposeSrc,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge,
							 Ti nedge,
							 Ti ncoeff)
    {
      if (btransposeDest) {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) = 0.0;
      }
      else {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) = 0.0;
      }
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing scalar
   * artificial dissipation proportional to the spectral radius
   * (largest eigenvector) of the cumulative Roe matrix.
   ****************************************************************************/

  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_SCALAR>
  {
    template <int nedgesimDest, int nedgesimSrc,
			  bool btransposeDest, bool btransposeSrc,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge,
							 Ti nedge,
							 Ti ncoeff)
    {
      Td ri,rj,hi,hj;
      if (btransposeSrc) {
		// Compute densities
		ri = DENSITY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc);
		rj = DENSITY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc);
	
		// Compute enthalpies
		hi = (TOTALENERGY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		hj = (TOTALENERGY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
      }
      else {
		// Compute densities
		ri = DENSITY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc);
		rj = DENSITY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc);
	
		// Compute enthalpies
		hi = (TOTALENERGY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		hj = (TOTALENERGY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
      }

      // Compute Roe mean values
      Td aux  = ROE_MEAN_RATIO(ri,rj);
      Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
      Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);

      // Compute skew-symmetric coefficient
      Td a[3];
      a[0] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge));
      a[1] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge));
      a[2] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge));
      
      // Compute auxiliary variables
      Td q_ij   = DCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
      Td vel_ij = u_ij * a[0] + v_ij * a[1] + w_ij * a[2];
          
      // Compute the speed of sound
      Td c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), DBL_EPSILON));
    
      // Compute scalar dissipation
      Td d_ij = abs(vel_ij) + sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])*c_ij;
    
      // Multiply the solution difference by the scalar dissipation
      if (btransposeDest) {
		if (btransposeSrc) {
		  // Both source and destination vector are transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Destination vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
      else {
		if (btransposeSrc) {
		  // Source vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Both vectors are not transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing scalar
   * artificial dissipation proportional to the spectral radius
   * (largest eigenvector) of the dimensional-split Roe matrix.
   ****************************************************************************/

  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_SCALAR_DSPLIT>
  {
    template <int nedgesimDest, int nedgesimSrc,
			  bool btransposeDest, bool btransposeSrc,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge,
							 Ti nedge,
							 Ti ncoeff)
    {
      Td ri,rj,hi,hj;
      if (btransposeSrc) {
		// Compute densities
		ri = DENSITY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc);
		rj = DENSITY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc);
	
		// Compute enthalpies
		hi = (TOTALENERGY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		hj = (TOTALENERGY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
      }
      else {
		// Compute densities
		ri = DENSITY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc);
		rj = DENSITY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc);
	
		// Compute enthalpies
		hi = (TOTALENERGY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		hj = (TOTALENERGY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
      }
    
      // Compute Roe mean values
      Td aux  = ROE_MEAN_RATIO(ri,rj);
      Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
      Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);

      // Compute auxiliary variables
      Td q_ij = DCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
    
      // Compute the speed of sound
      Td c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), DBL_EPSILON));

      // Compute skew-symmetric coefficient
      Td a[3];
      a[0] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge));
      a[1] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge));
      a[2] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge));
        
      // Compute scalar dissipation
      Td d_ij = ( abs(a[0]*u_ij) + abs(a[0])*c_ij +
				  abs(a[1]*v_ij) + abs(a[1])*c_ij +
				  abs(a[2]*w_ij) + abs(a[2])*c_ij );
    
      // Multiply the solution difference by the scalar dissipation
      if (btransposeDest) {
		if (btransposeSrc) {
		  // Both source and destination vector are transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Destination vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
      else {
		if (btransposeSrc) {
		  // Source vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Both vectors are not transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing
   * tensorial artificial dissipation of Roe-type.
   ****************************************************************************/

  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_ROE>
  {
    template <int nedgesimDest, int nedgesimSrc,
			  bool btransposeDest, bool btransposeSrc,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge, 
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi, 
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge, 
							 Ti nedge,
							 Ti ncoeff)
    {
      // Compute skew-symmetric coefficient
      Td a[3];
      a[0] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge));
      a[1] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge));
      a[2] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge));
      Td anorm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    
      if (anorm > DBL_EPSILON) {
      
		// Normalise the skew-symmetric coefficient
		a[0] = a[0]/anorm;
		a[1] = a[1]/anorm;
		a[2] = a[2]/anorm;
      
		Td ri,rj,hi,hj;
		if (btransposeSrc) {
		  // Compute densities
		  ri = DENSITY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc);
		  rj = DENSITY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc);
	  
		  // Compute enthalpies
		  hi = (TOTALENERGY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		  hj = (TOTALENERGY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
		}
		else {
		  // Compute densities
		  ri = DENSITY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc);
		  rj = DENSITY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc);
	  
		  // Compute enthalpies
		  hi = (TOTALENERGY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		  hj = (TOTALENERGY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
		}
      
		// Compute Roe mean values
		Td aux  = ROE_MEAN_RATIO(ri,rj);
		Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
		Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
		Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
		Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
		// Compute auxiliary variables
		Td vel_ij = u_ij * a[0] + v_ij * a[1] + w_ij * a[2];
		Td q_ij   = DCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
	
		// Compute the speed of sound
		Td c2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), DBL_EPSILON);
		Td c_ij  = sqrt(c2_ij);
	
		// Compute eigenvalues
		Td l1 = abs(vel_ij-c_ij);
		Td l2 = abs(vel_ij);
		Td l3 = abs(vel_ij+c_ij);
		Td l4 = abs(vel_ij);
		Td l5 = abs(vel_ij);
	
		// Compute solution difference U_j-U_i
		Td Diff[NVAR3D];
		if (btransposeSrc) {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}
		else {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}
		// Compute auxiliary quantities for characteristic variables
		Td aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*Diff[0]
											   -u_ij*Diff[1]
											   -v_ij*Diff[2]
											   -w_ij*Diff[3]
											   +Diff[4])/DCONST(2.0)/c2_ij;
		Td aux2 = (vel_ij*Diff[0]
				   -a[0]*Diff[1]
				   -a[1]*Diff[2]
				   -a[2]*Diff[3])/DCONST(2.0)/c_ij;
      
		// Get the dimension with largest coefficient
		if (a[0] >= a[1] && a[0] >= a[2]) {
	
		  // Compute characteristic variables multiplied by the corresponding eigenvalue
		  Td w1 = l1 * (aux1 + aux2);
		  Td w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff[0]
						+((HYDRO_GAMMA)-DCONST(1.0))*( u_ij*Diff[1]
													   +v_ij*Diff[2]
													   +w_ij*Diff[3]
													   -Diff[4])/c2_ij);
		  Td w3 = l3 * (aux1 - aux2);
		  Td w4 = l4 * ( (v_ij-vel_ij*a[1])/a[0]*Diff[0]
						 +a[1]*Diff[1]
						 +(a[1]*a[1]-DCONST(1.0))/a[0]*Diff[2]
						 +a[1]*a[2]/a[0]*Diff[3]);
		  Td w5 = l5 * ( (vel_ij*a[2]-w_ij)/a[0]*Diff[0]
						 -a[2]*Diff[1]
						 -a[1]*a[2]/a[0]*Diff[2]
						 +(DCONST(1.0)-a[2]*a[2])/a[0]*Diff[3]);
	
		  // Compute "R_ij * |Lbd_ij| * L_ij * dU"
		  if (btransposeDest) {
			IDX2T(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
			IDX2T(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
																		   (u_ij+c_ij*a[0])*w3 + a[1]*w4 - a[2]*w5 );
			IDX2T(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
																		   (v_ij+c_ij*a[1])*w3 - a[0]*w4 );
			IDX2T(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
																		   (w_ij+c_ij*a[2])*w3 + a[0]*w5 );
			IDX2T(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
																		   + (u_ij*a[1]-v_ij*a[0])*w4 + (w_ij*a[0]-u_ij*a[2])*w5 );
		  }
		  else {
			IDX2(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
			IDX2(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
																		  (u_ij+c_ij*a[0])*w3 + a[1]*w4 - a[2]*w5 );
			IDX2(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
																		  (v_ij+c_ij*a[1])*w3 - a[0]*w4 );
			IDX2(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
																		  (w_ij+c_ij*a[2])*w3 + a[0]*w5 );
			IDX2(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
																		  + (u_ij*a[1]-v_ij*a[0])*w4 + (w_ij*a[0]-u_ij*a[2])*w5 );
		  }
	  
		} else if (a[1] >= a[0] && a[1] >= a[2]) {
		  // Compute characteristic variables multiplied by the corresponding eigenvalue
		  Td w1 = l1 * (aux1 + aux2);
		  Td w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff[0]
						+((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*Diff[1]
													  +v_ij*Diff[2]
													  +w_ij*Diff[3]
													  -Diff[4])/c2_ij);
		  Td w3 = l3 * (aux1 - aux2);
		  Td w4 = l4 * ( (vel_ij*a[0]-u_ij)/a[1]*Diff[0]
						 +(DCONST(1.0)-a[0]*a[0])/a[1]*Diff[1]
						 -a[0]*Diff[2]
						 -a[0]*a[2]/a[1]*Diff[3]);
		  Td w5 = l5 * ( (w_ij-vel_ij*a[2])/a[1]*Diff[0]
						 +a[0]*a[2]/a[1]*Diff[1]
						 +a[2]*Diff[2]
						 +(a[2]*a[2]-DCONST(1.0))/a[1]*Diff[3]);
	  
		  // Compute "R_ij * |Lbd_ij| * L_ij * dU"
		  if (btransposeDest) {
			IDX2T(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
			IDX2T(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
																		   (u_ij+c_ij*a[0])*w3 + a[1]*w4 );
			IDX2T(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
																		   (v_ij+c_ij*a[1])*w3 - a[0]*w4 + a[2]*w5 );
			IDX2T(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
																		   (w_ij+c_ij*a[2])*w3 - a[1]*w5 );
			IDX2T(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
																		   + (u_ij*a[1]-v_ij*a[0])*w4 + (v_ij*a[2]-w_ij*a[1])*w5 );
		  }
		  else {
			IDX2(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
			IDX2(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
																		  (u_ij+c_ij*a[0])*w3 + a[1]*w4 );
			IDX2(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
																		  (v_ij+c_ij*a[1])*w3 - a[0]*w4 + a[2]*w5 );
			IDX2(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
																		  (w_ij+c_ij*a[2])*w3 - a[1]*w5 );
			IDX2(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
																		  + (u_ij*a[1]-v_ij*a[0])*w4 + (v_ij*a[2]-w_ij*a[1])*w5 );
		  }
	
		} else {
		  // Compute characteristic variables multiplied by the corresponding eigenvalue
		  Td w1 = l1 * (aux1 + aux2);
		  Td w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff[0]
						+((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*Diff[1]
													  +v_ij*Diff[2]
													  +w_ij*Diff[3]
													  -Diff[4])/c2_ij);
		  Td w3 = l3 * (aux1 - aux2);
		  Td w4 = l4 * ( (u_ij-vel_ij*a[0])/a[2]*Diff[0]
						 +(a[0]*a[0]-DCONST(1.0))/a[2]*Diff[1]
						 +a[0]*a[1]/a[2]*Diff[2]
						 +a[0]*Diff[3]);
		  Td w5 = l5 * ( (vel_ij*a[1]-v_ij)/a[2]*Diff[0]
						 -a[0]*a[1]/a[2]*Diff[1]
						 +(DCONST(1.0)-a[1]*a[1])/a[2]*Diff[2]
						 -a[1]*Diff[3]);
	  
		  // Compute "R_ij * |Lbd_ij| * L_ij * dU"
		  if (btransposeDest) {
			IDX2T(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
			IDX2T(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
																		   (u_ij+c_ij*a[0])*w3 - a[2]*w4 );
			IDX2T(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
																		   (v_ij+c_ij*a[1])*w3 + a[2]*w5 );
			IDX2T(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
																		   (w_ij+c_ij*a[2])*w3 + a[0]*w4 - a[1]*w5);
			IDX2T(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
																		   + (w_ij*a[0]-u_ij*a[2])*w4 + (v_ij*a[2]-w_ij*a[1])*w5 );
		  }
		  else {
			IDX2(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
			IDX2(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
																		  (u_ij+c_ij*a[0])*w3 - a[2]*w4 );
			IDX2(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
																		  (v_ij+c_ij*a[1])*w3 + a[2]*w5 );
			IDX2(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
																		  (w_ij+c_ij*a[2])*w3 + a[0]*w4 - a[1]*w5);
			IDX2(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
																		  + (w_ij*a[0]-u_ij*a[2])*w4 + (v_ij*a[2]-w_ij*a[1])*w5 );
		  }
		}
      } else {
		if (btransposeDest) {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) = 0.0;
		}
		else {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) = 0.0;
		}
      }
    }
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing
   * tensorial artificial dissipation of Roe-type using dimensional splitting.
   ****************************************************************************/
  
  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_ROE_DSPLIT>
  {
    template <int nedgesimDest, int nedgesimSrc,
			  bool btransposeDest, bool btransposeSrc,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge, 
							 Ti nedge,
							 Ti ncoeff)
    {
      // Compute skew-symmetric coefficient
      Td a[3];
      a[0] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge));
      a[1] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge));
      a[2] = DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
						  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge));
      Td anorm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
      
      if (anorm > DBL_EPSILON) {
      
		// Compute the absolute value
		a[0] = abs(a[0]);
		a[1] = abs(a[1]);
		a[2] = abs(a[2]);
      
		Td ri,rj,hi,hj;
		if (btransposeSrc) {
		  // Compute densities
		  ri = DENSITY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc);
		  rj = DENSITY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc);
	  
		  // Compute enthalpies
		  hi = (TOTALENERGY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		  hj = (TOTALENERGY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
		}
		else {
		  // Compute densities
		  ri = DENSITY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc);
		  rj = DENSITY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc);
	  
		  // Compute enthalpies
		  hi = (TOTALENERGY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc)+pi)/ri;
		  hj = (TOTALENERGY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc)+pj)/rj;
		}
      
		// Compute Roe mean values
		Td aux  = ROE_MEAN_RATIO(ri,rj);
		Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
		Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
		Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
		Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
		// Compute auxiliary variable
		Td q_ij = DCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
      
		// Compute the speed of sound
		Td c2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), DBL_EPSILON);
		Td c_ij  = sqrt(c2_ij);
      
		//----------------------------------------------------------------------
		// Dimensional splitting: x-direction
		//----------------------------------------------------------------------
      
		// Compute eigenvalues
		Td l1 = abs(u_ij-c_ij);
		Td l2 = abs(u_ij);
		Td l3 = abs(u_ij+c_ij);
		Td l4 = abs(u_ij);
		Td l5 = abs(u_ij);
      
		// Compute solution difference U_j-U_i
		Td Diff[NVAR3D];
		if (btransposeSrc) {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}
		else {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}

		// Compute auxiliary quantities for characteristic variables
		Td aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*Diff[0]
											   -u_ij*Diff[1]
											   -v_ij*Diff[2]
											   -w_ij*Diff[3]
											   +Diff[4])/DCONST(2.0)/c2_ij;
		Td aux2 = (u_ij*Diff[0]-Diff[1])/DCONST(2.0)/c_ij;
      
		// Compute characteristic variables multiplied by the corresponding eigenvalue
		Td w1 = l1 * (aux1 + aux2);
		Td w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff[0]
					  +((HYDRO_GAMMA)-DCONST(1.0))*( u_ij*Diff[1]
													 +v_ij*Diff[2]
													 +w_ij*Diff[3]
													 -Diff[4])/c2_ij);
		Td w3 = l3 * (aux1 - aux2);
		Td w4 = l4 * ( v_ij*Diff[0]-Diff[2]);
		Td w5 = l5 * (-w_ij*Diff[0]+Diff[3]);
      
		// Compute "R_ij * |Lbd_ij| * L_ij * dU"
		if (btransposeDest) {
		  IDX2T(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = a[0] * ( w1 + w2 + w3 );
		  IDX2T(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = a[0] * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 );
		  IDX2T(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = a[0] * ( v_ij*(w1 + w2 + w3) - w4 );
		  IDX2T(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = a[0] * ( w_ij*(w1 + w2 + w3) + w5 );
		  IDX2T(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = a[0] * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 + (H_ij+c_ij*u_ij)*w3
																		-v_ij*w4 + w_ij*w5 );
		}
		else {
		  IDX2(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) = a[0] * ( w1 + w2 + w3 );
		  IDX2(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) = a[0] * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 );
		  IDX2(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) = a[0] * ( v_ij*(w1 + w2 + w3) - w4 );
		  IDX2(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) = a[0] * ( w_ij*(w1 + w2 + w3) + w5 );
		  IDX2(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) = a[0] * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 + (H_ij+c_ij*u_ij)*w3
																	   -v_ij*w4 + w_ij*w5 );
		}
	
		//----------------------------------------------------------------------
		// Dimensional splitting: y-direction
		//----------------------------------------------------------------------
      
		// Compute eigenvalues
		l1 = abs(v_ij-c_ij);
		l2 = abs(v_ij);
		l3 = abs(v_ij+c_ij);
		l4 = abs(v_ij);
		l5 = abs(v_ij);
      
		// Compute solution difference U_j-U_i
		if (btransposeSrc) {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}
		else {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}
	      
		// Compute auxiliary quantities for characteristic variables
		aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*Diff[0]
											-u_ij*Diff[1]
											-v_ij*Diff[2]
											-w_ij*Diff[3]
											+Diff[4])/DCONST(2.0)/c2_ij;
		aux2 = (v_ij*Diff[0]-Diff[2])/DCONST(2.0)/c_ij;
      
		// Compute characteristic variables multiplied by the corresponding eigenvalue
		w1 = l1 * (aux1 + aux2);
		w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff[0]
				   +((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*Diff[1]
												 +v_ij*Diff[2]
												 +w_ij*Diff[3]
												 -Diff[4])/c2_ij);
		w3 = l3 * (aux1 - aux2);
		w4 = l4 * (-u_ij*Diff[0]+Diff[1]);
		w5 = l5 * ( w_ij*Diff[0]-Diff[3]);
      
		// Compute "R_ij * |Lbd_ij| * L_ij * dU"
		if (btransposeDest) {
		  IDX2T(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) += a[1] * ( w1 + w2 + w3 );
		  IDX2T(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) += a[1] * ( u_ij*(w1 + w2 + w3) + w4 );
		  IDX2T(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) += a[1] * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 );
		  IDX2T(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) += a[1] * ( w_ij*(w1 + w2 + w3) - w5 );
		  IDX2T(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) += a[1] * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 + (H_ij+c_ij*v_ij)*w3
																		 + u_ij*w4 -w_ij*w5 );
		}
		else {
		  IDX2(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) += a[1] * ( w1 + w2 + w3 );
		  IDX2(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) += a[1] * ( u_ij*(w1 + w2 + w3) + w4 );
		  IDX2(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) += a[1] * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 );
		  IDX2(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) += a[1] * ( w_ij*(w1 + w2 + w3) - w5 );
		  IDX2(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) += a[1] * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 + (H_ij+c_ij*v_ij)*w3
																		+ u_ij*w4 -w_ij*w5 ); 
		}
	
		//----------------------------------------------------------------------
		// Dimensional splitting: z-direction
		//----------------------------------------------------------------------
	
		// Compute eigenvalues
		l1 = abs(w_ij-c_ij);
		l2 = abs(w_ij);
		l3 = abs(w_ij+c_ij);
		l4 = abs(w_ij);
		l5 = abs(w_ij);

		// Compute solution difference U_j-U_i
		if (btransposeSrc) {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}
		else {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			Diff[i-1] = IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
			  -IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc);
		}

		// Compute auxiliary quantities for characteristic variables
		aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*Diff[0]
											-u_ij*Diff[1]
											-v_ij*Diff[2]
											-w_ij*Diff[3]
											+Diff[4])/DCONST(2.0)/c2_ij;
		aux2 = (w_ij*Diff[0]-Diff[2])/DCONST(2.0)/c_ij;
      
		// Compute characteristic variables multiplied by the corresponding eigenvalue
		w1 = l1 * (aux1 + aux2);
		w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff[0]
				   +((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*Diff[1]
												 +v_ij*Diff[2]
												 +w_ij*Diff[3]
												 -Diff[4])/c2_ij);
		w3 = l3 * (aux1 - aux2);
		w4 = l4 * ( u_ij*Diff[0]-Diff[1]);
		w5 = l5 * (-v_ij*Diff[0]+Diff[2]);
      
		// Compute "R_ij * |Lbd_ij| * L_ij * dU"
		if (btransposeDest) {
		  IDX2T(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) += a[2] * ( w1 + w2 + w3 );
		  IDX2T(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) += a[2] * ( u_ij*(w1 + w2 + w3) - w4 );
		  IDX2T(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) += a[2] * ( v_ij*(w1 + w2 + w3) + w5 );
		  IDX2T(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) += a[2] * ( (w_ij-c_ij)*w1 + w_ij*w2 + (w_ij+c_ij)*w3 );
		  IDX2T(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) += a[2] * ( (H_ij-c_ij*w_ij)*w1 + q_ij*w2 + (H_ij+c_ij*w_ij)*w3
																		 -u_ij*w4 + v_ij*w5 );
		}
		else {
		  IDX2(VectorAtEdge,1,iposDest,NVAR3D,nedgesimDest) += a[2] * ( w1 + w2 + w3 );
		  IDX2(VectorAtEdge,2,iposDest,NVAR3D,nedgesimDest) += a[2] * ( u_ij*(w1 + w2 + w3) - w4 );
		  IDX2(VectorAtEdge,3,iposDest,NVAR3D,nedgesimDest) += a[2] * ( v_ij*(w1 + w2 + w3) + w5 );
		  IDX2(VectorAtEdge,4,iposDest,NVAR3D,nedgesimDest) += a[2] * ( (w_ij-c_ij)*w1 + w_ij*w2 + (w_ij+c_ij)*w3 );
		  IDX2(VectorAtEdge,5,iposDest,NVAR3D,nedgesimDest) += a[2] * ( (H_ij-c_ij*w_ij)*w1 + q_ij*w2 + (H_ij+c_ij*w_ij)*w3
																		-u_ij*w4 + v_ij*w5 );
		}
      } else {
		if (btransposeDest) {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) = 0.0;
		}
		else {
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) = 0.0;
		} 
      }
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing
   * scalar artificial dissipation of Rusanov-type.
   ****************************************************************************/
  
  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_RUSANOV>
  {
    template <int nedgesimDest, int nedgesimSrc,
			  bool btransposeDest, bool btransposeSrc,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui, 
							 Td uj,
							 Td vi,
							 Td vj, 
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge, 
							 Ti nedge,
							 Ti ncoeff)
    {
      Td Ei,Ej;
      if (btransposeSrc) {
		// Compute specific energies
		Ei = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc);
		Ej = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc);
      }
      else {
		// Compute specific energies
		Ei = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc);
		Ej = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc);
      }
    
      // Compute the speed of sound
      Td ci = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*
					   (HYDRO_GAMMA)*(Ei-DCONST(0.5)*(ui*ui+vi*vi+wi*wi)), DBL_EPSILON));
      Td cj = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*
					   (HYDRO_GAMMA)*(Ej-DCONST(0.5)*(uj*uj+vj*vj+wj*wj)), DBL_EPSILON));
    
#ifdef HYDRO_USE_IBP
      // Compute scalar dissipation based on the skew-symmetric part
      // which does not include the symmetric boundary contribution
      Td d_ij = max( abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge))*uj+
						 DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge))*vj-
						 DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge))*wj)+
					 DCONST(0.5)*sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
										  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge),2)+
									  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
										  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge),2)+
									  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
										  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge),2))*cj,
					 abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge))*ui+
						 DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge))*vi+
						 DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge))*wi)+
					 DCONST(0.5)*sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)-
										  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge),2)+
									  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)-
										  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge),2)+
									  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)-
										  IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge),2))*ci );
#else
      // Compute scalar dissipation
      Td d_ij = max( abs(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*uj+
						 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*vj+
						 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*wj)+
					 sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge),2)+
						  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge),2)+
						  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge),2))*cj,
					 abs(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*ui+
						 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*vi+
						 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*wi)+
					 sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge),2)+
						  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge),2)+
						  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge),2))*ci );
#endif
    
      // Multiply the solution difference by the scalar dissipation
      if (btransposeDest) {
		if (btransposeSrc) {
		  // Both source and destination vector are transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Destination vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
      else {
		if (btransposeSrc) {
		  // Source vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Both vectors are not transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing
   * scalar artificial dissipation of Rusanov-type using dimensional splitting.
   ****************************************************************************/

  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_RUSANOV_DSPLIT>
  {
    template <int nedgesimDest, int nedgesimSrc, 
			  bool btransposeDest, bool btransposeSrc,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge, 
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi, 
							 Td pj,
							 Ti iposDest,
							 Ti iposSrc,
							 Ti iedge, 
							 Ti nedge,
							 Ti ncoeff)
    {
      Td Ei,Ej;
      if (btransposeSrc) {
		// Compute specific energies
		Ei = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3T,1,iposSrc,NVAR3D,2,nedgesimSrc);
		Ej = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3T,2,iposSrc,NVAR3D,2,nedgesimSrc);
      }
      else {
		// Compute specific energies
		Ei = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3,1,iposSrc,NVAR3D,2,nedgesimSrc);
		Ej = SPECIFICTOTALENERGY3_3D(DataAtEdge,IDX3,2,iposSrc,NVAR3D,2,nedgesimSrc);
      }
    
      // Compute the speed of sound
      Td ci = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*
					   (HYDRO_GAMMA)*(Ei-DCONST(0.5)*(ui*ui+vi*vi+wi*wi)), DBL_EPSILON));
      Td cj = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*
					   (HYDRO_GAMMA)*(Ej-DCONST(0.5)*(uj*uj+vj*vj+wj*wj)), DBL_EPSILON));
    
#ifdef HYDRO_USE_IBP
      // Compute scalar dissipation with dimensional splitting based on
      // the skew-symmetric part which does not include the symmetric
      // boundary contribution
      Td d_ij = max( abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge))*uj)+
					 abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)))*cj,
					 abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge))*ui)+
					 abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)-
									  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)))*ci )
		+ max( abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge))*vj)+
			   abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)))*cj,
			   abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge))*vi)+
			   abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)))*ci )
		+ max( abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge))*wj)+
			   abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)))*cj,
			   abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge))*wi)+
			   abs(DCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)-
								IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)))*ci );
#else
      // Compute scalar dissipation with dimensional splitting
      Td d_ij = max( abs(IDX3(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*uj)+
					 abs(IDX3(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge))*cj,
					 abs(IDX3(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*ui)+
					 abs(IDX3(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge))*ci )
		+ max( abs(IDX3(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*vj)+
			   abs(IDX3(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge))*cj,
			   abs(IDX3(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*vi)+
			   abs(IDX3(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge))*ci )
		+ max( abs(IDX3(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*wj)+
			   abs(IDX3(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge))*cj,
			   abs(IDX3(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*wi)+
			   abs(IDX3(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge))*ci );
#endif
    
      // Multiply the solution difference by the scalar dissipation
	  if (btransposeDest) {
		if (btransposeSrc) {
		  // Both source and destination vector are transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Destination vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2T(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
      else {
		if (btransposeSrc) {
		  // Source vector is transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3T(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
		else {
		  // Both vectors are not transposed
#pragma unroll
		  for (int i=1; i<=NVAR3D; i++)
			IDX2(VectorAtEdge,i,iposDest,NVAR3D,nedgesimDest) =
			  d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR3D,2,nedgesimSrc)
					-IDX3(DataAtEdge,i,1,iposSrc,NVAR3D,2,nedgesimSrc));
		}
      }
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipation: Artificial dissipation
   ****************************************************************************/
  
  template <int idissipationtype>
  struct InviscidFluxDissipation : public InviscidFluxDissipationBase<idissipationtype>
  {   
    // Enable use of inherited functions
    using InviscidFluxDissipationBase<idissipationtype>::calcEdgeData;

    /***************************************************************************
     * Wrapper routine for processing a single edge
     **************************************************************************/
    template <typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *VectorAtEdge,
							 Tc *CoeffsAtEdge,
							 Td *DataAtEdge,
							 Td ui,
							 Td uj,
							 Td vi,
							 Td vj,
							 Td wi,
							 Td wj,
							 Td pi,
							 Td pj,
							 Ti iedge, 
							 Ti nedge,
							 Ti ncoeff)
    {
      InviscidFluxDissipationBase<idissipationtype>::calcEdgeData<1,1,false,false>
		(VectorAtEdge,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,1,1,iedge,nedge,ncoeff);
    }
  };

  /*****************************************************************************
   * FluxBase
   ****************************************************************************/

  struct FluxBase
  {
    /*
     * Combine inviscid fluxes (not skew-symmetric) and artificial diffusion
     */
    template <int nedgesimDest, int nedgesimSrc, bool boverwrite,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
								Tc *CoeffsAtEdge,
								Td *Fxi,
								Td *Fxj, 
								Td *Fyi, 
								Td *Fyj,
								Td *Fzi, 
								Td *Fzj,
								Td *Diff,
								Td scale,
								Ti iposDest,
								Ti iposSrc,
								Ti iedge, 
								Ti nedge,
								Ti ncoeff)
    {
      if (boverwrite) {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*IDX2(Fxj,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*IDX2(Fyj,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*IDX2(Fzj,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*IDX2(Fxi,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*IDX2(Fyi,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*IDX2(Fzi,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc));

#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,2,iposDest,NVAR3D,2,nedgesimDest) = -IDX3(FluxesAtEdge,i,1,iposDest,NVAR3D,2,nedgesimDest);
      }
      else {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*IDX2(Fxj,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*IDX2(Fyj,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*IDX2(Fzj,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*IDX2(Fxi,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*IDX2(Fyi,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*IDX2(Fzi,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc));
	
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*IDX2(Fxj,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*IDX2(Fyj,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*IDX2(Fzj,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*IDX2(Fxi,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*IDX2(Fyi,i,iposSrc,NVAR3D,nedgesimSrc)-
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*IDX2(Fzi,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc));
      }
    }
    
    /*
     * Combine inviscid fluxes (skew-symmetric) and artificial diffusion
     */
    template <int nedgesimDest, int nedgesimSrc, bool boverwrite,
			  typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
								Tc *CoeffsAtEdge,
								Td *Fx_ij,
								Td *Fy_ij,
								Td *Fz_ij,
								Td *Diff,
								Td scale,
								Ti iposDest,
								Ti iposSrc,
								Ti iedge,
								Ti nedge,
								Ti ncoeff)
    {
      if (boverwrite) {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,1,iposDest,NVAR3D,2,nedgesimDest) = scale *
			(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*IDX2(Fz_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc));
    
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,2,iposDest,NVAR3D,2,nedgesimDest) = -scale *
			(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*IDX2(Fz_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc));
      }
      else {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,1,iposDest,NVAR3D,2,nedgesimDest) += scale *
			(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,3,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,3,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,3,ncoeff,nedge)*IDX2(Fz_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc));
    
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,2,iposDest,NVAR3D,2,nedgesimDest) -= scale *
			(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,3,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,3,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,3,ncoeff,nedge)*IDX2(Fz_ij,i,iposSrc,NVAR3D,nedgesimSrc)+
			 IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc));
      }
    }
    
    /*
     * Combine inviscid fluxes with artificial diffusion
     */
    template <int nedgesimDest, int nedgesimSrc, bool boverwrite,
			  typename Td, typename Ti>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
								Td *Diff,
								Td scale,
								Ti iposDest,
								Ti iposSrc)
    {
      if (boverwrite) {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,1,iposDest,NVAR3D,2,nedgesimDest) = scale * IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc);
    
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,2,iposDest,NVAR3D,2,nedgesimDest) = -scale * IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc);
      }
      else {
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,1,iposDest,NVAR3D,2,nedgesimDest) += scale * IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc);
    
#pragma unroll
		for (int i=1; i<=NVAR3D; i++)
		  IDX3(FluxesAtEdge,i,2,iposDest,NVAR3D,2,nedgesimDest) -= scale * IDX2(Diff,i,iposSrc,NVAR3D,nedgesimSrc);
      }
    }
  };

  /*****************************************************************************
   * Flux
   ****************************************************************************/

  struct Flux : public FluxBase
  {
    // Enable use of inherited functions
    using FluxBase::combineEdgeData;

    /***************************************************************************
     * Wrapper routines for processing a single edge
     **************************************************************************/

    /*
     * Combine inviscid fluxes (not skew-symmetric) and artificial diffusion
     */
    template <bool boverwrite, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
								Tc *CoeffsAtEdge,
								Td *Fxi,
								Td *Fxj, 
								Td *Fyi, 
								Td *Fyj,
								Td *Fzi, 
								Td *Fzj,
								Td *Diff,
								Td scale,       
								Ti iedge, 
								Ti nedge,
								Ti ncoeff)
    {
      FluxBase::combineEdgeData<1,1,boverwrite>
		(FluxesAtEdge,CoeffsAtEdge,Fxi,Fxj,Fyi,Fyj,Fzi,Fzj,Diff,scale,1,1,iedge,nedge,ncoeff);
    }

    /*
     * Combine inviscid fluxes (skew-symmetric) and artificial diffusion
     */
    template <bool boverwrite, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
								Tc *CoeffsAtEdge,
								Td *Fx_ij,
								Td *Fy_ij, 
								Td *Fz_ij, 
								Td *Diff,
								Td scale,       
								Ti iedge, 
								Ti nedge,
								Ti ncoeff)
    {
      FluxBase::combineEdgeData<1,1,boverwrite>
		(FluxesAtEdge,CoeffsAtEdge,Fx_ij,Fy_ij,Fz_ij,Diff,scale,1,1,iedge,nedge,ncoeff);
    }

    /*
     * Combine inviscid fluxes with artificial diffusion
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
								Td *Diff,
								Td scale)
    {
      FluxBase::combineEdgeData<1,1,boverwrite>
		(FluxesAtEdge,Diff,scale,1,1);
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (baseline implementation).
   ****************************************************************************/
  
  template <typename Tc,
			typename TdSrc,
			typename TdDest,
			typename Ti,
			int isystemformat,
			int idissipationtype,
			int threads_per_cta>
  __launch_bounds__(threads_per_cta)
  __global__ void hydro_calcFlux3d_baseline(Tc *CoeffsAtEdge,
											Ti *IedgeList,
											TdSrc *vecSrc,
											TdDest *vecDest,
											TdDest scale,
											Ti neq,
											Ti nedge,
											Ti ncoeff,
											Ti nedge_last,
											Ti nedge_per_thread=1,
											Ti nedge_offset=0)
  {
    // Loop over all items per thread
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {
      
      // Global edge ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*threads_per_cta+nedge_offset+threadIdx.x;
      
      if (threadIdx.x<threads_per_cta && idx<nedge_last)
		{
		  // Get positions of edge endpoints (idx starts at zero)
		  Ti i = IDX2_EDGELIST(IedgeList,1,idx+1,6,nedge);
		  Ti j = IDX2_EDGELIST(IedgeList,2,idx+1,6,nedge);
	  
		  // Local variables
		  TdDest DataAtEdge[2*NVAR3D];
	  
		  // Get solution values at edge endpoints
		  Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
			gatherEdgeData<true>(DataAtEdge,vecSrc,i,j,neq);
	  
		  // Compute velocities
		  TdDest ui = XVELOCITY2_3D(DataAtEdge,IDX2,1,NVAR3D,2);
		  TdDest vi = YVELOCITY2_3D(DataAtEdge,IDX2,1,NVAR3D,2);
		  TdDest wi = ZVELOCITY2_3D(DataAtEdge,IDX2,1,NVAR3D,2);
	  
		  TdDest uj = XVELOCITY2_3D(DataAtEdge,IDX2,2,NVAR3D,2);
		  TdDest vj = YVELOCITY2_3D(DataAtEdge,IDX2,2,NVAR3D,2);
		  TdDest wj = ZVELOCITY2_3D(DataAtEdge,IDX2,2,NVAR3D,2);
	  
		  // Compute pressures
		  TdDest pi = PRESSURE2_3D(DataAtEdge,IDX2,1,NVAR3D,2);
		  TdDest pj = PRESSURE2_3D(DataAtEdge,IDX2,2,NVAR3D,2);
	  
		  // Local variables
		  TdDest FluxAtEdge[2*NVAR3D];

		  // Compute the artificial viscosities
		  InviscidFluxDissipation<idissipationtype>::calcEdgeData
			(&FluxAtEdge[NVAR3D],CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,idx+1,nedge,ncoeff);
		  Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR3D],scale);

		  // Compute inviscid fluxes
		  InviscidFlux::calcEdgeData<false>
			(FluxAtEdge,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,scale,idx+1,nedge,ncoeff);

		  // Build fluxes into nodal vector
		  Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
			scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);
		}
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (shared memory implementation)
   ****************************************************************************/
  
  template <typename Tc,
			typename TdSrc,
			typename TdDest,
			typename Ti,
			int isystemformat,
			int idissipationtype,
			int threads_per_cta>
  __launch_bounds__(threads_per_cta)
  __global__ void hydro_calcFlux3d_shmem(Tc *CoeffsAtEdge,
										 Ti *IedgeList,
										 TdSrc *vecSrc,
										 TdDest *vecDest,
										 TdDest scale,
										 Ti neq,
										 Ti nedge,
										 Ti ncoeff,
										 Ti nedge_last,
										 Ti nedge_per_thread=1,
										 Ti nedge_offset=0)
  {
    // Shared memory
    __shared__ TdSrc s_DataAtEdge[2*NVAR3D*threads_per_cta];
    
    // Loop over all items per thread
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {
      
      // Global edge ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*threads_per_cta+nedge_offset+threadIdx.x;
      
      if (threadIdx.x<threads_per_cta && idx<nedge_last)
		{
		  // Get positions of edge endpoints (idx starts at zero)
		  Ti i = IDX2_EDGELIST(IedgeList,1,idx+1,6,nedge);
		  Ti j = IDX2_EDGELIST(IedgeList,2,idx+1,6,nedge);
	  
		  // Get solution values at edge endpoints
		  Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
			gatherEdgeData<threads_per_cta,SHMEM_DATA_TRANSPOSE,true>
			(s_DataAtEdge,vecSrc,(int)threadIdx.x+1,i,j,neq);
	  
		  // Compute velocities
		  TdDest ui = XVELOCITY3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,1,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
		  TdDest vi = YVELOCITY3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,1,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
		  TdDest wi = ZVELOCITY3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,1,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
	  
		  TdDest uj = XVELOCITY3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,2,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
		  TdDest vj = YVELOCITY3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,2,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
		  TdDest wj = ZVELOCITY3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,2,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
	  
		  // Compute pressures
		  TdDest pi = PRESSURE3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,1,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
		  TdDest pj = PRESSURE3_3D(s_DataAtEdge,SHMEM_DATA_IDX3,2,(int)threadIdx.x+1,NVAR3D,2,threads_per_cta);
	  
		  // Local variables
		  TdDest FluxAtEdge[2*NVAR3D];
	  
		  // Compute the artificial viscosities
		  InviscidFluxDissipation<idissipationtype>::
			calcEdgeData<1,threads_per_cta,false,SHMEM_DATA_TRANSPOSE>
			(&FluxAtEdge[NVAR3D],CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,
			 1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
		  Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR3D],scale);
	  
		  // Compute inviscid fluxes
		  InviscidFlux::calcEdgeData<1,threads_per_cta,false,SHMEM_DATA_TRANSPOSE,false>
			(FluxAtEdge,CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,
			 scale,1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	  
		  // Build fluxes into nodal vector
		  Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
			scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);
		}
    }
  };

  /*****************************************************************************
   * Internal C++ functions which invoke the CUDA kernels
   ****************************************************************************/

  template <typename Tc,
			typename TdSrc,
			typename TdDest,
			typename Ti,
			int idissipationtype>
  inline
  int hydro_calcFlux3d_cuda(__SIZET *d_CoeffsAtEdge,
							__SIZET *d_IedgeList,
							__SIZET *d_vecSrc,
							__SIZET *d_vecDest,
							TdDest scale,
							Ti nblocks,
							Ti neq,
							Ti nedge, 
							Ti ncoeff,
							Ti nedgeset,
							Ti iedgeset,
							cudaStream_t stream=0)
  {
    const cudaDeviceProp *devProp = coproc_getCurrentDeviceProp();
    
    // Strategy: run the largest possible number of blocks with a
    // predefined number of compute/dma threads per block and let each
    // compute thread process the minimal number of edges
    const int compute_threads_per_cta  = CUDADMA_COMPUTE_THREADS_PER_CTA;
    const int dma_threads_per_ld       = CUDADMA_THREADS_PER_LD;
    const int dma_lds                  = CUDADMA_DMA_LDS;
    int nedge_per_thread_cudaDMA       = CUDADMA_NEDGE_PER_THREAD;

    const int threads_per_cta_baseline = BASELINE_THREADS_PER_CTA;
    int nedge_per_thread_baseline      = BASELINE_NEDGE_PER_THREAD;
    
    int blocks, threads, nedge_cudaDMA, nedge_baseline;
    prepare_cudaDMA(devProp, nedgeset,
					&nedge_per_thread_cudaDMA,
					compute_threads_per_cta, dma_threads_per_ld,
					dma_lds, &blocks, &threads, &nedge_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(devProp, nedgeset-nedge_cudaDMA,
					 &nedge_per_thread_baseline, threads_per_cta_baseline,
					 &blocks, &threads, &nedge_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);

    TdSrc  *vecSrc = (TdSrc*)(*d_vecSrc);
    TdDest *vecDest = (TdDest*)(*d_vecDest);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
    
    if (nblocks == 1) {
#ifdef CUDADMA_KERNEL    
      if (grid_cudaDMA.x>0) {
      	// CudaDMA implementation
		CUDADMA_KERNEL
      	  <Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,idissipationtype,
          MAX(32,compute_threads_per_cta),MAX(32,dma_threads_per_ld)>
      	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
													   IedgeList,
													   vecSrc, vecDest, scale,
													   neq, nedge, ncoeff,
													   nedge_cudaDMA+iedgeset-1, 
													   nedge_per_thread_cudaDMA,
													   iedgeset-1);
      }
#endif

#ifdef BASELINE_KERNEL
      if (grid_baseline.x>0) {
		// Baseline implementation
		BASELINE_KERNEL
		  <Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,idissipationtype,
          threads_per_cta_baseline>
		  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
														 IedgeList,
														 vecSrc, vecDest, scale,
														 neq, nedge, ncoeff,
														 nedgeset+iedgeset-1, 
														 nedge_per_thread_baseline,
														 nedge_cudaDMA+iedgeset-1);
      }
#endif
    } else {
#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0) {
      	// CudaDMA implementation
      	CUDADMA_KERNEL
      	  <Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,idissipationtype,
          MAX(32,compute_threads_per_cta),MAX(32,dma_threads_per_ld)>
      	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
													   IedgeList,
													   vecSrc, vecDest, scale, 
													   neq, nedge, ncoeff,
													   nedge_cudaDMA+iedgeset-1, 
													   nedge_per_thread_cudaDMA,
													   iedgeset-1);
      }
#endif

#ifdef BASELINE_KERNEL
      if (grid_baseline.x>0) {
      	// Baseline implementation
		BASELINE_KERNEL
		  <Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,idissipationtype,
          threads_per_cta_baseline>
		  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
														 IedgeList,
														 vecSrc, vecDest, scale, 
														 neq, nedge, ncoeff,
														 nedgeset+iedgeset-1, 
														 nedge_per_thread_baseline,
														 nedge_cudaDMA+iedgeset-1);
      }
#endif
    }

    coproc_checkError("hydro_calcFlux3d_cuda");
    return 0;
  }; 
  
  /*****************************************************************************
   * External C functions which can be called from the Fortran code
   ****************************************************************************/

  extern "C"
  {
    __INT FNAME(hydro_calcfluxgalerkin3d_cuda)(__SIZET *d_CoeffsAtEdge,
											   __SIZET *d_IedgeList,
											   __SIZET *d_vecSrc,
											   __SIZET *d_vecDest,
											   __DP *scale,
											   __INT *nblocks,
											   __INT *neq,
											   __INT *nedge,
											   __INT *ncoeff,
											   __INT *nedges,
											   __INT *iedgeset,
											   __I64 *stream)
    {
      return (__INT) hydro_calcFlux3d_cuda
		<__DP,__DP,__DP,__INT,DISSIPATION_ZERO>
		(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
		 *scale, *nblocks, *neq, *nedge,
		 *ncoeff, *nedges, *iedgeset, 
		 (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcfluxscdiss3d_cuda)(__SIZET *d_CoeffsAtEdge,
											 __SIZET *d_IedgeList,
											 __SIZET *d_vecSrc,
											 __SIZET *d_vecDest,
											 __DP *scale,
											 __INT *nblocks,
											 __INT *neq,
											 __INT *nedge,
											 __INT *ncoeff,
											 __INT *nedges,
											 __INT *iedgeset,
											 __I64 *stream)
    {
      return (__INT) hydro_calcFlux3d_cuda
		<__DP,__DP,__DP,__INT,DISSIPATION_SCALAR>
		(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
		 *scale, *nblocks, *neq, *nedge,
		 *ncoeff, *nedges, *iedgeset, 
		 (cudaStream_t)(*stream));
    }

    /**************************************************************************/

    __INT FNAME(hydro_calcfluxscdissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
												 __SIZET *d_IedgeList,
												 __SIZET *d_vecSrc,
												 __SIZET *d_vecDest,
												 __DP *scale,
												 __INT *nblocks, 
												 __INT *neq, 
												 __INT *nedge, 
												 __INT *ncoeff,
												 __INT *nedges, 
												 __INT *iedgeset,
												 __I64 *stream)
    {
      return (__INT) hydro_calcFlux3d_cuda
		<__DP,__DP,__DP,__INT,DISSIPATION_SCALAR_DSPLIT>
		(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
		 *scale, *nblocks, *neq, *nedge,
		 *ncoeff, *nedges, *iedgeset,
		 (cudaStream_t)(*stream));
    }

    /**************************************************************************/

    __INT FNAME(hydro_calcfluxroediss3d_cuda)(__SIZET *d_CoeffsAtEdge,
											  __SIZET *d_IedgeList,
											  __SIZET *d_vecSrc,
											  __SIZET *d_vecDest,
											  __DP *scale,
											  __INT *nblocks, 
											  __INT *neq, 
											  __INT *nedge, 
											  __INT *ncoeff,
											  __INT *nedges, 
											  __INT *iedgeset,
											  __I64 *stream)
    {
      return (__INT) hydro_calcFlux3d_cuda
		<__DP,__DP,__DP,__INT,DISSIPATION_ROE>
		(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
		 *scale, *nblocks, *neq, *nedge,
		 *ncoeff, *nedges, *iedgeset,
		 (cudaStream_t)(*stream));
    }

	/***************************************************************************/

    __INT FNAME(hydro_calcfluxroedissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
												  __SIZET *d_IedgeList,
												  __SIZET *d_vecSrc,
												  __SIZET *d_vecDest,
												  __DP *scale,
												  __INT *nblocks, 
												  __INT *neq, 
												  __INT *nedge, 
												  __INT *ncoeff,
												  __INT *nedges, 
												  __INT *iedgeset,
												  __I64 *stream)
    {
      return (__INT) hydro_calcFlux3d_cuda
		<__DP,__DP,__DP,__INT,DISSIPATION_ROE_DSPLIT>
		(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
		 *scale, *nblocks, *neq, *nedge,
		 *ncoeff, *nedges, *iedgeset,
		 (cudaStream_t)*stream);
    }

    /**************************************************************************/

    __INT FNAME(hydro_calcfluxrusdiss3d_cuda)(__SIZET *d_CoeffsAtEdge,
											  __SIZET *d_IedgeList,
											  __SIZET *d_vecSrc,
											  __SIZET *d_vecDest,
											  __DP *scale,
											  __INT *nblocks, 
											  __INT *neq, 
											  __INT *nedge, 
											  __INT *ncoeff,
											  __INT *nedges, 
											  __INT *iedgeset,
											  __I64 *stream)
    {
      return (__INT)hydro_calcFlux3d_cuda
		<__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV>
		(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
		 *scale, *nblocks, *neq, *nedge,
		 *ncoeff, *nedges, *iedgeset,
		 (cudaStream_t)*stream);
    }
    
    /**************************************************************************/

    __INT FNAME(hydro_calcfluxrusdissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
												  __SIZET *d_IedgeList,
												  __SIZET *d_vecSrc,
												  __SIZET *d_vecDest,
												  __DP *scale,
												  __INT *nblocks, 
												  __INT *neq, 
												  __INT *nedge, 
												  __INT *ncoeff,
												  __INT *nedges, 
												  __INT *iedgeset,
												  __I64 *stream)
    {
      return (__INT) hydro_calcFlux3d_cuda
		<__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV_DSPLIT>
		(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
		 *scale, *nblocks, *neq, *nedge,
		 *ncoeff, *nedges, *iedgeset,
		 (cudaStream_t)*stream);
    }
  };
}
