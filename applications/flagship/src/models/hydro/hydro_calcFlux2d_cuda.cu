/*#############################################################################
 **************************************<****************************************
 * <name> hydro_calcFlux2d_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This file provides CUDA kernels to compute the fluxes for the low-order
 * scheme in 2D using different types if artificial viscosities.
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
#include "../../cudaDMA.h"
#include "../../cudaGatherScatter.h"

#define LANGUAGE LANGUAGE_C
#include "../../flagship.h"
#include "../../cudaMacros.h"

#define HYDRO_NDIM 2
#include "hydro.h"

// Define CUDA kernel which does not make use of the CUDADMA library
// and is applied to the remaining edges which are not processed in groups
// #define BASELINE_KERNEL hydro_calcFlux2d_shmem
#define BASELINE_KERNEL hydro_calcFlux2d_baseline

// Defines for baseline implementation
#define BASELINE_THREADS_PER_CTA  32*2
#define BASELINE_NEDGE_PER_THREAD 1

// Defines for shared memory implementation
#define SHMEM_DATA_TRANSPOSE   true
#define SHMEM_DATA_IDX3        IDX3T
#define SHMEM_NEDGE_PER_THREAD BASELINE_NEDGE_PER_THREAD

// Define CUDA kernel which makes use of the CUDADMA library to achive
// higher throughput between global and shared memory on the device
// #define CUDADMA_PREFETCH_SINGLE

// Defines for cudaDMA implementation without warp specialisation
#ifdef CUDADMA_NOSPEC
#define CUDADMA_KERNEL                  hydro_calcFlux2d_cudaDMA_nospec
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
#define CUDADMA_KERNEL                  hydro_calcFlux2d_cudaDMA_prefetch_single
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*4
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        1*1
#define CUDADMA_DMA_LDS_IND             0
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           1
#define CUDADMA_DMA_LDS                 (CUDADMA_DMA_LDS_IND + \
                                         CUDADMA_DMA_LDS_SRC + \
                                         CUDADMA_DMA_LDS_DEST +\
					 CUDADMA_DMA_LDS_COEFF)
#endif

// Defines for cudaDMA double buffer implementation with prefetching of indices
#ifdef CUDADMA_PREFETCH_DOUBLE
#define CUDADMA_KERNEL                  hydro_calcFlux2d_cudaDMA_prefetch_double
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*4
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        4*1
#define CUDADMA_DMA_LDS_IND             0
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           1
#define CUDADMA_DMA_LDS                 (CUDADMA_DMA_LDS_IND + \
                                         CUDADMA_DMA_LDS_SRC + \
                                         CUDADMA_DMA_LDS_DEST +\
					 CUDADMA_DMA_LDS_COEFF)
#endif

// Defines for cudaDMA double buffer implementation
#ifdef CUDADMA_DOUBLE
#define CUDADMA_KERNEL                  hydro_calcFlux2d_cudaDMA_double
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*2
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        6*1
#define CUDADMA_DMA_LDS_IND             1
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           0
#define CUDADMA_DMA_LDS                 (3*CUDADMA_DMA_LDS_IND + \
                                         2*CUDADMA_DMA_LDS_SRC + \
                                         2*CUDADMA_DMA_LDS_DEST)
#endif

// Defines for cudaDMA manual buffer implementation
#ifdef CUDADMA_MANUAL
#define CUDADMA_KERNEL                  hydro_calcFlux2d_cudaDMA_manual
#define CUDADMA_COMPUTE_THREADS_PER_CTA 32*2
#define CUDADMA_THREADS_PER_LD          32*1
#define CUDADMA_NEDGE_PER_THREAD        6*1
#define CUDADMA_DMA_LDS_IND             1
#define CUDADMA_DMA_LDS_SRC             1
#define CUDADMA_DMA_LDS_DEST            1
#define CUDADMA_DMA_LDS_COEFF           1
#define CUDADMA_DMA_LDS                 (CUDADMA_DMA_LDS_IND + \
                                         CUDADMA_DMA_LDS_SRC + \
                                         CUDADMA_DMA_LDS_DEST +\
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

namespace hydro2d_cuda
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
	    IDX2T(Fxi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2T(Fxj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fxi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2T(Fxj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fxi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2(Fxj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fxi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2(Fxj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	}
      }
      else {
	// Keep content of destination vector
	if (btransposeDest) {
	  if (btransposeSrc) {
	    // Both source and destination vector are transposed
	    IDX2T(Fxi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2T(Fxj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fxi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2T(Fxi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2T(Fxj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fxj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fxi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2(Fxj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fxi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    IDX2(Fxi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi);
	    
	    IDX2(Fxj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fxj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
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
	    IDX2T(Fyi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2T(Fyj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fyi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2T(Fyj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fyi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2(Fyj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fyi,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2(Fyj,1,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,2,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,3,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,4,iposDest,NVAR2D,nedgesimDest) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	}
      }
      else {
	// Keep content of destination vector
	if (btransposeDest) {
	  if (btransposeSrc) {
	    // Both source and destination vector are transposed
	    IDX2T(Fyi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2T(Fyj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fyi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2T(Fyi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2T(Fyj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fyj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fyi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2(Fyj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fyi,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	    IDX2(Fyi,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi);
	
	    IDX2(Fyj,1,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,2,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,3,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fyj,4,iposDest,NVAR2D,nedgesimDest) += INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
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
	    IDX2T(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3t,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	}
      }
      else {
	// Keep content of destination vector
	if (btransposeDest) {
	  if (btransposeSrc) {
	    // Both source and destination vector are transposed
	    IDX2T(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3t,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2T(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fx_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
	    IDX2(Fx_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
	      INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj);
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
	    IDX2T(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) =
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	}
      }
      else {
	// Keep content of destination vector
	if (btransposeDest) {
	  if (btransposeSrc) {
	    // Both source and destination vector are transposed
	    IDX2T(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Destination vector is transposed
	    IDX2T(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2T(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX2(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	  }
	  else {
	    // Both vectors are not transposed
	    IDX2(Fy_ij,1,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,2,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,3,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
	    IDX2(Fy_ij,4,iposDest,NVAR2D,nedgesimDest) +=
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
	      INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj);
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
    static void calcEdgeData(Td *Fxi,
			     Td *Fxj,
			     Td *Fyi,
			     Td *Fyj,
			     Td *DataAtEdge,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td pi,
			     Td pj,
			     Ti iposDest,
			     Ti iposSrc)
    {
      // Compute the Galerkin fluxes for x-direction
      InviscidFluxBase::
	calcFluxXdir<nedgesimDest,nedgesimSrc,btransposeDest,btransposeSrc,boverwrite>
	(Fxi,Fxj,DataAtEdge,ui,uj,pi,pj,iposDest,iposSrc);
      
      // Compute Galerkin fluxes for y-direction
      InviscidFluxBase::
	calcFluxYdir<nedgesimDest,nedgesimSrc,btransposeDest,btransposeSrc,boverwrite>
	(Fyi,Fyj,DataAtEdge,vi,vj,pi,pj,iposDest,iposSrc);
    }
    
    /*
     * Calculate the inviscid fluxes in all directions (skew-symmetric)
     */
    template <int nedgesimDest, int nedgesimSrc, 
	      bool btransposeDest, bool btransposeSrc,
	      bool boverwrite, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *Fx_ij,
			     Td *Fy_ij,
			     Td *DataAtEdge,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
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
    }

    /*
     * Calculate the inviscid fluxes in all directions (not skew-symmetric)
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
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }
      else {
	aux = scale *
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }

      if (boverwrite) {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3T(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
	else {
	  IDX3(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
      }
      else {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3T(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
	else {
	  IDX3(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
      }

      // Flux component 2
      if (btransposeSrc) {
	aux = scale *
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }
      else {
	aux = scale *
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }
      
      if (boverwrite) {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3T(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
	else {
	  IDX3(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
      }
      else {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3T(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
	else {
	  IDX3(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
      }

      // Flux component 3
      if (btransposeSrc) {
	aux = scale *
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }
      else {
	aux = scale *
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }
      
      if (boverwrite) {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3T(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
	else {
	  IDX3(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
      }
      else {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3T(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
	else {
	  IDX3(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
      }

      // Flux component 4
      if (btransposeSrc) {
	aux = scale *
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }
      else {
	aux = scale *
	  (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj)
	  +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)
	  -IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	   INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi));
      }

      if (boverwrite) {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3T(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
	else {
	  IDX3(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) =  aux;
	  IDX3(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) = -aux;
	}
      }
      else {
	if (btransposeDest) {
	  IDX3T(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3T(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
	else {
	  IDX3(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) += aux;
	  IDX3(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) -= aux;
	}
      }

#else
      // Calculate inviscid fluxes (not skew-symmetric)

      if (boverwrite) {
	// Overwrite destination vector
	if (btransposeDest) {
	  if (btransposeSrc) {
	    // Both source and destination vector are transposed
	    IDX3T(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	  }
	  else {
	    // Destination vector is transposed
	    IDX3T(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX3(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	  }
	  else {
	    // Both vectors are not transposed
	    IDX3(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	  }
	}
      }
      else {
	// Keep content of destination vector
	if (btransposeDest) {
	  if (btransposeSrc) {
	    // Both source and destination vector are transposed
	    IDX3T(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	  }
	  else {
	    // Destination vector is transposed
	    IDX3T(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3T(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3T(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	  }
	}
	else {
	  if (btransposeSrc) {
	    // Source vector is transposed
	    IDX3(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	  }
	  else {
	    // Both vectors are not transposed
	    IDX3(FluxesAtEdge,1,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,1,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,2,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,2,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,3,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,3,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    
	    IDX3(FluxesAtEdge,4,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
	    
	    IDX3(FluxesAtEdge,4,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	      (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,ui,pi)-
		INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,uj,pj))
	       +IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*
	       (INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc,vi,pi)-
		INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc,vj,pj)));
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
     * Calculate the inviscid fluxes in all directions (not skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcEdgeData(Td *Fxi,
			     Td *Fxj,
			     Td *Fyi,
			     Td *Fyj,
			     Td *DataAtEdge,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td pi,
			     Td pj)
    {
      InviscidFluxBase::calcEdgeData<1,1,false,false,boverwrite>
	(Fxi,Fxj,Fyi,Fyj,DataAtEdge,ui,uj,vi,vj,pi,pj,1,1);
    }

    /*
     * Calculate the inviscid fluxes in all directions (skew-symmetric)
     */
    template <bool boverwrite, typename Td>
    __device__ __forceinline__
    static void calcEdgeData(Td *Fx_ij,
			     Td *Fy_ij,
			     Td *DataAtEdge,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td pi,
			     Td pj)
    {
      InviscidFluxBase::calcEdgeData<1,1,false,false,boverwrite>
	(Fx_ij,Fy_ij,DataAtEdge,ui,uj,vi,vj,pi,pj,1,1);
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
			     Td pi,
			     Td pj,
			     Td scale,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      InviscidFluxBase::calcEdgeData<1,1,false,false,boverwrite>
	(FluxesAtEdge,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,pi,pj,
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
	for (int i=1; i<=NVAR2D; i++)
	  IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) = 0.0;
      }
      else {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) = 0.0;
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
	ri = DENSITY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc);
	rj = DENSITY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc);
	
	// Compute enthalpies
	hi = (TOTALENERGY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	hj = (TOTALENERGY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
      }
      else {
	// Compute densities
	ri = DENSITY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc);
	rj = DENSITY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc);
	
	// Compute enthalpies
	hi = (TOTALENERGY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	hj = (TOTALENERGY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
      }
    
      // Compute Roe mean values
      Td aux  = ROE_MEAN_RATIO(ri,rj);
      Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
    
      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));

      // Compute auxiliary variables
      Td q_ij   = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij);
      Td vel_ij = u_ij * a[0] + v_ij * a[1];
    
      // Compute the speed of sound
      Td c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON));
    
      // Compute scalar dissipation
      Td d_ij = abs(vel_ij) + sqrt(a[0] * a[0] + a[1] * a[1])*c_ij;
    
      // Multiply the solution difference by the scalar dissipation
      if (btransposeDest) {
	if (btransposeSrc) {
	  // Both source and destination vector are transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Destination vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
      }
      else {
	if (btransposeSrc) {
	  // Source vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Both vectors are not transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
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
	ri = DENSITY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc);
	rj = DENSITY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc);
	
	// Compute enthalpies
	hi = (TOTALENERGY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	hj = (TOTALENERGY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
      }
      else {
	// Compute densities
	ri = DENSITY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc);
	rj = DENSITY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc);
	
	// Compute enthalpies
	hi = (TOTALENERGY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	hj = (TOTALENERGY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
      }
    
      // Compute Roe mean values
      Td aux  = ROE_MEAN_RATIO(ri,rj);
      Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
    
      // Compute auxiliary variables
      Td q_ij = RCONST(0.5) *(u_ij * u_ij + v_ij * v_ij);
    
      // Compute the speed of sound
      Td c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON));
    
      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));

      // Compute scalar dissipation
      Td d_ij = ( abs(a[0]*u_ij) + abs(a[0])*c_ij +
		  abs(a[1]*v_ij) + abs(a[1])*c_ij );
    
      // Multiply the solution difference by the scalar dissipation
      if (btransposeDest) {
	if (btransposeSrc) {
	  // Both source and destination vector are transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Destination vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
      }
      else {
	if (btransposeSrc) {
	  // Source vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Both vectors are not transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
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
			     Td pi,
			     Td pj,
			     Ti iposDest,
			     Ti iposSrc,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      Td anorm = sqrt(a[0] * a[0] + a[1] * a[1]);
    
      if (anorm > DBL_EPSILON) {
      
	// Normalise the skew-symmetric coefficient
	a[0] = a[0]/anorm;
	a[1] = a[1]/anorm;
      
	Td ri,rj,hi,hj;
	if (btransposeSrc) {
	  // Compute densities
	  ri = DENSITY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc);
	  rj = DENSITY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc);
	  
	  // Compute enthalpies
	  hi = (TOTALENERGY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	  hj = (TOTALENERGY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
	}
	else {
	  // Compute densities
	  ri = DENSITY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc);
	  rj = DENSITY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc);
	  
	  // Compute enthalpies
	  hi = (TOTALENERGY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	  hj = (TOTALENERGY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
	}
      
	// Compute Roe mean values
	Td aux  = ROE_MEAN_RATIO(ri,rj);
	Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
	Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
	Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
	// Compute auxiliary variables
	Td vel_ij = u_ij * a[0] + v_ij * a[1];
	Td q_ij   = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij);
      
	// Compute the speed of sound
	Td c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON);
	Td c_ij  = sqrt(c2_ij);
      
	// Compute eigenvalues
	Td l1 = abs(vel_ij-c_ij);
	Td l2 = abs(vel_ij);
	Td l3 = abs(vel_ij+c_ij);
	Td l4 = abs(vel_ij);
      
	// Compute solution difference U_j-U_i
	Td Diff[NVAR2D];
	if (btransposeSrc) {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    Diff[i-1] = IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
	               -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc);
	}
	else {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    Diff[i-1] = IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
	               -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc);
	}
	      
	// Compute auxiliary quantities for characteristic variables
	Td aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*Diff[0]
					      -u_ij*Diff[1]
					      -v_ij*Diff[2]
					           +Diff[3])/RCONST(2.0)/c2_ij;
	Td aux2 = (vel_ij*Diff[0]
		    -a[0]*Diff[1]
		    -a[1]*Diff[2])/RCONST(2.0)/c_ij;
      
	// Compute characteristic variables multiplied by the corresponding eigenvalue
	Td w1 = l1 * (aux1 + aux2);
	Td w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff[0]
		      +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff[1]
						   +v_ij*Diff[2]
						        -Diff[3])/c2_ij);
	Td w3 = l3 * (aux1 - aux2);
	Td w4 = l4 * ((a[0]*v_ij-a[1]*u_ij)*Diff[0]
		                      +a[1]*Diff[1]
		                      -a[0]*Diff[2]);

	// Compute "R_ij * |Lbd_ij| * L_ij * dU"
	if (btransposeDest) {
	  IDX2T(VectorAtEdge,1,iposDest,NVAR2D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
	  IDX2T(VectorAtEdge,2,iposDest,NVAR2D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
									 (u_ij+c_ij*a[0])*w3 + a[1]*w4 );
	  IDX2T(VectorAtEdge,3,iposDest,NVAR2D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
									 (v_ij+c_ij*a[1])*w3 - a[0]*w4 );
	  IDX2T(VectorAtEdge,4,iposDest,NVAR2D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 +
									 (H_ij+c_ij*vel_ij)*w3 + (u_ij*a[1]-v_ij*a[0])*w4 );
	}
	else {
	  IDX2(VectorAtEdge,1,iposDest,NVAR2D,nedgesimDest) = anorm * ( w1 + w2 + w3 );
	  IDX2(VectorAtEdge,2,iposDest,NVAR2D,nedgesimDest) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
									(u_ij+c_ij*a[0])*w3 + a[1]*w4 );
	  IDX2(VectorAtEdge,3,iposDest,NVAR2D,nedgesimDest) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
									(v_ij+c_ij*a[1])*w3 - a[0]*w4 );
	  IDX2(VectorAtEdge,4,iposDest,NVAR2D,nedgesimDest) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 +
									(H_ij+c_ij*vel_ij)*w3 + (u_ij*a[1]-v_ij*a[0])*w4 );
	}
      } else {
	if (btransposeDest) {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) = 0.0;
	}
	else {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) = 0.0;
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
			     Td pi,
			     Td pj,
			     Ti iposDest,
			     Ti iposSrc,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      Td anorm = sqrt(a[0] * a[0] + a[1] * a[1]);
    
      if (anorm > DBL_EPSILON) {
      
	// Compute the absolute value
	a[0] = abs(a[0]);
	a[1] = abs(a[1]);
      
	Td ri,rj,hi,hj;
	if (btransposeSrc) {
	  // Compute densities
	  ri = DENSITY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc);
	  rj = DENSITY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc);
	  
	  // Compute enthalpies
	  hi = (TOTALENERGY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	  hj = (TOTALENERGY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
	}
	else {
	  // Compute densities
	  ri = DENSITY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc);
	  rj = DENSITY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc);
	  
	  // Compute enthalpies
	  hi = (TOTALENERGY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc)+pi)/ri;
	  hj = (TOTALENERGY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc)+pj)/rj;
	}
      
	// Compute Roe mean values
	Td aux  = ROE_MEAN_RATIO(ri,rj);
	Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
	Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
	Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
	// Compute auxiliary variable
	Td q_ij   = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij);
      
	// Compute the speed of sound
	Td c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON);
	Td c_ij  = sqrt(c2_ij);
      
	//----------------------------------------------------------------------
	// Dimensional splitting: x-direction
	//----------------------------------------------------------------------
      
	// Compute eigenvalues
	Td l1 = abs(u_ij-c_ij);
	Td l2 = abs(u_ij);
	Td l3 = abs(u_ij+c_ij);
	Td l4 = abs(u_ij);
      
	// Compute solution difference U_j-U_i
	Td Diff[NVAR2D];
	if (btransposeSrc) {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    Diff[i-1] = IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
	               -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc);
	}
	else {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    Diff[i-1] = IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
	               -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc);
	}
      
	// Compute auxiliary quantities for characteristic variables
	Td aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*Diff[0]
					      -u_ij*Diff[1]
					      -v_ij*Diff[2]
					           +Diff[3])/RCONST(2.0)/c2_ij;
	Td aux2 = (u_ij*Diff[0]
		   -Diff[1])/RCONST(2.0)/c_ij;
      
	// Compute characteristic variables multiplied by the corresponding eigenvalue
	Td w1 = l1 * (aux1 + aux2);
	Td w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff[0]
		      +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff[1]
						   +v_ij*Diff[2]
						        -Diff[3])/c2_ij);
	Td w3 = l3 * (aux1 - aux2);
	Td w4 = l4 * (v_ij*Diff[0]
		      -Diff[2]);
        
	// Compute "R_ij * |Lbd_ij| * L_ij * dU"
	if (btransposeDest) {
	  IDX2T(VectorAtEdge,1,iposDest,NVAR2D,nedgesimDest) = a[0] * ( w1 + w2 + w3 );
	  IDX2T(VectorAtEdge,2,iposDest,NVAR2D,nedgesimDest) = a[0] * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 );
	  IDX2T(VectorAtEdge,3,iposDest,NVAR2D,nedgesimDest) = a[0] * (        v_ij*w1 + v_ij*w2 +        v_ij*w3 - w4 );
	  IDX2T(VectorAtEdge,4,iposDest,NVAR2D,nedgesimDest) = a[0] * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 +
									(H_ij+c_ij*u_ij)*w3 - v_ij*w4 );
	}
	else {
	  IDX2(VectorAtEdge,1,iposDest,NVAR2D,nedgesimDest) = a[0] * ( w1 + w2 + w3 );
	  IDX2(VectorAtEdge,2,iposDest,NVAR2D,nedgesimDest) = a[0] * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 );
	  IDX2(VectorAtEdge,3,iposDest,NVAR2D,nedgesimDest) = a[0] * (        v_ij*w1 + v_ij*w2 +        v_ij*w3 - w4 );
	  IDX2(VectorAtEdge,4,iposDest,NVAR2D,nedgesimDest) = a[0] * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 +
								       (H_ij+c_ij*u_ij)*w3 - v_ij*w4 );
	}
      
	//----------------------------------------------------------------------
	// Dimensional splitting: y-direction
	//----------------------------------------------------------------------
      
	// Compute eigenvalues
	l1 = abs(v_ij-c_ij);
	l2 = abs(v_ij);
	l3 = abs(v_ij+c_ij);
	l4 = abs(v_ij);
      
	// Compute solution difference U_j-U_i
	if (btransposeSrc) {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    Diff[i-1] = IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
	               -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc);
	}
	else {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    Diff[i-1] = IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
	               -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc);
	}
       
	// Compute auxiliary quantities for characteristic variables
	aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*Diff[0]
					   -u_ij*Diff[1]
					   -v_ij*Diff[2]
					        +Diff[3])/RCONST(2.0)/c2_ij;
	aux2 = (v_ij*Diff[0]
		-Diff[2])/RCONST(2.0)/c_ij;
      
	// Compute characteristic variables multiplied by the corresponding eigenvalue
	w1 = l1 * (aux1 + aux2);
	w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff[0]
		   +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff[1]
						+v_ij*Diff[2]
						     -Diff[3])/c2_ij);
	w3 = l3 * (aux1 - aux2);
	w4 = l4 * (-u_ij*Diff[0]
		   +Diff[1]);
      
	// Compute "R_ij * |Lbd_ij| * L_ij * dU"
	if (btransposeDest) {
	  IDX2T(VectorAtEdge,1,iposDest,NVAR2D,nedgesimDest) += a[1] * ( w1 + w2 + w3 );
	  IDX2T(VectorAtEdge,2,iposDest,NVAR2D,nedgesimDest) += a[1] * (        u_ij*w1 + u_ij*w2 +        u_ij*w3 + w4 );
	  IDX2T(VectorAtEdge,3,iposDest,NVAR2D,nedgesimDest) += a[1] * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 );
	  IDX2T(VectorAtEdge,4,iposDest,NVAR2D,nedgesimDest) += a[1] * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 +
									 (H_ij+c_ij*v_ij)*w3 + u_ij*w4 );
	}
	else {
	  IDX2(VectorAtEdge,1,iposDest,NVAR2D,nedgesimDest) += a[1] * ( w1 + w2 + w3 );
	  IDX2(VectorAtEdge,2,iposDest,NVAR2D,nedgesimDest) += a[1] * (        u_ij*w1 + u_ij*w2 +        u_ij*w3 + w4 );
	  IDX2(VectorAtEdge,3,iposDest,NVAR2D,nedgesimDest) += a[1] * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 );
	  IDX2(VectorAtEdge,4,iposDest,NVAR2D,nedgesimDest) += a[1] * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 +
									(H_ij+c_ij*v_ij)*w3 + u_ij*w4 );
	}
      } else {
	if (btransposeDest) {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) = 0.0;
	}
	else {
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) = 0.0;
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
	Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc);
	Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc);
      }
      else {
	// Compute specific energies
	Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc);
	Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc);
      }
    
      // Compute the speed of sound
      Td ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ei-RCONST(0.5)*(ui*ui+vi*vi)), DBL_EPSILON));
      Td cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ej-RCONST(0.5)*(uj*uj+vj*vj)), DBL_EPSILON));
    
#ifdef HYDRO_USE_IBP
      // Compute scalar dissipation based on the skew-symmetric part
      // which does not include the symmetric boundary contribution
      Td d_ij = max( abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge))*uj+
			 RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge))*vj)+
		     RCONST(0.5)*sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
					  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
				      POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
					  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),2))*cj,
		     abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge))*ui+
			 RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge))*vi)+
		     RCONST(0.5)*sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
					  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
				      POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
					  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),2))*ci );
#else
      // Compute scalar dissipation
      Td d_ij = max( abs(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*uj+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*vj)+
		     sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
			  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),2))*cj,
		     abs(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*ui+
			 IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*vi)+
		     sqrt(POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
			  POW(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),2))*ci );
#endif
    
      // Multiply the solution difference by the scalar dissipation
      if (btransposeDest) {
	if (btransposeSrc) {
	  // Both source and destination vector are transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Destination vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
      }
      else {
	if (btransposeSrc) {
	  // Source vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Both vectors are not transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
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
	Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3T,1,iposSrc,NVAR2D,2,nedgesimSrc);
	Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3T,2,iposSrc,NVAR2D,2,nedgesimSrc);
      }
      else {
	// Compute specific energies
	Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,iposSrc,NVAR2D,2,nedgesimSrc);
	Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,iposSrc,NVAR2D,2,nedgesimSrc);
      }
    
      // Compute the speed of sound
      Td ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ei-RCONST(0.5)*(ui*ui+vi*vi)), DBL_EPSILON));
      Td cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ej-RCONST(0.5)*(uj*uj+vj*vj)), DBL_EPSILON));
    
#ifdef HYDRO_USE_IBP
      // Compute scalar dissipation with dimensional splitting based on
      // the skew-symmetric part which does not include the symmetric
      // boundary contribution
      Td d_ij = max( abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge))*uj)+
		     abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)))*cj,
		     abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge))*ui)+
		     abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)))*ci )
 	      + max( abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge))*vj)+
		     abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)))*cj,
		     abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge))*vi)+
		     abs(RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				      IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)))*ci );
#else
      // Compute scalar dissipation with dimensional splitting
      Td d_ij = max( abs(IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*uj)+
		     abs(IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge))*cj,
		     abs(IDX3(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*ui)+
		     abs(IDX3(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge))*ci )
	      + max( abs(IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*vj)+
		     abs(IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge))*cj,
		     abs(IDX3(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*vi)+
		     abs(IDX3(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge))*ci );
#endif
    
      // Multiply the solution difference by the scalar dissipation
if (btransposeDest) {
	if (btransposeSrc) {
	  // Both source and destination vector are transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Destination vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2T(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
      }
      else {
	if (btransposeSrc) {
	  // Source vector is transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3T(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3T(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
	}
	else {
	  // Both vectors are not transposed
#pragma unroll
	  for (int i=1; i<=NVAR2D; i++)
	    IDX2(VectorAtEdge,i,iposDest,NVAR2D,nedgesimDest) =
	      d_ij*(IDX3(DataAtEdge,i,2,iposSrc,NVAR2D,2,nedgesimSrc)
		   -IDX3(DataAtEdge,i,1,iposSrc,NVAR2D,2,nedgesimSrc));
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
			     Td pi,
			     Td pj,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      InviscidFluxDissipationBase<idissipationtype>::calcEdgeData<1,1,false,false>
	(VectorAtEdge,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,pi,pj,1,1,iedge,nedge,ncoeff);
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
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fxj,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fyj,i,iposSrc,NVAR2D,nedgesimSrc)-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fxi,i,iposSrc,NVAR2D,nedgesimSrc)-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fyi,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc));

#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,iposDest,NVAR2D,2,nedgesimDest) = -IDX3(FluxesAtEdge,i,1,iposDest,NVAR2D,2,nedgesimDest);
      }
      else {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fxj,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fyj,i,iposSrc,NVAR2D,nedgesimSrc)-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fxi,i,iposSrc,NVAR2D,nedgesimSrc)-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fyi,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc));
	
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fxj,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fyj,i,iposSrc,NVAR2D,nedgesimSrc)-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fxi,i,iposSrc,NVAR2D,nedgesimSrc)-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fyi,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc));
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
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,iposDest,NVAR2D,2,nedgesimDest) = scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc));
    
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,iposDest,NVAR2D,2,nedgesimDest) = -scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc));
      }
      else {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,iposDest,NVAR2D,2,nedgesimDest) += scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc));
    
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,iposDest,NVAR2D,2,nedgesimDest) -= scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fx_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*IDX2(Fy_ij,i,iposSrc,NVAR2D,nedgesimSrc)+
	     IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc));
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
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,iposDest,NVAR2D,2,nedgesimDest) = scale * IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc);
    
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,iposDest,NVAR2D,2,nedgesimDest) = -scale * IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc);
      }
      else {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,iposDest,NVAR2D,2,nedgesimDest) += scale * IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc);
    
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,iposDest,NVAR2D,2,nedgesimDest) -= scale * IDX2(Diff,i,iposSrc,NVAR2D,nedgesimSrc);
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
				Td *Diff,
				Td scale,       
				Ti iedge, 
				Ti nedge,
				Ti ncoeff)
    {
      FluxBase::combineEdgeData<1,1,boverwrite>
	(FluxesAtEdge,CoeffsAtEdge,Fxi,Fxj,Fyi,Fyj,Diff,scale,1,1,iedge,nedge,ncoeff);
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
				Td *Diff,
				Td scale,       
				Ti iedge, 
				Ti nedge,
				Ti ncoeff)
    {
      FluxBase::combineEdgeData<1,1,boverwrite>
	(FluxesAtEdge,CoeffsAtEdge,Fx_ij,Fy_ij,Diff,scale,1,1,iedge,nedge,ncoeff);
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
  __global__ void hydro_calcFlux2d_baseline(Tc *CoeffsAtEdge,
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
	  TdDest DataAtEdge[2*NVAR2D];

	  // Get solution values at edge endpoints
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    gatherEdgeData<true>(DataAtEdge,vecSrc,i,j,neq);
	  
	  // Compute velocities
	  TdDest ui = XVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
	  TdDest vi = YVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
	  
	  TdDest uj = XVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);
	  TdDest vj = YVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);
	  
	  // Compute pressures
	  TdDest pi = PRESSURE2(DataAtEdge,IDX2,1,NVAR2D,2);
	  TdDest pj = PRESSURE2(DataAtEdge,IDX2,2,NVAR2D,2);
	  
	  // Local variables
	  TdDest FluxAtEdge[2*NVAR2D];

	  // Compute the artificial viscosities
	  InviscidFluxDissipation<idissipationtype>::calcEdgeData
	    (&FluxAtEdge[NVAR2D],CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,pi,pj,idx+1,nedge,ncoeff);
	  Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);

	  // Compute inviscid fluxes
	  InviscidFlux::calcEdgeData<false>
	    (FluxAtEdge,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,pi,pj,scale,idx+1,nedge,ncoeff);

	  // Build fluxes into nodal vector
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
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
  __global__ void hydro_calcFlux2d_shmem(Tc *CoeffsAtEdge,
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
    __shared__ TdSrc s_DataAtEdge[2*NVAR2D*threads_per_cta];
    
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
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    gatherEdgeData<threads_per_cta,SHMEM_DATA_TRANSPOSE,true>
	    (s_DataAtEdge,vecSrc,(int)threadIdx.x+1,i,j,neq);
	  
	  // Compute velocities
	  TdDest ui = XVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,1,(int)threadIdx.x+1,NVAR2D,2,threads_per_cta);
	  TdDest vi = YVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,1,(int)threadIdx.x+1,NVAR2D,2,threads_per_cta);
	  
	  TdDest uj = XVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,2,(int)threadIdx.x+1,NVAR2D,2,threads_per_cta);
	  TdDest vj = YVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,2,(int)threadIdx.x+1,NVAR2D,2,threads_per_cta);
	  
	  // Compute pressures
	  TdDest pi = PRESSURE3(s_DataAtEdge,SHMEM_DATA_IDX3,1,(int)threadIdx.x+1,NVAR2D,2,threads_per_cta);
	  TdDest pj = PRESSURE3(s_DataAtEdge,SHMEM_DATA_IDX3,2,(int)threadIdx.x+1,NVAR2D,2,threads_per_cta);
	  
	  // Local variables
	  TdDest FluxAtEdge[2*NVAR2D];
	  
	  // Compute the artificial viscosities
	  InviscidFluxDissipation<idissipationtype>::
	    calcEdgeData<1,threads_per_cta,false,SHMEM_DATA_TRANSPOSE>
	    (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,pi,pj,
	     1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	  Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	  
	  // Compute inviscid fluxes
	  InviscidFlux::calcEdgeData<1,threads_per_cta,false,SHMEM_DATA_TRANSPOSE,false>
	    (FluxAtEdge,CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,pi,pj,
	     scale,1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	  
	  // Build fluxes into nodal vector
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);
	}
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (cudaDMA implementation
   * without warp specialisation).
   ****************************************************************************/

#define TOTAL_THREADS_PER_CTA  compute_threads_per_cta+dma_threads_per_ld* \
  (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)
  
  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipationtype,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __launch_bounds__(TOTAL_THREADS_PER_CTA)
  __global__ void hydro_calcFlux2d_cudaDMA_nospec(Tc *CoeffsAtEdge,
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
    __shared__ Ti s_IedgeList[2*compute_threads_per_cta];
    __shared__ TdSrc s_DataAtEdge[NVAR2D*2*compute_threads_per_cta];
    __shared__ Tc s_CoeffsAtEdge[HYDRO_NDIM*2*compute_threads_per_cta];
    
    //--------------------------------------------------------------------------

#if EDGELIST_DEVICE == SOA
    // List of edges is stored as structure of arrays, that is, we
    // have 6 integer subarrays of length nedge which store:
    //
    // 0-subarray: first end point i, 
    // 1-subarray: second end point j,
    // 2-subarray: matrix entry ij,
    // 3-subarray: matrix entry ji,
    // 4-subarray: matrix entry ii,
    // 5-subarray: matrix entry jj.
    //
    // For the flux assembly, only the two endpoints (i,j) are
    // required. Therefore, only subarrays 0 and 1 are transfered.

    // Sequential cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMASequential<false, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      TOTAL_THREADS_PER_CTA>dma_ind;
#else
    // List of edges is stored as array of structures, that is, we
    // have nedge integer subarrays of length 6 which store:
    //
    // (i,j,ij,jj,ii,jj) for each edge iedge
    //
    // For the flux assembly, only the two entpoins (i,j) are
    // required. Therefore, only the first two entries of each edge
    // are transfered using strided DMA.

    // Strided cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMAStrided<false, 2*sizeof(Ti), 2*sizeof(Ti), 
                   TOTAL_THREADS_PER_CTA,
		   compute_threads_per_cta>dma_ind(6*sizeof(Ti));
#endif

    //--------------------------------------------------------------------------

    // Indirect cudaDMA thread to transfer nodal data from vecSrc into
    // shared memory s_DataAtEdge, we need to distinguish between vecSrc
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, false,
      MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
               (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
      TOTAL_THREADS_PER_CTA, 2*compute_threads_per_cta>dma_vec;
    
    //--------------------------------------------------------------------------

#if COEFFSATEDGE_DEVICE == SOA
    // Coefficients at edges are stored as structure of arrays, that
    // is, we have 2*ncoeff subarrays of length nedge which store:
    //
    // 0-subarray: ij-coefficients for x-direction, 
    // 1-subarray: ji-coefficients for x-direction, 
    // 2-subarray: ij-coefficients for y-direction,
    // 3-subarray: ji-coefficients for y-direction,
    // ...
    // n-subarray: further coefficients not required here

    // Strided cudaDMA thread to transfer precomputed coefficients
    // CoeffsAtEdge into shared memory s_CoeffsAtEdge
    cudaDMAStrided<false, sizeof(Tc),
                   compute_threads_per_cta*sizeof(Tc),
                   TOTAL_THREADS_PER_CTA,
                   2*HYDRO_NDIM>dma_coeff(nedge*sizeof(Tc));
#else
    // Coefficients at edges are stored as array of structure, that
    // is, we have nedge real-valued subarray of length 2*ncoeff
    cudaDMAStrided<false, sizeof(Tc), HYDRO_NDIM*sizeof(Tc),
                   TOTAL_THREADS_PER_CTA,
                   2*compute_threads_per_cta>dma_coeff(ncoeff*sizeof(Tc));
#endif

    //--------------------------------------------------------------------------

    // Loop over all edge-groups to be processed by this block
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {

      if (nedge_per_thread>1)
	__syncthreads();

      //------------------------------------------------------------------------
      // Load the indices with all threads - no warp specialisation
      //------------------------------------------------------------------------
      dma_ind.execute_dma(&IedgeList[ ((ipt*gridDim.x+blockIdx.x)*
				       compute_threads_per_cta+nedge_offset)*
      				      (EDGELIST_DEVICE == SOA ? 2 : 6)], s_IedgeList);
      __syncthreads();

      dma_vec.execute_dma(s_IedgeList, vecSrc-NVAR2D, s_DataAtEdge);
      dma_coeff.execute_dma(&CoeffsAtEdge[ ((ipt*gridDim.x+blockIdx.x)*
					    compute_threads_per_cta+nedge_offset)*
					   (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
			    s_CoeffsAtEdge);
      __syncthreads();
      
      //--------------------------------------------------------------------------
      
      // Compute velocities
      TdDest ui = XVELOCITY3(s_DataAtEdge,
			     IDX3,1,(int)threadIdx.x+1,
			     NVAR2D,2,compute_threads_per_cta);
      TdDest vi = YVELOCITY3(s_DataAtEdge,
			     IDX3,1,(int)threadIdx.x+1,
			     NVAR2D,2,compute_threads_per_cta);
      
      TdDest uj = XVELOCITY3(s_DataAtEdge,
			     IDX3,2,(int)threadIdx.x+1,
			     NVAR2D,2,compute_threads_per_cta);
      TdDest vj = YVELOCITY3(s_DataAtEdge,
			     IDX3,2,(int)threadIdx.x+1,
			     NVAR2D,2,compute_threads_per_cta);
      
      // Compute pressures
      TdDest pi = PRESSURE3(s_DataAtEdge,
			    IDX3,1,(int)threadIdx.x+1,
			    NVAR2D,2,compute_threads_per_cta);
      TdDest pj = PRESSURE3(s_DataAtEdge,
			    IDX3,2,(int)threadIdx.x+1,
			    NVAR2D,2,compute_threads_per_cta);

      // Local variables
      TdDest FluxAtEdge[2*NVAR2D];
      
      // Compute the artificial viscosities
      InviscidFluxDissipation<idissipationtype>::
	calcEdgeData<1,compute_threads_per_cta,false,false>
	(&FluxAtEdge[NVAR2D],s_CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,pi,pj,
	 1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,2);
      Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
      
      // Compute inviscid fluxes
      InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	(FluxAtEdge,s_CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,pi,pj,scale,
	 1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,2);
      
      dma_vec.execute_dma(s_IedgeList, vecDest-NVAR2D, s_DataAtEdge);
      __syncthreads();
      
      // Get positions of edge endpoints (idx starts at zero)
      Ti i = IDX2(s_IedgeList,1,(int)threadIdx.x+1,2,compute_threads_per_cta);
      
#pragma unroll
      for (int ivar=1; ivar<=NVAR2D; ++ivar)
	IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	  IDX3(s_DataAtEdge, ivar, 1, (int)threadIdx.x+1,
	       NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
      
      // Get positions of edge endpoints (idx starts at zero)
      Ti j = IDX2(s_IedgeList,2,(int)threadIdx.x+1,2,compute_threads_per_cta);
      
#pragma unroll
      for (int ivar=1; ivar<=NVAR2D; ++ivar)
	IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	  IDX3(s_DataAtEdge, ivar, 2, (int)threadIdx.x+1,
	       NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
    }
  };

#undef TOTAL_THREADS_PER_CTA

 /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (cudaDMA implementation with
   * manual single buffering strategy with prefetching of indices).
   ****************************************************************************/

#define TOTAL_THREADS_PER_CTA compute_threads_per_cta+dma_threads_per_ld* \
  (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST+CUDADMA_DMA_LDS_COEFF)

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipationtype,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __launch_bounds__(TOTAL_THREADS_PER_CTA)
  __global__ void hydro_calcFlux2d_cudaDMA_prefetch_single(Tc *CoeffsAtEdge,
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
    __shared__ Ti s_IedgeList[1][2*compute_threads_per_cta];
    __shared__ TdSrc s_VecSrc[1][NVAR2D*2*compute_threads_per_cta];
    __shared__ TdDest s_VecDest[1][NVAR2D*2*compute_threads_per_cta];
    __shared__ Tc s_CoeffsAtEdge[1][HYDRO_NDIM*2*compute_threads_per_cta];

    //--------------------------------------------------------------------------

#if EDGELIST_DEVICE == SOA
    // List of edges is stored as structure of arrays, that is, we
    // have 6 integer subarrays of length nedge which store:
    //
    // 0-subarray: first end point i, 
    // 1-subarray: second end point j,
    // 2-subarray: matrix entry ij,
    // 3-subarray: matrix entry ji,
    // 4-subarray: matrix entry ii,
    // 5-subarray: matrix entry jj.
    //
    // For the flux assembly, only the two endpoints (i,j) are
    // required. Therefore, only subarrays 0 and 1 are transfered.

    // Sequential cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMASequential<false, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      TOTAL_THREADS_PER_CTA>dma_ind;
#else
    // List of edges is stored as array of structures, that is, we
    // have nedge integer subarrays of length 6 which store:
    //
    // (i,j,ij,jj,ii,jj) for each edge iedge
    //
    // For the flux assembly, only the two entpoins (i,j) are
    // required. Therefore, only the first two entries of each edge
    // are transfered using strided DMA.

    // Strided cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMAStrided<false, 2*sizeof(Ti), 2*sizeof(Ti), 
                   TOTAL_THREADS_PER_CTA,
		   compute_threads_per_cta>dma_ind(6*sizeof(Ti));
#endif

    //--------------------------------------------------------------------------

    // Indirect cudaDMA thread to transfer nodal data from vecSrc into
    // shared memory s_VecSrc, we need to distinguish between vecSrc
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_SRC*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecSrc0(0, compute_threads_per_cta, compute_threads_per_cta);

    // Indirect cudaDMA thread to transfer nodal data from vecDest into
    // shared memory s_VecDest, we need to distinguish between vecDest
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_DEST*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecDest0(1, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC)*dma_threads_per_ld);

    //--------------------------------------------------------------------------

#if COEFFSATEDGE_DEVICE == SOA
    // Coefficients at edges are stored as structure of arrays, that
    // is, we have 2*ncoeff subarrays of length nedge which store:
    //
    // 0-subarray: ij-coefficients for x-direction, 
    // 1-subarray: ji-coefficients for x-direction, 
    // 2-subarray: ij-coefficients for y-direction,
    // 3-subarray: ji-coefficients for y-direction,
    // ...
    // n-subarray: further coefficients not required here

    // Strided cudaDMA thread to transfer precomputed coefficients
    // CoeffsAtEdge into shared memory s_CoeffsAtEdge
    cudaDMAStrided<true, sizeof(Tc),
                   compute_threads_per_cta*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*HYDRO_NDIM>
      dma_coeff0(2, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
		 nedge*sizeof(Tc));
#else
    // Coefficients at edges are stored as array of structure, that
    // is, we have nedge real-valued subarray of length 2*ncoeff
    cudaDMAStrided<true, sizeof(Tc), HYDRO_NDIM*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*compute_threads_per_cta>
      dma_coeff0(2, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC+0*CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
		 ncoeff*sizeof(Tc));
#endif

    //--------------------------------------------------------------------------

    // Loop over all edge-groups to be processed by this block
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {

      //------------------------------------------------------------------------
      // Load the indices with all threads - no warp specialisation
      //------------------------------------------------------------------------
      if (nedge_per_thread>1)
	ptx_cudaDMA_barrier_blocking(5, TOTAL_THREADS_PER_CTA);

      dma_ind.execute_dma(&IedgeList[ ((ipt*gridDim.x+blockIdx.x)*
				       compute_threads_per_cta+nedge_offset)*
				      (EDGELIST_DEVICE == SOA ? 2 : 6)],
			  s_IedgeList[0]);

      ptx_cudaDMA_barrier_blocking(5, TOTAL_THREADS_PER_CTA);
      
      //------------------------------------------------------------------------
      // Warp specialisation
      //------------------------------------------------------------------------
      if (threadIdx.x<compute_threads_per_cta) {
      
	// Start DMA transfer of coefficients
	dma_coeff0.start_async_dma();
	
	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();
	
#define IBUF 0
#define DBUF 0
#define IOFF 0
	
	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	TdDest ui = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vi = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	TdDest uj = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vj = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	TdDest pi = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,1,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	TdDest pj = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,2,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	
	// Local variables
	TdDest FluxAtEdge[2*NVAR2D];
	
	// Wait for coefficients to be ready
	dma_coeff0.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],s_CoeffsAtEdge[0],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,s_CoeffsAtEdge[0],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	
	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
      
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);

#undef IBUF
#undef DBUF
#undef IOFF
      }

      //------------------------------------------------------------------------
      // DMA transfer warps
      //------------------------------------------------------------------------

      else if(dma_vecSrc0.owns_this_thread()) {
	dma_vecSrc0.execute_dma(s_IedgeList[0], vecSrc-NVAR2D, s_VecSrc[0]);
      }
      
      else if(dma_vecDest0.owns_this_thread()) {
      	dma_vecDest0.execute_dma(s_IedgeList[0], vecDest-NVAR2D, s_VecDest[0]);
      }

      else if(dma_coeff0.owns_this_thread()) {
	dma_coeff0.execute_dma(&CoeffsAtEdge[((ipt*gridDim.x+blockIdx.x)*
					      compute_threads_per_cta +
					      nedge_offset)*
					     (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
      			       s_CoeffsAtEdge[0]);
      }
    }
  };

#undef TOTAL_THREADS_PER_CTA

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (cudaDMA implementation with
   * manual double buffering strategy with prefetching of indices).
   ****************************************************************************/

#define TOTAL_THREADS_PER_CTA compute_threads_per_cta+dma_threads_per_ld* \
  (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST+CUDADMA_DMA_LDS_COEFF)

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipationtype,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __launch_bounds__(TOTAL_THREADS_PER_CTA)
  __global__ void hydro_calcFlux2d_cudaDMA_prefetch_double(Tc *CoeffsAtEdge,
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
    __shared__ Ti s_IedgeList[4][2*compute_threads_per_cta];
    __shared__ TdSrc s_VecSrc[2][NVAR2D*2*compute_threads_per_cta];
    __shared__ TdDest s_VecDest[2][NVAR2D*2*compute_threads_per_cta];
    __shared__ Tc s_CoeffsAtEdge[2][HYDRO_NDIM*2*compute_threads_per_cta];

    //--------------------------------------------------------------------------

#if EDGELIST_DEVICE == SOA
    // List of edges is stored as structure of arrays, that is, we
    // have 6 integer subarrays of length nedge which store:
    //
    // 0-subarray: first end point i, 
    // 1-subarray: second end point j,
    // 2-subarray: matrix entry ij,
    // 3-subarray: matrix entry ji,
    // 4-subarray: matrix entry ii,
    // 5-subarray: matrix entry jj.
    //
    // For the flux assembly, only the two endpoints (i,j) are
    // required. Therefore, only subarrays 0 and 1 are transfered.

    // Sequential cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMASequential<false, 2*sizeof(Ti),
                      4*2*compute_threads_per_cta*sizeof(Ti),
                      TOTAL_THREADS_PER_CTA>dma_ind;
#else
    // List of edges is stored as array of structures, that is, we
    // have nedge integer subarrays of length 6 which store:
    //
    // (i,j,ij,jj,ii,jj) for each edge iedge
    //
    // For the flux assembly, only the two entpoins (i,j) are
    // required. Therefore, only the first two entries of each edge
    // are transfered using strided DMA.

    // Strided cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMAStrided<false, 2*sizeof(Ti), 2*sizeof(Ti), 
                   TOTAL_THREADS_PER_CTA,
		   4*compute_threads_per_cta>dma_ind(6*sizeof(Ti));
#endif

    //--------------------------------------------------------------------------

    // Indirect cudaDMA thread to transfer nodal data from vecSrc into
    // shared memory s_VecSrc, we need to distinguish between vecSrc
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_SRC*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecSrc0(0, compute_threads_per_cta, compute_threads_per_cta);

    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_SRC*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecSrc1(1, compute_threads_per_cta, compute_threads_per_cta);

    // Indirect cudaDMA thread to transfer nodal data from vecDest into
    // shared memory s_VecDest, we need to distinguish between vecDest
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_DEST*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecDest0(2, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC)*dma_threads_per_ld);

    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_DEST*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecDest1(3, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC)*dma_threads_per_ld);

    //--------------------------------------------------------------------------

#if COEFFSATEDGE_DEVICE == SOA
    // Coefficients at edges are stored as structure of arrays, that
    // is, we have 2*ncoeff subarrays of length nedge which store:
    //
    // 0-subarray: ij-coefficients for x-direction, 
    // 1-subarray: ji-coefficients for x-direction, 
    // 2-subarray: ij-coefficients for y-direction,
    // 3-subarray: ji-coefficients for y-direction,
    // ...
    // n-subarray: further coefficients not required here

    // Strided cudaDMA thread to transfer precomputed coefficients
    // CoeffsAtEdge into shared memory s_CoeffsAtEdge
    cudaDMAStrided<true, sizeof(Tc),
                   compute_threads_per_cta*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*HYDRO_NDIM>
      dma_coeff0(4, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
		 nedge*sizeof(Tc));

    cudaDMAStrided<true, sizeof(Tc),
                   compute_threads_per_cta*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*HYDRO_NDIM>
      dma_coeff1(5, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
		 nedge*sizeof(Tc));
#else
    // Coefficients at edges are stored as array of structure, that
    // is, we have nedge real-valued subarray of length 2*ncoeff
    cudaDMAStrided<true, sizeof(Tc), HYDRO_NDIM*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*compute_threads_per_cta>
      dma_coeff0(4, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
		 ncoeff*sizeof(Tc));
  
    cudaDMAStrided<true, sizeof(Tc), HYDRO_NDIM*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*compute_threads_per_cta>
      dma_coeff1(5, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
		 ncoeff*sizeof(Tc));
#endif

    //--------------------------------------------------------------------------

    // Loop over all edge-groups to be processed by this block
    for (int ipt=0; ipt<nedge_per_thread; ipt+=4) {

      //------------------------------------------------------------------------
      // Load the indices with all threads - no warp specialisation
      //------------------------------------------------------------------------
      ptx_cudaDMA_barrier_blocking(11, TOTAL_THREADS_PER_CTA);
      dma_ind.execute_dma(&IedgeList[ ((ipt*gridDim.x+blockIdx.x)*
				       4*compute_threads_per_cta+nedge_offset)*
				      (EDGELIST_DEVICE == SOA ? 2 : 6)],
			  s_IedgeList[0]);
      ptx_cudaDMA_barrier_blocking(11, TOTAL_THREADS_PER_CTA);
      
      //------------------------------------------------------------------------
      // Warp specialisation
      //------------------------------------------------------------------------
      if (threadIdx.x<compute_threads_per_cta) {
      
	// Start DMA transfer of coefficients
	dma_coeff0.start_async_dma();
	dma_coeff1.start_async_dma();
	
	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();
	
	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();
	
#define IBUF 0
#define DBUF 0
#define IOFF 0
	
	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	TdDest ui = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vi = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	TdDest uj = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vj = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	TdDest pi = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,1,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	TdDest pj = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,2,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	
	// Local variables
	TdDest FluxAtEdge[2*NVAR2D];
	
	// Wait for coefficients to be ready
	dma_coeff0.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],s_CoeffsAtEdge[0],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,s_CoeffsAtEdge[0],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	
	// Start DMA transfer of coefficients
	dma_coeff0.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
      
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);

#undef IBUF
#undef DBUF
#undef IOFF

// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();

#define IBUF 1
#define DBUF 1
#define IOFF 1

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
		
	// Wait for coefficients to be ready
	dma_coeff1.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],s_CoeffsAtEdge[1],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,s_CoeffsAtEdge[1],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	
	// Start DMA transfer of coefficients
	dma_coeff1.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();

#define IBUF 2
#define DBUF 0
#define IOFF 2

	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
		
	// Wait for coefficients to be ready
	dma_coeff0.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],s_CoeffsAtEdge[0],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,s_CoeffsAtEdge[0],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	
	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

#define IBUF 3
#define DBUF 1
#define IOFF 3

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
		
	// Wait for coefficients to be ready
	dma_coeff1.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],s_CoeffsAtEdge[1],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,s_CoeffsAtEdge[1],s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,(int)threadIdx.x+1,compute_threads_per_cta,HYDRO_NDIM);
	
	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

      }

      //------------------------------------------------------------------------
      // DMA transfer warps
      //------------------------------------------------------------------------

      else if(dma_vecSrc0.owns_this_thread()) {
	dma_vecSrc0.execute_dma(s_IedgeList[0], vecSrc-NVAR2D, s_VecSrc[0]);
	dma_vecSrc1.execute_dma(s_IedgeList[1], vecSrc-NVAR2D, s_VecSrc[1]);
	dma_vecSrc0.execute_dma(s_IedgeList[2], vecSrc-NVAR2D, s_VecSrc[0]);
	dma_vecSrc1.execute_dma(s_IedgeList[3], vecSrc-NVAR2D, s_VecSrc[1]);
      }
      
      else if(dma_vecDest0.owns_this_thread()) {
	dma_vecDest0.execute_dma(s_IedgeList[0], vecDest-NVAR2D, s_VecDest[0]);
	dma_vecDest1.execute_dma(s_IedgeList[1], vecDest-NVAR2D, s_VecDest[1]);
	dma_vecDest0.execute_dma(s_IedgeList[2], vecDest-NVAR2D, s_VecDest[0]);
	dma_vecDest1.execute_dma(s_IedgeList[3], vecDest-NVAR2D, s_VecDest[1]);
      }

      else if(dma_coeff0.owns_this_thread()) {
	dma_coeff0.execute_dma(&CoeffsAtEdge[((ipt*gridDim.x+blockIdx.x)*
					      4*compute_threads_per_cta +
					      0*compute_threads_per_cta +
					      nedge_offset)*
					     (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
      			       s_CoeffsAtEdge[0]);
	dma_coeff1.execute_dma(&CoeffsAtEdge[((ipt*gridDim.x+blockIdx.x)*
					      4*compute_threads_per_cta +
					      1*compute_threads_per_cta +
					      nedge_offset)*
					     (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
      			       s_CoeffsAtEdge[1]);
	dma_coeff0.execute_dma(&CoeffsAtEdge[((ipt*gridDim.x+blockIdx.x)*
					      4*compute_threads_per_cta +
					      2*compute_threads_per_cta +
					      nedge_offset)*
					     (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
      			       s_CoeffsAtEdge[0]);
	dma_coeff1.execute_dma(&CoeffsAtEdge[((ipt*gridDim.x+blockIdx.x)*
					      4*compute_threads_per_cta +
					      3*compute_threads_per_cta +
					      nedge_offset)*
					     (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
      			       s_CoeffsAtEdge[1]);
      }
    }
  };

#undef TOTAL_THREADS_PER_CTA

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (cudaDMA implementation with
   * double buffering strategy).
   ****************************************************************************/

#define TOTAL_THREADS_PER_CTA compute_threads_per_cta+dma_threads_per_ld* \
  (3*CUDADMA_DMA_LDS_IND+2*CUDADMA_DMA_LDS_SRC+2*CUDADMA_DMA_LDS_DEST)

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipationtype,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __launch_bounds__(TOTAL_THREADS_PER_CTA)
  __global__ void hydro_calcFlux2d_cudaDMA_double(Tc *CoeffsAtEdge,
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
    __shared__ Ti s_IedgeList[3][2*compute_threads_per_cta];
    __shared__ TdSrc s_VecSrc[2][NVAR2D*2*compute_threads_per_cta];
    __shared__ TdDest s_VecDest[2][NVAR2D*2*compute_threads_per_cta];
    
    //--------------------------------------------------------------------------

#if EDGELIST_DEVICE == SOA
    // List of edges is stored as structure of arrays, that is, we
    // have 6 integer subarrays of length nedge which store:
    //
    // 0-subarray: first end point i, 
    // 1-subarray: second end point j,
    // 2-subarray: matrix entry ij,
    // 3-subarray: matrix entry ji,
    // 4-subarray: matrix entry ii,
    // 5-subarray: matrix entry jj.
    //
    // For the flux assembly, only the two endpoints (i,j) are
    // required. Therefore, only subarrays 0 and 1 are transfered.

    // Sequential cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMASequential<true, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      CUDADMA_DMA_LDS_IND*dma_threads_per_ld>
    dma_ind0(0, compute_threads_per_cta, compute_threads_per_cta);

    cudaDMASequential<true, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      CUDADMA_DMA_LDS_IND*dma_threads_per_ld>
    dma_ind1(1, compute_threads_per_cta,
	     compute_threads_per_cta+CUDADMA_DMA_LDS_IND*dma_threads_per_ld);

    cudaDMASequential<true, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      CUDADMA_DMA_LDS_IND*dma_threads_per_ld>
    dma_ind2(2, compute_threads_per_cta,
	     compute_threads_per_cta+2*CUDADMA_DMA_LDS_IND*dma_threads_per_ld);
#else
    // List of edges is stored as array of structures, that is, we
    // have nedge integer subarrays of length 6 which store:
    //
    // (i,j,ij,jj,ii,jj) for each edge iedge
    //
    // For the flux assembly, only the two entpoins (i,j) are
    // required. Therefore, only the first two entries of each edge
    // are transfered using strided DMA.

    // Strided cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMAStrided<true, 2*sizeof(Ti), 2*sizeof(Ti), 
                   CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
                   compute_threads_per_cta>
    dma_ind0(0, compute_threads_per_cta,
	     compute_threads_per_cta,
	     6*sizeof(Ti));

    cudaDMAStrided<true, 2*sizeof(Ti), 2*sizeof(Ti), 
                   CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
                   compute_threads_per_cta>
    dma_ind1(1, compute_threads_per_cta,
	     compute_threads_per_cta+CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
	     6*sizeof(Ti));

    cudaDMAStrided<true, 2*sizeof(Ti), 2*sizeof(Ti), 
                   CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
                   compute_threads_per_cta>
    dma_ind2(2, compute_threads_per_cta,
	     compute_threads_per_cta+2*CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
	     6*sizeof(Ti));
#endif

    //--------------------------------------------------------------------------

    // Indirect cudaDMA thread to transfer nodal data from vecSrc into
    // shared memory s_VecSrc, we need to distinguish between vecSrc
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_SRC*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecSrc0(3, compute_threads_per_cta,
		compute_threads_per_cta+3*CUDADMA_DMA_LDS_IND*dma_threads_per_ld);

    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_SRC*dma_threads_per_ld, 2*compute_threads_per_cta>
      dma_vecSrc1(4, compute_threads_per_cta, compute_threads_per_cta+
		  (3*CUDADMA_DMA_LDS_IND+1*CUDADMA_DMA_LDS_SRC)*dma_threads_per_ld);

    // Indirect cudaDMA thread to transfer nodal data from vecDest into
    // shared memory s_VecDest, we need to distinguish between vecDest
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_DEST*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecDest0(5, compute_threads_per_cta, compute_threads_per_cta+
		 (3*CUDADMA_DMA_LDS_IND+2*CUDADMA_DMA_LDS_SRC)*dma_threads_per_ld);

    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_DEST*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecDest1(6, compute_threads_per_cta, compute_threads_per_cta+
		 (3*CUDADMA_DMA_LDS_IND+2*CUDADMA_DMA_LDS_SRC
		 +1*CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld);

    //--------------------------------------------------------------------------

#if COEFFSATEDGE_DEVICE == SOA
    // Coefficients at edges are stored as structure of arrays, that
    // is, we have 2*ncoeff subarrays of length nedge which store:
    //
    // 0-subarray: ij-coefficients for x-direction, 
    // 1-subarray: ji-coefficients for x-direction, 
    // 2-subarray: ij-coefficients for y-direction,
    // 3-subarray: ji-coefficients for y-direction,
    // ...
    // n-subarray: further coefficients not required here

    // Strided cudaDMA thread to transfer precomputed coefficients
    // CoeffsAtEdge into shared memory s_CoeffsAtEdge
    // cudaDMAStrided<false, sizeof(Tc),
    //                compute_threads_per_cta*sizeof(Tc),
    //                TOTAL_THREADS_PER_CTA,
    //                2*HYDRO_NDIM>dma_coeff(nedge*sizeof(Tc));
#else
    // Coefficients at edges are stored as array of structure, that
    // is, we have nedge real-valued subarray of length 2*ncoeff
    // cudaDMAStrided<false, sizeof(Tc), HYDRO_NDIM*sizeof(Tc),
    //                TOTAL_THREADS_PER_CTA,
    //                2*compute_threads_per_cta>dma_coeff(ncoeff*sizeof(Tc));
#endif

    //--------------------------------------------------------------------------
    // Warp specialisation
    //--------------------------------------------------------------------------
    if (threadIdx.x<compute_threads_per_cta) {
      
      // Start DMA transfer of indices
      dma_ind0.start_async_dma();
      dma_ind1.start_async_dma();
      dma_ind2.start_async_dma();
      
      // Wait for indices to be ready
      dma_ind0.wait_for_dma_finish();
    
      // Start DMA transfer of indirect data
      dma_vecSrc0.start_async_dma();
      dma_vecDest0.start_async_dma();
      
      // Wait for indices to be ready
      dma_ind1.wait_for_dma_finish();

      // Start DMA transfer of indirect data
      dma_vecSrc1.start_async_dma();
      dma_vecDest1.start_async_dma();

      // Loop over all edge-groups to be processed by this block
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {
	
#define IBUF 0
#define DBUF 0
#define IOFF 0

	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	TdDest ui = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vi = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	TdDest uj = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vj = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	TdDest pi = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,1,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	TdDest pj = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,2,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	Ti idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	       + nedge_offset + threadIdx.x;
	
	// Local variables
	TdDest FluxAtEdge[2*NVAR2D];
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
      
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);

#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind0.start_async_dma();

	// Wait for indices to be ready
	dma_ind2.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();

#define IBUF 1
#define DBUF 1
#define IOFF 1

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind1.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind0.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();

#define IBUF 2
#define DBUF 0
#define IOFF 2

	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind2.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind1.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();

#define IBUF 0
#define DBUF 1
#define IOFF 3

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind0.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind2.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();

#define IBUF 1
#define DBUF 0
#define IOFF 4

	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind1.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind0.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();

#define IBUF 2
#define DBUF 1
#define IOFF 5

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind2.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind1.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();
      }
    }
    
    //--------------------------------------------------------------------------
    // DMA transfer warps
    //--------------------------------------------------------------------------
    else if(dma_ind0.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread+1; ipt+=3) {
	dma_ind0.execute_dma(&IedgeList[ ((ipt*gridDim.x+blockIdx.x)*
					  compute_threads_per_cta+nedge_offset)*
					 (EDGELIST_DEVICE == SOA ? 2 : 6)],
			     s_IedgeList[0]);
	}
    }

    else if(dma_ind1.owns_this_thread()) {
      for (int ipt=1; ipt<nedge_per_thread+2; ipt+=3) {
	dma_ind1.execute_dma(&IedgeList[ ((ipt*gridDim.x+blockIdx.x)*
					  compute_threads_per_cta+nedge_offset)*
					 (EDGELIST_DEVICE == SOA ? 2 : 6)],
			     s_IedgeList[1]);
	}
    }

    else if(dma_ind2.owns_this_thread()) {
      for (int ipt=2; ipt<nedge_per_thread+3; ipt+=3) {
	dma_ind2.execute_dma(&IedgeList[ ((ipt*gridDim.x+blockIdx.x)*
					  compute_threads_per_cta+nedge_offset)*
					 (EDGELIST_DEVICE == SOA ? 2 : 6)],
			     s_IedgeList[2]);
      }
    }   
    
    else if(dma_vecSrc0.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {
	
	dma_vecSrc0.wait_for_dma_start();
	dma_vecSrc0.execute_dma_no_sync(s_IedgeList[0], vecSrc-NVAR2D, s_VecSrc[0]);
	dma_vecSrc0.finish_async_dma();
	
      	dma_vecSrc0.wait_for_dma_start();
      	dma_vecSrc0.execute_dma_no_sync(s_IedgeList[2], vecSrc-NVAR2D, s_VecSrc[0]);
      	dma_vecSrc0.finish_async_dma();

      	dma_vecSrc0.wait_for_dma_start();
      	dma_vecSrc0.execute_dma_no_sync(s_IedgeList[1], vecSrc-NVAR2D, s_VecSrc[0]);
      	dma_vecSrc0.finish_async_dma();
      }
    }

    else if(dma_vecSrc1.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {
	
      	dma_vecSrc1.wait_for_dma_start();
      	dma_vecSrc1.execute_dma_no_sync(s_IedgeList[1], vecSrc-NVAR2D, s_VecSrc[1]);
      	dma_vecSrc1.finish_async_dma();

      	dma_vecSrc1.wait_for_dma_start();
      	dma_vecSrc1.execute_dma_no_sync(s_IedgeList[0], vecSrc-NVAR2D, s_VecSrc[1]);
      	dma_vecSrc1.finish_async_dma();

      	dma_vecSrc1.wait_for_dma_start();
      	dma_vecSrc1.execute_dma_no_sync(s_IedgeList[2], vecSrc-NVAR2D, s_VecSrc[1]);
      	dma_vecSrc1.finish_async_dma();
      }
    }

    else if(dma_vecDest0.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {
	
	dma_vecDest0.wait_for_dma_start();
	dma_vecDest0.execute_dma_no_sync(s_IedgeList[0], vecDest-NVAR2D, s_VecDest[0]);
	dma_vecDest0.finish_async_dma();

      	dma_vecDest0.wait_for_dma_start();
      	dma_vecDest0.execute_dma_no_sync(s_IedgeList[2], vecDest-NVAR2D, s_VecDest[0]);
      	dma_vecDest0.finish_async_dma();

      	dma_vecDest0.wait_for_dma_start();
      	dma_vecDest0.execute_dma_no_sync(s_IedgeList[1], vecDest-NVAR2D, s_VecDest[0]);
      	dma_vecDest0.finish_async_dma();
      }
    }

    else if(dma_vecDest1.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {

      	dma_vecDest1.wait_for_dma_start();
      	dma_vecDest1.execute_dma_no_sync(s_IedgeList[1], vecDest-NVAR2D, s_VecDest[1]);
      	dma_vecDest1.finish_async_dma();

      	dma_vecDest1.wait_for_dma_start();
      	dma_vecDest1.execute_dma_no_sync(s_IedgeList[0], vecDest-NVAR2D, s_VecDest[1]);
      	dma_vecDest1.finish_async_dma();

      	dma_vecDest1.wait_for_dma_start();
      	dma_vecDest1.execute_dma_no_sync(s_IedgeList[2], vecDest-NVAR2D, s_VecDest[1]);
      	dma_vecDest1.finish_async_dma();
      }
    }
  };
  
#undef TOTAL_THREADS_PER_CTA

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (cudaDMA implementation with
   * manual buffering strategy).
   ****************************************************************************/

#define TOTAL_THREADS_PER_CTA compute_threads_per_cta+dma_threads_per_ld* \
  (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST+CUDADMA_DMA_LDS_COEFF)

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipationtype,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __launch_bounds__(TOTAL_THREADS_PER_CTA)
  __global__ void hydro_calcFlux2d_cudaDMA_manual(Tc *CoeffsAtEdge,
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
    __shared__ Ti s_IedgeList[3][2*compute_threads_per_cta];
    __shared__ TdSrc s_VecSrc[2][NVAR2D*2*compute_threads_per_cta];
    __shared__ TdDest s_VecDest[2][NVAR2D*2*compute_threads_per_cta];
    __shared__ Tc s_CoeffsAtEdge[2][HYDRO_NDIM*2*compute_threads_per_cta];
    
    //--------------------------------------------------------------------------

#if EDGELIST_DEVICE == SOA
    // List of edges is stored as structure of arrays, that is, we
    // have 6 integer subarrays of length nedge which store:
    //
    // 0-subarray: first end point i, 
    // 1-subarray: second end point j,
    // 2-subarray: matrix entry ij,
    // 3-subarray: matrix entry ji,
    // 4-subarray: matrix entry ii,
    // 5-subarray: matrix entry jj.
    //
    // For the flux assembly, only the two endpoints (i,j) are
    // required. Therefore, only subarrays 0 and 1 are transfered.

    // Sequential cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMASequential<true, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      CUDADMA_DMA_LDS_IND*dma_threads_per_ld>
    dma_ind0(0, compute_threads_per_cta, compute_threads_per_cta);

    cudaDMASequential<true, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      CUDADMA_DMA_LDS_IND*dma_threads_per_ld>
    dma_ind1(1, compute_threads_per_cta, compute_threads_per_cta);

    cudaDMASequential<true, 2*sizeof(Ti),
                      2*compute_threads_per_cta*sizeof(Ti),
                      CUDADMA_DMA_LDS_IND*dma_threads_per_ld>
    dma_ind2(2, compute_threads_per_cta, compute_threads_per_cta);
#else
    // List of edges is stored as array of structures, that is, we
    // have nedge integer subarrays of length 6 which store:
    //
    // (i,j,ij,jj,ii,jj) for each edge iedge
    //
    // For the flux assembly, only the two entpoins (i,j) are
    // required. Therefore, only the first two entries of each edge
    // are transfered using strided DMA.

    // Strided cudaDMA thread to transfer edge list from integer
    // array IedgeList into shared memory s_IedgeList
    cudaDMAStrided<true, 2*sizeof(Ti), 2*sizeof(Ti), 
                   CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
                   compute_threads_per_cta>
    dma_ind0(0, compute_threads_per_cta, compute_threads_per_cta, 6*sizeof(Ti));

    cudaDMAStrided<true, 2*sizeof(Ti), 2*sizeof(Ti), 
                   CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
                   compute_threads_per_cta>
    dma_ind1(1, compute_threads_per_cta, compute_threads_per_cta, 6*sizeof(Ti));

    cudaDMAStrided<true, 2*sizeof(Ti), 2*sizeof(Ti), 
                   CUDADMA_DMA_LDS_IND*dma_threads_per_ld,
                   compute_threads_per_cta>
    dma_ind2(2, compute_threads_per_cta, compute_threads_per_cta, 6*sizeof(Ti));
#endif

    //--------------------------------------------------------------------------

    // Indirect cudaDMA thread to transfer nodal data from vecSrc into
    // shared memory s_VecSrc, we need to distinguish between vecSrc
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_SRC*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecSrc0(3, compute_threads_per_cta,
		compute_threads_per_cta+CUDADMA_DMA_LDS_IND*dma_threads_per_ld);

    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_SRC*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecSrc1(4, compute_threads_per_cta,
		compute_threads_per_cta+CUDADMA_DMA_LDS_IND*dma_threads_per_ld);

    // Indirect cudaDMA thread to transfer nodal data from vecDest into
    // shared memory s_VecDest, we need to distinguish between vecDest
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_DEST*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecDest0(5, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC)*dma_threads_per_ld);

    cudaDMAIndirect<true, true,
                    MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
                             (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
                    CUDADMA_DMA_LDS_DEST*dma_threads_per_ld, 2*compute_threads_per_cta>
    dma_vecDest1(6, compute_threads_per_cta, compute_threads_per_cta+
		 (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC)*dma_threads_per_ld);

    //--------------------------------------------------------------------------

#if COEFFSATEDGE_DEVICE == SOA
    // Coefficients at edges are stored as structure of arrays, that
    // is, we have 2*ncoeff subarrays of length nedge which store:
    //
    // 0-subarray: ij-coefficients for x-direction, 
    // 1-subarray: ji-coefficients for x-direction, 
    // 2-subarray: ij-coefficients for y-direction,
    // 3-subarray: ji-coefficients for y-direction,
    // ...
    // n-subarray: further coefficients not required here

    // Strided cudaDMA thread to transfer precomputed coefficients
    // CoeffsAtEdge into shared memory s_CoeffsAtEdge
    cudaDMAStrided<true, sizeof(Tc),
                   compute_threads_per_cta*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*HYDRO_NDIM>
    dma_coeff0(7, compute_threads_per_cta, compute_threads_per_cta+
    	       (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
    	       nedge*sizeof(Tc));

    // cudaDMAStrided<true, sizeof(Tc),
    //                compute_threads_per_cta*sizeof(Tc),
    //                CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
    //                2*HYDRO_NDIM>
    // dma_coeff1(8, compute_threads_per_cta, compute_threads_per_cta+
    // 	       (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
    // 	       nedge*sizeof(Tc));
#else
    // Coefficients at edges are stored as array of structure, that
    // is, we have nedge real-valued subarray of length 2*ncoeff
    cudaDMAStrided<true, sizeof(Tc), HYDRO_NDIM*sizeof(Tc),
                   CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
                   2*compute_threads_per_cta>
    dma_coeff0(7, compute_threads_per_cta, compute_threads_per_cta+
    	       (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
    	       ncoeff*sizeof(Tc));
  
    // cudaDMAStrided<true, sizeof(Tc), HYDRO_NDIM*sizeof(Tc),
    //                CUDADMA_DMA_LDS_COEFF*dma_threads_per_ld,
    //                2*compute_threads_per_cta>
    // dma_coeff1(8, compute_threads_per_cta, compute_threads_per_cta+
    // 	       (CUDADMA_DMA_LDS_IND+CUDADMA_DMA_LDS_SRC+CUDADMA_DMA_LDS_DEST)*dma_threads_per_ld,
    // 	       ncoeff*sizeof(Tc));
#endif

    //--------------------------------------------------------------------------
    // Warp specialisation
    //--------------------------------------------------------------------------
    if (threadIdx.x<compute_threads_per_cta) {
      
      // Start DMA transfer of indices
      dma_ind0.start_async_dma();
      dma_ind1.start_async_dma();
      dma_ind2.start_async_dma();
      
      // Start DMA transfer of coefficients
      dma_coeff0.start_async_dma();
      //      dma_coeff1.start_async_dma();

      // Wait for indices to be ready
      dma_ind0.wait_for_dma_finish();
    
      // Start DMA transfer of indirect data
      dma_vecSrc0.start_async_dma();
      dma_vecDest0.start_async_dma();
    
      // Wait for indices to be ready
      dma_ind1.wait_for_dma_finish();
      
      // Start DMA transfer of indirect data
      dma_vecSrc1.start_async_dma();
      dma_vecDest1.start_async_dma();

      // Loop over all edge-groups to be processed by this block
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {
	
#define IBUF 0
#define DBUF 0
#define IOFF 0
	
	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	TdDest ui = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vi = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,1,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	TdDest uj = XVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	TdDest vj = YVELOCITY3(s_VecSrc[DBUF],
			       IDX3,2,(int)threadIdx.x+1,
			       NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	TdDest pi = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,1,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	TdDest pj = PRESSURE3(s_VecSrc[DBUF],
			      IDX3,2,(int)threadIdx.x+1,
			      NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	Ti idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	       + nedge_offset + threadIdx.x;
	
	// Local variables
	TdDest FluxAtEdge[2*NVAR2D];
	
	// Wait for coefficients to be ready
	dma_coeff0.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Start DMA transfer of coefficients
	dma_coeff0.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
      
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	Ti j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);

#undef IBUF
#undef DBUF
#undef IOFF
	
	// Start DMA transfer of indices
	dma_ind0.start_async_dma();

	// Wait for indices to be ready
	dma_ind2.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();

#define IBUF 1
#define DBUF 1
#define IOFF 1

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Wait for coefficients to be ready
	//	dma_coeff1.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Start DMA transfer of coefficients
	//	dma_coeff1.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind1.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind0.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();

#define IBUF 2
#define DBUF 0
#define IOFF 2

	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Wait for coefficients to be ready
	dma_coeff0.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Start DMA transfer of coefficients
	dma_coeff0.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind2.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind1.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();

#define IBUF 0
#define DBUF 1
#define IOFF 3

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Wait for coefficients to be ready
	//	dma_coeff1.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Start DMA transfer of coefficients
	//	dma_coeff1.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind0.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind2.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();

#define IBUF 1
#define DBUF 0
#define IOFF 4

	// Wait for source vector to be ready
	dma_vecSrc0.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Wait for coefficients to be ready
	dma_coeff0.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Start DMA transfer of coefficients
	dma_coeff0.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest0.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind1.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind0.wait_for_dma_finish();

	// Start DMA transfer of indirect data
	dma_vecSrc0.start_async_dma();
	dma_vecDest0.start_async_dma();

#define IBUF 2
#define DBUF 1
#define IOFF 5

	// Wait for source vector to be ready
	dma_vecSrc1.wait_for_dma_finish();

	// Compute velocities
	ui = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,1,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_VecSrc[DBUF],
			IDX3,2,(int)threadIdx.x+1,
			NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,1,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_VecSrc[DBUF],
		       IDX3,2,(int)threadIdx.x+1,
		       NVAR2D,2,compute_threads_per_cta);
	
	// Global edge ID
	idx = ((ipt+IOFF)*gridDim.x+blockIdx.x)*compute_threads_per_cta
	  + nedge_offset + threadIdx.x;
	
	// Wait for coefficients to be ready
	//	dma_coeff1.wait_for_dma_finish();

	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,false>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,false,false>
	  (FluxAtEdge,CoeffsAtEdge,s_VecSrc[DBUF],ui,uj,vi,vj,pi,pj,scale,
	   1,(int)threadIdx.x+1,idx+1,nedge,ncoeff);
	
	// Start DMA transfer of coefficients
	//	dma_coeff1.start_async_dma();

	// Wait for destination vector to be ready
	dma_vecDest1.wait_for_dma_finish();
	
	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList[IBUF],1,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,i,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 1, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,1,NVAR2D,2);
	
	// Get positions of edge endpoints (idx starts at zero)
	j = IDX2(s_IedgeList[IBUF],2,(int)threadIdx.x+1,2,compute_threads_per_cta);
	
#pragma unroll
        for (int ivar=1; ivar<=NVAR2D; ++ivar)
	  IDX2_REVERSE(vecDest,ivar,j,NVAR2D,neq) =
	    IDX3(s_VecDest[DBUF], ivar, 2, (int)threadIdx.x+1,
		 NVAR2D, 2, compute_threads_per_cta) + IDX2(FluxAtEdge,ivar,2,NVAR2D,2);
	
#undef IBUF
#undef DBUF
#undef IOFF

	// Start DMA transfer of indices
	dma_ind2.start_async_dma();
	
	// Wait for indices to be ready
	dma_ind1.wait_for_dma_finish();
	
	// Start DMA transfer of indirect data
	dma_vecSrc1.start_async_dma();
	dma_vecDest1.start_async_dma();

      }
    }
    
    //--------------------------------------------------------------------------
    // DMA transfer warps
    //--------------------------------------------------------------------------
    else if(dma_ind0.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=3) {

	dma_ind0.execute_dma(&IedgeList[ (((ipt+0)*gridDim.x+blockIdx.x)*
					  compute_threads_per_cta+nedge_offset)*
					 (EDGELIST_DEVICE == SOA ? 2 : 6)],
			     s_IedgeList[0]);

	dma_ind1.execute_dma(&IedgeList[ (((ipt+1)*gridDim.x+blockIdx.x)*
					  compute_threads_per_cta+nedge_offset)*
					 (EDGELIST_DEVICE == SOA ? 2 : 6)],
			     s_IedgeList[1]);
	
	dma_ind2.execute_dma(&IedgeList[ (((ipt+2)*gridDim.x+blockIdx.x)*
					  compute_threads_per_cta+nedge_offset)*
					 (EDGELIST_DEVICE == SOA ? 2 : 6)],
					 s_IedgeList[2]);
      }
      dma_ind0.finish_async_dma();
      dma_ind1.finish_async_dma();
    }
    
    else if(dma_vecSrc0.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {
	
	dma_vecSrc0.wait_for_dma_start();
	dma_vecSrc0.execute_dma_no_sync(s_IedgeList[0], vecSrc-NVAR2D, s_VecSrc[0]);
	dma_vecSrc0.finish_async_dma();

	dma_vecSrc1.wait_for_dma_start();
	dma_vecSrc1.execute_dma_no_sync(s_IedgeList[1], vecSrc-NVAR2D, s_VecSrc[1]);
	dma_vecSrc1.finish_async_dma();

	dma_vecSrc0.wait_for_dma_start();
	dma_vecSrc0.execute_dma_no_sync(s_IedgeList[2], vecSrc-NVAR2D, s_VecSrc[0]);
	dma_vecSrc0.finish_async_dma();

	dma_vecSrc1.wait_for_dma_start();
	dma_vecSrc1.execute_dma_no_sync(s_IedgeList[0], vecSrc-NVAR2D, s_VecSrc[1]);
	dma_vecSrc1.finish_async_dma();

	dma_vecSrc0.wait_for_dma_start();
	dma_vecSrc0.execute_dma_no_sync(s_IedgeList[1], vecSrc-NVAR2D, s_VecSrc[0]);
	dma_vecSrc0.finish_async_dma();

	dma_vecSrc1.wait_for_dma_start();
	dma_vecSrc1.execute_dma_no_sync(s_IedgeList[2], vecSrc-NVAR2D, s_VecSrc[1]);
	dma_vecSrc1.finish_async_dma();
      }
    }

    else if(dma_vecDest0.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=6) {
	
	dma_vecDest0.wait_for_dma_start();
	dma_vecDest0.execute_dma_no_sync(s_IedgeList[0], vecDest-NVAR2D, s_VecDest[0]);
	dma_vecDest0.finish_async_dma();

	dma_vecDest1.wait_for_dma_start();
	dma_vecDest1.execute_dma_no_sync(s_IedgeList[1], vecDest-NVAR2D, s_VecDest[1]);
	dma_vecDest1.finish_async_dma();

	dma_vecDest0.wait_for_dma_start();
	dma_vecDest0.execute_dma_no_sync(s_IedgeList[2], vecDest-NVAR2D, s_VecDest[0]);
	dma_vecDest0.finish_async_dma();

	dma_vecDest1.wait_for_dma_start();
	dma_vecDest1.execute_dma_no_sync(s_IedgeList[0], vecDest-NVAR2D, s_VecDest[1]);
	dma_vecDest1.finish_async_dma();

	dma_vecDest0.wait_for_dma_start();
	dma_vecDest0.execute_dma_no_sync(s_IedgeList[1], vecDest-NVAR2D, s_VecDest[0]);
	dma_vecDest0.finish_async_dma();

	dma_vecDest1.wait_for_dma_start();
	dma_vecDest1.execute_dma_no_sync(s_IedgeList[2], vecDest-NVAR2D, s_VecDest[1]);
	dma_vecDest1.finish_async_dma();
      }
    }
    
    else if(dma_coeff0.owns_this_thread()) {
      for (int ipt=0; ipt<nedge_per_thread; ipt+=2) {
	
      	dma_coeff0.execute_dma(&CoeffsAtEdge[ ((ipt*gridDim.x+blockIdx.x)*
      					       compute_threads_per_cta+nedge_offset)*
      					      (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
      			       s_CoeffsAtEdge[0]);
	
      	// dma_coeff1.execute_dma(&CoeffsAtEdge[ (((ipt+1)*gridDim.x+blockIdx.x)*
      	// 				       compute_threads_per_cta+nedge_offset)*
      	// 				      (COEFFSATEDGE_DEVICE == SOA ? 1 : 2*ncoeff)],
      	// 		       s_CoeffsAtEdge[1]);
	
      }
    }
  };
  
#undef TOTAL_THREADS_PER_CTA
  
  /*****************************************************************************
   * Internal C++ functions which invoke the CUDA kernels
   ****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int idissipationtype>
  inline
  int hydro_calcFlux2d_cuda(__SIZET *d_CoeffsAtEdge,
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
    const int nedge_per_thread_cudaDMA = CUDADMA_NEDGE_PER_THREAD;

    const int threads_per_cta_baseline  = BASELINE_THREADS_PER_CTA;
    const int nedge_per_thread_baseline = BASELINE_NEDGE_PER_THREAD;
    
    int blocks, threads, nedge_cudaDMA, nedge_baseline;
    prepare_cudaDMA(devProp, nedgeset,
		    nedge_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &nedge_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(devProp, nedgeset-nedge_cudaDMA,
		     nedge_per_thread_baseline, threads_per_cta_baseline,
		     &blocks, &threads, &nedge_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);

    TdSrc  *vecSrc = (TdSrc*)(*d_vecSrc);
    TdDest *vecDest = (TdDest*)(*d_vecDest);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
    
    cout << "hydro_calcFlux2d_cuda" 
	 << " nblocks=" << nblocks
	 << " neq=" << neq
	 << " nedgeset=" << nedgeset << endl;
    cout << "Memory NEDGE: " << NVAR2D*nedgeset*sizeof(TdSrc)/1000000.0f
	 << " MB" << " NEDGE=" << nedgeset << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    if (nblocks == 1) {
#ifdef CUDADMA_KERNEL    
      if (grid_cudaDMA.x>0) {
	cudaEventRecord(start,stream);

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
	cudaEventRecord(stop,stream);
	cudaEventSynchronize(stop);
	
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	
	cout << " CudaDMA:"
	     << " #blocks=" << grid_cudaDMA.x << ","
	     << grid_cudaDMA.y << "," << grid_cudaDMA.z 
	     << " #threads per block=" << block_cudaDMA.x 
	     << "," << block_cudaDMA.y << "," << block_cudaDMA.z << endl;
	cout << "Elapsed time: " << elapsedTime << " ms" << endl;
	cout << "Bandwidth:    " << (2*nedge_cudaDMA*sizeof(Ti)+            // get i,j
				     2*NVAR2D*nedge_cudaDMA*sizeof(TdSrc)+  // gather solution
				     4*NVAR2D*nedge_cudaDMA*sizeof(TdDest)+ // gather and scatter vector
				     2*HYDRO_NDIM*nedge_cudaDMA*sizeof(Tc)) // gather coefficients
	  /1000000000.0f/elapsedTime*1000.0f
	     << " GB/s" << endl;
      }
#endif

#ifdef BASELINE_KERNEL
      if (grid_baseline.x>0) {
	cudaEventRecord(start,stream);
	
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
	cudaEventRecord(stop,stream);
	cudaEventSynchronize(stop);
	
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	
	cout << " Baseline:"
	     << " #blocks=" << grid_baseline.x << "," 
	     << grid_baseline.y << "," << grid_baseline.z 
	     << " #threads per block=" << block_baseline.x 
	     << "," << block_baseline.y << "," << block_baseline.z << endl;
	cout << "Elapsed time: " << elapsedTime << " ms" << endl;
	cout << "Bandwidth:    " << (2*nedge_baseline*sizeof(Ti)+            // get i,j
				     2*NVAR2D*nedge_baseline*sizeof(TdSrc)+  // gather solution
				     4*NVAR2D*nedge_baseline*sizeof(TdDest)+ // gather and scatter vector
				     2*HYDRO_NDIM*nedge_baseline*sizeof(Tc)) // gather coefficients
	  /1000000000.0f/elapsedTime*1000.0f
	     << " GB/s" << endl;
      }
#endif
    } else {
#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0) {
	cudaEventRecord(start,stream);
	
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
      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
      
      float elapsedTime;
      cudaEventElapsedTime(&elapsedTime, start, stop);
      
      cout << " CudaDMA:"
	   << " #blocks=" << grid_cudaDMA.x << ","
	   << grid_cudaDMA.y << "," << grid_cudaDMA.z 
	   << " #threads per block=" << block_cudaDMA.x 
	   << "," << block_cudaDMA.y << "," << block_cudaDMA.z << endl;
      cout << "Elapsed time: " << elapsedTime << " ms" << endl;
      cout << "Bandwidth:    " << (2*nedge_cudaDMA*sizeof(Ti)+            // get i,j
				   2*NVAR2D*nedge_cudaDMA*sizeof(TdSrc)+  // gather solution
				   4*NVAR2D*nedge_cudaDMA*sizeof(TdDest)+ // gather and scatter vector
				   2*HYDRO_NDIM*nedge_cudaDMA*sizeof(Tc)) // gather coefficients
	/1000000000.0f/elapsedTime*1000.0f
	   << " GB/s" << endl;
      }
#endif

#ifdef BASELINE_KERNEL
      if (grid_baseline.x>0) {
	cudaEventRecord(start,stream);
	
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
	cudaEventRecord(stop,stream);
	cudaEventSynchronize(stop);

	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	
	cout << " Baseline:"
	     << " #blocks=" << grid_baseline.x << "," 
	     << grid_baseline.y << "," << grid_baseline.z 
	     << " #threads per block=" << block_baseline.x 
	     << "," << block_baseline.y << "," << block_baseline.z << endl;
	cout << "Elapsed time: " << elapsedTime << " ms" << endl;
	cout << "Bandwidth:    " << (2*nedge_baseline*sizeof(Ti)+            // get i,j
				     2*NVAR2D*nedge_baseline*sizeof(TdSrc)+  // gather solution
				     4*NVAR2D*nedge_baseline*sizeof(TdDest)+ // gather and scatter vector
				     2*HYDRO_NDIM*nedge_baseline*sizeof(Tc)) // gather coefficients
	  /1000000000.0f/elapsedTime*1000.0f
	     << " GB/s" << endl;
      }
#endif
    }
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    coproc_checkErrors("hydro_calcFlux2d_cuda");

    return 0;
  };
  
  /*****************************************************************************
   * External C functions which can be called from the Fortran code
   ****************************************************************************/

  extern "C"
  { 
    __INT FNAME(hydro_calcfluxgalerkin2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
      return (__INT) hydro_calcFlux2d_cuda
	<__DP,__DP,__DP,__INT,DISSIPATION_ZERO>
	(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
	 *scale, *nblocks, *neq, *nedge,
	 *ncoeff, *nedges, *iedgeset, 
	 (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcfluxscdiss2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
      return (__INT) hydro_calcFlux2d_cuda
	<__DP,__DP,__DP,__INT,DISSIPATION_SCALAR>
	(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
	 *scale, *nblocks, *neq, *nedge,
	 *ncoeff, *nedges, *iedgeset, 
	 (cudaStream_t)(*stream));
    }

    /**************************************************************************/

    __INT FNAME(hydro_calcfluxscdissdisp2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
      return (__INT) hydro_calcFlux2d_cuda
	<__DP,__DP,__DP,__INT,DISSIPATION_SCALAR_DSPLIT>
	(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
	 *scale, *nblocks, *neq, *nedge,
	 *ncoeff, *nedges, *iedgeset,
	 (cudaStream_t)(*stream));
    }

    /**************************************************************************/

    __INT FNAME(hydro_calcfluxroediss2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
      return (__INT) hydro_calcFlux2d_cuda
	<__DP,__DP,__DP,__INT,DISSIPATION_ROE>
	(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
	 *scale, *nblocks, *neq, *nedge,
	 *ncoeff, *nedges, *iedgeset,
	 (cudaStream_t)(*stream));
    }

  /***************************************************************************/

    __INT FNAME(hydro_calcfluxroedissdisp2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
      return (__INT) hydro_calcFlux2d_cuda
	<__DP,__DP,__DP,__INT,DISSIPATION_ROE_DSPLIT>
	(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
	 *scale, *nblocks, *neq, *nedge,
	 *ncoeff, *nedges, *iedgeset,
	 (cudaStream_t)*stream);
    }

    /**************************************************************************/

    __INT FNAME(hydro_calcfluxrusdiss2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
      return (__INT)hydro_calcFlux2d_cuda
	<__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV>
	(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
	 *scale, *nblocks, *neq, *nedge,
	 *ncoeff, *nedges, *iedgeset,
	 (cudaStream_t)*stream);
    }
    
    /**************************************************************************/

    __INT FNAME(hydro_calcfluxrusdissdisp2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
      return (__INT) hydro_calcFlux2d_cuda
	<__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV_DSPLIT>
	(d_CoeffsAtEdge, d_IedgeList, d_vecSrc, d_vecDest,
	 *scale, *nblocks, *neq, *nedge,
	 *ncoeff, *nedges, *iedgeset,
	 (cudaStream_t)*stream);
    }
  };
}
