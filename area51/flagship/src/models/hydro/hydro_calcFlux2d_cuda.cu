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
      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      Td anorm = sqrt(a[0] * a[0] + a[1] * a[1]);
    
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
      Td c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON));
    
      // Compute scalar dissipation
      Td d_ij = abs(vel_ij) + anorm*c_ij;
    
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
      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    
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
    // Local variables
    TdDest DataAtEdge[2*NVAR2D];
    TdDest FluxAtEdge[2*NVAR2D];
    TdDest ui,uj,vi,vj,pi,pj;
    Ti  i,j;

    // Loop over all items per thread
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {
      
      // Global edge ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*threads_per_cta+nedge_offset+threadIdx.x;
      
      if (threadIdx.x<threads_per_cta && idx<nedge_last)
	{
	  // Get positions of edge endpoints (idx starts at zero)
	  i = IDX2_EDGELIST(IedgeList,1,idx+1,6,nedge);
	  j = IDX2_EDGELIST(IedgeList,2,idx+1,6,nedge);
	  
	  // Get solution values at edge endpoints
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    gatherEdgeData<true>(DataAtEdge,vecSrc,i,j,neq);
	  
	  // Compute velocities
	  ui = XVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
	  vi = YVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
	  
	  uj = XVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);
	  vj = YVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);
	  
	  // Compute pressures
	  pi = PRESSURE2(DataAtEdge,IDX2,1,NVAR2D,2);
	  pj = PRESSURE2(DataAtEdge,IDX2,2,NVAR2D,2);
	  
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
#define SHMEM_DATA_TRANSPOSE false
#define SHMEM_DATA_IDX3      IDX3

    // Local memory
    const int tid = threadIdx.x;
    
    // Local variables
    TdDest FluxAtEdge[2*NVAR2D];
    TdDest ui,uj,vi,vj,pi,pj;
    Ti i,j;

    // Shared memory
    __shared__ TdSrc s_DataAtEdge[2*NVAR2D*threads_per_cta];
    
    // Loop over all items per thread
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {
      
      // Global edge ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*threads_per_cta+nedge_offset+threadIdx.x;
      
      if (tid<threads_per_cta && idx<nedge_last)
	{
	  // Get positions of edge endpoints (idx starts at zero)
	  i = IDX2_EDGELIST(IedgeList,1,idx+1,6,nedge);
	  j = IDX2_EDGELIST(IedgeList,2,idx+1,6,nedge);
	  
	  // Get solution values at edge endpoints
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    gatherEdgeData<threads_per_cta,SHMEM_DATA_TRANSPOSE,true>
	    (s_DataAtEdge,vecSrc,tid+1,i,j,neq);
	  
	  // Compute velocities
	  ui = XVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,1,tid+1,NVAR2D,2,threads_per_cta);
	  vi = YVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,1,tid+1,NVAR2D,2,threads_per_cta);
	  
	  uj = XVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,2,tid+1,NVAR2D,2,threads_per_cta);
	  vj = YVELOCITY3(s_DataAtEdge,SHMEM_DATA_IDX3,2,tid+1,NVAR2D,2,threads_per_cta);
	  
	  // Compute pressures
	  pi = PRESSURE3(s_DataAtEdge,SHMEM_DATA_IDX3,1,tid+1,NVAR2D,2,threads_per_cta);
	  pj = PRESSURE3(s_DataAtEdge,SHMEM_DATA_IDX3,2,tid+1,NVAR2D,2,threads_per_cta);
	  
	  // Compute the artificial viscosities
	  InviscidFluxDissipation<idissipationtype>::
	    calcEdgeData<1,threads_per_cta,false,SHMEM_DATA_TRANSPOSE>
	    (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,pi,pj,
	     1,tid+1,idx+1,nedge,ncoeff);
	  Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	  
	  // Compute inviscid fluxes
	  InviscidFlux::calcEdgeData<1,threads_per_cta,false,SHMEM_DATA_TRANSPOSE,false>
	    (FluxAtEdge,CoeffsAtEdge,s_DataAtEdge,ui,uj,vi,vj,pi,pj,
	     scale,1,tid+1,idx+1,nedge,ncoeff);
	  
	  // Build fluxes into nodal vector
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);
	}
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes and applies
   * artificial dissipation if required) (cudaDMA implementation).
   ****************************************************************************/
  
  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipationtype,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __global__ void hydro_calcFlux2d_cudaDMA(Tc *CoeffsAtEdge,
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
    // Local memory
    const int tid = threadIdx.x;
    const int total_threads_per_cta = compute_threads_per_cta+2*dma_threads_per_ld;

    // Shared memory
    __shared__ Ti s_IedgeList0[2*compute_threads_per_cta];
    __shared__ Ti s_IedgeList1[2*compute_threads_per_cta];
    __shared__ Ti s_IedgeList2[2*compute_threads_per_cta];
    __shared__ Ti s_IedgeList3[2*compute_threads_per_cta];
    // __shared__ Tc s_CoeffsAtEdge0[2*HYDRO_NDIM*compute_threads_per_cta];
    // __shared__ Tc s_CoeffsAtEdge1[2*HYDRO_NDIM*compute_threads_per_cta];
    __shared__ TdSrc s_DataAtEdge0[2*NVAR2D*compute_threads_per_cta];
    __shared__ TdSrc s_DataAtEdge1[2*NVAR2D*compute_threads_per_cta];
    
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
    // array IedgeList into shared memory s_IedgeList0 and s_IedgeList1
    cudaDMASequential<false, 2*sizeof(Ti), 2*sizeof(Ti)*compute_threads_per_cta,
      total_threads_per_cta>dma_ind;
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
    // array IedgeList into shared memory s_IedgeList0 and s_IedgeList1
    
    cudaDMAStrided<false, 2*sizeof(Ti), 2*sizeof(Ti), 
      total_threads_per_cta, compute_threads_per_cta>dma_ind(6*sizeof(Ti));
#endif

    //--------------------------------------------------------------------------

    // Indirect cudaDMA thread to transfer nodal data from vecSrc into
    // shared memory s_DataAtEdge0, we need to distinguish between vecSrc
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
      MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
               (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
      dma_threads_per_ld, compute_threads_per_cta*2>
      dma_vec0(0, compute_threads_per_cta, compute_threads_per_cta);

    // Indirect cudaDMA thread to transfer nodal data from vecSrc into
    // shared memory s_DataAtEdge1, we need to distinguish between vecSrc
    // stored in interleaved format and vecSrc stored in block format
    cudaDMAIndirect<true, true,
      MAXALIGN((isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc)),
               (isystemformat==SYSTEM_BLOCK ? 1 : NVAR2D)*sizeof(TdSrc),
      dma_threads_per_ld, compute_threads_per_cta*2>
      dma_vec1(1, compute_threads_per_cta, compute_threads_per_cta+dma_threads_per_ld);

    //--------------------------------------------------------------------------

#if COEFFSATEDGE_DEVICE == SOA
    // Coefficients at edges are stored as structure of arrays, that
    // is, we have ncoeff subarrays of length nedge which store:
    //
    // 0-subarray: coefficients for x-direction, 
    // 1-subarray: coefficients for y-direction,
    // n-subarray: further coefficients not required here

    // Strided cudaDMA thread to transfer precomputed coefficients
    // CoeffsAtEdge into shared memory s_CoeffsAtEdge0 and s_CoeffsAtEdge1
    // cudaDMAStrided<true, sizeof(Tc), compute_threads_per_cta*sizeof(Tc),
    //                dma_threads_per_ld, HYDRO_NDIM>
    //   dma_coeff0(1, compute_threads_per_cta,
    // 		 compute_threads_per_cta+1*dma_threads_per_ld,
    // 		 nedge*sizeof(Tc));
    /*    cudaDMAStrided<true, sizeof(Tc), compute_threads_per_cta*sizeof(Tc),
                   dma_threads_per_ld, HYDRO_NDIM>
      dma_coeff1(3, compute_threads_per_cta,
    		 compute_threads_per_cta+3*dma_threads_per_ld,
    		 nedge*sizeof(Tc));*/
#else
    // Coefficients at edges are stored as array of structure, that
    // is, we have nedge real-valued subarray of length ncoeff
    // cudaDMAStrided<true, sizeof(Tc), sizeof(Tc), dma_threads_per_ld, compute_threads_per_cta>
    //   dma_coeff0(1, compute_threads_per_cta,
    // 		 compute_threads_per_cta+1*dma_threads_per_ld,
    // 		 ncoeff*sizeof(Tc), HYDRO_NDIM*sizeof(Tc));
    /*   cudaDMAStrided<true, sizeof(Tc), sizeof(Tc), dma_threads_per_ld, compute_threads_per_cta>
      dma_coeff1(3, compute_threads_per_cta,
    		 compute_threads_per_cta+3*dma_threads_per_ld,
    		 ncoeff*sizeof(Tc), HYDRO_NDIM*sizeof(Tc));*/
#endif

    //--------------------------------------------------------------------------

    // Loop over all edge-groups to be processed by this block
    for (int ipt=0; ipt<nedge_per_thread; ipt+=4) {

      //------------------------------------------------------------------------
      // Load the indices with all threads - no warp specialisation
      //------------------------------------------------------------------------
      ptx_cudaDMA_barrier_blocking(9, total_threads_per_cta);
      
      // Buffer0 - part A
      dma_ind.execute_dma(&IedgeList[ (((ipt+0)*gridDim.x+blockIdx.x)*
				       compute_threads_per_cta+nedge_offset)*
				      (EDGELIST_DEVICE == SOA ? 2 : 6)], s_IedgeList0);
      // Buffer1 - part A
      dma_ind.execute_dma(&IedgeList[ (((ipt+1)*gridDim.x+blockIdx.x)*
      				       compute_threads_per_cta+nedge_offset)*
      				      (EDGELIST_DEVICE == SOA ? 2 : 6)], s_IedgeList1);

      // Buffer0 - part B
      dma_ind.execute_dma(&IedgeList[ (((ipt+2)*gridDim.x+blockIdx.x)*
				       compute_threads_per_cta+nedge_offset)*
				      (EDGELIST_DEVICE == SOA ? 2 : 6)], s_IedgeList2);

      // Buffer1 - part B
      dma_ind.execute_dma(&IedgeList[ (((ipt+3)*gridDim.x+blockIdx.x)*
      				       compute_threads_per_cta+nedge_offset)*
      				      (EDGELIST_DEVICE == SOA ? 2 : 6)], s_IedgeList3);
      
      ptx_cudaDMA_barrier_blocking(9, total_threads_per_cta);

      //------------------------------------------------------------------------
      // Start warp specialisation
      //------------------------------------------------------------------------
      if (tid<compute_threads_per_cta) {

	// Local variables
	TdDest FluxAtEdge[2*NVAR2D];
	TdDest ui,uj,vi,vj,pi,pj;
	Ti idx,i,j;
	
	// Get solution values at edge endpoints
	dma_vec0.start_async_dma();
	dma_vec1.start_async_dma();
	
	// Get precomputed coefficients at edges
	//	dma_coeff0.start_async_dma();
	//	dma_coeff1.start_async_dma();

	//----------------------------------------------------------------------
	// Buffer0 - part A
	//----------------------------------------------------------------------

	// Wait for solution values to be available
	dma_vec0.wait_for_dma_finish();

	// Global edge ID
	idx = (ipt*gridDim.x+blockIdx.x)*compute_threads_per_cta+nedge_offset+threadIdx.x;

	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList0,1,tid+1,2,compute_threads_per_cta);
	j = IDX2(s_IedgeList0,2,tid+1,2,compute_threads_per_cta);

	// Compute velocities
	ui = XVELOCITY3(s_DataAtEdge0,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_DataAtEdge0,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_DataAtEdge0,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_DataAtEdge0,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_DataAtEdge0,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_DataAtEdge0,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Wait for precomputed coefficients to be available
	//	dma_coeff0.wait_for_dma_finish();
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_DataAtEdge0,ui,uj,vi,vj,pi,pj,
	   1,tid+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK,false>
	  (FluxAtEdge,CoeffsAtEdge,s_DataAtEdge0,ui,uj,vi,vj,pi,pj,
	   scale,1,tid+1,idx+1,nedge,ncoeff);

	// Get solution values at edge endpoints for part B
	dma_vec0.start_async_dma();
	
	// Build fluxes into nodal vector
	Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);

	//----------------------------------------------------------------------
	// Buffer1 - part A
	//----------------------------------------------------------------------

	// Wait for solution values to be available
	dma_vec1.wait_for_dma_finish();

	// Global edge ID
	idx = ((ipt+1)*gridDim.x+blockIdx.x)*compute_threads_per_cta+nedge_offset+threadIdx.x;

	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList1,1,tid+1,2,compute_threads_per_cta);
	j = IDX2(s_IedgeList1,2,tid+1,2,compute_threads_per_cta);
	
	// Compute velocities
	ui = XVELOCITY3(s_DataAtEdge1,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_DataAtEdge1,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_DataAtEdge1,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_DataAtEdge1,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_DataAtEdge1,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_DataAtEdge1,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Wait for precomputed coefficients to be available
	//	dma_coeff0.wait_for_dma_finish();
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_DataAtEdge1,ui,uj,vi,vj,pi,pj,
	   1,tid+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK,false>
	  (FluxAtEdge,CoeffsAtEdge,s_DataAtEdge1,ui,uj,vi,vj,pi,pj,
	   scale,1,tid+1,idx+1,nedge,ncoeff);

	// Get solution values at edge endpoints for part B
	dma_vec1.start_async_dma();

	// Build fluxes into nodal vector
	Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);

	//----------------------------------------------------------------------
	// Buffer0 - part B
	//----------------------------------------------------------------------

	// Wait for solution values to be available
	dma_vec0.wait_for_dma_finish();

	// Global edge ID
	idx = ((ipt+2)*gridDim.x+blockIdx.x)*compute_threads_per_cta+nedge_offset+threadIdx.x;

	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList2,1,tid+1,2,compute_threads_per_cta);
	j = IDX2(s_IedgeList2,2,tid+1,2,compute_threads_per_cta);

	// Compute velocities
	ui = XVELOCITY3(s_DataAtEdge0,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_DataAtEdge0,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_DataAtEdge0,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_DataAtEdge0,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_DataAtEdge0,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_DataAtEdge0,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Wait for precomputed coefficients to be available
	//	dma_coeff0.wait_for_dma_finish();
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_DataAtEdge0,ui,uj,vi,vj,pi,pj,
	   1,tid+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK,false>
	  (FluxAtEdge,CoeffsAtEdge,s_DataAtEdge0,ui,uj,vi,vj,pi,pj,
	   scale,1,tid+1,idx+1,nedge,ncoeff);
	
	// Build fluxes into nodal vector
	Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);

	//----------------------------------------------------------------------
	// Buffer1 - part B
	//----------------------------------------------------------------------

	// Wait for solution values to be available
	dma_vec1.wait_for_dma_finish();

	// Global edge ID
	idx = ((ipt+3)*gridDim.x+blockIdx.x)*compute_threads_per_cta+nedge_offset+threadIdx.x;

	// Get positions of edge endpoints (idx starts at zero)
	i = IDX2(s_IedgeList3,1,tid+1,2,compute_threads_per_cta);
	j = IDX2(s_IedgeList3,2,tid+1,2,compute_threads_per_cta);
	
	// Compute velocities
	ui = XVELOCITY3(s_DataAtEdge1,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	vi = YVELOCITY3(s_DataAtEdge1,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	uj = XVELOCITY3(s_DataAtEdge1,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	vj = YVELOCITY3(s_DataAtEdge1,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Compute pressures
	pi = PRESSURE3(s_DataAtEdge1,IDX3,1,tid+1,NVAR2D,2,compute_threads_per_cta);
	pj = PRESSURE3(s_DataAtEdge1,IDX3,2,tid+1,NVAR2D,2,compute_threads_per_cta);
	
	// Wait for precomputed coefficients to be available
	//	dma_coeff0.wait_for_dma_finish();
	
	// Compute the artificial viscosities
	InviscidFluxDissipation<idissipationtype>::
	  calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK>
	  (&FluxAtEdge[NVAR2D],CoeffsAtEdge,s_DataAtEdge1,ui,uj,vi,vj,pi,pj,
	   1,tid+1,idx+1,nedge,ncoeff);
	Flux::combineEdgeData<true>(FluxAtEdge,&FluxAtEdge[NVAR2D],scale);
	
	// Compute inviscid fluxes
	InviscidFlux::calcEdgeData<1,compute_threads_per_cta,false,isystemformat==SYSTEM_BLOCK,false>
	  (FluxAtEdge,CoeffsAtEdge,s_DataAtEdge1,ui,uj,vi,vj,pi,pj,
	   scale,1,tid+1,idx+1,nedge,ncoeff);

	// Build fluxes into nodal vector
	Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  scatterEdgeData<false>(vecDest,FluxAtEdge,i,j,neq);
      }

      //------------------------------------------------------------------------
      // DMA warps
      //------------------------------------------------------------------------

      else if(dma_vec0.owns_this_thread()) {
      	// Indirect cudaDMA transfer of global vector into s_DataAtEdge0

      	if (isystemformat==SYSTEM_BLOCK) {
	  // Buffer 0 - part A
	  dma_vec0.wait_for_dma_start();
	  for (int ivar=0; ivar<NVAR2D; ++ivar)
	    dma_vec0.execute_dma_no_sync(s_IedgeList0, &vecSrc[ivar*neq]-1,
					 &s_DataAtEdge0[ivar*compute_threads_per_cta*2]);
	  dma_vec0.finish_async_dma();

	  // Buffer 0 - part B
	  dma_vec0.wait_for_dma_start();
      	  // Transfer each block separately (index array is 1-based)
	  for (int ivar=0; ivar<NVAR2D; ++ivar)
	    dma_vec0.execute_dma_no_sync(s_IedgeList2, &vecSrc[ivar*neq]-1,
					 &s_DataAtEdge0[ivar*compute_threads_per_cta*2]);
	  dma_vec0.finish_async_dma();
      	}
      	else {
	  // Buffer 0 - part A
      	  dma_vec0.execute_dma(s_IedgeList0, vecSrc-NVAR2D, s_DataAtEdge0);

	  // Buffer 0 - part B
      	  dma_vec0.execute_dma(s_IedgeList2, vecSrc-NVAR2D, s_DataAtEdge0);
      	}
      }
      
      else if(dma_vec1.owns_this_thread()) {
      	// Indirect cudaDMA transfer of global vector into s_DataAtEdge1

      	if (isystemformat==SYSTEM_BLOCK) {
	  // Buffer 1 - part A
      	  dma_vec1.wait_for_dma_start();
      	  for (int ivar=0; ivar<NVAR2D; ++ivar)
      	    dma_vec1.execute_dma_no_sync(s_IedgeList1, &vecSrc[ivar*neq]-1,
      					 &s_DataAtEdge1[ivar*compute_threads_per_cta*2]);
      	  dma_vec1.finish_async_dma();

	  // Buffer 1 - part B
      	  dma_vec1.wait_for_dma_start();
      	  for (int ivar=0; ivar<NVAR2D; ++ivar)
      	    dma_vec1.execute_dma_no_sync(s_IedgeList3, &vecSrc[ivar*neq]-1,
      					 &s_DataAtEdge1[ivar*compute_threads_per_cta*2]);
      	  dma_vec1.finish_async_dma();
      	}
      	else {
	  // Buffer 0 - part A
      	  dma_vec1.execute_dma(s_IedgeList1, vecSrc-NVAR2D, s_DataAtEdge1);

	  // Buffer 0 - part B
      	  dma_vec1.execute_dma(s_IedgeList3, vecSrc-NVAR2D, s_DataAtEdge1);
      	}
      }

     //  else if(dma_coeff0.owns_this_thread()) {
// #if COEFFSATEDGE_DEVICE == SOA
// 	// Strided cudaDMA transfer of precomputed coefficients into s_CoeffsAtEdge0
// 	dma_coeff0.execute_dma(&CoeffsAtEdge[ (((ipt+0)*gridDim.x+blockIdx.x)*
// 					       compute_threads_per_cta+nedge_offset)],
// 			       s_CoeffsAtEdge0);
// #else
// 	// Strided cudaDMA transfer of precomputed coefficients into s_CoeffsAtEdge0
// 	for (int idim=0; idim<HYDRO_NDIM; ++idim)
// 	  dma_coeff0.execute_dma(&CoeffsAtEdge[ (((ipt+0)*gridDim.x+blockIdx.x)*
// 						 compute_threads_per_cta+nedge_offset)+idim],
// 				 &s_CoeffsAtEdge0[idim]);
// #endif
//       }
      /*
      else if(dma_coeff1.owns_this_thread()) {
#if COEFFSATEDGE_DEVICE == SOA
	// Strided cudaDMA transfer of precomputed coefficients into s_CoeffsAtEdge1
	dma_coeff1.execute_dma(&CoeffsAtEdge[ (((ipt+1)*gridDim.x+blockIdx.x)*
					       compute_threads_per_cta+nedge_offset)],
			       s_CoeffsAtEdge1);
#else
	// Strided cudaDMA transfer of precomputed coefficients into s_CoeffsAtEdge1
	for (int idim=0; idim<HYDRO_NDIM; ++idim)
	  dma_coeff1.execute_dma(&CoeffsAtEdge[ (((ipt+1)*gridDim.x+blockIdx.x)*
						 compute_threads_per_cta+nedge_offset)+idim],
				 &s_CoeffsAtEdge1[idim]);
#endif
}*/
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
    const int compute_threads_per_cta  = 32*0;
    const int dma_threads_per_ld       = 32*1;
    const int dma_lds                  = 2;
    const int nedge_per_thread_cudaDMA = 4;

    const int threads_per_cta_baseline  = 32*2;
    const int nedge_per_thread_baseline = 1;
    
    int blocks, threads, nedge_cudaDMA, nedge_baseline;
    prepare_cudaDMA(devProp, nedgeset, nedge_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &nedge_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(devProp, nedgeset-nedge_cudaDMA, nedge_per_thread_baseline,
		     threads_per_cta_baseline, &blocks, &threads, &nedge_baseline);
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
    cout << " CudaDMA:"
	 << " #blocks=" << grid_cudaDMA.x << ","
	 << grid_cudaDMA.y << "," << grid_cudaDMA.z 
	 << " #threads per block=" << block_cudaDMA.x 
	 << "," << block_cudaDMA.y << "," << block_cudaDMA.z << endl;
    cout << " Baseline:"
	 << " #blocks=" << grid_baseline.x << "," 
	 << grid_baseline.y << "," << grid_baseline.z 
	 << " #threads per block=" << block_baseline.x 
	 << "," << block_baseline.y << "," << block_baseline.z << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,stream);
    
    if (nblocks == 1) {
      if (grid_cudaDMA.x>0)
      	// CudaDMA implementation
      	hydro_calcFlux2d_cudaDMA
      	  <Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,idissipationtype,
	MAX(32,compute_threads_per_cta),MAX(32,dma_threads_per_ld)>
      	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
      						       IedgeList,
      						       vecSrc, vecDest, scale,
      						       neq, nedge, ncoeff,
      						       nedge_cudaDMA+iedgeset-1, 
      						       nedge_per_thread_cudaDMA,
      						       iedgeset-1);
      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcFlux2d_baseline
	  <Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,idissipationtype,
	   threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vecSrc, vecDest, scale,
							 neq, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
    } else {
      if (grid_cudaDMA.x>0)
      	// CudaDMA implementation
      	hydro_calcFlux2d_cudaDMA
      	  <Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,idissipationtype,
   	   MAX(32,compute_threads_per_cta),MAX(32,dma_threads_per_ld)>
      	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
      						       IedgeList,
      						       vecSrc, vecDest, scale, 
      						       neq, nedge, ncoeff,
      						       nedge_cudaDMA+iedgeset-1, 
      						       nedge_per_thread_cudaDMA,
      						       iedgeset-1);
      if (grid_baseline.x>0)
      	// Baseline implementation
      	hydro_calcFlux2d_baseline
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
    
    cudaEventRecord(stop,stream);
    cudaEventSynchronize(stop);

    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "Memory NEDGE: " << NVAR2D*nedgeset*sizeof(TdSrc)/1000000.0f
	 << " MB" << " NEDGE=" << nedgeset << endl;
    cout << "Elapsed time: " << elapsedTime << " ms" << endl;
    cout << "Bandwidth:    " << (2*nedgeset*sizeof(Ti)+            // get i,j
				 2*NVAR2D*nedgeset*sizeof(TdSrc)+  // gather solution
				 4*NVAR2D*nedgeset*sizeof(TdDest)+ // gather and scatter vector
				 2*HYDRO_NDIM*nedgeset*sizeof(Tc)) // gather coefficients
      /1000000000.0f/elapsedTime*1000.0f
	 << " GB/s" << endl;
    
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
