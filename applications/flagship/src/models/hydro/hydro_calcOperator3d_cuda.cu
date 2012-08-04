/*#############################################################################
******************************************************************************
* <name> hydro_calcOperator3d_cuda </name>
******************************************************************************
*
* <purpose>
* This file provides CUDA kernels to compute the operator for the low-order
* scheme in 3D using different types if artificial viscosities.
* </purpose>
*
*#############################################################################/
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

#define HYDRO_NDIM 3
#include "hydro.h"

// Defines for baseline implementation
#define BASELINE_THREADS_PER_CTA  32*2
#define BASELINE_NEQ_PER_THREAD   1
#define BASELINE_NEDGE_PER_THREAD 1

// Defines for empty cudaDMA implementation
#ifndef CUDADMA_KERNEL
#define CUDADMA_COMPUTE_THREADS_PER_CTA 0
#define CUDADMA_THREADS_PER_LD          0
#define CUDADMA_NEQ_PER_THREAD          0
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
   * CUDA kernels for hydrodynamic model in 3D
   ****************************************************************************/
  
  // Memory pool in constant device memory
  __device__ __constant__ __SIZET constMemPool[NVAR3D*NVAR3D];
  
  /*****************************************************************************
   * InviscidFluxJacobiMatrixBase (basic functionality and specialisations)
   ****************************************************************************/

  template <int isystemcoupling>
  struct InviscidFluxJacobiMatrixBase
  {
  };

  /*****************************************************************************
   * InviscidFluxJacobiMatrixBase: Specialization for block-diagonal matrix
   ****************************************************************************/

  template<>
  struct InviscidFluxJacobiMatrixBase<SYSTEM_SEGREGATED>
  {
    /***************************************************************************
     * Compute flux Jacobian matrix for neqsim nodes
     **************************************************************************/
    template <int neqsim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtDiag,
			     Tc *CoeffsAtDiag,
			     Td scale,
			     Td ui,
			     Td vi,
			     Td wi,
			     Ti ipos,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,2,ipos,NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,3,ipos,NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,4,ipos,NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,5,ipos,NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
#else
      // Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,2,ipos,NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,3,ipos,NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,4,ipos,NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
      IDX2(MatrixAtDiag,5,ipos,NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,_);
#endif
    }

    /***************************************************************************
     * Compute flux Jacobian matrix for nedgesim edges
     **************************************************************************/
    template <int nedgesim, bool bstabilise, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Ti ipos,
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,5,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);

      // Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,5,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
#else
      // Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      IDX3(MatrixAtEdge,5,1,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,_);
      
      // Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
      IDX3(MatrixAtEdge,5,2,ipos,NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,_);
#endif
    }
  };

  /*****************************************************************************
   * InviscidFluxJacobiMatrixBase: Specialization for full matrix
   ****************************************************************************/

  template <>
  struct InviscidFluxJacobiMatrixBase<SYSTEM_ALLCOUPLED>
  {
    /***************************************************************************
     * Compute flux Jacobian matrix for neqsim nodes
     **************************************************************************/
    template <int neqsim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtDiag,
			     Tc *CoeffsAtDiag,
			     Td scale,
			     Td ui,
			     Td vi,
			     Td wi,
			     Td Ei,
			     Ti ipos,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,2,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX21(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,3,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX31(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,4,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX41(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,5,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX51(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);

      IDX2(MatrixAtDiag,6,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX12(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,7,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,8,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX32(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,9,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX42(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,10,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX52(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,11,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX13(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,12,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX23(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,13,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,14,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX43(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,15,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX53(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,16,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX14(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,17,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX24(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,18,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX34(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,19,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,20,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX54(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,21,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX15(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,22,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX25(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,23,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX35(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,24,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX45(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,25,ipos,NVAR3D*NVAR3D,neqsim) =
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
#else
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,2,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX21(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,3,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX31(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,4,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX41(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,5,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX51(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,6,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX12(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,7,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,8,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX32(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,9,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX42(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,10,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX52(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,11,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX13(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,12,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX23(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,13,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,14,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX43(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,15,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX53(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,16,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX14(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,17,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX24(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,18,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX34(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,19,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,20,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX54(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,21,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX15(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,22,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX25(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,23,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX35(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,24,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX45(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
      IDX2(MatrixAtDiag,25,ipos,NVAR3D*NVAR3D,neqsim) = -
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),
				   IDX2_COEFFSATDIAG(CoeffsAtDiag,3,ieq,ncoeff,neq),ui,vi,wi,Ei);
#endif
    }

    /***************************************************************************
     * Compute flux Jacobian matrix for nedgesim edges
     **************************************************************************/
    template <int nedgesim, bool bstabilise, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Td Ei,
			     Td Ej,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX21(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX31(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX41(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,5,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX51(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,6,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX12(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,7,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,8,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX32(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,9,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX42(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,10,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX52(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,11,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX13(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,12,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX23(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,13,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,14,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX43(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,15,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX53(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,16,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX14(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,17,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX24(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,18,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX34(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,19,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,20,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX54(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,21,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX15(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,22,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX25(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,23,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX35(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,24,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX45(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,25,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);

      // Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX21(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX31(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX41(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,5,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX51(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,6,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX12(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,7,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,8,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX32(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,9,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX42(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,10,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX52(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,11,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX13(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,12,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX23(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,13,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,14,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX43(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,15,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX53(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,16,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX14(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,17,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX24(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,18,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX34(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,19,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,20,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX54(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,21,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX15(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,22,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX25(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,23,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX35(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,24,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX45(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,25,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
#else
      // Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX21(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX31(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX41(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,5,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX51(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,6,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX12(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,7,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,8,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX32(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,9,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX42(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,10,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX52(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,11,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX13(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,12,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX23(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,13,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,14,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX43(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,15,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX53(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,16,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX14(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,17,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX24(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,18,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX34(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,19,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,20,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX54(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,21,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX15(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,22,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX25(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,23,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX35(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,24,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX45(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);
      IDX3(MatrixAtEdge,25,1,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,wj,Ej);

      // Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX11(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX21(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX31(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX41(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,5,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX51(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,6,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX12(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,7,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX22(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,8,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX32(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,9,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX42(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,10,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX52(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,11,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX13(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,12,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX23(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,13,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX33(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,14,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX43(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,15,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX53(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,16,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX14(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,17,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX24(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,18,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX34(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,19,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX44(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,20,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX54(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,21,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX15(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,22,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX25(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,23,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX35(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,24,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX45(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
      IDX3(MatrixAtEdge,25,2,ipos,NVAR3D*NVAR3D,(bstabilise ? 3 : 2),nedgesim) =
	INVISCIDFLUXJACOBIMATRIX55(-scale,
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),
				   IDX3_COEFFSATEDGE(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,wi,Ei);
#endif
    }
  };

  /*****************************************************************************
   * InviscidFluxJacobiMatrix
   ****************************************************************************/

  template <int isystemcoupling>
  struct InviscidFluxJacobiMatrix : public InviscidFluxJacobiMatrixBase<isystemcoupling>
  {
    // Enable use of inherited functions
    using InviscidFluxJacobiMatrixBase<isystemcoupling>::calcNodeData;
    using InviscidFluxJacobiMatrixBase<isystemcoupling>::calcEdgeData;

    /**************************************************************************
     * Wrapper routine for processing a single node
     *************************************************************************/
    template <typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtDiag,
			     Tc *CoeffsAtDiag,
			     Td scale,
			     Td ui,
			     Td vi, 
			     Td wi,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
      InviscidFluxJacobiMatrixBase<isystemcoupling>::calcNodeData<1>
	(MatrixAtDiag,CoeffsAtDiag,scale,ui,vi,wi,1,ieq,neq,ncoeff);
    }

    /**************************************************************************
     * Wrapper routine for processing a single node
     *************************************************************************/
    template <typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtDiag,
			     Tc *CoeffsAtDiag,
			     Td scale,
			     Td ui,
			     Td vi,
			     Td wi,
			     Td Ei,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
      InviscidFluxJacobiMatrixBase<isystemcoupling>::calcNodeData<1>
	(MatrixAtDiag,CoeffsAtDiag,scale,ui,vi,wi,Ei,1,ieq,neq,ncoeff);
    }
    
    /**************************************************************************
     * Wrapper routine for processing a single edge
     *************************************************************************/
    template <bool bstabilise, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff)
    {
      InviscidFluxJacobiMatrixBase<isystemcoupling>::calcEdgeData<1,bstabilise>
	(MatrixAtEdge,CoeffsAtEdge,scale,ui,uj,vi,vj,wi,wj,1,iedge,nedge,ncoeff);
    }

    /**************************************************************************
     * Wrapper routine for processing a single edge
     *************************************************************************/
    template <bool bstabilise, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Td Ei,
			     Td Ej,
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff)
    {
      InviscidFluxJacobiMatrixBase<isystemcoupling>::calcEdgeData<1,bstabilise>
	(MatrixAtEdge,CoeffsAtEdge,scale,ui,uj,vi,vj,wi,wj,Ei,Ej,1,iedge,nedge,ncoeff);
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase (basic functionality individual specialisations)
   ****************************************************************************/
  
  template <int isystemcoupling, int idissipationtype>
  struct InviscidFluxDissipationMatrixBase
  {
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for block-diagonal matrix
   * computing zero artificial dissipation, aka standard Galerkin
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_SEGREGATED,DISSIPATION_ZERO>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
#pragma unroll
      for (int i=1; i<=NVAR3D; i++)
	IDX3(MatrixAtEdge,i,1,ipos,NVAR3D,3,nedgesim) = 0.0;
    }
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for full matrix
   * computing zero artificial dissipation, aka standard Galerkin
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_ALLCOUPLED,DISSIPATION_ZERO>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Td Ei,
			     Td Ej,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
#pragma unroll
      for (int i=1; i<=NVAR3D*NVAR3D; i++)
	IDX3(MatrixAtEdge,i,1,ipos,NVAR3D*NVAR3D,3,nedgesim) = 0.0;
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for block-diagonal matrix
   * computing scalar artificial dissipation proportional to the
   * spectral radius of the Roe matrix
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_SEGREGATED,DISSIPATION_SCALAR>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for full matrix
   * computing scalar artificial dissipation proportional to the
   * spectral radius of the Roe matrix
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_ALLCOUPLED,DISSIPATION_SCALAR>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Td Ei,
			     Td Ej,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for block-diagonal matrix
   * computing tensorial artificial dissipation of Roe-type
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_SEGREGATED,DISSIPATION_ROE>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for full matrix
   * computing tensorial artificial dissipation of Roe-type
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_ALLCOUPLED,DISSIPATION_ROE>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Td Ei,
			     Td Ej,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for block-diagonal matrix
   * computing scalar artificial dissipation of Rusanov-type
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_SEGREGATED,DISSIPATION_RUSANOV>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationMatrixBase: Specialization for full matrix
   * computing scalar artificial dissipation of Rusanov-type
   ****************************************************************************/

  template <>  
  struct InviscidFluxDissipationMatrixBase<SYSTEM_ALLCOUPLED,DISSIPATION_RUSANOV>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Td Ei,
			     Td Ej,
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };
  
  /*****************************************************************************
   * InviscidFluxDissipationMatrix: Artificial dissipation
   ****************************************************************************/
  
  template <int isystemcoupling, int idissipationtype>
  struct InviscidFluxDissipationMatrix : public InviscidFluxDissipationMatrixBase<isystemcoupling,idissipationtype>
  {   
    // Enable use of inherited functions
    using InviscidFluxDissipationMatrixBase<isystemcoupling,idissipationtype>::calcEdgeData;
    
    /***************************************************************************
     * Wrapper routine for processing a single edge
     **************************************************************************/
    template <typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      InviscidFluxDissipationMatrixBase<isystemcoupling,idissipationtype>::calcEdgeData<1>
	(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,scale,ui,uj,vi,vj,wi,wj,1,iedge,nedge,ncoeff);
    }

    /***************************************************************************
     * Wrapper routine for processing a single edge
     **************************************************************************/
    template <typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcEdgeData(Td *MatrixAtEdge,
			     Tc *CoeffsAtEdge,
			     Td *DataAtEdge,
			     Td scale,
			     Td ui,
			     Td uj,
			     Td vi,
			     Td vj,
			     Td wi,
			     Td wj,
			     Td Ei,
			     Td Ej,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      InviscidFluxDissipationMatrixBase<isystemcoupling,idissipationtype>::calcEdgeData<1>
	(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,scale,ui,uj,vi,vj,wi,wj,Ei,Ej,1,iedge,nedge,ncoeff);
    }
  }; 
  
  /*****************************************************************************
   * This CUDA kernel calculates the diagonal entries of the
   * block-diagonal global operator (baseline implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int threads_per_cta>
  __launch_bounds__(threads_per_cta)
  __global__ void hydro_calcMatDiagMatD3d_baseline(Tc *CoeffsAtDiag,
						   Ti *IdiagList,
						   Tv *vec,
						   Tm *mat,
						   Tm scale,
						   Ti neq,
						   Ti na,
						   Ti ncoeff,
						   Ti neq_last,
						   Ti neq_per_thread=1,
						   Ti neq_offset=0)
  {  
    // Loop over all items per thread
    for (int ipt=0; ipt<neq_per_thread; ++ipt) {
      
      // Global node ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*blockDim.x+neq_offset+threadIdx.x;
      
      if (idx < neq_last) {
	// Get actual equation number
	Ti ieq = IDX2_DIAGLIST(IdiagList,1,idx+1,2,neq);
	
	// Local data at node from local memory
	Tm DataAtDiag[NVAR3D];
	
	// Get solution values at node
	Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
	  gatherNodeData<true>(DataAtDiag,vec,ieq,neq);
	
	// Compute velocities
	Tm ui = XVELOCITY2(DataAtDiag,IDX2,1,NVAR3D,1);
	Tm vi = YVELOCITY2(DataAtDiag,IDX2,1,NVAR3D,1);
	Tm wi = ZVELOCITY2(DataAtDiag,IDX2,1,NVAR3D,1);
	
	// Compute Galerkin coefficient $K_ii$
	InviscidFluxJacobiMatrix<SYSTEM_SEGREGATED>::
	  calcNodeData(DataAtDiag,CoeffsAtDiag,scale,ui,vi,wi,ieq,neq,ncoeff);
	
	// Get diagonal position in the global matrix
	Ti ia  = IDX2_DIAGLIST(IdiagList,2,idx+1,2,neq);
	
	// Build coefficients into global operator
	Matrix<NVAR3D,isystemformat==SYSTEM_BLOCK>::
	  scatterNodeData<true>(mat,DataAtDiag,ia,na);
      }
    }
  };
  
  /*****************************************************************************
   * This CUDA kernel calculates the diagonal entries of the
   * full global operator (baseline implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int threads_per_cta>
  __launch_bounds__(threads_per_cta)
  __global__ void hydro_calcMatDiag3d_baseline(Tc *CoeffsAtDiag,
					       Ti *IdiagList,
					       Tv *vec,
					       Tm *mat,
					       Tm scale,
					       Ti neq,
					       Ti na,
					       Ti ncoeff,
					       Ti neq_last,
					       Ti neq_per_thread=1,
					       Ti neq_offset=0)
  {
    // Loop over all items per thread
    for (int ipt=0; ipt<neq_per_thread; ++ipt) {
      
      // Global node ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*blockDim.x+neq_offset+threadIdx.x;
      
      if (idx < neq_last) {
	// Get actual equation
	Ti ieq = IDX2_DIAGLIST(IdiagList,1,idx+1,2,neq);
	
	// Local solution data at node from local memory
	Tm DataAtDiag[NVAR3D*NVAR3D];
	
	// Get solution values at node
	Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
	  gatherNodeData<true>(DataAtDiag,vec,ieq,neq);
	
	// Compute velocities and energy
	Tm ui = XVELOCITY2(DataAtDiag,IDX2,1,NVAR3D*NVAR3D,1);
	Tm vi = YVELOCITY2(DataAtDiag,IDX2,1,NVAR3D*NVAR3D,1);
	Tm wi = ZVELOCITY2(DataAtDiag,IDX2,1,NVAR3D*NVAR3D,1);
	Tm Ei = SPECIFICTOTALENERGY2(DataAtDiag,IDX2,1,NVAR3D*NVAR3D,1);
	
	// Compute Galerkin coefficient $K_ii$
	InviscidFluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	  calcNodeData(DataAtDiag,CoeffsAtDiag,scale,ui,vi,wi,Ei,ieq,neq,ncoeff);
	
	// Get diagonal position in the global matrix
	Ti ia  = IDX2_DIAGLIST(IdiagList,2,idx+1,2,neq);
	
	// Build coefficients into global operator
	Matrix<NVAR3D*NVAR3D,isystemformat==SYSTEM_BLOCK>::
	  scatterNodeData<true>(mat,DataAtDiag,ia,na);
      }
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the off-diagonal entries of the
   * block-diagonal global operator and assembles the artificial
   * dissipation tensor if required (baseline implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int idissipation,
	    bool blumping,
	    int threads_per_cta>
  __launch_bounds__(threads_per_cta)
  __global__ void hydro_calcMatrixMatD3d_baseline(Tc *CoeffsAtEdge,
						  Ti *IedgeList,
						  Tv *vec,
						  Tm *mat,
						  Tm scale,
						  Ti neq,
						  Ti na,
						  Ti nedge,
						  Ti ncoeff,
						  Ti nedge_last,
						  Ti nedge_per_thread=1,
						  Ti nedge_offset=0)
  {
    // Loop over all items per thread
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {
      
      // Global edge ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*blockDim.x+nedge_offset+threadIdx.x;
      
      if (idx < nedge_last) {
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2_EDGELIST(IedgeList,1,idx+1,6,nedge);
	Ti j = IDX2_EDGELIST(IedgeList,2,idx+1,6,nedge);
	
	// Local solution data at edge from local memory
	Tm DataAtEdge[2*NVAR3D];
	
	// Get solution values at edge endpoints
	Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
	  gatherEdgeData<true>(DataAtEdge,vec,i,j,neq);
	
	// Compute velocities
	Tm ui = XVELOCITY2(DataAtEdge,IDX2,1,NVAR3D,2);
	Tm vi = YVELOCITY2(DataAtEdge,IDX2,1,NVAR3D,2);
	Tm wi = ZVELOCITY2(DataAtEdge,IDX2,1,NVAR3D,2);
	
	Tm uj = XVELOCITY2(DataAtEdge,IDX2,2,NVAR3D,2);
	Tm vj = YVELOCITY2(DataAtEdge,IDX2,2,NVAR3D,2);
	Tm wj = ZVELOCITY2(DataAtEdge,IDX2,2,NVAR3D,2);
	
	if (idissipation == DISSIPATION_ZERO) {
	  
	  // Local matrix data at edge from local memory
	  Tm MatrixAtEdge[2*NVAR3D];
	  
	  // Compute Galerkin coefficient $K_ij$ and $K_ji$
	  InviscidFluxJacobiMatrix<SYSTEM_SEGREGATED>::
	    calcEdgeData<false>(MatrixAtEdge,CoeffsAtEdge,
				scale,ui,uj,vi,vj,wi,wj,idx+1,nedge,ncoeff);
	  
	  // Build matrix coefficients into global operator
	  Matrix<NVAR3D,isystemformat==SYSTEM_BLOCK>::
	    scatterEdgeData<true,false,blumping>(mat,MatrixAtEdge,
						 IedgeList,idx+1,na,nedge);
	  
	} else {
	  
	  // Local matrix data at edge from local memory
	  Tm MatrixAtEdge[3*NVAR3D];
	  
	  // Compute Galerkin coefficient $K_ij$ and $K_ji$
	  InviscidFluxJacobiMatrix<SYSTEM_SEGREGATED>::
	    calcEdgeData<true>(MatrixAtEdge,CoeffsAtEdge,
			       scale,ui,uj,vi,vj,wi,wj,idx+1,nedge,ncoeff);
	  
	  // Compute contribution of artificial diffusion
	  InviscidFluxDissipationMatrix<SYSTEM_SEGREGATED,idissipation>::
	    calcEdgeData(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
			 scale,ui,uj,vi,vj,wi,wj,idx+1,nedge,ncoeff);
	  
	  // Build matrix coefficients into global operator
	  Matrix<NVAR3D,isystemformat==SYSTEM_BLOCK>::
	    scatterEdgeData<true,true,blumping>(mat,MatrixAtEdge,
						IedgeList,idx+1,nedge,na);
	}
      }
    }
  };
  
  /*****************************************************************************
   * This CUDA kernel calculates the off-diagonal entries of the full
   * global operator and assembles the artificial dissipation tensor
   * if required (baseline implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int idissipation,
	    bool blumping,
	    int threads_per_cta>
  __launch_bounds__(threads_per_cta)
  __global__ void hydro_calcMatrix3d_baseline(Tc *CoeffsAtEdge,
					      Ti *IedgeList,
					      Tv *vec,
					      Tm *mat,
					      Tm scale,
					      Ti neq,
					      Ti na,
					      Ti nedge,
					      Ti ncoeff,
					      Ti nedge_last,
					      Ti nedge_per_thread=1,
					      Ti nedge_offset=0)

  {
    // Loop over all items per thread
    for (int ipt=0; ipt<nedge_per_thread; ++ipt) {
      
      // Global edge ID
      Ti idx = (ipt*gridDim.x+blockIdx.x)*blockDim.x+nedge_offset+threadIdx.x;

      if (idx < nedge_last)
	{
	  // Get positions of edge endpoints (idx starts at zero)
	  Ti i = IDX2_EDGELIST(IedgeList,1,idx+1,6,nedge);
	  Ti j = IDX2_EDGELIST(IedgeList,2,idx+1,6,nedge);
	  
	  // Local solution data at edge from local memory
	  Tm DataAtEdge[2*NVAR3D];
	  
	  // Get solution values at edge endpoints
	  Vector<NVAR3D,isystemformat==SYSTEM_BLOCK>::
	    gatherEdgeData<true>(DataAtEdge,vec,i,j,neq);
	  
	  // Compute velocities
	  Tm ui = XVELOCITY2(DataAtEdge,IDX2,1,NVAR3D,2);
	  Tm vi = YVELOCITY2(DataAtEdge,IDX2,1,NVAR3D,2);
	  Tm wi = ZVELOCITY2(DataAtEdge,IDX2,1,NVAR3D,2);
	  
	  Tm uj = XVELOCITY2(DataAtEdge,IDX2,2,NVAR3D,2);
	  Tm vj = YVELOCITY2(DataAtEdge,IDX2,2,NVAR3D,2);
	  Tm wj = ZVELOCITY2(DataAtEdge,IDX2,2,NVAR3D,2);
	  
	  // Compute specific energies
	  Tm Ei = SPECIFICTOTALENERGY2(DataAtEdge,IDX2,1,NVAR3D,2);
	  Tm Ej = SPECIFICTOTALENERGY2(DataAtEdge,IDX2,2,NVAR3D,2);
	  
	  if (idissipation == DISSIPATION_ZERO) {
	    
	    // Local matrix data at edge from local memory
	    Tm MatrixAtEdge[2*NVAR3D*NVAR3D];
	    
	    // Compute Galerkin coefficient $K_ij$ and $K_ji$
	    InviscidFluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	      calcEdgeData<false>(MatrixAtEdge,CoeffsAtEdge,
				  scale,ui,uj,vi,vj,wi,wj,Ei,Ej,idx+1,nedge,ncoeff);
	    
	    // // Build matrix coefficients into global operator
	    Matrix<NVAR3D*NVAR3D,isystemformat==SYSTEM_BLOCK>::
	      scatterEdgeData<true,false,blumping>(mat,MatrixAtEdge,
						   IedgeList,idx+1,na,nedge);
	  } else {
	    
	    // Local matrix data at edge from local memory
	    Tm MatrixAtEdge[3*NVAR3D*NVAR3D];
	    
	    // Compute Galerkin coefficient $K_ij$ and $K_ji$
	    InviscidFluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	      calcEdgeData<true>(MatrixAtEdge,CoeffsAtEdge,
				 scale,ui,uj,vi,vj,wi,wj,Ei,Ej,idx+1,nedge,ncoeff);
	    
	    // Compute contribution of artificial diffusion
	    InviscidFluxDissipationMatrix<SYSTEM_ALLCOUPLED,idissipation>::
	      calcEdgeData(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
			   scale,ui,uj,vi,vj,wi,wj,Ei,Ej,idx+1,nedge,ncoeff);
	    
	    // Build matrix coefficients into global operator
	    Matrix<NVAR3D*NVAR3D,isystemformat==SYSTEM_BLOCK>::
	      scatterEdgeData<true,true,blumping>(mat,MatrixAtEdge,
						  IedgeList,idx+1,na,nedge);
	  }
	}
    }
  };
  
  /*****************************************************************************
   * Internal C++ functions which invoke the CUDA kernels
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti>
  inline
  int hydro_calcMatDiagMatD3d_cuda(__SIZET *d_CoeffsAtDiag,
				   __SIZET *d_IdiagList,
				   __SIZET *d_vec,
				   __SIZET *d_mat,
				   Tm scale,
				   Ti nblocks,
				   Ti neq,
				   Ti na,
				   Ti ncoeff,
				   cudaStream_t stream=0)
  {
    const cudaDeviceProp *devProp = coproc_getCurrentDeviceProp();
    
    // Strategy: run the largest possible number of blocks with a
    // predefined number of compute/dma threads per block and let each
    // compute thread process the minimal number of edges
    const int compute_threads_per_cta = CUDADMA_COMPUTE_THREADS_PER_CTA;
    const int dma_threads_per_ld      = CUDADMA_THREADS_PER_LD;
    const int dma_lds                 = CUDADMA_DMA_LDS;
    int neq_per_thread_cudaDMA        = CUDADMA_NEQ_PER_THREAD;

    const int threads_per_cta_baseline = BASELINE_THREADS_PER_CTA;
    int neq_per_thread_baseline        = BASELINE_NEQ_PER_THREAD;
    
    int blocks, threads, neq_cudaDMA, neq_baseline;
    prepare_cudaDMA(devProp, neq,
		    &neq_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &neq_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(devProp, neq-neq_cudaDMA,
		     &neq_per_thread_baseline, threads_per_cta_baseline,
		     &blocks, &threads, &neq_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);

    Tv *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
        
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);
    
#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0)
      	// CudaDMA implementation
      	hydro_calcMatDiagMatD3d_cudaDMA
      	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,compute_threads_per_cta,dma_threads_per_ld>
      	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
      						       IdiagList,
      						       vec, mat, scale,
      						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
#endif

      if (grid_baseline.x>0)
      	// Baseline implementation
      	hydro_calcMatDiagMatD3d_baseline
      	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,threads_per_cta_baseline>
      	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
      							 IdiagList,
      							 vec, mat, scale,
      							 neq, na, ncoeff,
							 neq, 
							 neq_per_thread_baseline,
							 neq_cudaDMA);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      __SIZET cmemPool[NVAR3D];
#pragma unroll
      for (int i=0; i<NVAR3D; i++)
	cmemPool[i] = d_mat[i*(NVAR3D+1)];
      
      cudaMemcpyToSymbolAsync("constMemPool", cmemPool,
      			      sizeof(__SIZET)*NVAR3D, 0,
      			      cudaMemcpyHostToDevice,
      			      stream);
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");

#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatDiagMatD3d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
						       IdiagList,
						       vec, mat, scale,
						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
#endif

      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatDiagMatD3d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
							 IdiagList,
							 vec, mat, scale,
							 neq, na, ncoeff,
							 neq,
							 neq_per_thread_baseline,
							 neq_cudaDMA);
    }
    
    coproc_checkError("hydro_calcMatDiagMatD3d_cuda");
    return 0;
  };

  /****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti>
  inline
  int hydro_calcMatDiag3d_cuda(__SIZET *d_CoeffsAtDiag,
			       __SIZET *d_IdiagList,
			       __SIZET *d_vec,
			       __SIZET *d_mat,
			       Tm scale,
			       Ti nblocks,
			       Ti neq,
			       Ti na,
			       Ti ncoeff,
			       cudaStream_t stream=0)
  {
    const cudaDeviceProp *devProp = coproc_getCurrentDeviceProp();
    
    // Strategy: run the largest possible number of blocks with a
    // predefined number of compute/dma threads per block and let each
    // compute thread process the minimal number of edges
    const int compute_threads_per_cta = CUDADMA_COMPUTE_THREADS_PER_CTA;
    const int dma_threads_per_ld      = CUDADMA_THREADS_PER_LD;
    const int dma_lds                 = CUDADMA_DMA_LDS;
    int neq_per_thread_cudaDMA        = CUDADMA_NEQ_PER_THREAD;

    const int threads_per_cta_baseline = BASELINE_THREADS_PER_CTA;
    int neq_per_thread_baseline        = BASELINE_NEQ_PER_THREAD;
    
    int blocks, threads, neq_cudaDMA, neq_baseline;
    prepare_cudaDMA(devProp, neq,
		    &neq_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &neq_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(devProp, neq-neq_cudaDMA,
		     &neq_per_thread_baseline, threads_per_cta_baseline,
		     &blocks, &threads, &neq_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);

    Tv *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
    
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);
      
#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0)
	//CudaDMA implementation
	hydro_calcMatDiag3d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
						       IdiagList,
						       vec, mat, scale,
						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
#endif

      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatDiag3d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
							 IdiagList,
							 vec, mat, scale,
							 neq, na, ncoeff,
							 neq,
							 neq_per_thread_baseline,
							 neq_cudaDMA);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      cudaMemcpyToSymbolAsync("constMemPool", d_mat,
			      sizeof(__SIZET)*NVAR3D*NVAR3D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");

#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatDiag3d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
						       IdiagList,
						       vec, mat, scale,
						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
#endif

      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatDiag3d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
							 IdiagList,
							 vec, mat, scale,
							 neq, na, ncoeff,
							 neq,
							 neq_per_thread_baseline,
							 neq_cudaDMA);
    }
    
    coproc_checkError("hydro_calcMatDiag2d_cuda");
    return 0;
  };

  /*****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int idissipationtype,
	    bool blumping>
  inline
  int hydro_calcMatrixMatD3d_cuda(__SIZET *d_CoeffsAtEdge,
				  __SIZET *d_IedgeList,
				  __SIZET *d_vec,
				  __SIZET *d_mat,
				  Tm scale,
				  Ti nblocks,
				  Ti neq,
				  Ti na,
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
    
    Tv  *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);

    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);

#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrixMatD3d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
#endif

      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatrixMatD3d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping,threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      __SIZET cmemPool[NVAR3D];
#pragma unroll
      for (int i=0; i<NVAR3D; i++)
	cmemPool[i] = d_mat[i*(NVAR3D+1)];
      
      cudaMemcpyToSymbolAsync("constMemPool", cmemPool,
			      sizeof(__SIZET)*NVAR3D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");

#ifdef CUDADMA_KERNEL
      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrixMatD3d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,idissipationtype,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
#endif

      if (grid_baseline.x>0)
	hydro_calcMatrixMatD3d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,idissipationtype,blumping,threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1,
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
    }
    
    coproc_checkError("hydro_calcMatrixMatD3d_cuda");
    return 0;
  };

  /****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int idissipationtype,
	    bool blumping>
  inline
  int hydro_calcMatrix3d_cuda(__SIZET *d_CoeffsAtEdge,
			      __SIZET *d_IedgeList,
			      __SIZET *d_vec,
			      __SIZET *d_mat,
			      Tm scale,
			      Ti nblocks,
			      Ti neq,
			      Ti na,
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

    Tv  *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);

    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);

#ifdef CUDADMA_KERNEL    
      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrix3d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
#endif

      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatrix3d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping,threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      cudaMemcpyToSymbolAsync("constMemPool", d_mat,
			      sizeof(__SIZET)*NVAR3D*NVAR3D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");

#ifdef CUDADMA_KERNEL    
      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrix3d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,idissipationtype,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
#endif

      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatrix3d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,idissipationtype,blumping,threads_per_cta_baseline>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
    }
        
    coproc_checkError("hydro_calcMatrix3d_cuda");
    return 0;
  };
  
  /*****************************************************************************
   * External C functions which can be called from the Fortran code
   ****************************************************************************/

  extern "C" {
    __INT FNAME(hydro_calcmatdiagmatd3d_cuda)(__SIZET *d_CoeffsAtDiag,
					      __SIZET *d_IdiagList,
					      __SIZET *d_vec,
					      __SIZET *d_mat,
					      __DP *scale,
					      __INT *nblocks,
					      __INT *neq,
					      __INT *na,
					      __INT *ncoeff,
					      __I64 *stream)
    {
      return (__INT) hydro_calcMatDiagMatD3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtDiag, d_IdiagList, d_vec, d_mat,
			       *scale, *nblocks, *neq, *na, *ncoeff,
			       (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatdiag3d_cuda)(__SIZET *d_CoeffsAtDiag,
					  __SIZET *d_IdiagList,
					  __SIZET *d_vec,
					  __SIZET *d_mat,
					  __DP *scale,
					  __INT *nblocks,
					  __INT *neq,
					  __INT *na,
					  __INT *ncoeff,
					  __I64 *stream)
    {
      return (__INT) hydro_calcMatDiag3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtDiag, d_IdiagList, d_vec, d_mat,
			       *scale, *nblocks, *neq, *na, *ncoeff,
			       (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatgalmatd3d_cuda)(__SIZET *d_CoeffsAtEdge,
					     __SIZET *d_IedgeList,
					     __SIZET *d_vec,
					     __SIZET *d_mat,
					     __DP *scale,
					     __INT *nblocks,
					     __INT *neq,
					     __INT *na,
					     __INT *nedge,
					     __INT *ncoeff,
					     __INT *nedgeset,
					     __INT *iedgeset,
					     __INT *cconstrType,
					     __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }
    
    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatgalerkin3d_cuda)(__SIZET *d_CoeffsAtEdge,
					      __SIZET *d_IedgeList,
					      __SIZET *d_vec,
					      __SIZET *d_mat,
					      __DP *scale,
					      __INT *nblocks,
					      __INT *neq,
					      __INT *na,
					      __INT *nedge,
					      __INT *ncoeff,
					      __INT *nedgeset,
					      __INT *iedgeset,
					      __INT *cconstrType,
					      __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatscdissmatd3d_cuda)(__SIZET *d_CoeffsAtEdge,
					     __SIZET *d_IedgeList,
					     __SIZET *d_vec,
					     __SIZET *d_mat,
					     __DP *scale,
					     __INT *nblocks,
					     __INT *neq,
					     __INT *na,
					     __INT *nedge,
					     __INT *ncoeff,
					     __INT *nedgeset,
					     __INT *iedgeset,
					     __INT *cconstrType,
					     __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatscdiss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					    __SIZET *d_IedgeList,
					    __SIZET *d_vec,
					    __SIZET *d_mat,
					    __DP *scale,
					    __INT *nblocks,
					    __INT *neq,
					    __INT *na,
					    __INT *nedge,
					    __INT *ncoeff,
					    __INT *nedgeset,
					    __INT *iedgeset,
					    __INT *cconstrType,
					    __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatroedissmatd3d_cuda)(__SIZET *d_CoeffsAtEdge,
						 __SIZET *d_IedgeList,
						 __SIZET *d_vec,
						 __SIZET *d_mat,
						 __DP *scale,
						 __INT *nblocks,
						 __INT *neq,
						 __INT *na,
						 __INT *nedge,
						 __INT *ncoeff,
						 __INT *nedgeset,
						 __INT *iedgeset,
						 __INT *cconstrType,
						 __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatroediss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					     __SIZET *d_IedgeList,
					     __SIZET *d_vec,
					     __SIZET *d_mat,
					     __DP *scale,
					     __INT *nblocks,
					     __INT *neq,
					     __INT *na,
					     __INT *nedge,
					     __INT *ncoeff,
					     __INT *nedgeset,
					     __INT *iedgeset,
					     __INT *cconstrType,
					     __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatrusdissmatd3d_cuda)(__SIZET *d_CoeffsAtEdge,
						 __SIZET *d_IedgeList,
						 __SIZET *d_vec,
						 __SIZET *d_mat,
						 __DP *scale,
						 __INT *nblocks,
						 __INT *neq,
						 __INT *na,
						 __INT *nedge,
						 __INT *ncoeff,
						 __INT *nedgeset,
						 __INT *iedgeset,
						 __INT *cconstrType,
						 __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatrusdiss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					     __SIZET *d_IedgeList,
					     __SIZET *d_vec,
					     __SIZET *d_mat,
					     __DP *scale,
					     __INT *nblocks,
					     __INT *neq,
					     __INT *na,
					     __INT *nedge,
					     __INT *ncoeff,
					     __INT *nedgeset,
					     __INT *iedgeset,
					     __INT *cconstrType,
					     __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix3d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }
  };
}
