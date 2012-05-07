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
#include "../../cudaDMA.h"

#define LANGUAGE LANGUAGE_C
#include "../../../../../kernel/System/idxmanager.h"

#define HYDRO_NDIM 2
#include "hydro.h"

// Warp size
#define WARP_SIZE 32

// Number of compute threads per cooperative thread block (CTA)
#define COMPUTE_THREADS_PER_CTA (WARP_SIZE * 16)

// Number of DMA threads per load/store operation
#define DMA_THREADS_PER_LD      (WARP_SIZE * 0)

// Memory pool in constant device memory
__device__ __constant__ __SIZET constMemPool[NVAR2D*NVAR2D];

#include "../../cudaGatherScatter.h"

// Define kernels
#define hydro_calcMatDiagMatD2d_knl hydro_calcMatDiagMatD2d_baseline
#define hydro_calcMatDiag2d_knl     hydro_calcMatDiag2d_baseline
#define hydro_calcMatrixMatD2d_knl  hydro_calcMatrixMatD2d_baseline
#define hydro_calcMatrix2d_knl      hydro_calcMatrix2d_baseline

namespace hydro2d_cuda
{

  /*****************************************************************************
   * CUDA kernels for hydrodynamic model in 2D
   ****************************************************************************/
  
  using namespace std;
  
  /*****************************************************************************
   * FluxJacobiMatrixBase (basic functionality and specialisations)
   ****************************************************************************/

  template <int isystemcoupling>
  struct FluxJacobiMatrixBase
  {
  };

  /*****************************************************************************
   * FluxJacobiMatrixBase: Specialization for block-diagonal matrix
   ****************************************************************************/

  template <>
  struct FluxJacobiMatrixBase<SYSTEM_SEGREGATED>
  {
    /***************************************************************************
     * Compute flux Jacobian matrix for neqsim nodes
     **************************************************************************/
    template <int neqsim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtNode,
			     Tc *CoeffsAtNode,
			     Td scale,
			     Td ui,
			     Td vi, 
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff,
			     Ti tid)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX2(MatrixAtNode,1,tid,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,2,tid,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,3,tid,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,4,tid,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
#else
      // Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX2(MatrixAtNode,1,tid,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,2,tid,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,3,tid,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtNode,4,tid,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,_);
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
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff,
			     Ti tid)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(MatrixAtEdge,1,1,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,2,1,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,3,1,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,4,1,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);

      // Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(MatrixAtEdge,1,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,2,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,3,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,4,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
#else
      // Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(MatrixAtEdge,1,1,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,2,1,tid,NVAR2D,2,1,) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,3,1,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,4,1,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      
      // Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(MatrixAtEdge,1,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,2,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,3,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,4,2,tid,NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
#endif
    }
  };

  /*****************************************************************************
   * FluxJacobiMatrixBase: Specialization for full matrix
   ****************************************************************************/

  template <>
  struct FluxJacobiMatrixBase<SYSTEM_ALLCOUPLED>
  {
    /***************************************************************************
     * Compute flux Jacobian matrix for neqsim nodes
     **************************************************************************/
    template <int neqsim, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtNode,
			     Tc *CoeffsAtNode,
			     Td scale,
			     Td ui,
			     Td vi,
			     Td Ei, 
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff,
			     Ti tid)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtNode,1,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,2,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX21(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,3,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX31(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,4,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX41(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,5,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX12(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,6,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,7,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX32(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,8,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX42(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,9,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX13(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,10,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX23(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,11,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,12,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX43(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,13,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX14(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,14,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX24(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,15,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX34(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,16,tid,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
#else
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtNode,1,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX11(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,2,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX21(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,3,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX31(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,4,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX41(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,5,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX12(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,6,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX22(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,7,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX32(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,8,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX42(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,9,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX13(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,10,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX23(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,11,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX33(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,12,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX43(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,13,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX14(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,14,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX24(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,15,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX34(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtNode,16,tid,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX44(scale,
			   IDX2T(CoeffsAtNode,1,ieq,ncoeff,neq),
			   IDX2T(CoeffsAtNode,2,ieq,ncoeff,neq),ui,vi,Ei);
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
			     Td Ei,
			     Td Ej,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff,
			     Ti tid)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(MatrixAtEdge,1,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,2,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX21(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,3,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX31(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,4,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX41(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,5,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX12(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,6,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,7,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX32(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,8,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX42(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,9,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX13(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,10,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX23(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,11,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,12,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX43(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,13,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX14(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,14,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX24(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,15,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX34(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,16,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);

      // Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(MatrixAtEdge,1,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,2,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX21(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,3,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX31(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,4,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX41(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,5,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX12(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,6,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,7,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX32(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,8,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX42(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,9,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX13(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,10,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX23(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,11,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,12,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX43(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,13,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX14(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,14,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX24(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,15,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX34(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,16,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
#else
      // Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(MatrixAtEdge,1,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,2,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX21(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,3,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX31(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,4,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX41(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,5,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX12(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,6,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,7,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX32(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,8,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX42(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,9,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX13(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,10,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX23(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,11,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,12,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX43(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,13,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX14(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,14,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX24(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,15,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX34(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,16,1,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      
      // Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(MatrixAtEdge,1,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,2,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX21(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,3,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX31(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,4,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX41(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,5,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX12(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,6,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,7,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX32(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,8,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX42(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,9,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX13(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,10,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX23(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,11,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,12,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX43(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,13,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX14(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,14,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX24(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,15,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX34(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,16,2,tid,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),neqsim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
#endif
    }
  };

  /*****************************************************************************
   * FluxJacobiMatrix
   ****************************************************************************/

  template <int isystemcoupling>
  struct FluxJacobiMatrix : public FluxJacobiMatrixBase<isystemcoupling>
  {
    // Enable use of inherited functions
    using FluxJacobiMatrixBase<isystemcoupling>::calcNodeData;
    using FluxJacobiMatrixBase<isystemcoupling>::calcEdgeData;

    /**************************************************************************
     * Wrapper routine for processing a single node
     *************************************************************************/
    template <typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtNode,
			     Tc *CoeffsAtNode,
			     Td scale,
			     Td ui,
			     Td vi, 
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
      FluxJacobiMatrixBase<isystemcoupling>::calcNodeData<1>
	(MatrixAtNode,CoeffsAtNode,scale,ui,vi,ieq,neq,ncoeff,1);
    }

    /**************************************************************************
     * Wrapper routine for processing a single node
     *************************************************************************/
    template <typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void calcNodeData(Td *MatrixAtNode,
			     Tc *CoeffsAtNode,
			     Td scale,
			     Td ui,
			     Td vi,
			     Td Ei,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
      FluxJacobiMatrixBase<isystemcoupling>::calcNodeData<1>
	(MatrixAtNode,CoeffsAtNode,scale,ui,vi,Ei,ieq,neq,ncoeff,1);
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
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff)
    {
      FluxJacobiMatrixBase<isystemcoupling>::calcEdgeData<1,bstabilise>
	(MatrixAtEdge,CoeffsAtEdge,scale,ui,uj,vi,vj,iedge,nedge,ncoeff,1);
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
			     Td Ei,
			     Td Ej,
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff)
    {
      FluxJacobiMatrixBase<isystemcoupling>::calcEdgeData<1,bstabilise>
	(MatrixAtEdge,CoeffsAtEdge,scale,ui,uj,vi,vj,Ei,Ej,iedge,nedge,ncoeff,1);
    }
  };

  /*****************************************************************************
   * DissipationBase (basic functionality individual specialisations)
   ****************************************************************************/
  
  template <int isystemcoupling, int idissipationtype>
  struct DissipationBase
  {
  };
  
  /*****************************************************************************
   * DissipationBase: Specialization for block-diagonal matrix
   * computing zero artificial dissipation, aka standard Galerkin
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_SEGREGATED,DISSIPATION_ZERO>
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
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff,
			     Ti tid)
    {
#pragma unroll
      for (int i=1; i<=NVAR2D; i++)
	IDX3(MatrixAtEdge,i,1,tid,NVAR2D,3,nedgesim) = 0.0;
    }
  };
  
  /*****************************************************************************
   * DissipationBase: Specialization for full matrix
   * computing zero artificial dissipation, aka standard Galerkin
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_ALLCOUPLED,DISSIPATION_ZERO>
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
			     Td Ei,
			     Td Ej,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff,
			     Ti tid)
    {
#pragma unroll
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	IDX3(MatrixAtEdge,i,1,tid,NVAR2D*NVAR2D,3,nedgesim) = 0.0;
    }
  };
  
  /*****************************************************************************
   * Dissipation: Artificial dissipation
   ****************************************************************************/
  
  template <int isystemcoupling, int idissipationtype>
  struct Dissipation : public DissipationBase<isystemcoupling,idissipationtype>
  {   
    // Enable use of inherited functions
    using DissipationBase<isystemcoupling,idissipationtype>::calcEdgeData;
    
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
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      DissipationBase<isystemcoupling,idissipationtype>::calcEdgeData<1>
	(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,scale,ui,uj,vi,vj,iedge,nedge,ncoeff,1);
    }

    /*
     * Wrapper routine for processing a single edge
     */
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
			     Td Ei,
			     Td Ej,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      DissipationBase<isystemcoupling,idissipationtype>::calcEdgeData<1>
	(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,scale,ui,uj,vi,vj,Ei,Ej,iedge,nedge,ncoeff,1);
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
	    int isystemformat>
  __global__ void hydro_calcMatDiagMatD2d_baseline(Tc *CoeffsAtDiag,
						   Ti *IdiagList,
						   Tv *vec,
						   Tm *mat,
						   Tm scale,
						   Ti neq,
						   Ti na,
						   Ti ncoeff)
  {  
    // Global node ID
    Ti idx = blockIdx.x * COMPUTE_THREADS_PER_CTA + threadIdx.x;
    
    if (idx<neq) {
      // Get actual equation number
      Ti ieq = IDX2(IdiagList,1,idx+1,2,neq);
      
      // Local data at node from local memory
      Tm DataAtNode[NVAR2D];
      
      // Get solution values at node
      Vector<NVAR2D,isystemformat>::
	gatherNodeData(DataAtNode,vec,ieq,neq);
      
      // Compute velocities
      Tm ui = XVELOCITY2(DataAtNode,IDX2,1,NVAR2D,1);
      Tm vi = YVELOCITY2(DataAtNode,IDX2,1,NVAR2D,1);
      
      // Compute Galerkin coefficient $K_ii$
      FluxJacobiMatrix<SYSTEM_SEGREGATED>::
	calcNodeData(DataAtNode,CoeffsAtDiag,scale,ui,vi,ieq,neq,ncoeff);
      
      // Get diagonal position in the global matrix
      Ti ia  = IDX2(IdiagList,2,idx+1,2,neq);
      
      // Build coefficients into global operator
      Matrix<NVAR2D,isystemformat>::
	scatterNodeData(mat,DataAtNode,ia,na);
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the diagonal entries of the
   * block-diagonal global operator (shared memory implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat>
  __global__ void hydro_calcMatDiagMatD2d_shmem(Tc *CoeffsAtDiag,
						Ti *IdiagList,
						Tv *vec,
						Tm *mat,
						Tm scale,
						Ti neq,
						Ti na,
						Ti ncoeff)
  {
    // Shared memory
    Tm DataAtNode[NVAR2D*COMPUTE_THREADS_PER_CTA];
    
    // Global node ID
    Ti idx = blockIdx.x * COMPUTE_THREADS_PER_CTA + threadIdx.x;
    
    if (idx<neq) {
      // Local thread ID
      Ti tid = threadIdx.x;

      // Get actual equation number
      Ti ieq = IDX2(IdiagList,1,idx+1,2,neq);
      
      // Get solution values at node
      Vector<NVAR2D,isystemformat>::
	gatherNodeData<COMPUTE_THREADS_PER_CTA>(DataAtNode,vec,ieq,neq,tid+1);

      // Compute velocities
      Tm ui = XVELOCITY3(DataAtNode,IDX3,1,tid+1,NVAR2D,1,COMPUTE_THREADS_PER_CTA);
      Tm vi = YVELOCITY3(DataAtNode,IDX3,1,tid+1,NVAR2D,1,COMPUTE_THREADS_PER_CTA);
      
      // Compute Galerkin coefficient $K_ii$
      FluxJacobiMatrix<SYSTEM_SEGREGATED>::
	calcNodeData<COMPUTE_THREADS_PER_CTA>
	(DataAtNode,CoeffsAtDiag,scale,ui,vi,ieq,neq,ncoeff,tid+1);

      // Get diagonal position in the global matrix
      Ti ia  = IDX2(IdiagList,2,idx+1,2,neq);
      
      // Build coefficients into global operator
      Matrix<NVAR2D,isystemformat>::
	scatterNodeData<COMPUTE_THREADS_PER_CTA>(mat,DataAtNode,ia,na,tid+1);
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
	    int isystemformat>
  __global__ void hydro_calcMatDiag2d_baseline(Tc *CoeffsAtDiag,
					       Ti *IdiagList,
					       Tv *vec,
					       Tm *mat,
					       Tm scale,
					       Ti neq,
					       Ti na,
					       Ti ncoeff)
  {
    // Global node ID
    Ti idx = blockIdx.x * COMPUTE_THREADS_PER_CTA + threadIdx.x;
    
    if (idx<neq) {
      // Get actual equation
      Ti ieq = IDX2(IdiagList,1,idx+1,2,neq);
      
      // Local solution data at node from local memory
      Tm DataAtNode[NVAR2D*NVAR2D];
      
      // Get solution values at node
      Vector<NVAR2D,isystemformat>::
	gatherNodeData(DataAtNode,vec,ieq,neq);
      
      // Compute velocities and energy
      Tm ui = XVELOCITY2(DataAtNode,IDX2,1,NVAR2D*NVAR2D,1);
      Tm vi = YVELOCITY2(DataAtNode,IDX2,1,NVAR2D*NVAR2D,1);
      Tm Ei = SPECIFICTOTALENERGY2(DataAtNode,IDX2,1,NVAR2D*NVAR2D,1);
      
      // Compute Galerkin coefficient $K_ii$
      FluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	calcNodeData(DataAtNode,CoeffsAtDiag,scale,ui,vi,Ei,ieq,neq,ncoeff);
      
      // Get diagonal position in the global matrix
      Ti ia  = IDX2(IdiagList,2,idx+1,2,neq);
      
      // Build coefficients into global operator
      Matrix<NVAR2D*NVAR2D,isystemformat>::
	scatterNodeData(mat,DataAtNode,ia,na);
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
	    bool blumping>
  __global__ void hydro_calcMatrixMatD2d_baseline(Tc *CoeffsAtEdge,
						  Ti *IedgeList,
						  Tv *vec,
						  Tm *mat,
						  Tm scale,
						  Ti neq,
						  Ti na,
						  Ti nedge,
						  Ti ncoeff,
						  Ti nedges,
						  Ti iedgeset)
  {
    // Global node ID
    Ti idx = blockIdx.x * COMPUTE_THREADS_PER_CTA + threadIdx.x;

    if (idx<nedges) {
      // Get positions of edge endpoints (idx starts at zero)
      Ti i = IDX2(IedgeList,1,iedgeset+idx,6,nedge);
      Ti j = IDX2(IedgeList,2,iedgeset+idx,6,nedge);
      
      // Local solution data at edge from local memory
      Tm DataAtEdge[2*NVAR2D];
      
      // Get solution values at edge endpoints
      Vector<NVAR2D,isystemformat>::
	gatherEdgeData(DataAtEdge,vec,i,j,neq);
      
      // Compute velocities
      Tm ui = XVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
      Tm vi = YVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
      
      Tm uj = XVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);
      Tm vj = YVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);

      if (idissipation == DISSIPATION_ZERO) {

	// Local matrix data at edge from local memory
	Tm MatrixAtEdge[2*NVAR2D];

	// Compute Galerkin coefficient $K_ij$ and $K_ji$
	FluxJacobiMatrix<SYSTEM_SEGREGATED>::
	  calcEdgeData<false>(MatrixAtEdge,CoeffsAtEdge,
			      scale,ui,uj,vi,vj,iedgeset+idx,nedge,ncoeff);

	// Build matrix coefficients into global operator
	Matrix<NVAR2D,isystemformat>::
	  scatterEdgeData<false,blumping>(mat,MatrixAtEdge,
					  IedgeList,iedgeset+idx,na,nedge);
	
      } else {
	
	// Local matrix data at edge from local memory
	Tm MatrixAtEdge[3*NVAR2D];
	
	// Compute Galerkin coefficient $K_ij$ and $K_ji$
	FluxJacobiMatrix<SYSTEM_SEGREGATED>::
	  calcEdgeData<true>(MatrixAtEdge,CoeffsAtEdge,
			     scale,ui,uj,vi,vj,iedgeset+idx,nedge,ncoeff);
	
	// Compute contribution of artificial diffusion
	Dissipation<SYSTEM_SEGREGATED,idissipation>::
	  calcEdgeData(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
		       scale,ui,uj,vi,vj,iedgeset+idx,nedge,ncoeff);
	
	// Build matrix coefficients into global operator
	Matrix<NVAR2D,isystemformat>::
	  scatterEdgeData<true,blumping>(mat,MatrixAtEdge,
					 IedgeList,iedgeset+idx,na,nedge);
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
	    bool blumping>
  __global__ void hydro_calcMatrix2d_baseline(Tc *CoeffsAtEdge,
					      Ti *IedgeList,
					      Tv *vec,
					      Tm *mat,
					      Tm scale,
					      Ti neq,
					      Ti na,
					      Ti nedge,
					      Ti ncoeff,
					      Ti nedges,
					      Ti iedgeset)
  {
    // Global node ID
    Ti idx = blockIdx.x * COMPUTE_THREADS_PER_CTA + threadIdx.x;
    
    if (idx<nedges)
      {
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(IedgeList,1,iedgeset+idx,6,nedge);
	Ti j = IDX2(IedgeList,2,iedgeset+idx,6,nedge);
	
	// Local solution data at edge from local memory
	Tm DataAtEdge[2*NVAR2D];
	
	// Get solution values at edge endpoints
	Vector<NVAR2D,isystemformat>::
	  gatherEdgeData(DataAtEdge,vec,i,j,neq);
	
	// Compute velocities
	Tm ui = XVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
	Tm vi = YVELOCITY2(DataAtEdge,IDX2,1,NVAR2D,2);
	
	Tm uj = XVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);
	Tm vj = YVELOCITY2(DataAtEdge,IDX2,2,NVAR2D,2);
	
	// Compute specific energies
	Tm Ei = SPECIFICTOTALENERGY2(DataAtEdge,IDX2,1,NVAR2D,2);
	Tm Ej = SPECIFICTOTALENERGY2(DataAtEdge,IDX2,2,NVAR2D,2);
	
	if (idissipation == DISSIPATION_ZERO) {

	  // Local matrix data at edge from local memory
	  Tm MatrixAtEdge[2*NVAR2D*NVAR2D];
	  
	  // Compute Galerkin coefficient $K_ij$ and $K_ji$
	  FluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	    calcEdgeData<false>(MatrixAtEdge,CoeffsAtEdge,
				scale,ui,uj,vi,vj,Ei,Ej,iedgeset+idx,nedge,ncoeff);
	  
	  // // Build matrix coefficients into global operator
	  Matrix<NVAR2D*NVAR2D,isystemformat>::
	    scatterEdgeData<false,blumping>(mat,MatrixAtEdge,
	  				    IedgeList,iedgeset+idx,na,nedge);
	} else {
	  
	  // Local matrix data at edge from local memory
	  Tm MatrixAtEdge[3*NVAR2D*NVAR2D];
	  
	  // Compute Galerkin coefficient $K_ij$ and $K_ji$
	  FluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	    calcEdgeData<true>(MatrixAtEdge,CoeffsAtEdge,
			       scale,ui,uj,vi,vj,Ei,Ej,iedgeset+idx,nedge,ncoeff);
	  
	  // Compute contribution of artificial diffusion
	  Dissipation<SYSTEM_ALLCOUPLED,idissipation>::
	    calcEdgeData(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
			 scale,ui,uj,vi,vj,Ei,Ej,iedgeset+idx,nedge,ncoeff);
	  
	  // Build matrix coefficients into global operator
	  Matrix<NVAR2D*NVAR2D,isystemformat>::
	    scatterEdgeData<true,blumping>(mat,MatrixAtEdge,
					   IedgeList,iedgeset+idx,na,nedge);
	}
      }
  }
  
  /*****************************************************************************
   * Internal C++ functions which invoke the CUDA kernels
   *****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti>
  inline
  int hydro_calcMatDiagMatD2d_cuda(__SIZET *d_CoeffsAtDiag,
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
    Tv *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
    
    // The total number of blocks depends on the problem size NEQ and
    // on the number of compute threads per cooperative thread block
    int blocks = (neq+COMPUTE_THREADS_PER_CTA-1)/COMPUTE_THREADS_PER_CTA;
    dim3 grid(blocks, 1, 1);

    // The total number of threads per cooperative thread block
    // depends on the number of compute threads plus the number od DMA
    // threads multiplied by the number of load/store operations
    dim3 block( COMPUTE_THREADS_PER_CTA + 1*DMA_THREADS_PER_LD, 1, 1);
    
    cout << "hydro_calcMatDiagMatD2d_cuda" 
	 << " nblocks=" << nblocks
	 << " neq=" << neq
	 << " grid=" << grid.x << " block=" << block.x << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);
      
      cudaEventRecord(start,stream);
      
      hydro_calcMatDiagMatD2d_knl
       	<Tc,Tv,Tm,Ti,SYSTEM_SCALAR>
      	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
       				     IdiagList,
				     vec, mat, scale,
				     neq, na, ncoeff);
      
      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      __SIZET cmemPool[NVAR2D];
#pragma unroll
      for (int i=0; i<NVAR2D; i++)
	cmemPool[i] = d_mat[i*(NVAR2D+1)];
      
      cudaEventRecord(start,stream);
      
      cudaMemcpyToSymbolAsync("constMemPool", cmemPool,
			      sizeof(__SIZET)*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");
      
      hydro_calcMatDiagMatD2d_knl
	<Tc,Tv,Tm,Ti,SYSTEM_BLOCK>
	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
				     IdiagList,
				     vec, mat, scale,
				     neq, na, ncoeff);

      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    }

    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "Memory NEQ:   " << NVAR2D*neq*sizeof(Tv)/1000000.0f
	 << " MB" << " NEQ=" << neq << endl;
    cout << "Elapsed time: " << elapsedTime << " ms" << endl;
    cout << "Bandwidth:    " << (4*NVAR2D*neq*sizeof(Tv)+
				 2*neq*sizeof(Ti))/1000000000.0f/elapsedTime*1000.0f
	 << " GB/s" << endl;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    coproc_checkErrors("hydro_calcMatDiagMatD2d_cuda");
    return 0;
  };

  /*****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti>
  inline
  int hydro_calcMatDiag2d_cuda(__SIZET *d_CoeffsAtDiag,
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
    Tv *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
  
    // The total number of blocks depends on the problem size NEQ and
    // on the number of compute threads per cooperative thread block
    int blocks = (neq+COMPUTE_THREADS_PER_CTA-1)/COMPUTE_THREADS_PER_CTA;
    dim3 grid(blocks, 1, 1);

    // The total number of threads per cooperative thread block
    // depends on the number of compute threads plus the number od DMA
    // threads multiplied by the number of load/store operations
    dim3 block( COMPUTE_THREADS_PER_CTA + 1*DMA_THREADS_PER_LD, 1, 1);

    cout << "hydro_calcMatDiag2d_cuda" 
	 << " nblocks=" << nblocks
	 << " neq=" << neq
	 << " grid=" << grid.x << " block=" << block.x << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);
      
      cudaEventRecord(start,stream);

      hydro_calcMatDiag2d_knl
	<Tc,Tv,Tm,Ti,SYSTEM_SCALAR>
	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
				     IdiagList,
				     vec, mat, scale,
				     neq, na, ncoeff);

      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    } else {
      cudaEventRecord(start,stream);
      
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      cudaMemcpyToSymbolAsync("constMemPool", d_mat,
			      sizeof(__SIZET)*NVAR2D*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");
      
      hydro_calcMatDiag2d_knl
	<Tc,Tv,Tm,Ti,SYSTEM_BLOCK>
	<<<grid, block, 0, stream>>>(CoeffsAtDiag,
				     IdiagList,
				     vec, mat, scale,
				     neq, na, ncoeff);
      
      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    }
        
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "Memory NEQ:   " << NVAR2D*NVAR2D*neq*sizeof(Tv)/1000000.0f
	 << " MB" << " NEQ=" << neq << endl;
    cout << "Elapsed time: " << elapsedTime << " ms" << endl;
    cout << "Bandwidth:    " << (NVAR2D*neq*sizeof(Tv)+
				 3*NVAR2D*NVAR2D*neq*sizeof(Tv)+
				 2*neq*sizeof(Ti))/1000000000.0f/elapsedTime*1000.0f
	 << " GB/s" << endl;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    coproc_checkErrors("hydro_calcMatDiag2d_cuda");
    return 0;
  };

  /*****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    bool blumping>
  inline
  int hydro_calcMatGalMatD2d_cuda(__SIZET *d_CoeffsAtEdge,
				  __SIZET *d_IedgeList,
				  __SIZET *d_vec,
				  __SIZET *d_mat,
				  Tm scale,
				  Ti nblocks,
				  Ti neq,
				  Ti na,
				  Ti nedge,
				  Ti ncoeff,
				  Ti nedges,
				  Ti iedgeset,
				  cudaStream_t stream=0)
  {
    Tv  *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);

    // The total number of blocks depends on the problem size NEQ and
    // on the number of compute threads per cooperative thread block
    int blocks = (nedges+COMPUTE_THREADS_PER_CTA-1)/COMPUTE_THREADS_PER_CTA;
    dim3 grid(blocks, 1, 1);

    // The total number of threads per cooperative thread block
    // depends on the number of compute threads plus the number od DMA
    // threads multiplied by the number of load/store operations
    dim3 block( COMPUTE_THREADS_PER_CTA + 1*DMA_THREADS_PER_LD, 1, 1);
    
    cout << "hydro_calcMatGalMatD2d_cuda" 
	 << " nblocks=" << nblocks
	 << " nedges=" << nedges
	 << " grid=" << grid.x << " block=" << block.x << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);

      cudaEventRecord(start,stream);

      hydro_calcMatrixMatD2d_knl
	<Tc,Tv,Tm,Ti,SYSTEM_SCALAR,DISSIPATION_ZERO,blumping>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     vec, mat, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);
      
      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);  
    } else {
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      __SIZET cmemPool[NVAR2D];
#pragma unroll
      for (int i=0; i<NVAR2D; i++)
	cmemPool[i] = d_mat[i*(NVAR2D+1)];
      
      cudaEventRecord(start,stream);
      
      cudaMemcpyToSymbolAsync("constMemPool", cmemPool,
			      sizeof(__SIZET)*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");

      hydro_calcMatrixMatD2d_knl
	<Tc,Tv,Tm,Ti,SYSTEM_BLOCK,DISSIPATION_ZERO,blumping>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     vec, mat, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);

      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    }

    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "Memory NEDGE: " << NVAR2D*nedges*sizeof(Tv)/1000000.0f
	 << " MB" << " NEDGE=" << nedges << endl;
    cout << "Elapsed time: " << elapsedTime << " ms" << endl;
    cout << "Bandwidth:    " << (6*nedges*sizeof(Ti)+
				 3*NVAR2D*nedges*sizeof(Tv))/1000000000.0f/elapsedTime*1000.0f
	 << " GB/s" << endl;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    coproc_checkErrors("hydro_calcMatGalMatD2d_cuda");
    return 0;
  };

  /*****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    bool blumping>
  inline
  int hydro_calcMatGalerkin2d_cuda(__SIZET *d_CoeffsAtEdge,
				   __SIZET *d_IedgeList,
				   __SIZET *d_vec,
				   __SIZET *d_mat,
				   Tm scale,
				   Ti nblocks,
				   Ti neq,
				   Ti na,
				   Ti nedge,
				   Ti ncoeff,
				   Ti nedges,
				   Ti iedgeset,
				   cudaStream_t stream=0)
  {
    Tv  *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);

    // The total number of blocks depends on the problem size NEQ and
    // on the number of compute threads per cooperative thread block
    int blocks = (nedges+COMPUTE_THREADS_PER_CTA-1)/COMPUTE_THREADS_PER_CTA;
    dim3 grid(blocks, 1, 1);

    // The total number of threads per cooperative thread block
    // depends on the number of compute threads plus the number od DMA
    // threads multiplied by the number of load/store operations
    dim3 block( COMPUTE_THREADS_PER_CTA + 1*DMA_THREADS_PER_LD, 1, 1);
    
    cout << "hydro_calcMatGalerkin2d_cuda" 
	 << " nblocks=" << nblocks
	 << " nedges=" << nedges
	 << " grid=" << grid.x << " block=" << block.x << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);

      cudaEventRecord(start,stream);

      hydro_calcMatrix2d_knl
	<Tc,Tv,Tm,Ti,SYSTEM_SCALAR,DISSIPATION_ZERO,blumping>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     vec, mat, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);

      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    } else {
      cudaEventRecord(start,stream);
      
      // Matrix is stored in block format, that is, the data of each
      // scalar submatrix resides in an individual device memory
      // block; thus we transfer the starting addresses of each memory
      // block into constant device memory and pass a dummy argument
      cudaMemcpyToSymbolAsync("constMemPool", d_mat,
			      sizeof(__SIZET)*NVAR2D*NVAR2D, 0,
			      cudaMemcpyHostToDevice,
			      stream);
      
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");

      hydro_calcMatrix2d_knl
	<Tc,Tv,Tm,Ti,SYSTEM_BLOCK,DISSIPATION_ZERO,blumping>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     vec, mat, scale,
				     neq, na, nedge, ncoeff,
				     nedges, iedgeset);

      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    }
    
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "Memory NEDGE: " << NVAR2D*NVAR2D*nedges*sizeof(Tv)/1000000.0f
	 << " MB" << " NEDGE=" << nedges << endl;
    cout << "Elapsed time: " << elapsedTime << " ms" << endl;
    cout << "Bandwidth:    " << (6*nedges*sizeof(Ti)+
				 2*NVAR2D*nedges*sizeof(Tv)+
				 2*NVAR2D*NVAR2D*nedges*sizeof(Tv))/1000000000.0f/elapsedTime*1000.0f
	 << " GB/s" << endl;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    coproc_checkErrors("hydro_calcMatGalerkin2d_cuda");
    return 0;
  };
  
  /*****************************************************************************
   * External C functions which can be called from the Fortran code
   *****************************************************************************/

  extern "C" {
    __INT FNAME(hydro_calcmatdiagmatd2d_cuda)(__SIZET *d_CoeffsAtDiag,
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
      return (__INT) hydro_calcMatDiagMatD2d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtDiag, d_IdiagList, d_vec, d_mat,
			       *scale, *nblocks, *neq, *na, *ncoeff,
			       (cudaStream_t)(*stream));
    }

    /***************************************************************************/
    
    __INT FNAME(hydro_calcmatdiag2d_cuda)(__SIZET *d_CoeffsAtDiag,
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
      return (__INT) hydro_calcMatDiag2d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtDiag, d_IdiagList, d_vec, d_mat,
			       *scale, *nblocks, *neq, *na, *ncoeff,
			       (cudaStream_t)(*stream));
    }

    /***************************************************************************/
    
    __INT FNAME(hydro_calcmatgalmatd2d_cuda)(__SIZET *d_CoeffsAtEdge,
					     __SIZET *d_IedgeList,
					     __SIZET *d_vec,
					     __SIZET *d_mat,
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
	  <__DP,__DP,__DP,__INT,false>(d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
				       *scale, *nblocks, *neq, *na, *nedge,
				       *ncoeff, *nedges, *iedgeset,
				       (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatGalMatD2d_cuda
	  <__DP,__DP,__DP,__INT,true>(d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
				      *scale, *nblocks, *neq, *na, *nedge,
				      *ncoeff, *nedges, *iedgeset,
				      (cudaStream_t)(*stream));
    }
    
    /***************************************************************************/
    
    __INT FNAME(hydro_calcmatgalerkin2d_cuda)(__SIZET *d_CoeffsAtEdge,
					      __SIZET *d_IedgeList,
					      __SIZET *d_vec,
					      __SIZET *d_mat,
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
	  <__DP,__DP,__DP,__INT,false>(d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
				       *scale, *nblocks, *neq, *na, *nedge,
				       *ncoeff, *nedges, *iedgeset,
				       (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatGalerkin2d_cuda
	  <__DP,__DP,__DP,__INT,true>(d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
				      *scale, *nblocks, *neq, *na, *nedge,
				      *ncoeff, *nedges, *iedgeset,
				      (cudaStream_t)(*stream));
    }
  };    
}
