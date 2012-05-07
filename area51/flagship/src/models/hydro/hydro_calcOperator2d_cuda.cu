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
#include <coproc_core.h>
#include <coproc_storage_cuda.h>
#include "../../cudaDMA.h"
#include "../../cudaGatherScatter.h"

#define LANGUAGE LANGUAGE_C
#include "../../flagship.h"
#include "../../cudaMacros.h"

#define HYDRO_NDIM 2
#include "hydro.h"

// Memory pool in constant device memory
__device__ __constant__ __SIZET constMemPool[NVAR2D*NVAR2D];


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

  template<>
  struct FluxJacobiMatrixBase<SYSTEM_SEGREGATED>
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
			     Ti ipos,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtDiag,2,ipos,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtDiag,3,ipos,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtDiag,4,ipos,NVAR2D,neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
#else
      // Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX11(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtDiag,2,ipos,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX22(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtDiag,3,ipos,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX33(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
      IDX2(MatrixAtDiag,4,ipos,NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX44(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,_);
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
			     Ti ipos,
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);

      // Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
#else
      // Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR2D,2,1,) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,_);
      
      // Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,_);
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
    static void calcNodeData(Td *MatrixAtDiag,
			     Tc *CoeffsAtDiag,
			     Td scale,
			     Td ui,
			     Td vi,
			     Td Ei,
			     Ti ipos,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,2,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX21(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,3,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX31(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,4,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX41(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,5,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX12(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,6,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,7,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX32(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,8,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX42(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,9,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX13(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,10,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX23(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,11,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,12,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX43(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,13,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX14(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,14,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX24(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,15,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX34(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,16,ipos,NVAR2D*NVAR2D,neqsim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
#else
      // Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX2(MatrixAtDiag,1,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX11(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,2,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX21(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,3,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX31(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,4,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX41(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,5,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX12(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,6,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX22(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,7,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX32(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,8,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX42(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,9,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX13(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,10,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX23(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,11,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX33(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,12,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX43(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,13,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX14(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,14,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX24(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,15,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX34(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
      IDX2(MatrixAtDiag,16,ipos,NVAR2D*NVAR2D,neqsim) = -
	FLUXJACOBIMATRIX44(scale,
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,1,ieq,ncoeff,neq),
			   IDX2_COEFFSATDIAG(CoeffsAtDiag,2,ieq,ncoeff,neq),ui,vi,Ei);
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
#ifdef HYDRO_USE_IBP
      // Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX21(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX31(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX41(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,5,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX12(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,6,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,7,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX32(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,8,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX42(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,9,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX13(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,10,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX23(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,11,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,12,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX43(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,13,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX14(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,14,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX24(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,15,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX34(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,16,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);

      // Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX21(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX31(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX41(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,5,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX12(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,6,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX22(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,7,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX32(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,8,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX42(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,9,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX13(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,10,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX23(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,11,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,12,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX43(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,13,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX14(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,14,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX24(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,15,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX34(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,16,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
#else
      // Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(MatrixAtEdge,1,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,2,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX21(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,3,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX31(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,4,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX41(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,5,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX12(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,6,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,7,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX32(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,8,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX42(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,9,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX13(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,10,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX23(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,11,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,12,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX43(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,13,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX14(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,14,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX24(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,15,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX34(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      IDX3(MatrixAtEdge,16,1,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),uj,vj,Ej);
      
      // Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(MatrixAtEdge,1,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX11(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,2,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX21(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,3,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX31(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,4,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX41(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,5,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX12(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,6,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX22(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,7,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX32(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,8,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX42(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,9,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX13(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,10,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX23(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,11,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX33(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,12,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX43(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,13,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX14(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,14,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX24(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,15,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX34(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
      IDX3(MatrixAtEdge,16,2,ipos,NVAR2D*NVAR2D,(bstabilise ? 3 : 2),nedgesim) =
	FLUXJACOBIMATRIX44(-scale,
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),
			   IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),ui,vi,Ei);
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
    static void calcNodeData(Td *MatrixAtDiag,
			     Tc *CoeffsAtDiag,
			     Td scale,
			     Td ui,
			     Td vi, 
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
      FluxJacobiMatrixBase<isystemcoupling>::calcNodeData<1>
	(MatrixAtDiag,CoeffsAtDiag,scale,ui,vi,1,ieq,neq,ncoeff);
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
			     Td Ei,
			     Ti ieq,
			     Ti neq,
			     Ti ncoeff)
    {
      FluxJacobiMatrixBase<isystemcoupling>::calcNodeData<1>
	(MatrixAtDiag,CoeffsAtDiag,scale,ui,vi,Ei,1,ieq,neq,ncoeff);
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
	(MatrixAtEdge,CoeffsAtEdge,scale,ui,uj,vi,vj,1,iedge,nedge,ncoeff);
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
	(MatrixAtEdge,CoeffsAtEdge,scale,ui,uj,vi,vj,Ei,Ej,1,iedge,nedge,ncoeff);
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
#pragma unroll
      for (int i=1; i<=NVAR2D; i++)
	IDX3(MatrixAtEdge,i,1,ipos,NVAR2D,3,nedgesim) = 0.0;
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
#pragma unroll
      for (int i=1; i<=NVAR2D*NVAR2D; i++)
	IDX3(MatrixAtEdge,i,1,ipos,NVAR2D*NVAR2D,3,nedgesim) = 0.0;
    }
  };

  /*****************************************************************************
   * DissipationBase: Specialization for block-diagonal matrix
   * computing scalar artificial dissipation proportional to the
   * spectral radius of the Roe matrix
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_SEGREGATED,DISSIPATION_SCALAR>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };
  
  /*****************************************************************************
   * DissipationBase: Specialization for full matrix
   * computing scalar artificial dissipation proportional to the
   * spectral radius of the Roe matrix
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_ALLCOUPLED,DISSIPATION_SCALAR>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };

  /*****************************************************************************
   * DissipationBase: Specialization for block-diagonal matrix
   * computing tensorial artificial dissipation of Roe-type
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_SEGREGATED,DISSIPATION_ROE>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };
  
  /*****************************************************************************
   * DissipationBase: Specialization for full matrix
   * computing tensorial artificial dissipation of Roe-type
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_ALLCOUPLED,DISSIPATION_ROE>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };

  /*****************************************************************************
   * DissipationBase: Specialization for block-diagonal matrix
   * computing scalar artificial dissipation of Rusanov-type
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_SEGREGATED,DISSIPATION_RUSANOV>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

    }
  };
  
  /*****************************************************************************
   * DissipationBase: Specialization for full matrix
   * computing scalar artificial dissipation of Rusanov-type
   ****************************************************************************/

  template <>  
  struct DissipationBase<SYSTEM_ALLCOUPLED,DISSIPATION_RUSANOV>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {

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
	(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,scale,ui,uj,vi,vj,1,iedge,nedge,ncoeff);
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
	(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,scale,ui,uj,vi,vj,Ei,Ej,1,iedge,nedge,ncoeff);
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
	Tm DataAtDiag[NVAR2D];
	
	// Get solution values at node
	Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  gatherNodeData<true>(DataAtDiag,vec,ieq,neq);
	
	// Compute velocities
	Tm ui = XVELOCITY2(DataAtDiag,IDX2,1,NVAR2D,1);
	Tm vi = YVELOCITY2(DataAtDiag,IDX2,1,NVAR2D,1);
	
	// Compute Galerkin coefficient $K_ii$
	FluxJacobiMatrix<SYSTEM_SEGREGATED>::
	  calcNodeData(DataAtDiag,CoeffsAtDiag,scale,ui,vi,ieq,neq,ncoeff);
	
	// Get diagonal position in the global matrix
	Ti ia  = IDX2_DIAGLIST(IdiagList,2,idx+1,2,neq);
	
	// Build coefficients into global operator
	Matrix<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  scatterNodeData<true>(mat,DataAtDiag,ia,na);
      }
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the diagonal entries of the
   * block-diagonal global operator (cudaDMA implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __global__ void hydro_calcMatDiagMatD2d_cudaDMA(Tc *CoeffsAtDiag,
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
    // Not implemented yet
    printf("Not implemented\n");
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
	Tm DataAtDiag[NVAR2D*NVAR2D];
	
	// Get solution values at node
	Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  gatherNodeData<true>(DataAtDiag,vec,ieq,neq);
	
	// Compute velocities and energy
	Tm ui = XVELOCITY2(DataAtDiag,IDX2,1,NVAR2D*NVAR2D,1);
	Tm vi = YVELOCITY2(DataAtDiag,IDX2,1,NVAR2D*NVAR2D,1);
	Tm Ei = SPECIFICTOTALENERGY2(DataAtDiag,IDX2,1,NVAR2D*NVAR2D,1);
	
	// Compute Galerkin coefficient $K_ii$
	FluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	  calcNodeData(DataAtDiag,CoeffsAtDiag,scale,ui,vi,Ei,ieq,neq,ncoeff);
	
	// Get diagonal position in the global matrix
	Ti ia  = IDX2_DIAGLIST(IdiagList,2,idx+1,2,neq);
	
	// Build coefficients into global operator
	Matrix<NVAR2D*NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  scatterNodeData<true>(mat,DataAtDiag,ia,na);
      }
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the diagonal entries of the
   * full global operator (cudaDMA implementation).
   ****************************************************************************/
  
  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __global__ void hydro_calcMatDiag2d_cudaDMA(Tc *CoeffsAtDiag,
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
    // Not implemented yet
    printf("Not implemented\n");
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
	Tm DataAtEdge[2*NVAR2D];
	
	// Get solution values at edge endpoints
	Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	  gatherEdgeData<true>(DataAtEdge,vec,i,j,neq);
	
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
				scale,ui,uj,vi,vj,idx+1,nedge,ncoeff);
	  
	  // Build matrix coefficients into global operator
	  Matrix<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    scatterEdgeData<true,false,blumping>(mat,MatrixAtEdge,
						 IedgeList,idx+1,na,nedge);
	  
	} else {
	  
	  // Local matrix data at edge from local memory
	  Tm MatrixAtEdge[3*NVAR2D];
	  
	  // Compute Galerkin coefficient $K_ij$ and $K_ji$
	  FluxJacobiMatrix<SYSTEM_SEGREGATED>::
	    calcEdgeData<true>(MatrixAtEdge,CoeffsAtEdge,
			       scale,ui,uj,vi,vj,idx+1,nedge,ncoeff);
	  
	  // Compute contribution of artificial diffusion
	  Dissipation<SYSTEM_SEGREGATED,idissipation>::
	    calcEdgeData(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
			 scale,ui,uj,vi,vj,idx+1,nedge,ncoeff);
	  
	  // Build matrix coefficients into global operator
	  Matrix<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    scatterEdgeData<true,true,blumping>(mat,MatrixAtEdge,
						IedgeList,idx+1,na,nedge);
	}
      }
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the off-diagonal entries of the
   * block-diagonal global operator and assembles the artificial
   * dissipation tensor if required (cudaDMA implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int idissipation,
	    bool blumping,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __global__ void hydro_calcMatrixMatD2d_cudaDMA(Tc *CoeffsAtEdge,
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
    // Not implemented yet
    printf("Not implemented\n");
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
	  Tm DataAtEdge[2*NVAR2D];
	  
	  // Get solution values at edge endpoints
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    gatherEdgeData<true>(DataAtEdge,vec,i,j,neq);
	  
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
				  scale,ui,uj,vi,vj,Ei,Ej,idx+1,nedge,ncoeff);
	    
	    // // Build matrix coefficients into global operator
	    Matrix<NVAR2D*NVAR2D,isystemformat==SYSTEM_BLOCK>::
	      scatterEdgeData<true,false,blumping>(mat,MatrixAtEdge,
						   IedgeList,idx+1,na,nedge);
	  } else {
	    
	    // Local matrix data at edge from local memory
	    Tm MatrixAtEdge[3*NVAR2D*NVAR2D];
	    
	    // Compute Galerkin coefficient $K_ij$ and $K_ji$
	    FluxJacobiMatrix<SYSTEM_ALLCOUPLED>::
	      calcEdgeData<true>(MatrixAtEdge,CoeffsAtEdge,
				 scale,ui,uj,vi,vj,Ei,Ej,idx+1,nedge,ncoeff);
	    
	    // Compute contribution of artificial diffusion
	    Dissipation<SYSTEM_ALLCOUPLED,idissipation>::
	      calcEdgeData(MatrixAtEdge,CoeffsAtEdge,DataAtEdge,
			   scale,ui,uj,vi,vj,Ei,Ej,idx+1,nedge,ncoeff);
	    
	    // Build matrix coefficients into global operator
	    Matrix<NVAR2D*NVAR2D,isystemformat==SYSTEM_BLOCK>::
	      scatterEdgeData<true,true,blumping>(mat,MatrixAtEdge,
						  IedgeList,idx+1,na,nedge);
	  }
	}
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the off-diagonal entries of the full
   * global operator and assembles the artificial dissipation tensor
   * if required (cudaDMA implementation).
   ****************************************************************************/

  template <typename Tc,
	    typename Tv,
	    typename Tm,
	    typename Ti,
	    int isystemformat,
	    int idissipation,
	    bool blumping,
	    int compute_threads_per_cta,
	    int dma_threads_per_ld>
  __global__ void hydro_calcMatrix2d_cudaDMA(Tc *CoeffsAtEdge,
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
    // Not implemented yet
    printf("Not implemented\n");
  };
  
  /*****************************************************************************
   * Internal C++ functions which invoke the CUDA kernels
   ****************************************************************************/

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
    // Strategy: run the largest possible number of blocks with a
    // predefined number of compute/dma threads per block and let each
    // compute thread process the minimal number of equations
    const int compute_threads_per_cta  = 32*0;
    const int dma_threads_per_ld       = 32*1;
    const int dma_lds                  = 1;
    const int neq_per_thread_cudaDMA   = 1;

    const int threads_per_cta_baseline = 32*1;
    const int neq_per_thread_baseline  = 1;
    
    int blocks, threads, neq_cudaDMA, neq_baseline;
    prepare_cudaDMA(neq, neq_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &neq_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(neq-neq_cudaDMA, neq_per_thread_baseline,
		     threads_per_cta_baseline, &blocks, &threads, &neq_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);
    
    Tv *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
        
    cout << "hydro_calcMatDiagMatD2d_cuda" 
	 << " nblocks=" << nblocks
	 << " neq=" << neq
	 << " CudaDMA:"
	 << " #blocks=" << grid_cudaDMA.x << ","
	 << grid_cudaDMA.y << "," << grid_cudaDMA.z 
	 << " #threads per block=" << block_cudaDMA.x 
	 << "," << block_cudaDMA.y << "," << block_cudaDMA.z 
	 << " Baseline:"
	 << " #blocks=" << grid_baseline.x << "," 
	 << grid_baseline.y << "," << grid_baseline.z 
	 << " #threads per block=" << block_baseline.x 
	 << "," << block_baseline.y << "," << block_baseline.z 
	 << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);
    
      cudaThreadSynchronize();
      cudaEventRecord(start,stream);
      
      if (grid_cudaDMA.x>0)
      	// CudaDMA implementation
      	hydro_calcMatDiagMatD2d_cudaDMA
      	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,compute_threads_per_cta,dma_threads_per_ld>
      	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
      						       IdiagList,
      						       vec, mat, scale,
      						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
      if (grid_baseline.x>0)
      	// Baseline implementation
      	hydro_calcMatDiagMatD2d_baseline
      	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR>
      	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
      							 IdiagList,
      							 vec, mat, scale,
      							 neq, na, ncoeff,
							 neq, 
							 neq_per_thread_baseline,
							 neq_cudaDMA);
      
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
      
      cudaMemcpyToSymbolAsync("constMemPool", cmemPool,
      			      sizeof(__SIZET)*NVAR2D, 0,
      			      cudaMemcpyHostToDevice,
      			      stream);
      Tm *mat;
      cudaGetSymbolAddress(((void**)&mat), "constMemPool");
      
      cudaThreadSynchronize();
      cudaEventRecord(start,stream);

      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatDiagMatD2d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
						       IdiagList,
						       vec, mat, scale,
						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatDiagMatD2d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
							 IdiagList,
							 vec, mat, scale,
							 neq, na, ncoeff,
							 neq,
							 neq_per_thread_baseline,
							 neq_cudaDMA);
      
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

  /****************************************************************************/

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
    // Strategy: run the largest possible number of blocks with a
    // predefined number of compute/dma threads per block and let each
    // compute thread process the minimal number of equations
    const int compute_threads_per_cta  = 32*0;
    const int dma_threads_per_ld       = 32*1;
    const int dma_lds                  = 1;
    const int neq_per_thread_cudaDMA   = 1;

    const int threads_per_cta_baseline = 32*1;
    const int neq_per_thread_baseline  = 1;
    
    int blocks, threads, neq_cudaDMA, neq_baseline;
    prepare_cudaDMA(neq, neq_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &neq_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(neq-neq_cudaDMA, neq_per_thread_baseline,
		     threads_per_cta_baseline, &blocks, &threads, &neq_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);

    Tv *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtDiag = (Tc*)(*d_CoeffsAtDiag);
    Ti *IdiagList = (Ti*)(*d_IdiagList);
  
    cout << "hydro_calcMatDiag2d_cuda" 
	 << " nblocks=" << nblocks
	 << " neq=" << neq
	 << " CudaDMA:"
	 << " #blocks=" << grid_cudaDMA.x << ","
	 << grid_cudaDMA.y << "," << grid_cudaDMA.z 
	 << " #threads per block=" << block_cudaDMA.x 
	 << "," << block_cudaDMA.y << "," << block_cudaDMA.z 
	 << " Baseline:"
	 << " #blocks=" << grid_baseline.x << "," 
	 << grid_baseline.y << "," << grid_baseline.z 
	 << " #threads per block=" << block_baseline.x 
	 << "," << block_baseline.y << "," << block_baseline.z 
	 << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);
      
      cudaThreadSynchronize();
      cudaEventRecord(start,stream);

      if (grid_cudaDMA.x>0)
	//CudaDMA implementation
	hydro_calcMatDiag2d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
						       IdiagList,
						       vec, mat, scale,
						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatDiag2d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
							 IdiagList,
							 vec, mat, scale,
							 neq, na, ncoeff,
							 neq,
							 neq_per_thread_baseline,
							 neq_cudaDMA);
      
      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    } else {
      cudaThreadSynchronize();
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

      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatDiag2d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtDiag,
						       IdiagList,
						       vec, mat, scale,
						       neq, na, ncoeff,
						       neq_cudaDMA,
						       neq_per_thread_cudaDMA);
      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatDiag2d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtDiag,
							 IdiagList,
							 vec, mat, scale,
							 neq, na, ncoeff,
							 neq,
							 neq_per_thread_baseline,
							 neq_cudaDMA);
      
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
	    int idissipationtype,
	    bool blumping>
  inline
  int hydro_calcMatrixMatD2d_cuda(__SIZET *d_CoeffsAtEdge,
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
    // Strategy: run the largest possible number of blocks with a
    // predefined number of compute/dma threads per block and let each
    // compute thread process the minimal number of edges
    const int compute_threads_per_cta  = 32*0;
    const int dma_threads_per_ld       = 32*1;
    const int dma_lds                  = 1;
    const int nedge_per_thread_cudaDMA = 1;

    const int threads_per_cta_baseline  = 32*1;
    const int nedge_per_thread_baseline = 1;
    
    int blocks, threads, nedge_cudaDMA, nedge_baseline;
    prepare_cudaDMA(nedgeset, nedge_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &nedge_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(nedgeset-nedge_cudaDMA, nedge_per_thread_baseline,
		     threads_per_cta_baseline, &blocks, &threads, &nedge_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);
    
    Tv  *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);

    cout << "hydro_calcMatrixMatD2d_cuda" 
	 << " nblocks=" << nblocks
	 << " nedgeset=" << nedgeset
	 << " CudaDMA:"
	 << " #blocks=" << grid_cudaDMA.x << ","
	 << grid_cudaDMA.y << "," << grid_cudaDMA.z 
	 << " #threads per block=" << block_cudaDMA.x 
	 << "," << block_cudaDMA.y << "," << block_cudaDMA.z 
	 << " Baseline:"
	 << " #blocks=" << grid_baseline.x << "," 
	 << grid_baseline.y << "," << grid_baseline.z 
	 << " #threads per block=" << block_baseline.x 
	 << "," << block_baseline.y << "," << block_baseline.z 
	 << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  
    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);

      cudaEventRecord(start,stream);

      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrixMatD2d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatrixMatD2d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
      
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

      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrixMatD2d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,idissipationtype,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
      if (grid_baseline.x>0)
	hydro_calcMatrixMatD2d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,idissipationtype,blumping>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1,
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
      
      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    }

    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "Memory NEDGE: " << NVAR2D*nedgeset*sizeof(Tv)/1000000.0f
	 << " MB" << " NEDGE=" << nedgeset << endl;
    cout << "Elapsed time: " << elapsedTime << " ms" << endl;
    cout << "Bandwidth:    " << (6*nedgeset*sizeof(Ti)+
				 3*NVAR2D*nedgeset*sizeof(Tv))/1000000000.0f/elapsedTime*1000.0f
	 << " GB/s" << endl;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    coproc_checkErrors("hydro_calcMatrixMatD2d_cuda");
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
  int hydro_calcMatrix2d_cuda(__SIZET *d_CoeffsAtEdge,
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
    // Strategy: run the largest possible number of blocks with a
    // predefined number of compute/dma threads per block and let each
    // compute thread process the minimal number of edges
    const int compute_threads_per_cta  = 32*0;
    const int dma_threads_per_ld       = 32*1;
    const int dma_lds                  = 1;
    const int nedge_per_thread_cudaDMA = 1;

    const int threads_per_cta_baseline  = 32*1;
    const int nedge_per_thread_baseline = 1;
    
    int blocks, threads, nedge_cudaDMA, nedge_baseline;
    prepare_cudaDMA(nedgeset, nedge_per_thread_cudaDMA,
		    compute_threads_per_cta, dma_threads_per_ld,
		    dma_lds, &blocks, &threads, &nedge_cudaDMA);
    dim3 grid_cudaDMA(blocks, 1, 1);
    dim3 block_cudaDMA(threads, 1, 1);

    prepare_baseline(nedgeset-nedge_cudaDMA, nedge_per_thread_baseline,
		     threads_per_cta_baseline, &blocks, &threads, &nedge_baseline);
    dim3 grid_baseline(blocks, 1, 1);
    dim3 block_baseline(threads, 1, 1);

    Tv  *vec = (Tv*)(*d_vec);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);

    cout << "hydro_calcMatrix2d_cuda" 
	 << " nblocks=" << nblocks
	 << " nedgeset=" << nedgeset
	 << " CudaDMA:"
	 << " #blocks=" << grid_cudaDMA.x << ","
	 << grid_cudaDMA.y << "," << grid_cudaDMA.z 
	 << " #threads per block=" << block_cudaDMA.x 
	 << "," << block_cudaDMA.y << "," << block_cudaDMA.z 
	 << " Baseline:"
	 << " #blocks=" << grid_baseline.x << "," 
	 << grid_baseline.y << "," << grid_baseline.z 
	 << " #threads per block=" << block_baseline.x 
	 << "," << block_baseline.y << "," << block_baseline.z 
	 << endl;

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    if (nblocks == 1) {
      // Matrix is store in interleaved matrix so that all matrix data
      // are stored contiguously in one single device memory block
      Tm *mat = (Tm*)(*d_mat);

      cudaEventRecord(start,stream);

      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrix2d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatrix2d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_SCALAR,idissipationtype,blumping>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
      
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

      if (grid_cudaDMA.x>0)
	// CudaDMA implementation
	hydro_calcMatrix2d_cudaDMA
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,DISSIPATION_ZERO,blumping,
	   compute_threads_per_cta,dma_threads_per_ld>
	  <<<grid_cudaDMA, block_cudaDMA, 0, stream>>>(CoeffsAtEdge,
						       IedgeList,
						       vec, mat, scale,
						       neq, na, nedge, ncoeff,
						       nedge_cudaDMA+iedgeset-1, 
						       nedge_per_thread_cudaDMA,
						       iedgeset-1);
      
      if (grid_baseline.x>0)
	// Baseline implementation
	hydro_calcMatrix2d_baseline
	  <Tc,Tv,Tm,Ti,SYSTEM_BLOCK,DISSIPATION_ZERO,blumping>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vec, mat, scale,
							 neq, na, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
      
      cudaEventRecord(stop,stream);
      cudaEventSynchronize(stop);
    }
    
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "Memory NEDGE: " << NVAR2D*NVAR2D*nedgeset*sizeof(Tv)/1000000.0f
	 << " MB" << " NEDGE=" << nedgeset << endl;
    cout << "Elapsed time: " << elapsedTime << " ms" << endl;
    cout << "Bandwidth:    " << (6*nedgeset*sizeof(Ti)+
				 2*NVAR2D*nedgeset*sizeof(Tv)+
				 2*NVAR2D*NVAR2D*nedgeset*sizeof(Tv))/1000000000.0f/elapsedTime*1000.0f
	 << " GB/s" << endl;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    coproc_checkErrors("hydro_calcMatrix2d_cuda");
    return 0;
  };
  
  /*****************************************************************************
   * External C functions which can be called from the Fortran code
   ****************************************************************************/

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

    /**************************************************************************/
    
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

    /**************************************************************************/
    
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
					     __INT *nedgeset,
					     __INT *iedgeset,
					     __INT *cconstrType,
					     __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }
    
    /**************************************************************************/
    
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
					      __INT *nedgeset,
					      __INT *iedgeset,
					      __INT *cconstrType,
					      __I64 *stream)
    {
      if (*cconstrType == 0)
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ZERO,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatscdissmatd2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatscdiss2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_SCALAR,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatroedissmatd2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatroediss2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_ROE,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatrusdissmatd2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrixMatD2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }

    /**************************************************************************/
    
    __INT FNAME(hydro_calcmatrusdiss2d_cuda)(__SIZET *d_CoeffsAtEdge,
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
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,false>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
      else
	return (__INT) hydro_calcMatrix2d_cuda
	  <__DP,__DP,__DP,__INT,DISSIPATION_RUSANOV,true>
	  (d_CoeffsAtEdge, d_IedgeList, d_vec, d_mat,
	   *scale, *nblocks, *neq, *na, *nedge,
	   *ncoeff, *nedgeset, *iedgeset,
	   (cudaStream_t)(*stream));
    }
  };
}
