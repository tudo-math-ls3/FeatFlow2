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
    template <int nedgesim, typename Td, typename Ti>
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
			     Ti ipos)
    {
      // Compute the Galerkin fluxes for x-direction
      IDX2(Fxi,1,ipos,NVAR2D,nedgesim) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi);
      IDX2(Fxi,2,ipos,NVAR2D,nedgesim) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi);
      IDX2(Fxi,3,ipos,NVAR2D,nedgesim) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi);
      IDX2(Fxi,4,ipos,NVAR2D,nedgesim) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi);
      
      IDX2(Fxj,1,ipos,NVAR2D,nedgesim) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      IDX2(Fxj,2,ipos,NVAR2D,nedgesim) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      IDX2(Fxj,3,ipos,NVAR2D,nedgesim) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      IDX2(Fxj,4,ipos,NVAR2D,nedgesim) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      
      // Compute Galerkin fluxes for y-direction
      IDX2(Fyi,1,ipos,NVAR2D,nedgesim) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi);
      IDX2(Fyi,2,ipos,NVAR2D,nedgesim) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi);
      IDX2(Fyi,3,ipos,NVAR2D,nedgesim) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi);
      IDX2(Fyi,4,ipos,NVAR2D,nedgesim) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi);
      
      IDX2(Fyj,1,ipos,NVAR2D,nedgesim) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
      IDX2(Fyj,2,ipos,NVAR2D,nedgesim) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
      IDX2(Fyj,3,ipos,NVAR2D,nedgesim) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
      IDX2(Fyj,4,ipos,NVAR2D,nedgesim) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
    }
    
    template <int nedgesim, typename Td, typename Ti>
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
			     Ti ipos)
    {
      // Compute Galerkin flux difference for x-direction
      IDX2(Fx_ij,1,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi)-
	INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      IDX2(Fx_ij,2,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi)-
	INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      IDX2(Fx_ij,3,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi)-
	INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      IDX2(Fx_ij,4,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,ui,pi)-
	INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,uj,pj);
      
      // Compute Galerkin flux difference for y-direction
      IDX2(Fy_ij,1,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi)-
	INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
      IDX2(Fy_ij,2,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi)-
	INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
      IDX2(Fy_ij,3,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi)-
	INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
      IDX2(Fy_ij,4,ipos,NVAR2D,nedgesim) =
	INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim,vi,pi)-
	INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim,vj,pj);
    }
  };

  /*****************************************************************************
   * InviscidFlux
   ****************************************************************************/

  struct InviscidFlux : public InviscidFluxBase
  {
    // Enable use of inherited functions
    using InviscidFluxBase::calcEdgeData;

    /**************************************************************************
     * Wrapper routine for processing a single edge
     *************************************************************************/
    template <typename Td>
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
      InviscidFluxBase::calcEdgeData<1>
	(Fxi,Fxj,Fyi,Fyj,DataAtEdge,ui,uj,vi,vj,pi,pj,1);
    }

    template <typename Td>
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
      InviscidFluxBase::calcEdgeData<1>
	(Fx_ij,Fy_ij,DataAtEdge,ui,uj,vi,vj,pi,pj,1);
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
    template <int nedgesim, typename Tc, typename Td, typename Ti>
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
			     Ti ipos,
			     Ti iedge,
			     Ti nedge,
			     Ti ncoeff)
    {
#pragma unroll
      for (int i=1; i<=NVAR2D; i++)
	IDX2(VectorAtEdge,i,ipos,NVAR2D,nedgesim) = 0.0;
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
    template <int nedgesim, typename Tc, typename Td, typename Ti>
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
			     Ti ipos,
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
    
      // Compute densities
      Td ri = DENSITY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim);
      Td rj = DENSITY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim);
    
      // Compute enthalpies
      Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim)+pi)/ri;
      Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim)+pj)/rj;
    
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
#pragma unroll
      for (int i=1; i<=NVAR2D; i++)
	IDX2(VectorAtEdge,i,ipos,NVAR2D,nedgesim) =
	  d_ij*(IDX3(DataAtEdge,i,2,ipos,NVAR2D,2,nedgesim)
	       -IDX3(DataAtEdge,i,1,ipos,NVAR2D,2,nedgesim));
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
    template <int nedgesim, typename Tc, typename Td, typename Ti>
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
			     Ti ipos,
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
    
      // Compute densities
      Td ri = DENSITY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim);
      Td rj = DENSITY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim);
    
      // Compute enthalpies
      Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim)+pi)/ri;
      Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim)+pj)/rj;
    
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
#pragma unroll
      for (int i=1; i<=NVAR2D; i++)
	IDX2(VectorAtEdge,i,ipos,NVAR2D,nedgesim) =
	  d_ij*(IDX3(DataAtEdge,i,2,ipos,NVAR2D,2,nedgesim)
	       -IDX3(DataAtEdge,i,1,ipos,NVAR2D,2,nedgesim));
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing
   * tensorial artificial dissipation of Roe-type.
   ****************************************************************************/

  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_ROE>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
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
			     Ti ipos,
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
      
	// Compute densities
	Td ri = DENSITY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim);
	Td rj = DENSITY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim);
      
	// Compute enthalpies
	Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim)+pi)/ri;
	Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim)+pj)/rj;
      
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
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  Diff[i-1] = IDX3(DataAtEdge,i,2,ipos,NVAR2D,2,nedgesim)
	             -IDX3(DataAtEdge,i,1,ipos,NVAR2D,2,nedgesim);
	      
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
	IDX2(VectorAtEdge,1,ipos,NVAR2D,nedgesim) = anorm * ( w1 + w2 + w3 );
	IDX2(VectorAtEdge,2,ipos,NVAR2D,nedgesim) = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
							      (u_ij+c_ij*a[0])*w3 + a[1]*w4 );
	IDX2(VectorAtEdge,3,ipos,NVAR2D,nedgesim) = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
							      (v_ij+c_ij*a[1])*w3 - a[0]*w4 );
	IDX2(VectorAtEdge,4,ipos,NVAR2D,nedgesim) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 +
							      (H_ij+c_ij*vel_ij)*w3 + (u_ij*a[1]-v_ij*a[0])*w4 );
      } else {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX2(VectorAtEdge,i,ipos,NVAR2D,nedgesim) = 0.0;
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
    template <int nedgesim, typename Tc, typename Td, typename Ti>
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
			     Ti ipos,
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
      
	// Compute densities
	Td ri = DENSITY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim);
	Td rj = DENSITY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim);
      
	// Compute enthalpies
	Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim)+pi)/ri;
	Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim)+pj)/rj;
      
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
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  Diff[i-1] = IDX3(DataAtEdge,i,2,ipos,NVAR2D,2,nedgesim)
	             -IDX3(DataAtEdge,i,1,ipos,NVAR2D,2,nedgesim);
      
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
	IDX2(VectorAtEdge,1,ipos,NVAR2D,nedgesim) = a[0] * ( w1 + w2 + w3 );
	IDX2(VectorAtEdge,2,ipos,NVAR2D,nedgesim) = a[0] * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 );
	IDX2(VectorAtEdge,3,ipos,NVAR2D,nedgesim) = a[0] * (        v_ij*w1 + v_ij*w2 +        v_ij*w3 - w4 );
	IDX2(VectorAtEdge,4,ipos,NVAR2D,nedgesim) = a[0] * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 +
							     (H_ij+c_ij*u_ij)*w3 - v_ij*w4 );
      
	//----------------------------------------------------------------------
	// Dimensional splitting: y-direction
	//----------------------------------------------------------------------
      
	// Compute eigenvalues
	l1 = abs(v_ij-c_ij);
	l2 = abs(v_ij);
	l3 = abs(v_ij+c_ij);
	l4 = abs(v_ij);
      
	// Compute solution difference U_j-U_i
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  Diff[i-1] = IDX3(DataAtEdge,i,2,ipos,NVAR2D,2,nedgesim)
	             -IDX3(DataAtEdge,i,1,ipos,NVAR2D,2,nedgesim);
       
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
	IDX2(VectorAtEdge,1,ipos,NVAR2D,nedgesim) += a[1] * ( w1 + w2 + w3 );
	IDX2(VectorAtEdge,2,ipos,NVAR2D,nedgesim) += a[1] * (        u_ij*w1 + u_ij*w2 +        u_ij*w3 + w4 );
	IDX2(VectorAtEdge,3,ipos,NVAR2D,nedgesim) += a[1] * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 );
	IDX2(VectorAtEdge,4,ipos,NVAR2D,nedgesim) += a[1] * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 +
							      (H_ij+c_ij*v_ij)*w3 + u_ij*w4 );
      } else {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX2(VectorAtEdge,i,ipos,NVAR2D,nedgesim) = 0.0;
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
    template <int nedgesim, typename Tc, typename Td, typename Ti>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      // Compute specific energies
      Td Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim);
      Td Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim);
    
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
#pragma unroll
      for (int i=1; i<=NVAR2D; i++)
	IDX2(VectorAtEdge,i,ipos,NVAR2D,nedgesim) =
	  d_ij*(IDX3(DataAtEdge,i,2,ipos,NVAR2D,2,nedgesim)
		-IDX3(DataAtEdge,i,1,ipos,NVAR2D,2,nedgesim));
    }
  };

  /*****************************************************************************
   * InviscidFluxDissipationBase: Specialisation for computing
   * scalar artificial dissipation of Rusanov-type using dimensional splitting.
   ****************************************************************************/

  template <>
  struct InviscidFluxDissipationBase<DISSIPATION_RUSANOV_DSPLIT>
  {
    template <int nedgesim, typename Tc, typename Td, typename Ti>
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
			     Ti ipos,
			     Ti iedge, 
			     Ti nedge,
			     Ti ncoeff)
    {
      // Compute specific energies
      Td Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,ipos,NVAR2D,2,nedgesim);
      Td Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,ipos,NVAR2D,2,nedgesim);
    
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
#pragma unroll
      for (int i=1; i<=NVAR2D; i++)
	IDX2(VectorAtEdge,i,ipos,NVAR2D,nedgesim) =
	  d_ij*(IDX3(DataAtEdge,i,2,ipos,NVAR2D,2,nedgesim)
		-IDX3(DataAtEdge,i,1,ipos,NVAR2D,2,nedgesim));
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
      InviscidFluxDissipationBase<idissipationtype>::calcEdgeData<1>
	(VectorAtEdge,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,pi,pj,1,iedge,nedge,ncoeff);
    }
  };

  /*****************************************************************************
   * FluxBase
   ****************************************************************************/

  struct FluxBase
  {
    template <int nedgesim, bool boverwrite, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
				Tc *CoeffsAtEdge,
				Td *Fxi,
				Td *Fxj, 
				Td *Fyi, 
				Td *Fyj,
				Td *Diff,
				Td scale,
				Ti ipos,
				Ti iedge, 
				Ti nedge,
				Ti ncoeff)
    {
      if (boverwrite) {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,ipos,NVAR2D,2,nedgesim) = scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[i-1]+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[i-1]-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[i-1]-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[i-1] + Diff[i-1]);

#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,ipos,NVAR2D,2,nedgesim) = -IDX3(FluxesAtEdge,i,1,ipos,NVAR2D,2,nedgesim);
      }
      else {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,ipos,NVAR2D,2,nedgesim) += scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[i-1]+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[i-1]-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[i-1]-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[i-1] + Diff[i-1]);
	
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,ipos,NVAR2D,2,nedgesim) -= scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[i-1]+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[i-1]-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[i-1]-
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[i-1] + Diff[i-1]);
      }
    }
    
    template <bool boverwrite, typename Tc, typename Td, typename Ti>
    __device__ __forceinline__
    static void combineEdgeData(Td *FluxesAtEdge,
				Tc *CoeffsAtEdge,
				Td *Fx_ij,
				Td *Fy_ij,
				Td *Diff,
				Td scale,
				Ti ipos,
				Ti iedge,
				Ti nedge,
				Ti ncoeff)
    {
      if (boverwrite) {
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,1,ipos,NVAR2D,2,nedgesim) = scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[i-1]+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[i-1] + Diff[i-1]);
    
#pragma unroll
	for (int i=1; i<=NVAR2D; i++)
	  IDX3(FluxesAtEdge,i,2,ipos,NVAR2D,2,nedgesim) = -scale *
	    (IDX3_COEFFSATEDGE(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[i-1]+
	     IDX3_COEFFSATEDGE(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[i-1] + Diff[i-1]);
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
      Flux::combineEdgeData<1,boverwrite>
	(FluxesAtEdge,CoeffsAtEdge,Fxi,Fxj,Fyi,Fyj,Diff,scale,1,iedge,nedge,ncoeff);
    }

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
      Flux::combineEdgeData<1,boverwrite>
	(FluxesAtEdge,CoeffsAtEdge,Fx_ij,Fy_ij,Diff,scale,1,iedge,nedge,ncoeff);
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
	    int idissipationtype>
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
      Ti idx = (ipt*gridDim.x+blockIdx.x)*blockDim.x+nedge_offset+threadIdx.x;
      
      if (idx<nedge_last)
	{
	  // Get positions of edge endpoints (idx starts at zero)
	  Ti i = IDX2_EDGELIST(IedgeList,1,idx+1,6,nedge);
	  Ti j = IDX2_EDGELIST(IedgeList,2,idx+1,6,nedge);
	  
	  // Local data at edge from local memory
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
	  
	  TdDest Diff[NVAR2D];
	  // Compute the artificial viscosities
	  InviscidFluxDissipation<idissipationtype>::
	    calcEdgeData(Diff,CoeffsAtEdge,DataAtEdge,
			 ui,uj,vi,vj,pi,pj,idx+1,nedge,ncoeff);
	  
#ifdef HYDRO_USE_IBP
	  TdDest Fxi[NVAR2D];
	  TdDest Fxj[NVAR2D];
	  TdDest Fyi[NVAR2D];
	  TdDest Fyj[NVAR2D];
	  
	  // Compute the Galerkin fluxes
	  InviscidFlux::
	    calcEdgeData(Fxi,Fxj,Fyi,Fyj,DataAtEdge,ui,uj,vi,vj,pi,pj);

	  // Build both contributions into the fluxes
	  Flux::
	    combineEdgeData<true>(DataAtEdge,CoeffsAtEdge,Fxi,Fxj,Fyi,Fyj,Diff,
				  scale,idx+1,nedge,ncoeff);
#else
	  TdDest Fx_ij[NVAR2D];
	  TdDest Fy_ij[NVAR2D];
	  
	  // Compute the Galerkin fluxes
	  InviscidFlux::
	    calcEdgeData(Fx_ij,Fy_ij,DataAtEdge,ui,uj,vi,vj,pi,pj);
	  
	  // Build both contributions into the fluxes
	  Flux::
	    combineEdgeData<true>(DataAtEdge,CoeffsAtEdge,Fx_ij,Fy_ij,Diff,
				  scale,idx+1,nedge,ncoeff);
#endif

	  // Build fluxes into nodal vector
	  Vector<NVAR2D,isystemformat==SYSTEM_BLOCK>::
	    scatterEdgeData<false>(vecDest,DataAtEdge,i,j,neq);
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
    // Not implemented yet
    printf("Not implemented\n");
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
    const int dma_lds                  = 1;
    const int nedge_per_thread_cudaDMA = 1;

    const int threads_per_cta_baseline  = 32*1;
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
    
    if (nblocks == 1) {

      if (grid_cudaDMA.x>0)
      	// CudaDMA implementation
      	hydro_calcFlux2d_cudaDMA
	  <Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,idissipationtype,
   	   compute_threads_per_cta,dma_threads_per_ld>
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
	  <Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,idissipationtype>
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
	   compute_threads_per_cta,dma_threads_per_ld>
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
	  <Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,idissipationtype>
	  <<<grid_baseline, block_baseline, 0, stream>>>(CoeffsAtEdge,
							 IedgeList,
							 vecSrc, vecDest, scale, 
							 neq, nedge, ncoeff,
							 nedgeset+iedgeset-1, 
							 nedge_per_thread_baseline,
							 nedge_cudaDMA+iedgeset-1);
    }
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
