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
 *#############################################################################/
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "coproc_core.h"
#include "coproc_storage_cuda.h"

#define LANGUAGE LANGUAGE_C
#include "../../../../../kernel/System/idxmanager.h"

#define HYDRO_NDIM 3
#include "hydro.h"

/*******************************************************************************
 * CUDA kernels
 *******************************************************************************/

template <int isystemformat>
struct gather_DataAtEdge
{ 
};

/*******************************************************************************/

template <>
struct gather_DataAtEdge<SYSTEM_SEGREGATED>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval (Td *DataAtEdge,
		    Td *Dx,
		    Ti i,
		    Ti j,
		    Ti neq)
  {
    // Solution vector is stored in interleaved format
    IDX2(DataAtEdge,1,1,NVAR3D,2) = IDX2_REVERSE(Dx,1,i,NVAR3D,neq);
    IDX2(DataAtEdge,2,1,NVAR3D,2) = IDX2_REVERSE(Dx,2,i,NVAR3D,neq);
    IDX2(DataAtEdge,3,1,NVAR3D,2) = IDX2_REVERSE(Dx,3,i,NVAR3D,neq);
    IDX2(DataAtEdge,4,1,NVAR3D,2) = IDX2_REVERSE(Dx,4,i,NVAR3D,neq);
    IDX2(DataAtEdge,5,1,NVAR3D,2) = IDX2_REVERSE(Dx,5,i,NVAR3D,neq);
    
    IDX2(DataAtEdge,1,2,NVAR3D,2) = IDX2_REVERSE(Dx,1,j,NVAR3D,neq);
    IDX2(DataAtEdge,2,2,NVAR3D,2) = IDX2_REVERSE(Dx,2,j,NVAR3D,neq);
    IDX2(DataAtEdge,3,2,NVAR3D,2) = IDX2_REVERSE(Dx,3,j,NVAR3D,neq);
    IDX2(DataAtEdge,4,2,NVAR3D,2) = IDX2_REVERSE(Dx,4,j,NVAR3D,neq);
    IDX2(DataAtEdge,5,2,NVAR3D,2) = IDX2_REVERSE(Dx,5,j,NVAR3D,neq);
  }
};

/*******************************************************************************/

template <>
struct gather_DataAtEdge<SYSTEM_ALLCOUPLED>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval (Td *DataAtEdge,
		    Td *Dx,
		    Ti i,
		    Ti j,
		    Ti neq)
  {
    // Solution vector is stored in block format
    IDX2(DataAtEdge,1,1,NVAR3D,2) = IDX2_FORWARD(Dx,1,i,NVAR3D,neq);
    IDX2(DataAtEdge,2,1,NVAR3D,2) = IDX2_FORWARD(Dx,2,i,NVAR3D,neq);
    IDX2(DataAtEdge,3,1,NVAR3D,2) = IDX2_FORWARD(Dx,3,i,NVAR3D,neq);
    IDX2(DataAtEdge,4,1,NVAR3D,2) = IDX2_FORWARD(Dx,4,i,NVAR3D,neq);
    IDX2(DataAtEdge,5,1,NVAR3D,2) = IDX2_FORWARD(Dx,5,i,NVAR3D,neq);
    
    IDX2(DataAtEdge,1,2,NVAR3D,2) = IDX2_FORWARD(Dx,1,j,NVAR3D,neq);
    IDX2(DataAtEdge,2,2,NVAR3D,2) = IDX2_FORWARD(Dx,2,j,NVAR3D,neq);
    IDX2(DataAtEdge,3,2,NVAR3D,2) = IDX2_FORWARD(Dx,3,j,NVAR3D,neq);
    IDX2(DataAtEdge,4,2,NVAR3D,2) = IDX2_FORWARD(Dx,4,j,NVAR3D,neq);
    IDX2(DataAtEdge,5,2,NVAR3D,2) = IDX2_FORWARD(Dx,5,j,NVAR3D,neq);
  }
};

/*******************************************************************************/

struct calc_GalerkinFlux
{
  template <typename Td>
#ifdef HYDRO_USE_IBP
  __device__ inline
  static void eval(Td *Fxi,
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
    // Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      
      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);

      // Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      
      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);

      // Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      
      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
  }
#else
  __device__ inline
  static void eval(Td *Fx_ij,
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
    // Compute flux difference for x-direction
    IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                    INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
    IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                    INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
    IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                    INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
    IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                    INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
    IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                    INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
    
    // Compute flux difference for y-direction
    IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      
    // Compute flux difference for z-direction
    IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
    IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                    INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
  }
#endif
};

/*******************************************************************************/

template <int idissipationtype>
struct calc_Dissipation
{
};

/*******************************************************************************/

template <>
struct calc_Dissipation<DISSIPATION_ZERO>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval(Td *Diff,
		   Td *CoeffsAtEdge,
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
    Diff[0] = 0.0;
    Diff[1] = 0.0;
    Diff[2] = 0.0;
    Diff[3] = 0.0;
    Diff[4] = 0.0;
  }
};

/*******************************************************************************/

template <>
struct calc_Dissipation<DISSIPATION_SCALAR>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval(Td *Diff,
		   Td *CoeffsAtEdge,
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
    //------------------------------------------------------------------------
    // Evaluate the scalar dissipation proportional to the spectral
    // radius (largest eigenvalue) of the Roe-matrix
    //------------------------------------------------------------------------
    
    // Compute skew-symmetric coefficient
    Td a[HYDRO_NDIM];
    a[0] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[1] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[2] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    Td anorm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    
    // Compute densities
    Td ri = DENSITY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
    Td rj = DENSITY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
    
    // Compute enthalpies
    Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1)+pi)/ri;
    Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1)+pj)/rj;
    
    // Compute Roe mean values
    Td aux  = ROE_MEAN_RATIO(ri,rj);
    Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
    Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
    Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
    Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
    
    // Compute auxiliary variables
    Td vel_ij = u_ij * a[0] + v_ij * a[1] + w_ij * a[2];
    Td q_ij   = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
    
    // Compute the speed of sound
    //TODO echtes double epsilon einbauen
    Td c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), 1e-14));
    
    // Compute scalar dissipation
    Td d_ij = abs(vel_ij) + anorm*c_ij;
    
    // Multiply the solution difference by the scalar dissipation
    Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1));
    Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1));
    Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1));
    Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1));
    Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1));
  }
};

/*******************************************************************************/

template <>
struct calc_Dissipation<DISSIPATION_SCALAR_DSPLIT>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval(Td *Diff,
		   Td *CoeffsAtEdge,
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
    //------------------------------------------------------------------------
    // Evaluate the scalar dissipation proportional to the spectral
    // radius (largest eigenvalue) of the dimensional-split Roe-matrix
    //------------------------------------------------------------------------
    
    // Compute skew-symmetric coefficient
    Td a[HYDRO_NDIM];
    a[0] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[1] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[2] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    
    // Compute densities
    Td ri = DENSITY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
    Td rj = DENSITY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
    
    // Compute enthalpies
    Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1)+pi)/ri;
    Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1)+pj)/rj;
    
    // Compute Roe mean values
    Td aux  = ROE_MEAN_RATIO(ri,rj);
    Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
    Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
    Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
    Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
    
    // Compute auxiliary variables
    Td q_ij = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
    
    // Compute the speed of sound
    //TODO echtes double epsilon einbauen
    Td c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), 1e-14));
    
    // Compute scalar dissipation
    Td d_ij = ( abs(a[0]*u_ij) + abs(a[0])*c_ij +
 	        abs(a[1]*v_ij) + abs(a[1])*c_ij +
	        abs(a[2]*w_ij) + abs(a[2])*c_ij );
    
    // Multiply the solution difference by the scalar dissipation
    Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1));
    Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1));
    Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1));
    Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1));
    Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1));
  }
};

/*******************************************************************************/

template <>
struct calc_Dissipation<DISSIPATION_ROE>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval(Td *Diff, 
		   Td *CoeffsAtEdge,
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
    //------------------------------------------------------------------------
    // Evaluate the dissipation tensor of Roe-type
    //------------------------------------------------------------------------
    
    // Compute skew-symmetric coefficient
    Td a[HYDRO_NDIM];
    a[0] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[1] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[2] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    Td anorm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    
    if (anorm > 1e-14) {
      
      // Normalise the skew-symmetric coefficient
      a[0] = a[0]/anorm;
      a[1] = a[1]/anorm;
      a[2] = a[2]/anorm;
      
      // Compute densities
      Td ri = DENSITY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
      Td rj = DENSITY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
      
      // Compute enthalpies
      Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1)+pi)/ri;
      Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1)+pj)/rj;
      
      // Compute Roe mean values
      Td aux  = ROE_MEAN_RATIO(ri,rj);
      Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
      Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
      // Compute auxiliary variables
      Td vel_ij = u_ij * a[0] + v_ij * a[1] + w_ij * a[2];
      Td q_ij   = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
      
      // Compute the speed of sound
      //TODO echtes double epsilon einbauen
      Td c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), 1e-14);
      Td c_ij  = sqrt(c2_ij);
      
      // Compute eigenvalues
      Td l1 = abs(vel_ij-c_ij);
      Td l2 = abs(vel_ij);
      Td l3 = abs(vel_ij+c_ij);
      Td l4 = abs(vel_ij);
      Td l5 = abs(vel_ij);
      
      // Compute solution difference U_j-U_i
      Diff[0] = IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1);
      Diff[1] = IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1);
      Diff[2] = IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1);
      Diff[3] = IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1);
      Diff[4] = IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1);
      
      // Compute auxiliary quantities for characteristic variables
      Td aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*Diff[0]
					    -u_ij*Diff[1]
					    -v_ij*Diff[2]
					    -w_ij*Diff[3]
					         +Diff[4])/RCONST(2.0)/c2_ij;
      Td aux2 = (vel_ij*Diff[0]
		  -a[0]*Diff[1]
	          -a[1]*Diff[2]
	          -a[2]*Diff[3])/RCONST(2.0)/c_ij;
      
      // Get the dimension with largest coefficient
      if (a[0] >= a[1] && a[0] >= a[2]) {
	
	// Compute characteristic variables multiplied by the corresponding eigenvalue
	Td w1 = l1 * (aux1 + aux2);
	Td w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff[0]
                                       +((HYDRO_GAMMA)-RCONST(1.0))*( u_ij*Diff[1]
								     +v_ij*Diff[2]
                                                                     +w_ij*Diff[3]
								          -Diff[4])/c2_ij);
	Td w3 = l3 * (aux1 - aux2);
	Td w4 = l4 * ( (v_ij-vel_ij*a[1])/a[0]*Diff[0]
			                 +a[1]*Diff[1]
		 +(a[1]*a[1]-RCONST(1.0))/a[0]*Diff[2]
		               +a[1]*a[2]/a[0]*Diff[3]);
	Td w5 = l5 * ( (vel_ij*a[2]-w_ij)/a[0]*Diff[0]
                                         -a[2]*Diff[1]
                               -a[1]*a[2]/a[0]*Diff[2]
	         +(RCONST(1.0)-a[2]*a[2])/a[0]*Diff[3]);
	
	// Compute "R_ij * |Lbd_ij| * L_ij * dU"
	Diff[0] = anorm * ( w1 + w2 + w3 );
	Diff[1] = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
			    (u_ij+c_ij*a[0])*w3 + a[1]*w4 - a[2]*w5 );
	Diff[2] = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
			    (v_ij+c_ij*a[1])*w3 - a[0]*w4 );
	Diff[3] = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
			    (w_ij+c_ij*a[2])*w3 + a[0]*w5 );
	Diff[4] = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
			    + (u_ij*a[1]-v_ij*a[0])*w4 + (w_ij*a[0]-u_ij*a[2])*w5 );
	
      } else if (a[1] >= a[0] && a[1] >= a[2]) {
	// Compute characteristic variables multiplied by the corresponding eigenvalue
	Td w1 = l1 * (aux1 + aux2);
	Td w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff[0]
		                        +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff[1]
								     +v_ij*Diff[2]
                                                                     +w_ij*Diff[3]
								          -Diff[4])/c2_ij);
	Td w3 = l3 * (aux1 - aux2);
	Td w4 = l4 * ( (vel_ij*a[0]-u_ij)/a[1]*Diff[0]
	 	 +(RCONST(1.0)-a[0]*a[0])/a[1]*Diff[1]
		                         -a[0]*Diff[2]
		               -a[0]*a[2]/a[1]*Diff[3]);
	Td w5 = l5 * ( (w_ij-vel_ij*a[2])/a[1]*Diff[0]
		               +a[0]*a[2]/a[1]*Diff[1]
		                         +a[2]*Diff[2]
		 +(a[2]*a[2]-RCONST(1.0))/a[1]*Diff[3]);
	  
	// Compute "R_ij * |Lbd_ij| * L_ij * dU"
	Diff[0] = anorm * ( w1 + w2 + w3 );
	Diff[1] = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
			    (u_ij+c_ij*a[0])*w3 + a[1]*w4 );
	Diff[2] = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
			    (v_ij+c_ij*a[1])*w3 - a[0]*w4 + a[2]*w5 );
	Diff[3] = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
			    (w_ij+c_ij*a[2])*w3 - a[1]*w5 );
	Diff[4] = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
			    + (u_ij*a[1]-v_ij*a[0])*w4 + (v_ij*a[2]-w_ij*a[1])*w5 );
	
      } else {
	// Compute characteristic variables multiplied by the corresponding eigenvalue
	Td w1 = l1 * (aux1 + aux2);
	Td w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff[0]
		                        +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff[1]
							 	     +v_ij*Diff[2]
								     +w_ij*Diff[3]
								          -Diff[4])/c2_ij);
	Td w3 = l3 * (aux1 - aux2);
	Td w4 = l4 * ( (u_ij-vel_ij*a[0])/a[2]*Diff[0]
	 	 +(a[0]*a[0]-RCONST(1.0))/a[2]*Diff[1]
		               +a[0]*a[1]/a[2]*Diff[2]
		                         +a[0]*Diff[3]);
	Td w5 = l5 * ( (vel_ij*a[1]-v_ij)/a[2]*Diff[0]
		               -a[0]*a[1]/a[2]*Diff[1]
		 +(RCONST(1.0)-a[1]*a[1])/a[2]*Diff[2]
		                         -a[1]*Diff[3]);
	  
	// Compute "R_ij * |Lbd_ij| * L_ij * dU"
	Diff[0] = anorm * ( w1 + w2 + w3 );
	Diff[1] = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
			    (u_ij+c_ij*a[0])*w3 - a[2]*w4 );
	Diff[2] = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
			    (v_ij+c_ij*a[1])*w3 + a[2]*w5 );
	Diff[3] = anorm * ( (w_ij-c_ij*a[2])*w1 + w_ij*w2 +
			    (w_ij+c_ij*a[2])*w3 + a[0]*w4 - a[1]*w5);
	Diff[4] = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3
			    + (w_ij*a[0]-u_ij*a[2])*w4 + (v_ij*a[2]-w_ij*a[1])*w5 );
      }
    } else {
      Diff[0] = 0.0;
      Diff[1] = 0.0;
      Diff[2] = 0.0;
      Diff[3] = 0.0;
      Diff[4] = 0.0;
    }
  }
};

/*******************************************************************************/

template <>
struct calc_Dissipation<DISSIPATION_ROE_DSPLIT>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval(Td *Diff,
		   Td *CoeffsAtEdge,
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
    //------------------------------------------------------------------------
    // Evaluate the dissipation tensor of Roe-type with dimensional splitting
    //------------------------------------------------------------------------
    
  // Compute skew-symmetric coefficient
    Td a[HYDRO_NDIM];
    a[0] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[1] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    a[2] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    Td anorm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    
    Td DiffAux[NVAR3D];
    if (anorm > 1e-14) {
      
      // Compute the absolute value
      a[0] = abs(a[0]);
      a[1] = abs(a[1]);
      a[2] = abs(a[2]);
      
      // Compute densities
      Td ri = DENSITY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
      Td rj = DENSITY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
      
      // Compute enthalpies
      Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1)+pi)/ri;
      Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1)+pj)/rj;
      
      // Compute Roe mean values
      Td aux  = ROE_MEAN_RATIO(ri,rj);
      Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
      Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
      // Compute auxiliary variable
      Td q_ij = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
      
      // Compute the speed of sound
      //TODO echtes double epsilon einbauen
      Td c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), 1e-14);
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
      DiffAux[0] = IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1);
      DiffAux[1] = IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1);
      DiffAux[2] = IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1);
      DiffAux[3] = IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1);
      DiffAux[4] = IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1);
      
      // Compute auxiliary quantities for characteristic variables
      Td aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*DiffAux[0]
					    -u_ij*DiffAux[1]
					    -v_ij*DiffAux[2]
					    -w_ij*DiffAux[3]
					         +DiffAux[4])/RCONST(2.0)/c2_ij;
      Td aux2 = (u_ij*DiffAux[0]-DiffAux[1])/RCONST(2.0)/c_ij;
      
      // Compute characteristic variables multiplied by the corresponding eigenvalue
      Td w1 = l1 * (aux1 + aux2);
      Td w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffAux[0]
		                     +((HYDRO_GAMMA)-RCONST(1.0))*( u_ij*DiffAux[1]
								   +v_ij*DiffAux[2]
                                                                   +w_ij*DiffAux[3]
								        -DiffAux[4])/c2_ij);
      Td w3 = l3 * (aux1 - aux2);
      Td w4 = l4 * ( v_ij*DiffAux[0]-DiffAux[2]);
      Td w5 = l5 * (-w_ij*DiffAux[0]+DiffAux[3]);
      
      // Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff[0] = a[0] * ( w1 + w2 + w3 );
      Diff[1] = a[0] * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 );
      Diff[2] = a[0] * ( v_ij*(w1 + w2 + w3) - w4 );
      Diff[3] = a[0] * ( w_ij*(w1 + w2 + w3) + w5 );
      Diff[4] = a[0] * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 + (H_ij+c_ij*u_ij)*w3
			 -v_ij*w4 + w_ij*w5 );
      
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
      DiffAux[0] = IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1);
      DiffAux[1] = IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1);
      DiffAux[2] = IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1);
      DiffAux[3] = IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1);
      DiffAux[4] = IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1);
      
      // Compute auxiliary quantities for characteristic variables
      aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*DiffAux[0]
					 -u_ij*DiffAux[1]
					 -v_ij*DiffAux[2]
					 -w_ij*DiffAux[3]
					      +DiffAux[4])/RCONST(2.0)/c2_ij;
      aux2 = (v_ij*DiffAux[0]-DiffAux[2])/RCONST(2.0)/c_ij;
      
      // Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2);
      w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffAux[0]
		                   +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*DiffAux[1]
								+v_ij*DiffAux[2]
								+w_ij*DiffAux[3]
								     -DiffAux[4])/c2_ij);
      w3 = l3 * (aux1 - aux2);
      w4 = l4 * (-u_ij*DiffAux[0]+DiffAux[1]);
      w5 = l5 * ( w_ij*DiffAux[0]-DiffAux[3]);
      
      // Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff[0] += a[1] * ( w1 + w2 + w3 );
      Diff[1] += a[1] * ( u_ij*(w1 + w2 + w3) + w4 );
      Diff[2] += a[1] * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 );
      Diff[3] += a[1] * ( w_ij*(w1 + w2 + w3) - w5 );
      Diff[4] += a[1] * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 + (H_ij+c_ij*v_ij)*w3
			  + u_ij*w4 -w_ij*w5 );
      
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
      DiffAux[0] = IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1);
      DiffAux[1] = IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1);
      DiffAux[2] = IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1);
      DiffAux[3] = IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1);
      DiffAux[4] = IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1);
      
      // Compute auxiliary quantities for characteristic variables
      aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*DiffAux[0]
					 -u_ij*DiffAux[1]
					 -v_ij*DiffAux[2]
					 -w_ij*DiffAux[3]
					      +DiffAux[4])/RCONST(2.0)/c2_ij;
      aux2 = (w_ij*DiffAux[0]-DiffAux[2])/RCONST(2.0)/c_ij;
      
      // Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2);
      w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffAux[0]
		                   +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*DiffAux[1]
								+v_ij*DiffAux[2]
								+w_ij*DiffAux[3]
								     -DiffAux[4])/c2_ij);
      w3 = l3 * (aux1 - aux2);
      w4 = l4 * ( u_ij*DiffAux[0]-DiffAux[1]);
      w5 = l5 * (-v_ij*DiffAux[0]+DiffAux[2]);
      
      // Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff[0] += a[2] * ( w1 + w2 + w3 );
      Diff[1] += a[2] * ( u_ij*(w1 + w2 + w3) - w4 );
      Diff[2] += a[2] * ( v_ij*(w1 + w2 + w3) + w5 );
      Diff[3] += a[2] * ( (w_ij-c_ij)*w1 + w_ij*w2 + (w_ij+c_ij)*w3 );
      Diff[4] += a[2] * ( (H_ij-c_ij*w_ij)*w1 + q_ij*w2 + (H_ij+c_ij*w_ij)*w3
			  -u_ij*w4 + v_ij*w5 );
    } else {
      Diff[0] = 0.0;
      Diff[1] = 0.0;
      Diff[2] = 0.0;
      Diff[3] = 0.0;
      Diff[4] = 0.0;
      }
  }
};

/*******************************************************************************/

template <>
struct calc_Dissipation<DISSIPATION_RUSANOV>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval(Td *Diff,
		   Td *CoeffsAtEdge,
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
    //------------------------------------------------------------------------
    // Evaluate the scalar dissipation of Rusanov-type
    //------------------------------------------------------------------------
    
    // Compute specific energies
    Td Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
    Td Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
    
    // Compute the speed of sound
    Td ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*(Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), 1e-14));
    Td cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*(Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), 1e-14));
    
#ifdef HYDRO_USE_IBP
    // Compute scalar dissipation based on the skew-symmetric part
    // which does not include the symmetric boundary contribution
    Td d_ij = max( abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge))*uj+
		       RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			 	    IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge))*vj-
		       RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge))*wj)+
	       RCONST(0.5)*sqrt(POW(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
			        POW(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
		 	        POW(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),2))*cj,
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge))*ui+
		       RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge))*vi+
	 	       RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge))*wi)+
	       RCONST(0.5)*sqrt(POW(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
			        POW(IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
			        POW(IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),2))*ci );
#else
    // Compute scalar dissipation
    Td d_ij = max( abs(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*uj+
		       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*vj+
		       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*j)+
	      sqrt(POW(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
	 	   POW(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
		   POW(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge),2))*cj,
		   abs(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*ui+
		       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*vi+
		       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*wi)+
	      sqrt(POW(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
	 	   POW(IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge),2)+
		   POW(IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge),2))*ci );
#endif
    
    // Multiply the solution difference by the scalar dissipation
    Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1));
    Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1));
    Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1));
    Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1));
    Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1));
  }
};

/*******************************************************************************/

template <>
struct calc_Dissipation<DISSIPATION_RUSANOV_DSPLIT>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval(Td *Diff, 
		   Td *CoeffsAtEdge,
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
    //------------------------------------------------------------------------
    // Evaluate the scalar dissipation of Rusanov-type with dimensional splitting
    //------------------------------------------------------------------------
    
    // Compute specific energies
    Td Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
    Td Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
    
    // Compute the speed of sound
    Td ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*(Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), 1e-14));
    Td cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*(Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), 1e-14));
    
#ifdef HYDRO_USE_IBP
    // Compute scalar dissipation with dimensional splitting based on
    // the skew-symmetric part which does not include the symmetric
    // boundary contribution
    Td d_ij = max( abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge))*uj)+
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)))*cj,
	   	   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
		 	  	    IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge))*ui)+
	 	   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)))*ci )
            + max( abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge))*vj)+
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)))*cj,
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge))*vi)+
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)))*ci )
            + max( abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge))*wj)+
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)))*cj,
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge))*wi)+
		   abs(RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)-
				    IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)))*ci );
#else
    // Compute scalar dissipation with dimensional splitting
    Td d_ij = max( abs(IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*uj)+
		   abs(IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge))*cj,
		   abs(IDX3(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*ui)+
		   abs(IDX3(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge))*ci )
            + max( abs(IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*vj)+
		   abs(IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge))*cj,
		   abs(IDX3(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*vi)+
		   abs(IDX3(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge))*ci )
            + max( abs(IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*wj)+
		   abs(IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge))*cj,
		   abs(IDX3(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*wi)+
		   abs(IDX3(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge))*ci );
#endif
    
    // Multiply the solution difference by the scalar dissipation
    Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,1,1,1,NVAR3D,2,1));
    Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,2,1,1,NVAR3D,2,1));
    Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,3,1,1,NVAR3D,2,1));
    Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,4,1,1,NVAR3D,2,1));
    Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DataAtEdge,5,1,1,NVAR3D,2,1));
  }
};

template <int isystemformat>
struct scatter_FluxesAtEdge
{ 
};

/*******************************************************************************/

struct calc_FluxesAtEdge
{
  template <typename Td, typename Ti>
#ifdef HYDRO_USE_IBP
  __device__ inline
  static void eval(Td *FluxesAtEdge,
		   Td *CoeffsAtEdge,
		   Td *Fxi,
		   Td *Fxj, 
		   Td *Fyi, 
		   Td *Fyj,
		   Td *Fzi,
		   Td *Fzj, 
		   Td *Diff,
		   Td scale,
		   Ti iedge, 
		   Ti ncoeff,
		   Ti nedge)
  {
    IDX3(FluxesAtEdge,1,1,1,NVAR3D,2,1) = scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[0]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[0]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[0]-
       IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[0]-
       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[0]-
       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[0] + Diff[0]);
    
    IDX3(FluxesAtEdge,2,1,1,NVAR3D,2,1) = scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[1]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[1]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[1]-
       IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[1]-
       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[1]-
       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[1] + Diff[1]);
    
    IDX3(FluxesAtEdge,3,1,1,NVAR3D,2,1) = scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[2]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[2]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[2]-
       IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[2]-
       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[2]-
       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[2] + Diff[2]);
    
    IDX3(FluxesAtEdge,4,1,1,NVAR3D,2,1) = scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[3]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[3]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[3]-
       IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[3]-
       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[3]-
       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[3] + Diff[3]);
    
    IDX3(FluxesAtEdge,5,1,1,NVAR3D,2,1) = scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[4]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[4]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[4]-
       IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[4]-
       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[4]-
       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[4] + Diff[4]);
    
    
    IDX3(FluxesAtEdge,1,2,1,NVAR3D,2,1) = -IDX3(FluxesAtEdge,1,1,1,NVAR3D,2,1);
    IDX3(FluxesAtEdge,2,2,1,NVAR3D,2,1) = -IDX3(FluxesAtEdge,2,1,1,NVAR3D,2,1);
    IDX3(FluxesAtEdge,3,2,1,NVAR3D,2,1) = -IDX3(FluxesAtEdge,3,1,1,NVAR3D,2,1);
    IDX3(FluxesAtEdge,4,2,1,NVAR3D,2,1) = -IDX3(FluxesAtEdge,4,1,1,NVAR3D,2,1);
    IDX3(FluxesAtEdge,5,2,1,NVAR3D,2,1) = -IDX3(FluxesAtEdge,5,1,1,NVAR3D,2,1);
  }
#else
  __device__ inline
  static void eval(Td *FluxesAtEdge,
		   Td *CoeffsAtEdge,
		   Td *Fx_ij,
		   Td *Fy_ij,
		   Td *Fz_ij,
		   Td *Diff,
		   Td scale,
		   Ti iedge,
		   Ti ncoeff,
		   Ti nedge)
  {
    IDX3(FluxesAtEdge,1,1,1,NVAR3D,2,1) = scale *
      (IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[0]+
       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[0]+
       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[0] + Diff[0]);
    
    IDX3(FluxesAtEdge,2,1,1,NVAR3D,2,1) = scale *
      (IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[1]+
       IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[1]+
       IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[1] + Diff[1]);
    
    IDX3(FluxesAtEdge,3,1,1,NVAR3D,2,1) = scale *
      (IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[2]+
       IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[2]+
       IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[2] + Diff[2]);
    
    IDX3(FluxesAtEdge,4,1,1,NVAR3D,2,1) = scale *
      (IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[3]+
       IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[3]+
       IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[3] + Diff[3]);
    
    IDX3(FluxesAtEdge,5,1,1,NVAR3D,2,1) = scale *
      (IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[4]+
       IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[4]+
       IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[4] + Diff[4]);
    
    
    IDX3(FluxesAtEdge,1,2,1,NVAR3D,2,1) = -scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[0]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[0]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[0] + Diff[0]);
    
    IDX3(FluxesAtEdge,2,2,1,NVAR3D,2,1) = -scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[1]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[1]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[1] + Diff[1]);
    
    IDX3(FluxesAtEdge,3,2,1,NVAR3D,2,1) = -scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[2]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[2]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[2] + Diff[2]);
    
    IDX3(FluxesAtEdge,4,2,1,NVAR3D,2,1) = -scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[3]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[3]+
       IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[3] + Diff[3]);
    
    IDX3(FluxesAtEdge,5,2,1,NVAR3D,2,1) = -scale *
      (IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[4]+
       IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[4]+
       IDX3T(CoeffsAtEdge,3,2,iedge+idx,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[4] + Diff[4]);
  }
#endif
};

/*******************************************************************************/

template <>
struct scatter_FluxesAtEdge<SYSTEM_SEGREGATED>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval (Td *FluxesAtEdge,
		    Td *Dy,
		    Ti i,
		    Ti j,
		    Ti neq)
  {
    // Solution vector is stored in interleaved format
    IDX2_REVERSE(Dy,1,i,NVAR3D,neq) += IDX3(FluxesAtEdge,1,1,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,2,i,NVAR3D,neq) += IDX3(FluxesAtEdge,2,1,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,3,i,NVAR3D,neq) += IDX3(FluxesAtEdge,3,1,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,4,i,NVAR3D,neq) += IDX3(FluxesAtEdge,4,1,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,5,i,NVAR3D,neq) += IDX3(FluxesAtEdge,5,1,1,NVAR3D,2,1);
    
    IDX2_REVERSE(Dy,1,j,NVAR3D,neq) += IDX3(FluxesAtEdge,1,2,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,2,j,NVAR3D,neq) += IDX3(FluxesAtEdge,2,2,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,3,j,NVAR3D,neq) += IDX3(FluxesAtEdge,3,2,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,4,j,NVAR3D,neq) += IDX3(FluxesAtEdge,4,2,1,NVAR3D,2,1);
    IDX2_REVERSE(Dy,5,j,NVAR3D,neq) += IDX3(FluxesAtEdge,5,2,1,NVAR3D,2,1);
  }
};

/*******************************************************************************/

template <>
struct scatter_FluxesAtEdge<SYSTEM_ALLCOUPLED>
{
  template <typename Td, typename Ti>
  __device__ inline
  static void eval (Td *FluxesAtEdge,
		    Td *Dy,
		    Ti i,
		    Ti j,
		    Ti neq)
  {
    // Solution vector is stored in block format
    IDX2_FORWARD(Dy,1,i,NVAR3D,neq) += IDX3(FluxesAtEdge,1,1,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,2,i,NVAR3D,neq) += IDX3(FluxesAtEdge,2,1,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,3,i,NVAR3D,neq) += IDX3(FluxesAtEdge,3,1,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,4,i,NVAR3D,neq) += IDX3(FluxesAtEdge,4,1,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,5,i,NVAR3D,neq) += IDX3(FluxesAtEdge,5,1,1,NVAR3D,2,1);
    
    IDX2_FORWARD(Dy,1,j,NVAR3D,neq) += IDX3(FluxesAtEdge,1,2,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,2,j,NVAR3D,neq) += IDX3(FluxesAtEdge,2,2,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,3,j,NVAR3D,neq) += IDX3(FluxesAtEdge,3,2,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,4,j,NVAR3D,neq) += IDX3(FluxesAtEdge,4,2,1,NVAR3D,2,1);
    IDX2_FORWARD(Dy,5,j,NVAR3D,neq) += IDX3(FluxesAtEdge,5,2,1,NVAR3D,2,1);
  }
};

/*******************************************************************************/
  
template <typename Td, typename Ti, int isystemformat, int idissipationtype>
__global__ void hydro_calcFlux3d_knl(Td *CoeffsAtEdge,
				     Ti *IedgeList,
				     Td *Dx,
				     Td *Dy,
				     Td scale,
				     Ti neq,
				     Ti nvar,
				     Ti nedge,
				     Ti ncoeff,
				     Ti nedges,
				     Ti iedgeset)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx<nedges)
    {
      // Get positions of edge endpoints (idx starts at zero)
      Ti i = IDX2(IedgeList,1,iedgeset+idx,6,nedge);
      Ti j = IDX2(IedgeList,2,iedgeset+idx,6,nedge);
      
      // Local variables
      Td DataAtEdge[2*NVAR3D];
      
      // Get solution values at edge endpoints
      gather_DataAtEdge<isystemformat>::
	eval(DataAtEdge,Dx,i,j,neq);
      
      // Compute velocities
      Td ui = XVELOCITY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
      Td vi = YVELOCITY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
      Td wi = ZVELOCITY3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
      
      Td uj = XVELOCITY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
      Td vj = YVELOCITY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
      Td wj = ZVELOCITY3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);
      
      // Compute pressures
      Td pi = PRESSURE3(DataAtEdge,IDX3,1,1,NVAR3D,2,1);
      Td pj = PRESSURE3(DataAtEdge,IDX3,2,1,NVAR3D,2,1);

#ifdef HYDRO_USE_IBP
      Td Fxi[NVAR3D];
      Td Fxj[NVAR3D];
      Td Fyi[NVAR3D];
      Td Fyj[NVAR3D];
      Td Fzi[NVAR3D];
      Td Fzj[NVAR3D];
      
      // Compute the Galerkin fluxes
      calc_GalerkinFlux::
	eval(Fxi,Fxj,Fyi,Fyj,Fzi,Fzj,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj);
#else
      Td Fx_ij[NVAR3D];
      Td Fy_ij[NVAR3D];
      Td Fz_ij[NVAR3D];

      // Compute the Galerkin fluxes
      calc_GalerkinFlux::
	calc_eval(Fx_ij,Fy_ij,Fz_ij,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj);
#endif
      
      Td Diff[NVAR3D];
      // Compute the artificial viscosities
      calc_Dissipation<idissipationtype>::
	eval(Diff,CoeffsAtEdge,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,iedgeset+idx,nedge,ncoeff);
      
      // Build both contributions into the fluxes
#ifdef HYDRO_USE_IBP
      calc_FluxesAtEdge::
	eval(DataAtEdge,CoeffsAtEdge,Fxi,Fxj,Fyi,Fyj,Fzi,Fzj,Diff,scale,iedgeset+idx,ncoeff,nedge);
#else
      calc_FluxesAtEdge::
	eval(DataAtEdge,CoeffsAtEdge,Fx_ij,Fy_ij,Fz_ij,Diff,scale,iedgeset+idx,ncoeff,nedge);
#endif
        
      // Build fluxes into nodal vector
      scatter_FluxesAtEdge<isystemformat>::
	eval(DataAtEdge,Dy,i,j,neq);
    }
}

/*******************************************************************************
 * External C functions which can be called from the Fortran code
 *******************************************************************************/

template <typename Td, typename Ti>
inline
int hydro_calcFluxGalerkin3d_cuda(__I64 *d_CoeffsAtEdge,
				  __I64 *d_IedgeList,
				  __I64 *d_Dx,
				  __I64 *d_Dy,
				  Td scale,
				  Ti nblocks,
				  Ti neq,
				  Ti nvar,
				  Ti nedge, 
				  Ti ncoeff,
				  Ti nedges,
				  Ti iedgeset,
				  cudaStream_t stream=0)
{
  Td *Dx = (Td*)(*d_Dx);
  Td *Dy = (Td*)(*d_Dy);
  Td *CoeffsAtEdge = (Td*)(*d_CoeffsAtEdge);
  Ti *IedgeList = (Ti*)(*d_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
  if (nblocks == 1) {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_SEGREGATED,DISSIPATION_ZERO>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale,
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  } else {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_ALLCOUPLED,DISSIPATION_ZERO>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  }
  coproc_checkErrors("hydro_calcFluxGalerkin3d_cuda");
  return 0;
}

/*******************************************************************************/
extern "C" {
  __INT FNAME(hydro_calcfluxgalerkin3d_cuda)(__I64 *d_CoeffsAtEdge,
					     __I64 *d_IedgeList,
					     __I64 *d_Dx,
					     __I64 *d_Dy,
#ifdef HAS_CUDADOUBLEPREC
					     __DP *scale,
#else
					     __SP *scale,
#endif
					     __INT *nblocks,
					     __INT *neq,
					     __INT *nvar,
					     __INT *nedge,
					     __INT *ncoeff,
					     __INT *nedges,
					     __INT *iedgeset,
					     __I64 *stream)
  {
    return (__INT) hydro_calcFluxGalerkin3d_cuda(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
						 *scale, *nblocks, *neq, *nvar, *nedge,
						 *ncoeff, *nedges, *iedgeset,
						 (cudaStream_t)*stream);
  }
}

template <typename Td, typename Ti>
inline
int hydro_calcFluxScDiss3d_cuda(__I64 *d_CoeffsAtEdge,
				__I64 *d_IedgeList,
				__I64 *d_Dx,
				__I64 *d_Dy,
				Td scale,
				Ti nblocks,
				Ti neq,
				Ti nvar,
				Ti nedge, 
				Ti ncoeff,
				Ti nedges,
				Ti iedgeset,
				cudaStream_t stream=0)
{
  Td *Dx = (Td*)(*d_Dx);
  Td *Dy = (Td*)(*d_Dy);
  Td *CoeffsAtEdge = (Td*)(*d_CoeffsAtEdge);
  Ti *IedgeList = (Ti*)(*d_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
  if (nblocks == 1) {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_SEGREGATED,DISSIPATION_SCALAR>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale,
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  } else {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_ALLCOUPLED,DISSIPATION_SCALAR>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  }
  coproc_checkErrors("hydro_calcFluxScDiss3d_cuda");
  return 0;
}

/*******************************************************************************/
extern "C" {
  __INT FNAME(hydro_calcfluxscdiss3d_cuda)(__I64 *d_CoeffsAtEdge,
					   __I64 *d_IedgeList,
					   __I64 *d_Dx,
					   __I64 *d_Dy,
#ifdef HAS_CUDADOUBLEPREC
					   __DP *scale,
#else
					   __SP *scale,
#endif
					   __INT *nblocks,
					   __INT *neq,
					   __INT *nvar,
					   __INT *nedge,
					   __INT *ncoeff,
					   __INT *nedges,
					   __INT *iedgeset,
					   __I64 *stream)
  {
    return (__INT) hydro_calcFluxScDiss3d_cuda(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
					       *scale, *nblocks, *neq, *nvar, *nedge,
					       *ncoeff, *nedges, *iedgeset,
					       (cudaStream_t)*stream);
  }
}

/*******************************************************************************/

template <typename Td, typename Ti>
inline
int hydro_calcFluxScDissDiSp3d_cuda(__I64 *d_CoeffsAtEdge,
				    __I64 *d_IedgeList,
				    __I64 *d_Dx,
				    __I64 *d_Dy,
				    Td scale,
				    Ti nblocks,
				    Ti neq, 
				    Ti nvar,
				    Ti nedge,
				    Ti ncoeff,
				    Ti nedges,
				    Ti iedgeset,
				    cudaStream_t stream=0)
{
  Td *Dx = (Td*)(*d_Dx);
  Td *Dy = (Td*)(*d_Dy);
  Td *CoeffsAtEdge = (Td*)(*d_CoeffsAtEdge);
  Ti *IedgeList = (Ti*)(*d_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
  if (nblocks == 1) {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_SEGREGATED,DISSIPATION_SCALAR_DSPLIT>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  } else {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_ALLCOUPLED,DISSIPATION_SCALAR_DSPLIT>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  }
  coproc_checkErrors("hydro_calcFluxScDissDiSp3d_cuda");
  return 0;
}

/*******************************************************************************/
extern "C" {
  __INT FNAME(hydro_calcfluxscdissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
					       __SIZET *d_IedgeList,
					       __SIZET *d_Dx,
					       __SIZET *d_Dy,
#ifdef HAS_CUDADOUBLEPREC
					       __DP *scale,
#else
					       __SP *scale,
#endif
					       __INT *nblocks, 
					       __INT *neq, 
					       __INT *nvar,
					       __INT *nedge, 
					       __INT *ncoeff,
					       __INT *nedges, 
					       __INT *iedgeset,
					       __I64 *stream=0)
  {
    return (__INT) hydro_calcFluxScDissDiSp3d_cuda(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
						   *scale, *nblocks, *neq, *nvar, *nedge,
						   *ncoeff, *nedges, *iedgeset,
						   (cudaStream_t)(*stream));
  }
}

/*******************************************************************************/

template <typename Td, typename Ti>
inline
int hydro_calcFluxRoeDiss3d_cuda(__I64 *d_CoeffsAtEdge,
				 __I64 *d_IedgeList,
				 __I64 *d_Dx,
				 __I64 *d_Dy,
				 Td scale,
				 Ti nblocks, 
				 Ti neq, 
				 Ti nvar,
				 Ti nedge,
				 Ti ncoeff,
				 Ti nedges, 
				 Ti iedgeset,
				 cudaStream_t stream=0)
{
  Td *Dx = (Td*)(*d_Dx);
  Td *Dy = (Td*)(*d_Dy);
  Td *CoeffsAtEdge = (Td*)(*d_CoeffsAtEdge);
  Ti *IedgeList = (Ti*)(*d_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
  if (nblocks == 1) {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_SEGREGATED,DISSIPATION_ROE>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  } else {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_ALLCOUPLED,DISSIPATION_ROE>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  }
  coproc_checkErrors("hydro_calcFluxRoeDiss3d_cuda");
  return 0;
}

/*******************************************************************************/
extern "C" {
  __INT FNAME(hydro_calcfluxroediss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					    __SIZET *d_IedgeList,
					    __SIZET *d_Dx,
					    __SIZET *d_Dy,
#ifdef HAS_CUDADOUBLEPREC
					    __DP *scale,
#else
					    __SP *scale,
#endif
					    __INT *nblocks, 
					    __INT *neq, 
					    __INT *nvar,
					    __INT *nedge, 
					    __INT *ncoeff,
					    __INT *nedges, 
					    __INT *iedgeset,
					    __I64 *stream=0)
  {
    return (__INT) hydro_calcFluxRoeDiss3d_cuda(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
						*scale, *nblocks, *neq, *nvar, *nedge,
						*ncoeff, *nedges, *iedgeset,
						(cudaStream_t)(*stream));
  }
}

/*******************************************************************************/

template <typename Td, typename Ti>
inline
int hydro_calcFluxRoeDissDiSp3d_cuda(__I64 *d_CoeffsAtEdge,
				     __I64 *d_IedgeList,
				     __I64 *d_Dx,
				     __I64 *d_Dy,
				     Td scale,
				     Ti nblocks, 
				     Ti neq, 
				     Ti nvar,
				     Ti nedge,
				     Ti ncoeff,
				     Ti nedges, 
				     Ti iedgeset,
				     cudaStream_t stream=0)
{
  Td *Dx = (Td*)(*d_Dx);
  Td *Dy = (Td*)(*d_Dy);
  Td *CoeffsAtEdge = (Td*)(*d_CoeffsAtEdge);
  Ti *IedgeList = (Ti*)(*d_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
  if (nblocks == 1) {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_SEGREGATED,DISSIPATION_ROE_DSPLIT>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  } else {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_ALLCOUPLED,DISSIPATION_ROE_DSPLIT>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  }
  coproc_checkErrors("hydro_calcFluxRoeDissDiSp3d_cuda");
  return 0;
}

/*******************************************************************************/
extern "C" {
  __INT FNAME(hydro_calcfluxroedissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
						__SIZET *d_IedgeList,
						__SIZET *d_Dx,
						__SIZET *d_Dy,
#ifdef HAS_CUDADOUBLEPREC
						__DP *scale,
#else
						__SP *scale,
#endif
						__INT *nblocks, 
						__INT *neq, 
						__INT *nvar,
						__INT *nedge, 
						__INT *ncoeff,
						__INT *nedges, 
						__INT *iedgeset,
						__I64 *stream=0)
  {
    return (__INT) hydro_calcFluxRoeDissDiSp3d_cuda(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
						    *scale, *nblocks, *neq, *nvar, *nedge,
						    *ncoeff, *nedges, *iedgeset,
						    (cudaStream_t)*stream);
  }
}

/*******************************************************************************/

template <typename Td, typename Ti>
inline
int hydro_calcFluxRusDiss3d_cuda(__I64 *d_CoeffsAtEdge,
				 __I64 *d_IedgeList,
				 __I64 *d_Dx,
				 __I64 *d_Dy,
				 Td scale,
				 Ti nblocks, 
				 Ti neq, 
				 Ti nvar,
				 Ti nedge, 
				 Ti ncoeff,
				 Ti nedges, 
				 Ti iedgeset,
				 cudaStream_t stream=0)
{
  Td *Dx = (Td*)(*d_Dx);
  Td *Dy = (Td*)(*d_Dy);
  Td *CoeffsAtEdge = (Td*)(*d_CoeffsAtEdge);
  Ti *IedgeList = (Ti*)(*d_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
  if (nblocks == 1) {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_SEGREGATED,DISSIPATION_RUSANOV>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  } else {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_ALLCOUPLED,DISSIPATION_RUSANOV>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  }
  coproc_checkErrors("hydro_calcFluxRusDiss3d_cuda");
  return 0;
}

/*******************************************************************************/
extern "C" {
  __INT FNAME(hydro_calcfluxrusdiss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					    __SIZET *d_IedgeList,
					    __SIZET *d_Dx,
					    __SIZET *d_Dy,
#ifdef HAS_CUDADOUBLEPREC
					    __DP *scale,
#else
					    __SP *scale,
#endif
					    __INT *nblocks, 
					    __INT *neq, 
					    __INT *nvar,
					    __INT *nedge, 
					    __INT *ncoeff,
					    __INT *nedges, 
					    __INT *iedgeset,
					    __I64 *stream=0)
  {
    return (__INT)hydro_calcFluxRusDiss3d_cuda(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
					       *scale, *nblocks, *neq, *nvar, *nedge,
					       *ncoeff, *nedges, *iedgeset,
					       (cudaStream_t)*stream);
  }
}

/*******************************************************************************/

template <typename Td, typename Ti>
inline
int hydro_calcFluxRusDissDiSp3d_cuda(__I64 *d_CoeffsAtEdge,
				     __I64 *d_IedgeList,
				     __I64 *d_Dx,
				     __I64 *d_Dy,
				     Td scale,
				     Ti nblocks, 
				     Ti neq,
				     Ti nvar,
				     Ti nedge, 
				     Ti ncoeff,
				     Ti nedges, 
				     Ti iedgeset,
				     cudaStream_t stream=0)
{
  Td *Dx = (Td*)(*d_Dx);
  Td *Dy = (Td*)(*d_Dy);
  Td *CoeffsAtEdge = (Td*)(*d_CoeffsAtEdge);
  Ti *IedgeList = (Ti*)(*d_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
  if (nblocks == 1) {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_SEGREGATED,DISSIPATION_RUSANOV_DSPLIT>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  } else {
    hydro_calcFlux3d_knl
      <Td,Ti,SYSTEM_ALLCOUPLED,DISSIPATION_RUSANOV_DSPLIT>
      <<<grid, block, 0, stream>>>(CoeffsAtEdge,
				   IedgeList,
				   Dx, Dy, scale, 
				   neq, nvar,
				   nedge, ncoeff,
				   nedges, iedgeset);
  }
  coproc_checkErrors("hydro_calcFluxRusDissDiSp3d_cuda");
  return 0;
}

/*******************************************************************************/
extern "C" {
  __INT FNAME(hydro_calcfluxrusdissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
						__SIZET *d_IedgeList,
						__SIZET *d_Dx,
						__SIZET *d_Dy,
#ifdef HAS_CUDADOUBLEPREC
						__DP *scale,
#else
						__SP *scale,
#endif
						__INT *nblocks, 
						__INT *neq, 
						__INT *nvar,
						__INT *nedge, 
						__INT *ncoeff,
						__INT *nedges, 
						__INT *iedgeset,
						__I64 *stream)
  {
    return (__INT) hydro_calcFluxRusDissDiSp3d_cuda(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
						    *scale, *nblocks, *neq, *nvar, *nedge,
						    *ncoeff, *nedges, *iedgeset,
						    (cudaStream_t)*stream);
  }
}
