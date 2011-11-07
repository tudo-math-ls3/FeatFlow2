/*#############################################################################
 ******************************************************************************
 * <name> hydro_calcFluxScDiss2d_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This CUDA kernel computes the fluxes for the low-order scheme in 2D
 * using tensorial artificial viscosities of Roe-type.
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

#define HYDRO_NDIM 2
#include "hydro.h"

extern "C"
{
  int hydro_calcFluxRoeDiss2d_cuda(unsigned long * h_DmatrixCoeffsAtEdge,
				   unsigned long * h_IverticesAtEdge,
				   unsigned long * h_Dx,
				   unsigned long * h_Dy,
				   double * dscale,
				   int * nblocks,
				   int * neq,
				   int * nvar,
				   int * nedge,
				   int * nmatcoeff,
				   int * nedges,
				   int * iedgeset);
  int FNAME(hydro_calcfluxroediss2d_cuda)(unsigned long * h_DmatrixCoeffsAtEdge,
					  unsigned long * h_IverticesAtEdge,
					  unsigned long * h_Dx,
					  unsigned long * h_Dy,
					  double * dscale,
					  int * nblocks,
					  int * neq,
					  int * nvar,
					  int * nedge,
					  int * nmatcoeff,
					  int * nedges,
					  int * iedgeset);
}

/*******************************************************************************/
template <int isystemformat>
__global__ void hydro_calcFluxRoeDiss2d_knl(double * DmatrixCoeffsAtEdge,
					    int * IverticesAtEdge,
					    double * Dx,
					    double * Dy,
					    double dscale,
					    int neq,
					    int nvar,
					    int nedge,
					    int nmatcoeff,
					    int nedges,
					    int iedgeset)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx<nedges)
  {
    // Get positions of edge endpoints (idx starts at zero)
    int i = IDX2T(IverticesAtEdge,1,iedgeset+idx,4,nedge);
    int j = IDX2T(IverticesAtEdge,2,iedgeset+idx,4,nedge);

    // Get solution values at edge endpoints
    double DdataAtEdge[2*NVAR2D];

    if (isystemformat == 0) {
      // Solution vector is stored in interleaved format
      IDX2(DdataAtEdge,1,1,NVAR2D,2) = IDX2_REVERSE(Dx,1,i,NVAR2D,neq);
      IDX2(DdataAtEdge,2,1,NVAR2D,2) = IDX2_REVERSE(Dx,2,i,NVAR2D,neq);
      IDX2(DdataAtEdge,3,1,NVAR2D,2) = IDX2_REVERSE(Dx,3,i,NVAR2D,neq);
      IDX2(DdataAtEdge,4,1,NVAR2D,2) = IDX2_REVERSE(Dx,4,i,NVAR2D,neq);
      
      IDX2(DdataAtEdge,1,2,NVAR2D,2) = IDX2_REVERSE(Dx,1,j,NVAR2D,neq);
      IDX2(DdataAtEdge,2,2,NVAR2D,2) = IDX2_REVERSE(Dx,2,j,NVAR2D,neq);
      IDX2(DdataAtEdge,3,2,NVAR2D,2) = IDX2_REVERSE(Dx,3,j,NVAR2D,neq);
      IDX2(DdataAtEdge,4,2,NVAR2D,2) = IDX2_REVERSE(Dx,4,j,NVAR2D,neq);

    } else {
      // Solution vector is stored in block format
      IDX2(DdataAtEdge,1,1,NVAR2D,2) = IDX2_FORWARD(Dx,1,i,NVAR2D,neq);
      IDX2(DdataAtEdge,2,1,NVAR2D,2) = IDX2_FORWARD(Dx,2,i,NVAR2D,neq);
      IDX2(DdataAtEdge,3,1,NVAR2D,2) = IDX2_FORWARD(Dx,3,i,NVAR2D,neq);
      IDX2(DdataAtEdge,4,1,NVAR2D,2) = IDX2_FORWARD(Dx,4,i,NVAR2D,neq);

      IDX2(DdataAtEdge,1,2,NVAR2D,2) = IDX2_FORWARD(Dx,1,j,NVAR2D,neq);
      IDX2(DdataAtEdge,2,2,NVAR2D,2) = IDX2_FORWARD(Dx,2,j,NVAR2D,neq);
      IDX2(DdataAtEdge,3,2,NVAR2D,2) = IDX2_FORWARD(Dx,3,j,NVAR2D,neq);
      IDX2(DdataAtEdge,4,2,NVAR2D,2) = IDX2_FORWARD(Dx,4,j,NVAR2D,neq);
    }

    //--------------------------------------------------------------------------
    // Evaluate the Galerkin fluxes
    //--------------------------------------------------------------------------

    // Compute velocities
    double ui = XVELOCITY3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1);
    double vi = YVELOCITY3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1);

    double uj = XVELOCITY3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1);
    double vj = YVELOCITY3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1);

    // Compute pressures
    double pi = PRESSURE3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1);
    double pj = PRESSURE3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1);

#ifdef HYDRO_USE_IBP
    double Fxi[NVAR2D];
    double Fxj[NVAR2D];

    // Compute fluxes for x-direction
    IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi);
    IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi);
    IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi);
    IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi);

    IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);
    IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);
    IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);
    IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);

    double Fyi[NVAR2D];
    double Fyj[NVAR2D];

    // Compute fluxes for x-direction
    IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi);
    IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi);
    IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi);
    IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi);
      
    IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
    IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
    IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
    IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
#else
    double Fx_ij[NVAR2D];
    
    // Compute flux difference for x-direction
    IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi)-
                    INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);
    IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi)-
                    INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);
    IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi)-
                    INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);
    IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,ui,pi)-
                    INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,uj,pj);

    double Fy_ij[NVAR2D];

    // Compute flux difference for y-direction
    IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi)-
                    INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
    IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi)-
                    INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
    IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi)-
                    INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
    IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1,vi,pi)-
                    INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1,vj,pj);
#endif

    //--------------------------------------------------------------------------
    // Evaluate the scalar dissipation tensor of Roe-type
    //--------------------------------------------------------------------------

    // Compute skew-symmetric coefficient
    double a[HYDRO_NDIM];
    a[0] = RCONST(0.5)*(IDX3T(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
			IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge));
    a[1] = RCONST(0.5)*(IDX3T(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
			IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge));
    double anorm = sqrt(a[0] * a[0] + a[1] * a[1]);

    double Diff[NVAR2D];
    if (anorm > 1e-14) {

      // Normalize the skew-symmetric coefficient
      a[0] = a[0]/anorm;
      a[1] = a[1]/anorm;
      
      // Compute densities
      double ri = DENSITY3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1);
      double rj = DENSITY3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1);
      
      // Compute enthalpies
      double hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,1,NVAR2D,2,1)+pi)/ri;
      double hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,1,NVAR2D,2,1)+pj)/rj;
      
      //! Compute Roe mean values
      double aux  = ROE_MEAN_RATIO(ri,rj);
      double u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      double v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      double H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
      // Compute auxiliary variables
      double vel_ij = u_ij * a[0] + v_ij * a[1];
      double q_ij   = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij);
      
      // Compute the speed of sound
      //TODO echtes double epsilon einbauen
      double c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), 1e-14);
      double c_ij  = sqrt(c2_ij);
      
      // Compute eigenvalues
      double l1 = abs(vel_ij-c_ij);
      double l2 = abs(vel_ij);
      double l3 = abs(vel_ij+c_ij);
      double l4 = abs(vel_ij);
      
      // Compute solution difference U_j-U_i
      Diff[0] = IDX3(DdataAtEdge,1,2,1,NVAR2D,2,1)-IDX3(DdataAtEdge,1,1,1,NVAR2D,2,1);
      Diff[1] = IDX3(DdataAtEdge,2,2,1,NVAR2D,2,1)-IDX3(DdataAtEdge,2,1,1,NVAR2D,2,1);
      Diff[2] = IDX3(DdataAtEdge,3,2,1,NVAR2D,2,1)-IDX3(DdataAtEdge,3,1,1,NVAR2D,2,1);
      Diff[3] = IDX3(DdataAtEdge,4,2,1,NVAR2D,2,1)-IDX3(DdataAtEdge,4,1,1,NVAR2D,2,1);
      
      // Compute auxiliary quantities for characteristic variables
      double aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*Diff[0]
					        -u_ij*Diff[1]
					        -v_ij*Diff[2]
					             +Diff[3])/RCONST(2.0)/c2_ij;
      double aux2 = (vel_ij*Diff[0]
		      -a[0]*Diff[1]
		      -a[1]*Diff[2])/RCONST(2.0)/c_ij;

      // Compute characteristic variables multiplied by the corresponding eigenvalue
      double w1 = l1 * (aux1 + aux2);
      double w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff[0]
                                          +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff[1]
                                                                       +v_ij*Diff[2]
								            -Diff[3])/c2_ij);
      double w3 = l3 * (aux1 - aux2);
      double w4 = l4 * ((a[0]*v_ij-a[1]*u_ij)*Diff[0]
                                        +a[1]*Diff[1]
		                        -a[0]*Diff[2]);

      // Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff[0] = anorm * ( w1 + w2 + w3 );
      Diff[1] = anorm * ( (u_ij-c_ij*a[0])*w1 + u_ij*w2 +
			  (u_ij+c_ij*a[0])*w3 + a[1]*w4 );
      Diff[2] = anorm * ( (v_ij-c_ij*a[1])*w1 + v_ij*w2 +
			  (v_ij+c_ij*a[1])*w3 - a[0]*w4 );
      Diff[3] = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 +
			  (H_ij+c_ij*vel_ij)*w3 + (u_ij*a[1]-v_ij*a[0])*w4 );
      
    } else {
      Diff[0] = 0.0;
      Diff[1] = 0.0;
      Diff[2] = 0.0;
      Diff[3] = 0.0;
    }

    //--------------------------------------------------------------------------
    // Build both contributions into the fluxes
    //--------------------------------------------------------------------------
    
#ifdef HYDRO_USE_IBP
    double DfluxesAtEdge[2*NVAR2D];
    IDX3(DfluxesAtEdge,1,1,1,NVAR2D,2,1) = dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[0]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[0]-
       IDX3T(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[0]-
       IDX3T(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[0] + Diff[0]);
    
    IDX3(DfluxesAtEdge,2,1,1,NVAR2D,2,1) = dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[1]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[1]-
       IDX3T(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[1]-
       IDX3T(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[1] + Diff[1]);
    
    IDX3(DfluxesAtEdge,3,1,1,NVAR2D,2,1) = dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[2]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[2]-
       IDX3T(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[2]-
       IDX3T(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[2] + Diff[2]);
    
    IDX3(DfluxesAtEdge,4,1,1,NVAR2D,2,1) = dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[3]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[3]-
       IDX3T(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[3]-
       IDX3T(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[3] + Diff[3]);
    
    
    IDX3(DfluxesAtEdge,1,2,1,NVAR2D,2,1) = -IDX3(DfluxesAtEdge,1,1,1,NVAR2D,2,1);
    IDX3(DfluxesAtEdge,2,2,1,NVAR2D,2,1) = -IDX3(DfluxesAtEdge,2,1,1,NVAR2D,2,1);
    IDX3(DfluxesAtEdge,3,2,1,NVAR2D,2,1) = -IDX3(DfluxesAtEdge,3,1,1,NVAR2D,2,1);
    IDX3(DfluxesAtEdge,4,2,1,NVAR2D,2,1) = -IDX3(DfluxesAtEdge,4,1,1,NVAR2D,2,1);
#else
    double DfluxesAtEdge[2*NVAR2D];    
    IDX3(DfluxesAtEdge,1,1,1,NVAR2D,2,1) = dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[0]+
       IDX3T(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[0] + Diff[0]);
    
    IDX3(DfluxesAtEdge,2,1,1,NVAR2D,2,1) = dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[1]+
       IDX3T(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[1] + Diff[1]);
    
    IDX3(DfluxesAtEdge,3,1,1,NVAR2D,2,1) = dscale *
      (IDX3(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[2]+
       IDX3(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[2] + Diff[2]);
    
    IDX3(DfluxesAtEdge,4,1,1,NVAR2D,2,1) = dscale *
      (IDX3(DmatrixCoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[3]+
       IDX3(DmatrixCoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[3] + Diff[3]);
    
    
    IDX3(DfluxesAtEdge,1,2,1,NVAR2D,2,1) = -dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[0]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[0] + Diff[0]);
    
    IDX3(DfluxesAtEdge,2,2,1,NVAR2D,2,1) = -dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[1]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[1] + Diff[1]);
    
    IDX3(DfluxesAtEdge,3,2,1,NVAR2D,2,1) = -dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[2]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[2] + Diff[2]);
    
    IDX3(DfluxesAtEdge,4,2,1,NVAR2D,2,1) = -dscale *
      (IDX3T(DmatrixCoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[3]+
       IDX3T(DmatrixCoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[3] + Diff[3]);
#endif
    
    //--------------------------------------------------------------------------
    // Build fluxes into nodal vector
    //--------------------------------------------------------------------------

    if (isystemformat == 0) {
      // Solution vector is stored in interleaved format
      IDX2_REVERSE(Dy,1,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,1,1,1,NVAR2D,2,1);
      IDX2_REVERSE(Dy,2,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,2,1,1,NVAR2D,2,1);
      IDX2_REVERSE(Dy,3,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,3,1,1,NVAR2D,2,1);
      IDX2_REVERSE(Dy,4,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,4,1,1,NVAR2D,2,1);
      
      IDX2_REVERSE(Dy,1,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,1,2,1,NVAR2D,2,1);
      IDX2_REVERSE(Dy,2,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,2,2,1,NVAR2D,2,1);
      IDX2_REVERSE(Dy,3,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,3,2,1,NVAR2D,2,1);
      IDX2_REVERSE(Dy,4,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,4,2,1,NVAR2D,2,1);

    } else {
      // Solution vector is stored in block format
      IDX2_FORWARD(Dy,1,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,1,1,1,NVAR2D,2,1);
      IDX2_FORWARD(Dy,2,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,2,1,1,NVAR2D,2,1);
      IDX2_FORWARD(Dy,3,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,3,1,1,NVAR2D,2,1);
      IDX2_FORWARD(Dy,4,i,NVAR2D,neq) += IDX3(DfluxesAtEdge,4,1,1,NVAR2D,2,1);

      IDX2_FORWARD(Dy,1,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,1,2,1,NVAR2D,2,1);
      IDX2_FORWARD(Dy,2,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,2,2,1,NVAR2D,2,1);
      IDX2_FORWARD(Dy,3,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,3,2,1,NVAR2D,2,1);
      IDX2_FORWARD(Dy,4,j,NVAR2D,neq) += IDX3(DfluxesAtEdge,4,2,1,NVAR2D,2,1);
    }
  }
}

/*******************************************************************************/

int hydro_calcFluxRoeDiss2d_cuda(unsigned long * h_DmatrixCoeffsAtEdge,
				 unsigned long * h_IverticesAtEdge,
				 unsigned long * h_Dx,
				 unsigned long * h_Dy,
				 double * dscale,
				 int * nblocks, int * neq, int * nvar,
				 int * nedge, int * nmatcoeff,
				 int * nedges, int * iedgeset)
{
  double * d_Dx = (double*)(*h_Dx);
  double * d_Dy = (double*)(*h_Dy);
  double * d_DmatrixCoeffsAtEdge = (double*)(*h_DmatrixCoeffsAtEdge);
  int * d_IverticesAtEdge = (int*)(*h_IverticesAtEdge);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((*nedges)/(double)(block.x));

  if (*nblocks == 1) {
    hydro_calcFluxRoeDiss2d_knl<0><<<grid, block>>>(d_DmatrixCoeffsAtEdge,
						   d_IverticesAtEdge,
						   d_Dx, d_Dy, (*dscale), 
						   (*neq), (*nvar),
						   (*nedge), (*nmatcoeff),
						   (*nedges), (*iedgeset));
  } else {
    hydro_calcFluxRoeDiss2d_knl<1><<<grid, block>>>(d_DmatrixCoeffsAtEdge,
						   d_IverticesAtEdge,
						   d_Dx, d_Dy, (*dscale), 
						   (*neq), (*nvar),
						   (*nedge), (*nmatcoeff),
						   (*nedges), (*iedgeset));
  }
  coproc_checkErrors("hydro_calcFluxRoeDiss2d_cuda");
  return 0;
}

int FNAME(hydro_calcfluxroediss2d_cuda)(unsigned long * h_DmatrixCoeffsAtEdge,
					unsigned long * h_IverticesAtEdge,
					unsigned long * h_Dx,
					unsigned long * h_Dy,
					double * dscale,
					int * nblocks, int * neq, int * nvar,
					int * nedge,   int * nmatcoeff,
					int * nedges,  int * iedgeset)
{
  return hydro_calcFluxRoeDiss2d_cuda(h_DmatrixCoeffsAtEdge, h_IverticesAtEdge,
				      h_Dx, h_Dy, dscale, nblocks, neq, nvar,
				      nedge, nmatcoeff, nedges, iedgeset);
}
