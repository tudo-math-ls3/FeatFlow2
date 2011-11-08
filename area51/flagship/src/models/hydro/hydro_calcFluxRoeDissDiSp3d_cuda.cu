/*#############################################################################
 ******************************************************************************
 * <name>  hydro_calcFluxRoeDissDiSp3d_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This CUDA kernel computes the fluxes for the low-order scheme in 3D
 * using tensorial artificial viscosities of Roe-type, whereby dimensional
 * splitting is employed.
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

extern "C"
{
  int hydro_calcFluxRoeDissDiSp3d_cuda(unsigned long * h_DcoeffsAtEdge,
				       unsigned long * h_IedgeList,
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
  int FNAME(hydro_calcfluxroedissdisp3d_cuda)(unsigned long * h_DcoeffsAtEdge,
					      unsigned long * h_IedgeList,
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
__global__ void hydro_calcFluxRoeDissDiSp3d_knl(double * DcoeffsAtEdge,
						int * IedgeList,
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
      int i = IDX2T(IedgeList,1,iedgeset+idx,6,nedge);
      int j = IDX2T(IedgeList,2,iedgeset+idx,6,nedge);
      
      // Get solution values at edge endpoints
      double DdataAtEdge[2*NVAR3D];
      
      if (isystemformat == 0) {
	// Solution vector is stored in interleaved format
	IDX2(DdataAtEdge,1,1,NVAR3D,2) = IDX2_REVERSE(Dx,1,i,NVAR3D,neq);
	IDX2(DdataAtEdge,2,1,NVAR3D,2) = IDX2_REVERSE(Dx,2,i,NVAR3D,neq);
	IDX2(DdataAtEdge,3,1,NVAR3D,2) = IDX2_REVERSE(Dx,3,i,NVAR3D,neq);
	IDX2(DdataAtEdge,4,1,NVAR3D,2) = IDX2_REVERSE(Dx,4,i,NVAR3D,neq);
	IDX2(DdataAtEdge,5,1,NVAR3D,2) = IDX2_REVERSE(Dx,5,i,NVAR3D,neq);
	
	IDX2(DdataAtEdge,1,2,NVAR3D,2) = IDX2_REVERSE(Dx,1,j,NVAR3D,neq);
	IDX2(DdataAtEdge,2,2,NVAR3D,2) = IDX2_REVERSE(Dx,2,j,NVAR3D,neq);
	IDX2(DdataAtEdge,3,2,NVAR3D,2) = IDX2_REVERSE(Dx,3,j,NVAR3D,neq);
	IDX2(DdataAtEdge,4,2,NVAR3D,2) = IDX2_REVERSE(Dx,4,j,NVAR3D,neq);
	IDX2(DdataAtEdge,5,2,NVAR3D,2) = IDX2_REVERSE(Dx,5,j,NVAR3D,neq);
	
      } else {
	// Solution vector is stored in block format
	IDX2(DdataAtEdge,1,1,NVAR3D,2) = IDX2_FORWARD(Dx,1,i,NVAR3D,neq);
	IDX2(DdataAtEdge,2,1,NVAR3D,2) = IDX2_FORWARD(Dx,2,i,NVAR3D,neq);
	IDX2(DdataAtEdge,3,1,NVAR3D,2) = IDX2_FORWARD(Dx,3,i,NVAR3D,neq);
	IDX2(DdataAtEdge,4,1,NVAR3D,2) = IDX2_FORWARD(Dx,4,i,NVAR3D,neq);
	IDX2(DdataAtEdge,5,1,NVAR3D,2) = IDX2_FORWARD(Dx,5,i,NVAR3D,neq);
	
	IDX2(DdataAtEdge,1,2,NVAR3D,2) = IDX2_FORWARD(Dx,1,j,NVAR3D,neq);
	IDX2(DdataAtEdge,2,2,NVAR3D,2) = IDX2_FORWARD(Dx,2,j,NVAR3D,neq);
	IDX2(DdataAtEdge,3,2,NVAR3D,2) = IDX2_FORWARD(Dx,3,j,NVAR3D,neq);
	IDX2(DdataAtEdge,4,2,NVAR3D,2) = IDX2_FORWARD(Dx,4,j,NVAR3D,neq);
	IDX2(DdataAtEdge,5,2,NVAR3D,2) = IDX2_FORWARD(Dx,5,j,NVAR3D,neq);
      }

      //------------------------------------------------------------------------
      // Evaluate the Galerkin fluxes
      //------------------------------------------------------------------------
      
      // Compute velocities
      double ui = XVELOCITY3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1);
      double vi = YVELOCITY3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1);
      double wi = ZVELOCITY3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1);
      
      double uj = XVELOCITY3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1);
      double vj = YVELOCITY3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1);
      double wj = ZVELOCITY3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1);
      
      // Compute pressures
      double pi = PRESSURE3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1);
      double pj = PRESSURE3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1);
      
#ifdef HYDRO_USE_IBP
      double Fxi[NVAR3D];
      double Fxj[NVAR3D];
      
      // Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi);
      
      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      
      double Fyi[NVAR3D];
      double Fyj[NVAR3D];
      
      // Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      
      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      
      double Fzi[NVAR3D];
      double Fzj[NVAR3D];
      
      // Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi);
      
      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
#else
      double Fx_ij[NVAR3D];
      
      // Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                      INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                      INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                      INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                      INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,ui,pi)-
                      INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,uj,pj);
      
      double Fy_ij[NVAR3D];
      
      // Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);

      double Fz_ij[NVAR3D];
      
      // Compute flux difference for y-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
	              INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1,vi,pi)-
                      INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1,vj,pj);
#endif

      //------------------------------------------------------------------------
      // Evaluate the dissipation tensor of Roe-type
      //------------------------------------------------------------------------
      
      // Compute skew-symmetric coefficient
      double a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
			  IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
			  IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge));
      a[2] = RCONST(0.5)*(IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
			  IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge));
      double anorm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
      
      double DiffX[NVAR3D];
      double DiffY[NVAR3D];
      double DiffZ[NVAR3D];
      if (anorm > 1e-14) {
	
	// Compute the absolute value
	a[0] = abs(a[0]);
	a[1] = abs(a[1]);
	a[2] = abs(a[2]);
	
	// Compute densities
	double ri = DENSITY3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1);
	double rj = DENSITY3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1);
	
	// Compute enthalpies
	double hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1)+pi)/ri;
	double hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1)+pj)/rj;
	
	// Compute Roe mean values
	double aux  = ROE_MEAN_RATIO(ri,rj);
	double u_ij = ROE_MEAN_VALUE(ui,uj,aux);
	double v_ij = ROE_MEAN_VALUE(vi,vj,aux);
	double w_ij = ROE_MEAN_VALUE(wi,wj,aux);
	double H_ij = ROE_MEAN_VALUE(hi,hj,aux);
	
	// Compute auxiliary variable
	double q_ij   = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
	
	// Compute the speed of sound
	//TODO echtes double epsilon einbauen
	double c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), 1e-14);
	double c_ij  = sqrt(c2_ij);
	
	//----------------------------------------------------------------------
	// Dimensional splitting: x-direction
	//----------------------------------------------------------------------
        
	// Compute eigenvalues
	double l1 = abs(u_ij-c_ij);
	double l2 = abs(u_ij);
	double l3 = abs(u_ij+c_ij);
	double l4 = abs(u_ij);
	double l5 = abs(u_ij);
        
	// Compute solution difference U_j-U_i
	DiffX[0] = IDX3(DdataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,1,1,1,NVAR3D,2,1);
	DiffX[1] = IDX3(DdataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,2,1,1,NVAR3D,2,1);
	DiffX[2] = IDX3(DdataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,3,1,1,NVAR3D,2,1);
	DiffX[3] = IDX3(DdataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,4,1,1,NVAR3D,2,1);
	DiffX[4] = IDX3(DdataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,5,1,1,NVAR3D,2,1);

	// Compute auxiliary quantities for characteristic variables
        double aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*DiffX[0]
					          -u_ij*DiffX[1]
						  -v_ij*DiffX[2]
						  -w_ij*DiffX[3]
					               +DiffX[4])/RCONST(2.0)/c2_ij;
        double aux2 = (u_ij*DiffX[0]-DiffX[1])/RCONST(2.0)/c_ij;

        // Compute characteristic variables multiplied by the corresponding eigenvalue
        double w1 = l1 * (aux1 + aux2);
        double w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffX[0]
                                           +((HYDRO_GAMMA)-RCONST(1.0))*( u_ij*DiffX[1]
                                                                         +v_ij*DiffX[2]
                                                                         +w_ij*DiffX[3]
									      -DiffX[4])/c2_ij);
        double w3 = l3 * (aux1 - aux2);
        double w4 = l4 * ( v_ij*DiffX[0]-DiffX[2]);
        double w5 = l5 * (-w_ij*DiffX[0]+DiffX[3]);

        // Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffX[0] = a[0] * ( w1 + w2 + w3 );
        DiffX[1] = a[0] * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 );
        DiffX[2] = a[0] * ( v_ij*(w1 + w2 + w3) - w4 );
        DiffX[3] = a[0] * ( w_ij*(w1 + w2 + w3) + w5 );
        DiffX[4] = a[0] * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 + (H_ij+c_ij*u_ij)*w3
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
	DiffY[0] = IDX3(DdataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,1,1,1,NVAR3D,2,1);
	DiffY[1] = IDX3(DdataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,2,1,1,NVAR3D,2,1);
	DiffY[2] = IDX3(DdataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,3,1,1,NVAR3D,2,1);
	DiffY[3] = IDX3(DdataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,4,1,1,NVAR3D,2,1);
	DiffY[4] = IDX3(DdataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,5,1,1,NVAR3D,2,1);
	
	// Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*DiffY[0]
                                           -u_ij*DiffY[1]
                                           -v_ij*DiffY[2]
                                           -w_ij*DiffY[3]
					        +DiffY[4])/RCONST(2.0)/c2_ij;
        aux2 = (v_ij*DiffY[0]-DiffY[2])/RCONST(2.0)/c_ij;

        // Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2);
        w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffY[0]
                                     +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*DiffY[1]
                                                                  +v_ij*DiffY[2]
                                                                  +w_ij*DiffY[3]
								       -DiffY[4])/c2_ij);
        w3 = l3 * (aux1 - aux2);
        w4 = l4 * (-u_ij*DiffY[0]+DiffY[1]);
        w5 = l5 * ( w_ij*DiffY[0]-DiffY[3]);

        // Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffY[0] = a[1] * ( w1 + w2 + w3 );
        DiffY[1] = a[1] * ( u_ij*(w1 + w2 + w3) + w4 );
        DiffY[2] = a[1] * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 );
        DiffY[3] = a[1] * ( w_ij*(w1 + w2 + w3) - w5 );
        DiffY[4] = a[1] * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 + (H_ij+c_ij*v_ij)*w3
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
	DiffZ[0] = IDX3(DdataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,1,1,1,NVAR3D,2,1);
	DiffZ[1] = IDX3(DdataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,2,1,1,NVAR3D,2,1);
	DiffZ[2] = IDX3(DdataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,3,1,1,NVAR3D,2,1);
	DiffZ[3] = IDX3(DdataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,4,1,1,NVAR3D,2,1);
	DiffZ[4] = IDX3(DdataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,5,1,1,NVAR3D,2,1);
	
	// Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-RCONST(1.0))*(q_ij*DiffZ[0]
                                           -u_ij*DiffZ[1]
                                           -v_ij*DiffZ[2]
                                           -w_ij*DiffZ[3]
					        +DiffZ[4])/RCONST(2.0)/c2_ij;
        aux2 = (w_ij*DiffZ[0]-DiffZ[2])/RCONST(2.0)/c_ij;

        // Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2);
        w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffZ[0]
                                     +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*DiffZ[1]
                                                                  +v_ij*DiffZ[2]
                                                                  +w_ij*DiffZ[3]
								       -DiffZ[4])/c2_ij);
        w3 = l3 * (aux1 - aux2);
        w4 = l4 * ( u_ij*DiffZ[0]-DiffZ[1]);
        w5 = l5 * (-v_ij*DiffZ[0]+DiffZ[2]);
        
        // Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffZ[0] = a[2] * ( w1 + w2 + w3 );
        DiffZ[1] = a[2] * ( u_ij*(w1 + w2 + w3) - w4 );
        DiffZ[2] = a[2] * ( v_ij*(w1 + w2 + w3) + w5 );
        DiffZ[3] = a[2] * ( (w_ij-c_ij)*w1 + w_ij*w2 + (w_ij+c_ij)*w3 );
        DiffZ[4] = a[2] * ( (H_ij-c_ij*w_ij)*w1 + q_ij*w2 + (H_ij+c_ij*w_ij)*w3
                            -u_ij*w4 + v_ij*w5 );
      } else {
	DiffX[0] = 0.0;
	DiffX[1] = 0.0;
	DiffX[2] = 0.0;
	DiffX[3] = 0.0;
	DiffX[4] = 0.0;
	
	DiffY[0] = 0.0;
	DiffY[1] = 0.0;
	DiffY[2] = 0.0;
	DiffY[3] = 0.0;
	DiffY[4] = 0.0;

	DiffZ[0] = 0.0;
	DiffZ[1] = 0.0;
	DiffZ[2] = 0.0;
	DiffZ[3] = 0.0;
	DiffZ[4] = 0.0;
      }
      
      //------------------------------------------------------------------------
      // Build both contributions into the fluxes
      //------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
      double DfluxesAtEdge[2*NVAR3D];
      IDX3(DfluxesAtEdge,1,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[0]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[0]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[0]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[0]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[0]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[0]+
	 DiffX[0]+DiffY[0]+DiffZ[0]);
      
      IDX3(DfluxesAtEdge,2,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[1]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[1]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[1]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[1]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[1]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[1]+
	 DiffX[1]+DiffY[1]+DiffZ[1]);
      
      IDX3(DfluxesAtEdge,3,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[2]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[2]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[2]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[2]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[2]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[2]+
	 DiffX[2]+DiffY[2]+DiffZ[2]);
      
      IDX3(DfluxesAtEdge,4,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[3]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[3]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[3]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[3]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[3]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[3]+
	 DiffX[3]+DiffY[3]+DiffZ[3]);

      IDX3(DfluxesAtEdge,5,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[4]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[4]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[4]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[4]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[4]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[4]+
	 DiffX[4]+DiffY[4]+DiffZ[4]);
      
      
      IDX3(DfluxesAtEdge,1,2,1,NVAR3D,2,1) = -IDX3(DfluxesAtEdge,1,1,1,NVAR3D,2,1);
      IDX3(DfluxesAtEdge,2,2,1,NVAR3D,2,1) = -IDX3(DfluxesAtEdge,2,1,1,NVAR3D,2,1);
      IDX3(DfluxesAtEdge,3,2,1,NVAR3D,2,1) = -IDX3(DfluxesAtEdge,3,1,1,NVAR3D,2,1);
      IDX3(DfluxesAtEdge,4,2,1,NVAR3D,2,1) = -IDX3(DfluxesAtEdge,4,1,1,NVAR3D,2,1);
      IDX3(DfluxesAtEdge,5,2,1,NVAR3D,2,1) = -IDX3(DfluxesAtEdge,5,1,1,NVAR3D,2,1);
#else
      double DfluxesAtEdge[2*NVAR3D];    
      IDX3(DfluxesAtEdge,1,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[0]+
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[0]+
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[0]+
	 DiffX[0]+DiffY[0]+DiffZ[0]);
      
      IDX3(DfluxesAtEdge,2,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[1]+
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[1]+
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[1]+
	 DiffX[1]+DiffY[1]+DiffZ[1]);
      
      IDX3(DfluxesAtEdge,3,1,1,NVAR3D,2,1) = dscale *
	(IDX3(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[2]+
	 IDX3(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[2]+
	 IDX3(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[2]+
	 DiffX[2]+DiffY[2]+DiffZ[2]);
      
      IDX3(DfluxesAtEdge,4,1,1,NVAR3D,2,1) = dscale *
	(IDX3(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[3]+
	 IDX3(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[3]+
	 IDX3(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[3]+
	 DiffX[3]+DiffY[3]+DiffZ[3]);

      IDX3(DfluxesAtEdge,5,1,1,NVAR3D,2,1) = dscale *
	(IDX3(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[4]+
	 IDX3(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[4]+
	 IDX3(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[4]+
	 DiffX[4]+DiffY[4]+DiffZ[4]);
      
      
      IDX3(DfluxesAtEdge,1,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[0]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[0]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[0]+
	 DiffX[0]+DiffY[0]+DiffZ[0]);
      
      IDX3(DfluxesAtEdge,2,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[1]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[1]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[1]+
	 DiffX[1]+DiffY[1]+DiffZ[1]);
      
      IDX3(DfluxesAtEdge,3,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[2]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[2]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[2]+
	 DiffX[2]+DiffY[2]+DiffZ[2]);
      
      IDX3(DfluxesAtEdge,4,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[3]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[3]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[3]+
	 DiffX[3]+DiffY[3]+DiffZ[3]);
	 
      IDX3(DfluxesAtEdge,5,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[4]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[4]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[4]+
	 DiffX[4]+DiffY[4]+DiffZ[4]);
#endif
      
      //--------------------------------------------------------------------------
      // Build fluxes into nodal vector
      //--------------------------------------------------------------------------
      
       if (isystemformat == 0) {
	// Solution vector is stored in interleaved format
	IDX2_REVERSE(Dy,1,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,1,1,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,2,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,2,1,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,3,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,3,1,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,4,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,4,1,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,5,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,5,1,1,NVAR3D,2,1);
	
	IDX2_REVERSE(Dy,1,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,1,2,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,2,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,2,2,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,3,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,3,2,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,4,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,4,2,1,NVAR3D,2,1);
	IDX2_REVERSE(Dy,5,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,5,2,1,NVAR3D,2,1);
      } else {
	// Solution vector is stored in block format
	IDX2_FORWARD(Dy,1,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,1,1,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,2,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,2,1,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,3,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,3,1,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,4,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,4,1,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,5,i,NVAR3D,neq) += IDX3(DfluxesAtEdge,5,1,1,NVAR3D,2,1);
	
	IDX2_FORWARD(Dy,1,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,1,2,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,2,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,2,2,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,3,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,3,2,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,4,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,4,2,1,NVAR3D,2,1);
	IDX2_FORWARD(Dy,5,j,NVAR3D,neq) += IDX3(DfluxesAtEdge,5,2,1,NVAR3D,2,1);
      }
    }
}

/*******************************************************************************/

int hydro_calcFluxRoeDissDiSp3d_cuda(unsigned long * h_DcoeffsAtEdge,
				     unsigned long * h_IedgeList,
				     unsigned long * h_Dx,
				     unsigned long * h_Dy,
				     double * dscale,
				     int * nblocks, int * neq, int * nvar,
				     int * nedge, int * nmatcoeff,
				     int * nedges, int * iedgeset)
{
  double * d_Dx = (double*)(*h_Dx);
  double * d_Dy = (double*)(*h_Dy);
  double * d_DcoeffsAtEdge = (double*)(*h_DcoeffsAtEdge);
  int * d_IedgeList = (int*)(*h_IedgeList);
  
  // Define number of threads per block
  int blocksize = 128;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((*nedges)/(double)(block.x));
  
  if (*nblocks == 1) {
    hydro_calcFluxRoeDissDiSp3d_knl<0><<<grid, block>>>(d_DcoeffsAtEdge,
							d_IedgeList,
							d_Dx, d_Dy, (*dscale), 
							(*neq), (*nvar),
							(*nedge), (*nmatcoeff),
							(*nedges), (*iedgeset));
  } else {
    hydro_calcFluxRoeDissDiSp3d_knl<1><<<grid, block>>>(d_DcoeffsAtEdge,
							d_IedgeList,
							d_Dx, d_Dy, (*dscale), 
							(*neq), (*nvar),
							(*nedge), (*nmatcoeff),
							(*nedges), (*iedgeset));
  }
  coproc_checkErrors("hydro_calcFluxRoeDissDiSp3d_cuda");
  return 0;
}

int FNAME(hydro_calcfluxroedissdisp3d_cuda)(unsigned long * h_DcoeffsAtEdge,
					    unsigned long * h_IedgeList,
					    unsigned long * h_Dx,
					    unsigned long * h_Dy,
					    double * dscale,
					    int * nblocks, int * neq, int * nvar,
					    int * nedge,   int * nmatcoeff,
					    int * nedges,  int * iedgeset)
{
  return hydro_calcFluxRoeDissDiSp3d_cuda(h_DcoeffsAtEdge, h_IedgeList,
					  h_Dx, h_Dy, dscale, nblocks, neq, nvar,
					  nedge, nmatcoeff, nedges, iedgeset);
}
