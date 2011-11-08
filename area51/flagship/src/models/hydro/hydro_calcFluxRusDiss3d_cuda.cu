/*#############################################################################
 ******************************************************************************
 * <name> hydro_calcFluxRusDiss3d_cuda </name>
 ******************************************************************************
 *
 * <purpose>
 * This CUDA kernel computes the fluxes for the low-order scheme in 3D
 * using scalar artificial viscosities of Rusanov-type.
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
  int hydro_calcFluxRusDiss3d_cuda(unsigned long * h_DcoeffsAtEdge,
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
  int FNAME(hydro_calcfluxrusdiss3d_cuda)(unsigned long * h_DcoeffsAtEdge,
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
__global__ void hydro_calcFluxRusDiss3d_knl(double * DcoeffsAtEdge,
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
      
      // Compute flux difference for z-direction
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
      // Evaluate the scalar dissipation of Rusanov-type
      //------------------------------------------------------------------------
      
      // Compute specific energies
      double Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,1,NVAR3D,2,1);
      double Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,1,NVAR3D,2,1);
      
      // Compute the speed of sound
      double ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*(Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), 1e-14));
      double cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*(Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), 1e-14));
      
#ifdef HYDRO_USE_IBP
      // Compute scalar dissipation based on the skew-symmetric part
      // which does not include the symmetric boundary contribution
      double d_ij = max( abs(RCONST(0.5)*(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge))*uj+
			     RCONST(0.5)*(IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge))*vj-
			     RCONST(0.5)*(IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge))*wj)+
		     RCONST(0.5)*sqrt(POW(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
				      POW(IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
				      POW(IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2))*cj,
			 abs(RCONST(0.5)*(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge))*ui+
			     RCONST(0.5)*(IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge))*vi+
			     RCONST(0.5)*(IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge))*wi)+
		     RCONST(0.5)*sqrt(POW(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
				      POW(IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
				      POW(IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)-
					  IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2))*ci );
#else
      // Compute scalar dissipation
      double d_ij = max( abs(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*uj+
			     IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*vj+
			     IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*j)+
			 sqrt(POW(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
			      POW(IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
			      POW(IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2))*cj,
			 abs(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*ui+
			     IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*vi+
			     IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*wi)+
			 sqrt(POW(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
			      POW(IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2)+
			      POW(IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge),2))*ci );
#endif
      
      // Multiply the solution difference by the scalar dissipation
      double Diff[NVAR3D];
      Diff[0] = d_ij*(IDX3(DdataAtEdge,1,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,1,1,1,NVAR3D,2,1));
      Diff[1] = d_ij*(IDX3(DdataAtEdge,2,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,2,1,1,NVAR3D,2,1));
      Diff[2] = d_ij*(IDX3(DdataAtEdge,3,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,3,1,1,NVAR3D,2,1));
      Diff[3] = d_ij*(IDX3(DdataAtEdge,4,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,4,1,1,NVAR3D,2,1));
      Diff[4] = d_ij*(IDX3(DdataAtEdge,5,2,1,NVAR3D,2,1)-IDX3(DdataAtEdge,5,1,1,NVAR3D,2,1));
      
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
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[0] + Diff[0]);
      
      IDX3(DfluxesAtEdge,2,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[1]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[1]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[1]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[1]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[1]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[1] + Diff[1]);
      
      IDX3(DfluxesAtEdge,3,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[2]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[2]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[2]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[2]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[2]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[2] + Diff[2]);
      
      IDX3(DfluxesAtEdge,4,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[3]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[3]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[3]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[3]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[3]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[3] + Diff[3]);
      
      IDX3(DfluxesAtEdge,5,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxj[4]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyj[4]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzj[4]-
	 IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fxi[4]-
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fyi[4]-
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fzi[4] + Diff[4]);
      
      
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
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[0] + Diff[0]);
      
      IDX3(DfluxesAtEdge,2,1,1,NVAR3D,2,1) = dscale *
	(IDX3T(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[1]+
	 IDX3T(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[1]+
	 IDX3T(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[1] + Diff[1]);
      
      IDX3(DfluxesAtEdge,3,1,1,NVAR3D,2,1) = dscale *
	(IDX3(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[2]+
	 IDX3(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[2]+
	 IDX3(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[2] + Diff[2]);
      
      IDX3(DfluxesAtEdge,4,1,1,NVAR3D,2,1) = dscale *
	(IDX3(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[3]+
	 IDX3(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[3]+
	 IDX3(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[3] + Diff[3]);
      
      IDX3(DfluxesAtEdge,5,1,1,NVAR3D,2,1) = dscale *
	(IDX3(DcoeffsAtEdge,1,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[4]+
	 IDX3(DcoeffsAtEdge,2,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[4]+
	 IDX3(DcoeffsAtEdge,3,1,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[4] + Diff[4]);
      
      
      IDX3(DfluxesAtEdge,1,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[0]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[0]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[0] + Diff[0]);
      
      IDX3(DfluxesAtEdge,2,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[1]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[1]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[1] + Diff[1]);
      
      IDX3(DfluxesAtEdge,3,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[2]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[2]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[2] + Diff[2]);
      
      IDX3(DfluxesAtEdge,4,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[3]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[3]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[3] + Diff[3]);
      
      IDX3(DfluxesAtEdge,5,2,1,NVAR3D,2,1) = -dscale *
	(IDX3T(DcoeffsAtEdge,1,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fx_ij[4]+
	 IDX3T(DcoeffsAtEdge,2,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fy_ij[4]+
	 IDX3T(DcoeffsAtEdge,3,2,iedgeset+idx,HYDRO_NDIM,nmatcoeff,nedge)*Fz_ij[4] + Diff[4]);
#endif
      
      //------------------------------------------------------------------------
      // Build fluxes into nodal vector
      //------------------------------------------------------------------------
      
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

int hydro_calcFluxRusDiss3d_cuda(unsigned long * h_DcoeffsAtEdge,
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
    hydro_calcFluxRusDiss3d_knl<0><<<grid, block>>>(d_DcoeffsAtEdge,
						    d_IedgeList,
						    d_Dx, d_Dy, (*dscale), 
						    (*neq), (*nvar),
						    (*nedge), (*nmatcoeff),
						    (*nedges), (*iedgeset));
  } else {
    hydro_calcFluxRusDiss3d_knl<1><<<grid, block>>>(d_DcoeffsAtEdge,
						    d_IedgeList,
						    d_Dx, d_Dy, (*dscale), 
						    (*neq), (*nvar),
						    (*nedge), (*nmatcoeff),
						    (*nedges), (*iedgeset));
  }
  coproc_checkErrors("hydro_calcFluxRusDiss3d_cuda");
  return 0;
}

int FNAME(hydro_calcfluxrusdiss3d_cuda)(unsigned long * h_DcoeffsAtEdge,
					unsigned long * h_IedgeList,
					unsigned long * h_Dx,
					unsigned long * h_Dy,
					double * dscale,
					int * nblocks, int * neq, int * nvar,
					int * nedge,   int * nmatcoeff,
					int * nedges,  int * iedgeset)
{
  return hydro_calcFluxRusDiss3d_cuda(h_DcoeffsAtEdge, h_IedgeList,
				      h_Dx, h_Dy, dscale, nblocks, neq, nvar,
				      nedge, nmatcoeff, nedges, iedgeset);
}
