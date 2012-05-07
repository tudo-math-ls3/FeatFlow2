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

#include <cmath>
#include <cfloat>
#include <iostream>
#include "coproc_core.h"
#include "coproc_storage_cuda.h"

#define LANGUAGE LANGUAGE_C
#include "../../../../../kernel/System/idxmanager.h"

#define HYDRO_NDIM 3
#include "hydro.h"

#ifdef ENABLE_COPROC_SHMEM
#define SHMEM_IDX idx
#define SHMEM_BLOCKSIZE blockDim.x
#else
#define SHMEM_IDX 1
#define SHMEM_BLOCKSIZE 1
#endif

namespace hydro3d_cuda
{
  
  /*****************************************************************************
   * CUDA kernels for hydrodynamic model in 3D
   ****************************************************************************/

  using namespace std;

  /*****************************************************************************
   * This CUDA kernel collects the nodal solution data at the two
   * endpoints of the given edge from the global solution vector.
   ****************************************************************************/
  
  template <int isystemformat>
  struct gather_DataAtEdge
  { 
  };
  
  /*****************************************************************************
   * Input:  solution vector Dx stored in interleaved format
   * Output: DataAtEdge vector
   ****************************************************************************/

  template <>
  struct gather_DataAtEdge<SYSTEM_SCALAR>
  {
    template <typename TdSrc,
	      typename TdDest,
	      typename Ti>
    __device__ inline
    static void eval (TdDest *DataAtEdge,
		      TdSrc *Dx,
		      Ti i,
		      Ti j,
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in interleaved format

      // Gather solution data at first end point i
      IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,1,i,NVAR3D,neq);
      IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,2,i,NVAR3D,neq);
      IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,3,i,NVAR3D,neq);
      IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,4,i,NVAR3D,neq);
      IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,5,i,NVAR3D,neq);
    
      // Gather solution data at second end point j
      IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,1,j,NVAR3D,neq);
      IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,2,j,NVAR3D,neq);
      IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,3,j,NVAR3D,neq);
      IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,4,j,NVAR3D,neq);
      IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_REVERSE(Dx,5,j,NVAR3D,neq);
    }
  };
  
  /*****************************************************************************
   * Input:  solution vector Dx stored in block format
   * Output: DataAtEdge vector
   ****************************************************************************/

  template <>
  struct gather_DataAtEdge<SYSTEM_BLOCK>
  {
    template <typename TdSrc,
	      typename TdDest,
	      typename Ti>
    __device__ inline
    static void eval (TdDest *DataAtEdge,
		      TdSrc *Dx,
		      Ti i,
		      Ti j,
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in block format

      // Gather solution data at first end point i
      IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,1,i,NVAR3D,neq);
      IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,2,i,NVAR3D,neq);
      IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,3,i,NVAR3D,neq);
      IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,4,i,NVAR3D,neq);
      IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,5,i,NVAR3D,neq);
    
      // Gather solution data at second end point j
      IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,1,j,NVAR3D,neq);
      IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,2,j,NVAR3D,neq);
      IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,3,j,NVAR3D,neq);
      IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,4,j,NVAR3D,neq);
      IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = IDX2_FORWARD(Dx,5,j,NVAR3D,neq);
    }
  };

  /*****************************************************************************
   * This CUDA kernel scatters the fluxes into the global solution vector.
   ****************************************************************************/

  template <int isystemformat>
  struct scatter_FluxesAtEdge
  { 
  };

  /*****************************************************************************
   * Input:  FluxesAtEdge
   * Output: right-hand side vector Dy stored in interleaved format
   ****************************************************************************/

  template <>
  struct scatter_FluxesAtEdge<SYSTEM_SCALAR>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval (Td *FluxesAtEdge,
		      Td *Dy,
		      Ti i,
		      Ti j,
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in interleaved format

      // Scatter flux to first node i
      IDX2_REVERSE(Dy,1,i,NVAR3D,neq) += IDX3(FluxesAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,2,i,NVAR3D,neq) += IDX3(FluxesAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,3,i,NVAR3D,neq) += IDX3(FluxesAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,4,i,NVAR3D,neq) += IDX3(FluxesAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,5,i,NVAR3D,neq) += IDX3(FluxesAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    
      // Scatter flux to first node j
      IDX2_REVERSE(Dy,1,j,NVAR3D,neq) += IDX3(FluxesAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,2,j,NVAR3D,neq) += IDX3(FluxesAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,3,j,NVAR3D,neq) += IDX3(FluxesAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,4,j,NVAR3D,neq) += IDX3(FluxesAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_REVERSE(Dy,5,j,NVAR3D,neq) += IDX3(FluxesAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    }
  };

  /*****************************************************************************
   * Input:  FluxesAtEdge
   * Output: right-hand side vector Dy stored in block format
   ****************************************************************************/

  template <>
  struct scatter_FluxesAtEdge<SYSTEM_BLOCK>
  {
    template <typename Td,
	      typename Ti>
    __device__ inline
    static void eval (Td *FluxesAtEdge,
		      Td *Dy,
		      Ti i,
		      Ti j,
		      Ti neq,
		      Ti idx)
    {
      // Solution vector is stored in block format

      // Scatter flux to first node i
      IDX2_FORWARD(Dy,1,i,NVAR3D,neq) += IDX3(FluxesAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,2,i,NVAR3D,neq) += IDX3(FluxesAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,3,i,NVAR3D,neq) += IDX3(FluxesAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,4,i,NVAR3D,neq) += IDX3(FluxesAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,5,i,NVAR3D,neq) += IDX3(FluxesAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    
      // Scatter flux to first node j
      IDX2_FORWARD(Dy,1,j,NVAR3D,neq) += IDX3(FluxesAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,2,j,NVAR3D,neq) += IDX3(FluxesAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,3,j,NVAR3D,neq) += IDX3(FluxesAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,4,j,NVAR3D,neq) += IDX3(FluxesAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX2_FORWARD(Dy,5,j,NVAR3D,neq) += IDX3(FluxesAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the Galerkin fluxes at the given edge.
   ****************************************************************************/

  struct calc_GalerkinFluxAtEdge
  {
    template <typename Td,
	      typename Ti>
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
		     Td pj,
		     Ti idx)
    {
      // Compute the Galerkin fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi);
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi);
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi);
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi);
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi);
      
      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);

      // Compute Galerkin fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      
      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);

      // Compute Galerkin fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi);
      
      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
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
		     Td pj,
		     Ti idx)
    {
      // Compute Galerkin flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi)-
                      INVISCIDFLUX1_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi)-
	              INVISCIDFLUX2_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi)-
	              INVISCIDFLUX3_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi)-
	              INVISCIDFLUX4_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,ui,pi)-
	              INVISCIDFLUX5_XDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,uj,pj);
    
      // Compute Galerkin flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX1_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX2_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX3_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX4_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX5_YDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      
      // Compute Galerkin flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX1_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX2_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX3_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
	              INVISCIDFLUX4_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vi,pi)-
  	              INVISCIDFLUX5_ZDIR3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE,vj,pj);
    }
#endif
  };

  /*****************************************************************************
   * This CUDA kernel calculates the artificial dissipation at the given edge
   ****************************************************************************/

  template <int idissipationtype>
  struct calc_DissipationAtEdge
  {
  };

  /*****************************************************************************
   * Zero artificial dissipation, aka standard Galerkin approach
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<DISSIPATION_ZERO>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *Diff,
		     Tc *CoeffsAtEdge,
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
		     Ti ncoeff,
		     Ti idx)
    {
      Diff[0] = 0.0;
      Diff[1] = 0.0;
      Diff[2] = 0.0;
      Diff[3] = 0.0;
      Diff[4] = 0.0;
    }
  };

  /*****************************************************************************
   * Scalar artificial dissipation proportional to the spectral radius
   * (largest eigenvector) of the cumulative Roe matrix.
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<DISSIPATION_SCALAR>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *Diff,
		     Tc *CoeffsAtEdge,
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
		     Ti ncoeff,
		     Ti idx)
    {
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
      Td ri = DENSITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      Td rj = DENSITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    
      // Compute enthalpies
      Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pi)/ri;
      Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pj)/rj;
    
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
      Td c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON));
    
      // Compute scalar dissipation
      Td d_ij = abs(vel_ij) + anorm*c_ij;
    
      // Multiply the solution difference by the scalar dissipation
      Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
    }
  };

  /*****************************************************************************
   * Scalar artificial dissipation proportional to the spectral radius
   * (largest eigenvector) of the dimensional-split Roe matrix.
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<DISSIPATION_SCALAR_DSPLIT>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *Diff,
		     Tc *CoeffsAtEdge,
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
		     Ti ncoeff,
		     Ti idx)
    {
      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[2] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge));
    
      // Compute densities
      Td ri = DENSITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      Td rj = DENSITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    
      // Compute enthalpies
      Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pi)/ri;
      Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pj)/rj;
    
      // Compute Roe mean values
      Td aux  = ROE_MEAN_RATIO(ri,rj);
      Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
      Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
      Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
      Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
    
      // Compute auxiliary variables
      Td q_ij = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
    
      // Compute the speed of sound
      Td c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON));
    
      // Compute scalar dissipation
      Td d_ij = ( abs(a[0]*u_ij) + abs(a[0])*c_ij +
		  abs(a[1]*v_ij) + abs(a[1])*c_ij +
		  abs(a[2]*w_ij) + abs(a[2])*c_ij );
    
      // Multiply the solution difference by the scalar dissipation
      Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
    }
  };

  /*****************************************************************************
   * Tensorial artificial dissipation of Roe-type.
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<DISSIPATION_ROE>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *Diff, 
		     Tc *CoeffsAtEdge,
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
		     Ti ncoeff,
		     Ti idx)
    {

      // Compute skew-symmetric coefficient
      Td a[HYDRO_NDIM];
      a[0] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[1] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      a[2] = RCONST(0.5)*(IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)-
			  IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge));
      Td anorm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    
      if (anorm > DBL_EPSILON) {
      
	// Normalise the skew-symmetric coefficient
	a[0] = a[0]/anorm;
	a[1] = a[1]/anorm;
	a[2] = a[2]/anorm;
      
	// Compute densities
	Td ri = DENSITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	Td rj = DENSITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
	// Compute enthalpies
	Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pi)/ri;
	Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pj)/rj;
      
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
	Td c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), DBL_EPSILON);
	Td c_ij  = sqrt(c2_ij);
      
	// Compute eigenvalues
	Td l1 = abs(vel_ij-c_ij);
	Td l2 = abs(vel_ij);
	Td l3 = abs(vel_ij+c_ij);
	Td l4 = abs(vel_ij);
	Td l5 = abs(vel_ij);
      
	// Compute solution difference U_j-U_i
	Diff[0] = IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	         -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	Diff[1] = IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		 -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	Diff[2] = IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		 -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	Diff[3] = IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		 -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	Diff[4] = IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		 -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
      
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

  /*****************************************************************************
   * Tensorial artificial dissipation of Roe-type using dimensional splitting.
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<DISSIPATION_ROE_DSPLIT>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *Diff,
		     Tc *CoeffsAtEdge,
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
		     Ti ncoeff,
		     Ti idx)
    {
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
      if (anorm > DBL_EPSILON) {
      
	// Compute the absolute value
	a[0] = abs(a[0]);
	a[1] = abs(a[1]);
	a[2] = abs(a[2]);
      
	// Compute densities
	Td ri = DENSITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	Td rj = DENSITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
	// Compute enthalpies
	Td hi = (TOTALENERGY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pi)/ri;
	Td hj = (TOTALENERGY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)+pj)/rj;
      
	// Compute Roe mean values
	Td aux  = ROE_MEAN_RATIO(ri,rj);
	Td u_ij = ROE_MEAN_VALUE(ui,uj,aux);
	Td v_ij = ROE_MEAN_VALUE(vi,vj,aux);
	Td w_ij = ROE_MEAN_VALUE(wi,wj,aux);
	Td H_ij = ROE_MEAN_VALUE(hi,hj,aux);
      
	// Compute auxiliary variable
	Td q_ij = RCONST(0.5) * (u_ij * u_ij + v_ij * v_ij + w_ij * w_ij);
      
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
	Td l5 = abs(u_ij);
      
	// Compute solution difference U_j-U_i
	DiffAux[0] = IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[1] = IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[2] = IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[3] = IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[4] = IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
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
	DiffAux[0] = IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[1] = IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[2] = IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[3] = IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[4] = IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
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
	DiffAux[0] = IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[1] = IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[2] = IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[3] = IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	DiffAux[4] = IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
	            -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
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

  /*****************************************************************************
   * Scalar artificial dissipation of Rusanov-type.
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<DISSIPATION_RUSANOV>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *Diff,
		     Tc *CoeffsAtEdge,
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
		     Ti ncoeff,
		     Ti idx)
    {
      // Compute specific energies
      Td Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      Td Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    
      // Compute the speed of sound
      Td ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ei-RCONST(0.5)*(ui*ui+vi*vi)), DBL_EPSILON));
      Td cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ej-RCONST(0.5)*(uj*uj+vj*vj)), DBL_EPSILON));
    
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
      Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
    }
  };

  /*****************************************************************************
   * Scalar artificial dissipation of Rusanov-type using dimensional splitting.
   ****************************************************************************/

  template <>
  struct calc_DissipationAtEdge<DISSIPATION_RUSANOV_DSPLIT>
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
    __device__ inline
    static void eval(Td *Diff, 
		     Tc *CoeffsAtEdge,
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
		     Ti ncoeff,
		     Ti idx)
    {
      // Compute specific energies
      Td Ei = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      Td Ej = SPECIFICTOTALENERGY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    
      // Compute the speed of sound
      Td ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ei-RCONST(0.5)*(ui*ui+vi*vi)), DBL_EPSILON));
      Td cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*
		       (HYDRO_GAMMA)*(Ej-RCONST(0.5)*(uj*uj+vj*vj)), DBL_EPSILON));
    
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
      Diff[0] = d_ij*(IDX3(DataAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[1] = d_ij*(IDX3(DataAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[2] = d_ij*(IDX3(DataAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[3] = d_ij*(IDX3(DataAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
      Diff[4] = d_ij*(IDX3(DataAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE)
		     -IDX3(DataAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE));
    }
  };

  /*****************************************************************************
   * This CUDA kernel calculates the fluxes at the given edge.
   ****************************************************************************/

  struct calc_FluxesAtEdge
  {
    template <typename Tc,
	      typename Td,
	      typename Ti>
#ifdef HYDRO_USE_IBP
    __device__ inline
    static void eval(Td *FluxesAtEdge,
		     Tc *CoeffsAtEdge,
		     Td *Fxi,
		     Td *Fxj, 
		     Td *Fyi, 
		     Td *Fyj,
		     Td *Fzi,
		     Td *Fzj, 
		     Td *Diff,
		     Td scale,
		     Ti iedge, 
		     Ti nedge,
		     Ti ncoeff,
		     Ti idx)
    {
      IDX3(FluxesAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[0]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[0]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[0]-
	 IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[0]-
	 IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[0]-
	 IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[0] + Diff[0]);
    
      IDX3(FluxesAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[1]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[1]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[1]-
	 IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[1]-
	 IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[1]-
	 IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[1] + Diff[1]);
    
      IDX3(FluxesAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[2]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[2]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[2]-
	 IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[2]-
	 IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[2]-
	 IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[2] + Diff[2]);
    
      IDX3(FluxesAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[3]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[3]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[3]-
	 IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[3]-
	 IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[3]-
	 IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[3] + Diff[3]);
    
      IDX3(FluxesAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxj[4]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyj[4]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzj[4]-
	 IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fxi[4]-
	 IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fyi[4]-
	 IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fzi[4] + Diff[4]);
    
      IDX3(FluxesAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) =
     -IDX3(FluxesAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX3(FluxesAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) =
     -IDX3(FluxesAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX3(FluxesAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) =
     -IDX3(FluxesAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX3(FluxesAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) =
     -IDX3(FluxesAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      IDX3(FluxesAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) =
     -IDX3(FluxesAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
    }
#else
    __device__ inline
    static void eval(Td *FluxesAtEdge,
		     Tc *CoeffsAtEdge,
		     Td *Fx_ij,
		     Td *Fy_ij,
		     Td *Fz_ij,
		     Td *Diff,
		     Td scale,
		     Ti iedge,
		     Ti nedge,
		     Ti ncoeff,
		     Ti idx)
    {
      IDX3(FluxesAtEdge,1,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[0]+
	 IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[0]+
	 IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[0] + Diff[0]);
    
      IDX3(FluxesAtEdge,2,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3T(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[1]+
	 IDX3T(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[1]+
	 IDX3T(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[1] + Diff[1]);
    
      IDX3(FluxesAtEdge,3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[2]+
	 IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[2]+
	 IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[2] + Diff[2]);
    
      IDX3(FluxesAtEdge,4,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[3]+
	 IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[3]+
	 IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[3] + Diff[3]);
    
      IDX3(FluxesAtEdge,5,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = scale *
	(IDX3(CoeffsAtEdge,1,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[4]+
	 IDX3(CoeffsAtEdge,2,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[4]+
	 IDX3(CoeffsAtEdge,3,1,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[4] + Diff[4]);
    
    
      IDX3(FluxesAtEdge,1,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = -scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[0]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[0]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[0] + Diff[0]);
    
      IDX3(FluxesAtEdge,2,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = -scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[1]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[1]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[1] + Diff[1]);
    
      IDX3(FluxesAtEdge,3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = -scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[2]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[2]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[2] + Diff[2]);
    
      IDX3(FluxesAtEdge,4,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = -scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[3]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[3]+
	 IDX3T(CoeffsAtEdge,3,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[3] + Diff[3]);
    
      IDX3(FluxesAtEdge,5,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE) = -scale *
	(IDX3T(CoeffsAtEdge,1,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fx_ij[4]+
	 IDX3T(CoeffsAtEdge,2,2,iedge,HYDRO_NDIM,ncoeff,nedge)*Fy_ij[4]+
	 IDX3T(CoeffsAtEdge,3,2,iedge+idx,HYDRO_NDIM,ncoeff,nedge)*Fz_ij[4] + Diff[4]);
    }
#endif
  };

  /*****************************************************************************
   * This CUDA kernel calculates the inviscid fluxes.
   ****************************************************************************/
  
  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti,
	    int isystemformat,
	    int idissipationtype>
  __global__ void hydro_calcFlux3d_knl(Tc *CoeffsAtEdge,
				       Ti *IedgeList,
				       TdSrc *Dx,
				       TdDest *Dy,
				       TdDest scale,
				       Ti neq,
				       Ti nedge,
				       Ti ncoeff,
				       Ti nedges,
				       Ti iedgeset)
  {
#ifdef ENABLE_COPROC_SHMEM
    // Use shared memory
    extern __shared__ TdDest shmemData[];
#endif

    Ti idx = blockIdx.x * blockDim.x + threadIdx.x;
  
    if (idx<nedges)
      {
	// Get positions of edge endpoints (idx starts at zero)
	Ti i = IDX2(IedgeList,1,iedgeset+idx,6,nedge);
	Ti j = IDX2(IedgeList,2,iedgeset+idx,6,nedge);
      
#ifdef ENABLE_COPROC_SHMEM
	// Local data at edge from shared memory
	TdDest *DataAtEdge = shmemData;
#else
	// Local data at edge from local memory
	TdDest DataAtEdge[2*NVAR3D];
#endif
      
	// Get solution values at edge endpoints
	gather_DataAtEdge<isystemformat>::
	  eval(DataAtEdge,Dx,i,j,neq,idx);
      
	// Compute velocities
	TdDest ui = XVELOCITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	TdDest vi = YVELOCITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	TdDest wi = ZVELOCITY3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
	TdDest uj = XVELOCITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	TdDest vj = YVELOCITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	TdDest wj = ZVELOCITY3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
      
	// Compute pressures
	TdDest pi = PRESSURE3(DataAtEdge,IDX3,1,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);
	TdDest pj = PRESSURE3(DataAtEdge,IDX3,2,SHMEM_IDX,NVAR3D,2,SHMEM_BLOCKSIZE);

#ifdef HYDRO_USE_IBP
	TdDest Fxi[NVAR3D];
	TdDest Fxj[NVAR3D];
	TdDest Fyi[NVAR3D];
	TdDest Fyj[NVAR3D];
	TdDest Fzi[NVAR3D];
	TdDest Fzj[NVAR3D];
      
	// Compute the Galerkin fluxes
	calc_GalerkinFluxAtEdge::
	  eval(Fxi,Fxj,Fyi,Fyj,Fzi,Fzj,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,idx);
#else
	TdDest Fx_ij[NVAR3D];
	TdDest Fy_ij[NVAR3D];
	TdDest Fz_ij[NVAR3D];

	// Compute the Galerkin fluxes
	calc_GalerkinFluxAtEdge::
	  eval(Fx_ij,Fy_ij,Fz_ij,DataAtEdge,ui,uj,vi,vj,wi,wj,pi,pj,idx);
#endif
      
	TdDest Diff[NVAR3D];
	// Compute the artificial viscosities
	calc_DissipationAtEdge<idissipationtype>::
	  eval(Diff,CoeffsAtEdge,DataAtEdge,
	       ui,uj,vi,vj,wi,wj,pi,pj,iedgeset+idx,nedge,ncoeff,idx);
      
	// Build both contributions into the fluxes
#ifdef HYDRO_USE_IBP
	calc_FluxesAtEdge::
	  eval(DataAtEdge,CoeffsAtEdge,Fxi,Fxj,Fyi,Fyj,Fzi,Fzj,Diff,
	       scale,iedgeset+idx,nedge,ncoeff,idx);
#else
	calc_FluxesAtEdge::
	  eval(DataAtEdge,CoeffsAtEdge,Fx_ij,Fy_ij,Fz_ij,Diff,
	       scale,iedgeset+idx,nedge,ncoeff,idx);
#endif

	// Build fluxes into nodal vector
	scatter_FluxesAtEdge<isystemformat>::
	  eval(DataAtEdge,Dy,i,j,neq,idx);
      }
  };

  /*****************************************************************************
   * Internal C++ functions which invoke the CUDA kernels
   ****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcFluxGalerkin3d_cuda(__SIZET *d_CoeffsAtEdge,
				    __SIZET *d_IedgeList,
				    __SIZET *d_Dx,
				    __SIZET *d_Dy,
				    TdDest scale,
				    Ti nblocks,
				    Ti neq,
				    Ti nedge, 
				    Ti ncoeff,
				    Ti nedges,
				    Ti iedgeset,
				    cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    TdDest *Dy = (TdDest*)(*d_Dy);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_ZERO>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale,
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_ZERO>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcFluxGalerkin3d_cuda");
    return 0;
  }; 

  /****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcFluxScDiss3d_cuda(__SIZET *d_CoeffsAtEdge,
				  __SIZET *d_IedgeList,
				  __SIZET *d_Dx,
				  __SIZET *d_Dy,
				  TdDest scale,
				  Ti nblocks,
				  Ti neq,
				  Ti nedge, 
				  Ti ncoeff,
				  Ti nedges,
				  Ti iedgeset,
				  cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    TdDest *Dy = (TdDest*)(*d_Dy);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_SCALAR>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale,
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_SCALAR>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcFluxScDiss3d_cuda");
    return 0;
  };

  /****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcFluxScDissDiSp3d_cuda(__SIZET *d_CoeffsAtEdge,
				      __SIZET *d_IedgeList,
				      __SIZET *d_Dx,
				      __SIZET *d_Dy,
				      TdDest scale,
				      Ti nblocks,
				      Ti neq, 
				      Ti nedge,
				      Ti ncoeff,
				      Ti nedges,
				      Ti iedgeset,
				      cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    TdDest *Dy = (TdDest*)(*d_Dy);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_SCALAR_DSPLIT>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_SCALAR_DSPLIT>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcFluxScDissDiSp3d_cuda");
    return 0;
  };

  /****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcFluxRoeDiss3d_cuda(__SIZET *d_CoeffsAtEdge,
				   __SIZET *d_IedgeList,
				   __SIZET *d_Dx,
				   __SIZET *d_Dy,
				   TdDest scale,
				   Ti nblocks, 
				   Ti neq, 
				   Ti nedge,
				   Ti ncoeff,
				   Ti nedges, 
				   Ti iedgeset,
				   cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    TdDest *Dy = (TdDest*)(*d_Dy);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_ROE>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_ROE>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcFluxRoeDiss3d_cuda");
    return 0;
  };

  /****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcFluxRoeDissDiSp3d_cuda(__SIZET *d_CoeffsAtEdge,
				       __SIZET *d_IedgeList,
				       __SIZET *d_Dx,
				       __SIZET *d_Dy,
				       TdDest scale,
				       Ti nblocks, 
				       Ti neq, 
				       Ti nedge,
				       Ti ncoeff,
				       Ti nedges, 
				       Ti iedgeset,
				       cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    TdDest *Dy = (TdDest*)(*d_Dy);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_ROE_DSPLIT>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_ROE_DSPLIT>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcFluxRoeDissDiSp3d_cuda");
    return 0;
  };

  /****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcFluxRusDiss3d_cuda(__SIZET *d_CoeffsAtEdge,
				   __SIZET *d_IedgeList,
				   __SIZET *d_Dx,
				   __SIZET *d_Dy,
				   TdDest scale,
				   Ti nblocks, 
				   Ti neq, 
				   Ti nedge, 
				   Ti ncoeff,
				   Ti nedges, 
				   Ti iedgeset,
				   cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    TdDest *Dy = (TdDest*)(*d_Dy);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_RUSANOV>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_RUSANOV>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcFluxRusDiss3d_cuda");
    return 0;
  }

  /****************************************************************************/

  template <typename Tc,
	    typename TdSrc,
	    typename TdDest,
	    typename Ti>
  inline
  int hydro_calcFluxRusDissDiSp3d_cuda(__SIZET *d_CoeffsAtEdge,
				       __SIZET *d_IedgeList,
				       __SIZET *d_Dx,
				       __SIZET *d_Dy,
				       TdDest scale,
				       Ti nblocks, 
				       Ti neq,
				       Ti nedge, 
				       Ti ncoeff,
				       Ti nedges, 
				       Ti iedgeset,
				       cudaStream_t stream=0)
  {
    TdSrc  *Dx = (TdSrc*)(*d_Dx);
    TdDest *Dy = (TdDest*)(*d_Dy);
    Tc *CoeffsAtEdge = (Tc*)(*d_CoeffsAtEdge);
    Ti *IedgeList = (Ti*)(*d_IedgeList);
  
    // Define number of threads per block
    int blocksize = 128;
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((nedges)/(double)(block.x));
  
    if (nblocks == 1) {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_SCALAR,DISSIPATION_RUSANOV_DSPLIT>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    } else {
      hydro_calcFlux3d_knl
	<Tc,TdSrc,TdDest,Ti,SYSTEM_BLOCK,DISSIPATION_RUSANOV_DSPLIT>
	<<<grid, block, 0, stream>>>(CoeffsAtEdge,
				     IedgeList,
				     Dx, Dy, scale, 
				     neq, nedge, ncoeff,
				     nedges, iedgeset);
    }
    coproc_checkErrors("hydro_calcFluxRusDissDiSp3d_cuda");
    return 0;
  };

  /*****************************************************************************
   * External C functions which can be called from the Fortran code
   ****************************************************************************/

  extern "C"
  {
    __INT FNAME(hydro_calcfluxgalerkin3d_cuda)(__SIZET *d_CoeffsAtEdge,
					       __SIZET *d_IedgeList,
					       __SIZET *d_Dx,
					       __SIZET *d_Dy,
					       __DP *scale,
					       __INT *nblocks,
					       __INT *neq,
					       __INT *nedge,
					       __INT *ncoeff,
					       __INT *nedges,
					       __INT *iedgeset,
					       __I64 *stream)
    {
      return (__INT) hydro_calcFluxGalerkin3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
			       *scale, *nblocks, *neq, *nedge,
			       *ncoeff, *nedges, *iedgeset,
			       (cudaStream_t)(*stream));
    }
    
    /**************************************************************************/
    
    __INT FNAME(hydro_calcfluxscdiss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					     __SIZET *d_IedgeList,
					     __SIZET *d_Dx,
					     __SIZET *d_Dy,
					     __DP *scale,
					     __INT *nblocks,
					     __INT *neq,
					     __INT *nedge,
					     __INT *ncoeff,
					     __INT *nedges,
					     __INT *iedgeset,
					     __I64 *stream)
    {
      return (__INT) hydro_calcFluxScDiss3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
			       *scale, *nblocks, *neq, *nedge,
			       *ncoeff, *nedges, *iedgeset,
			       (cudaStream_t)*stream);
    }

    /**************************************************************************/

    __INT FNAME(hydro_calcfluxscdissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
						 __SIZET *d_IedgeList,
						 __SIZET *d_Dx,
						 __SIZET *d_Dy,
						 __DP *scale,
						 __INT *nblocks, 
						 __INT *neq, 
						 __INT *nedge, 
						 __INT *ncoeff,
						 __INT *nedges, 
						 __INT *iedgeset,
						 __I64 *stream=0)
    {
      return (__INT) hydro_calcFluxScDissDiSp3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
			       *scale, *nblocks, *neq, *nedge,
			       *ncoeff, *nedges, *iedgeset,
			       (cudaStream_t)(*stream));
    }

  /***************************************************************************/

    __INT FNAME(hydro_calcfluxroediss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					      __SIZET *d_IedgeList,
					      __SIZET *d_Dx,
					      __SIZET *d_Dy,
					      __DP *scale,
					      __INT *nblocks, 
					      __INT *neq, 
					      __INT *nedge, 
					      __INT *ncoeff,
					      __INT *nedges, 
					      __INT *iedgeset,
					      __I64 *stream)
    {
      return (__INT) hydro_calcFluxRoeDiss3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
			       *scale, *nblocks, *neq, *nedge,
			       *ncoeff, *nedges, *iedgeset,
			       (cudaStream_t)(*stream));
    }

  /***************************************************************************/

    __INT FNAME(hydro_calcfluxroedissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
						  __SIZET *d_IedgeList,
						  __SIZET *d_Dx,
						  __SIZET *d_Dy,
						  __DP *scale,
						  __INT *nblocks, 
						  __INT *neq, 
						  __INT *nedge, 
						  __INT *ncoeff,
						  __INT *nedges, 
						  __INT *iedgeset,
						  __I64 *stream)
    {
      return (__INT) hydro_calcFluxRoeDissDiSp3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
			       *scale, *nblocks, *neq, *nedge,
			       *ncoeff, *nedges, *iedgeset,
			       (cudaStream_t)*stream);
    }

  /***************************************************************************/

    __INT FNAME(hydro_calcfluxrusdiss3d_cuda)(__SIZET *d_CoeffsAtEdge,
					      __SIZET *d_IedgeList,
					      __SIZET *d_Dx,
					      __SIZET *d_Dy,
					      __DP *scale,
					      __INT *nblocks, 
					      __INT *neq, 
					      __INT *nedge, 
					      __INT *ncoeff,
					      __INT *nedges, 
					      __INT *iedgeset,
					      __I64 *stream)
    {
      return (__INT)hydro_calcFluxRusDiss3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
			       *scale, *nblocks, *neq, *nedge,
			       *ncoeff, *nedges, *iedgeset,
			       (cudaStream_t)*stream);
    }

  /***************************************************************************/

    __INT FNAME(hydro_calcfluxrusdissdisp3d_cuda)(__SIZET *d_CoeffsAtEdge,
						  __SIZET *d_IedgeList,
						  __SIZET *d_Dx,
						  __SIZET *d_Dy,
						  __DP *scale,
						  __INT *nblocks, 
						  __INT *neq, 
						  __INT *nedge, 
						  __INT *ncoeff,
						  __INT *nedges, 
						  __INT *iedgeset,
						  __I64 *stream)
    {
      return (__INT) hydro_calcFluxRusDissDiSp3d_cuda
	<__DP,__DP,__DP,__INT>(d_CoeffsAtEdge, d_IedgeList, d_Dx, d_Dy,
			       *scale, *nblocks, *neq, *nedge,
			       *ncoeff, *nedges, *iedgeset,
			       (cudaStream_t)*stream);
    }
  }
}
