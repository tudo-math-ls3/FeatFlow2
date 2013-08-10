!##############################################################################
!# ****************************************************************************
!# <name> zpinch_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the simplified MHD equations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) zpinch_calcMatDiagConvIntlP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#        for hydrodynamic systems stored in interleaved format
!#
!# 2.) zpinch_calcMatRusConvIntlP2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#        and applies scalar artificial viscosities of Rusanov-type
!#        for hydrodynamic systems stored in interleaved format
!#
!# 3.) zpinch_calcMatDiagConvIntlD2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#        for hydrodynamic systems stored in interleaved format
!#
!# 4.) zpinch_calcMatRusConvIntlD2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#        and applies scalar artificial viscosities of Rusanov-type
!#        for hydrodynamic systems stored in interleaved format
!#
!# 5.) zpinch_calcMatDiagConvBlockP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#        for hydrodynamic systems stored in block format
!#
!# 6.) zpinch_calcMatRusConvBlockP2d_sim
!#      -> Calculates the off-diagonal Galerkin transport coefficients
!#         for linear convection in 2D (primal formulation)
!#         and applies scalar artificial viscosities of Rusanov-type
!#         for hydrodynamic systems stored in block format
!#
!# 7.) zpinch_calcMatDiagConvBlockD2d_sim
!#      -> Calculates the diagonal Galerkin transport coefficients
!#         for linear convection in 2D (dual formulation)
!#         for hydrodynamic systems stored in block format
!#
!# 8.) zpinch_calcMatRusConvBlockD2d_sim
!#      -> Calculates the off-diagonal Galerkin transport coefficients
!#         for linear convection in 2D (dual formulation)
!#         and applies scalar artificial viscosities of Rusanov-type
!#         for hydrodynamic systems stored in block format
!#
!# </purpose>
!##############################################################################

module zpinch_callback2d

#include "flagship.h"
#include "models/hydro/hydro.h"

!$use omp_lib
  use collection
  use hydro_basic
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use storage

  implicit none

  private

  public :: zpinch_calcMatDiagConvIntlP2d_sim
  public :: zpinch_calcMatRusConvIntlP2d_sim
  public :: zpinch_calcMatDiagConvIntlD2d_sim
  public :: zpinch_calcMatRusConvIntlD2d_sim
  public :: zpinch_calcMatDiagConvBlockP2d_sim
  public :: zpinch_calcMatRusConvBlockP2d_sim
  public :: zpinch_calcMatDiagConvBlockD2d_sim
  public :: zpinch_calcMatRusConvBlockD2d_sim
  
contains

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvIntlP2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)
    
!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for the velocity vector $v=v(x,y,t)$
    ! for the primal problem in 2D.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in interleaved format.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    
    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))

#else
      ! Compute convective coefficient  $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvIntlP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvIntlP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the velocity vector $v=v(x,y,t)$ for
    ! the primal problem in 2D.  Moreover, scalar artificial viscosity
    ! of Rusanov-type is applied.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in interleaved format.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DsolutionHydro
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP) :: Ei,Ej,ui,uj,vi,vj,ci,cj
    integer :: iedge,idx,jdx

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    ! This subroutine assumes that a second collection structure is
    ! attached to rcollection, and moreover, the first quick access
    ! vector of this collection structure points to the solution of
    ! the hydrodynamic mode
    call lsysbl_getbase_double(rcollection%p_rnextCollection%&
        p_rvectorQuickAccess1, p_DsolutionHydro)
    
    do iedge = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------
      
#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ij} = v_j*C_{ji}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient  $k_{ji} = v_i*Cx_{ij}$
      DmatrixAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient  $k_{ij} = -v_j*C_{ij}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient  $k_{ji} = -v_i*Cx_{ji}$
      DmatrixAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Compute base indices
      idx = 4*(IedgeList(1,iedge)-1)
      jdx = 4*(IedgeList(2,iedge)-1)
      
      ! Compute auxiliary variables
      ui = p_DsolutionHydro(idx+2)/p_DsolutionHydro(idx+1)
      vi = p_DsolutionHydro(idx+3)/p_DsolutionHydro(idx+1)
      Ei = p_DsolutionHydro(idx+4)/p_DsolutionHydro(idx+1)
      uj = p_DsolutionHydro(jdx+2)/p_DsolutionHydro(jdx+1)
      vj = p_DsolutionHydro(jdx+3)/p_DsolutionHydro(jdx+1)
      Ej = p_DsolutionHydro(jdx+4)/p_DsolutionHydro(jdx+1)
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))**2+&
                          (DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))**2+&
                          (DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))**2)*ci )
#else
      ! Compute dissipation tensor D_ij
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(DcoeffsAtEdge(1,1,iedge)*uj +&
                   DcoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DcoeffsAtEdge(1,1,iedge)**2 +&
                        DcoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DcoeffsAtEdge(1,2,iedge)*ui +&
                   DcoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DcoeffsAtEdge(1,2,iedge)**2 +&
                        DcoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvIntlP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvIntlD2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)
    
!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for the velocity vector $v=v(x,y,t)$
    ! for the dual problem in 2D.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in interleaved format.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    
    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient  $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvIntlD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvIntlD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the velocity vector $v=v(x,y,t)$ for
    ! the dual problem in 2D.  Moreover, scalar artificial viscosity
    ! of Rusanov-type is applied.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in interleaved format.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DsolutionHydro
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP) :: Ei,Ej,ui,uj,vi,vj,ci,cj
    integer :: iedge,idx,jdx

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    
    ! This subroutine assumes that a second collection structure is
    ! attached to rcollection, and moreover, the first quick access
    ! vector of this collection structure points to the solution of
    ! the hydrodynamic mode
    call lsysbl_getbase_double(rcollection%p_rnextCollection%&
        p_rvectorQuickAccess1, p_DsolutionHydro)

    do iedge = 1, nedges

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji] = -v_i*C_{ij}$
      DmatrixAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji] = v_i*C_{ji}$
      DmatrixAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

      ! Compute base indices
      idx = 4*(IedgeList(1,iedge)-1)
      jdx = 4*(IedgeList(2,iedge)-1)
      
      ! Compute auxiliary variables
      ui = p_DsolutionHydro(idx+2)/p_DsolutionHydro(idx+1)
      vi = p_DsolutionHydro(idx+3)/p_DsolutionHydro(idx+1)
      Ei = p_DsolutionHydro(idx+4)/p_DsolutionHydro(idx+1)
      uj = p_DsolutionHydro(jdx+2)/p_DsolutionHydro(jdx+1)
      vj = p_DsolutionHydro(jdx+3)/p_DsolutionHydro(jdx+1)
      Ej = p_DsolutionHydro(jdx+4)/p_DsolutionHydro(jdx+1)
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))**2+&
                          (DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))**2+&
                          (DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))**2)*ci )
#else
      ! Compute dissipation tensor D_ij
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(DcoeffsAtEdge(1,1,iedge)*uj +&
                   DcoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DcoeffsAtEdge(1,1,iedge)**2 +&
                        DcoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DcoeffsAtEdge(1,2,iedge)*ui +&
                   DcoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DcoeffsAtEdge(1,2,iedge)**2 +&
                        DcoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvIntlD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvBlockP2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)
    
!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for the velocity vector $v=v(x,y,t)$
    ! for the primal problem in 2D.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in block format.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    
    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvBlockP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvBlockP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the velocity vector $v=v(x,y,t)$ for
    ! the dual problem in 2D.  Moreover, scalar artificial viscosity
    ! of Rusanov-type is applied.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in interleaved format.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DsolutionHydro
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP) :: Ei,Ej,ui,uj,vi,vj,ci,cj
    integer :: iedge,neq,i,j

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    ! This subroutine assumes that a second collection structure is
    ! attached to rcollection, and moreover, the first quick access
    ! vector of this collection structure points to the solution of
    ! the hydrodynamic mode
    call lsysbl_getbase_double(rcollection%p_rnextCollection%&
        p_rvectorQuickAccess1, p_DsolutionHydro)

    do iedge = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DmatrixAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DmatrixAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Get number of equations for single variable
      neq = size(p_DvelocityX)

      ! Get equations i and j
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)
      
      ! Compute auxiliary variables
      ui = p_DsolutionHydro(neq  +i)/p_DsolutionHydro(i)
      vi = p_DsolutionHydro(neq*2+i)/p_DsolutionHydro(i)
      Ei = p_DsolutionHydro(neq*3+i)/p_DsolutionHydro(i)
      uj = p_DsolutionHydro(neq  +j)/p_DsolutionHydro(j)
      vj = p_DsolutionHydro(neq*2+j)/p_DsolutionHydro(j)
      Ej = p_DsolutionHydro(neq*3+j)/p_DsolutionHydro(j)
            
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))

#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))**2+&
                          (DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))**2+&
                          (DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))**2)*ci )
#else
      ! Compute dissipation tensor D_ij
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(DcoeffsAtEdge(1,1,iedge)*uj +&
                   DcoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DcoeffsAtEdge(1,1,iedge)**2 +&
                        DcoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DcoeffsAtEdge(1,2,iedge)*ui +&
                   DcoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DcoeffsAtEdge(1,2,iedge)**2 +&
                        DcoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvBlockP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvBlockD2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)
    
!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for the velocity vector $v=v(x,y,t)$
    ! for the dual problem in 2D.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in block format.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    
    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvBlockD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvBlockD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the velocity vector $v=v(x,y,t)$ for
    ! the dual problem in 2D.  Moreover, scalar artificial viscosity
    ! of Rusanov-type is applied.
    !
    ! This subroutine assumes that the conservative variables of the
    ! Hydrodynamic system are stored in block format.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DsolutionHydro
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP) :: Ei,Ej,ui,uj,vi,vj,ci,cj
    integer :: iedge,neq,i,j

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    ! This subroutine assumes that a second collection structure is
    ! attached to rcollection, and moreover, the first quick access
    ! vector of this collection structure points to the solution of
    ! the hydrodynamic mode
    call lsysbl_getbase_double(rcollection%p_rnextCollection%&
        p_rvectorQuickAccess1, p_DsolutionHydro)

    do iedge = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji] = -v_i*C_{ij}$
      DmatrixAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji] = v_i*C_{ji}$
      DmatrixAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Get number of equations for single variable
      neq = size(p_DvelocityX)
      
      ! Get equations i and j
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)

      ! Compute auxiliary variables
      ui = p_DsolutionHydro(neq  +i)/p_DsolutionHydro(i)
      vi = p_DsolutionHydro(neq*2+i)/p_DsolutionHydro(i)
      Ei = p_DsolutionHydro(neq*3+i)/p_DsolutionHydro(i)
      uj = p_DsolutionHydro(neq  +j)/p_DsolutionHydro(j)
      vj = p_DsolutionHydro(neq*2+j)/p_DsolutionHydro(j)
      Ej = p_DsolutionHydro(neq*3+j)/p_DsolutionHydro(j)
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,1,iedge)-&
                           DcoeffsAtEdge(1,2,iedge))**2+&
                          (DcoeffsAtEdge(2,1,iedge)-&
                           DcoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DcoeffsAtEdge(1,2,iedge)-&
                           DcoeffsAtEdge(1,1,iedge))**2+&
                          (DcoeffsAtEdge(2,2,iedge)-&
                           DcoeffsAtEdge(2,1,iedge))**2)*ci )
#else
      ! Compute dissipation tensor D_ij
      DmatrixAtEdge(1,iedge) = dscale*&
          max( abs(DcoeffsAtEdge(1,1,iedge)*uj +&
                   DcoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DcoeffsAtEdge(1,1,iedge)**2 +&
                        DcoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DcoeffsAtEdge(1,2,iedge)*ui +&
                   DcoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DcoeffsAtEdge(1,2,iedge)**2 +&
                        DcoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvBlockD2d_sim

end module zpinch_callback2d
