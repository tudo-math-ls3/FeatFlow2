!##############################################################################
!# ****************************************************************************
!# <name> transport_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in 3D.
!#
!# ****************************************************************************
!#
!# The following routines for linear velocity case are available:
!#
!# 1.) transp_calcMatDiagConvP3d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 3D (primal formulation)
!#
!# 2.) transp_calcMatGalConvP3d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 3D (primal formulation)
!#
!# 3.) transp_calcMatUpwConvP3d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 3D (primal formulation)
!#
!# 4.) transp_calcMatDiagConvD3d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 3D (dual formulation)
!#
!# 5.) transp_calcMatGalConvD3d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 3D (dual formulation)
!#
!# 6.) transp_calcMatUpwConvD3d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 3D (dual formulation)
!#
!# 7.) transp_calcVecBdrConvP3d_sim
!#     -> Calculates the group finite element coefficients for the
!#        linear form in 3D (primal formulation)
!#
!# 8.) transp_calcMatBdrConvP3d_sim
!#     -> Calculates the group finite element coefficients for the
!#        bilinear form in 3D (primal formulation)
!#
!# 9.) transp_calcVecBdrConvD3d_sim
!#     -> Calculates the group finite element coefficients for the
!#        linear form in 3D (dual formulation)
!#
!# 10.) transp_calcMatBdrConvD3d_sim
!#      -> Calculates the group finite element coefficients for the
!#         bilinear form in 3D (dual formulation)
!#
!#!# </purpose>
!##############################################################################

module transport_callback3d

  use basicgeometry
  use boundarycondaux
  use collection
  use fparser
  use fsystem
  use genoutput
  use linearsystemscalar
  use linearsystemblock
  use storage

  ! Modules from transport model
  use transport_basic

  implicit none

  private

  ! linear velocity in 3D - primal formulation
  public :: transp_calcMatDiagConvP3d_sim
  public :: transp_calcMatGalConvP3d_sim
  public :: transp_calcMatUpwConvP3d_sim
  public :: transp_calcMatBdrConvP3d_sim
  public :: transp_calcVecBdrConvP3d_sim

  ! linear velocity in 3D - dual formulation
  public :: transp_calcMatDiagConvD3d_sim
  public :: transp_calcMatGalConvD3d_sim
  public :: transp_calcMatUpwConvD3d_sim
  public :: transp_calcMatBdrConvD3d_sim
  public :: transp_calcVecBdrConvD3d_sim

contains
  
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvP3d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D.
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
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(3), p_DvelocityZ)

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode)&
          +p_DvelocityZ(InodeList(1,inode))*DcoeffsAtNode(3,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode)&
          +p_DvelocityZ(InodeList(1,inode))*DcoeffsAtNode(3,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvP3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvP3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D.
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
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(3), p_DvelocityZ)

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
      DmatrixAtEdge(1,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,2,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DmatrixAtEdge(1,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,2,iedge))
#endif
    end do

  end subroutine transp_calcMatGalConvP3d_sim
  
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvP3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(3), p_DvelocityZ)

    if (dscale .gt. 0.0_DP) then

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,2,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,1,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(-DmatrixAtEdge(2,iedge), 0.0_DP,&
                -DmatrixAtEdge(3,iedge))
      end do

    else

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,2,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,1,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(DmatrixAtEdge(2,iedge), 0.0_DP,&
                DmatrixAtEdge(3,iedge))
      end do

    end if

  end subroutine transp_calcMatUpwConvP3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvD3d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D.
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
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(3), p_DvelocityZ)
    
    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode)&
          +p_DvelocityZ(InodeList(1,inode))*DcoeffsAtNode(3,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode)&
          +p_DvelocityZ(InodeList(1,inode))*DcoeffsAtNode(3,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvD3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvD3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D.
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
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>
    
    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(3), p_DvelocityZ)
    
    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
      DmatrixAtEdge(1,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,2,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DmatrixAtEdge(1,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,2,iedge))
#endif
    end do
    
  end subroutine transp_calcMatGalConvD3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvD3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(3), p_DvelocityZ)
    
    if (dscale .gt. 0.0_DP) then

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,2,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,1,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(-DmatrixAtEdge(2,iedge), 0.0_DP,&
                -DmatrixAtEdge(3,iedge))
      end do

    else

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,2,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge)&
            +p_DvelocityZ(IedgeList(2,iedge))*DcoeffsAtEdge(3,1,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge)&
            +p_DvelocityZ(IedgeList(1,iedge))*DcoeffsAtEdge(3,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(DmatrixAtEdge(2,iedge), 0.0_DP,&
                DmatrixAtEdge(3,iedge))
      end do

    end if

  end subroutine transp_calcMatUpwConvD3d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcMatBdrConvP3d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local matrix entries
    ! DmatrixAtNode for the node $i$.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation in 3D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    !   DIMENSION(ncoeffs,nnodes)
    ! with ncoeffs the number of matrix coefficients at the node
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP) :: dtime,dnv
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcMatBdrConvP3d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(3), p_DvelocityZ)
    end if
    
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral (if any)

      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in node
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)

          ! Scale normal velocity by scaling parameter
          DmatrixAtNode(1,inode) = -dscale*dnv
        end do
        
      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

      ! Loop over all noodes at the boundary
      do inode =  1, nnodes

        ! Impose Dirichlet boundary conditions via penalty method
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode))+&
            abs(DcoeffsAtNode(3,inode)) .gt. SYS_EPSREAL_DP) then
          DmatrixAtNode(1,inode) = -dscale*BDRC_DIRICHLET_PENALTY
        else
          DmatrixAtNode(1,inode) = 0.0_DP
        end if
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form

      DmatrixAtNode = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcMatBdrConvP3d_sim')


    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions
      
      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes

          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)

          ! Check if node i is at the primal outflow boundary
          if (dnv .gt. SYS_EPSREAL_DP) then
            DmatrixAtNode(1,inode) = -dscale*dnv
          else
            DmatrixAtNode(1,inode) = 0.0_DP
          end if
        end do

      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if
      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcMatBdrConvP3d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_calcMatBdrConvP3d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcVecBdrConvP3d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DvectorAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local vector entries
    ! DvectorAtNode for the node $i$.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation in 3D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the vector for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(out) :: DvectorAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime,dnv,dval
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcVecBdrConvP3d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(3), p_DvelocityZ)
    end if
        
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)
    
    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))
      
    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      
      DvectorAtNode = 0.0_DP
      
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcVecBdrConvP3d_sim')
      
      
    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Loop over all noodes at the boundary
      do inode =  1,  nnodes
        
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode))+&
            abs(DcoeffsAtNode(3,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
          Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Multiply by scaling coefficient
          DvectorAtNode(inode) = dscale*dval
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow ({\bf v}u)\cdot{\bf n}=({\bf v}g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Loop over all noodes at the boundary
      do inode = 1, nnodes

        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode))+&
            abs(DcoeffsAtNode(3,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
          Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Impose Dirichlet value via penalty method
          DvectorAtNode(inode) = dscale*dval*BDRC_DIRICHLET_PENALTY
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      if (associated(p_rvelocity)) then
        
        ! Loop over all noodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)
          
          if (abs(dnv) .gt. SYS_EPSREAL_DP) then
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
            Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at Robin boundary
            DvectorAtNode(inode) = -dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      
      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      if (associated(p_rvelocity)) then
        
        ! Loop over all nodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)

          ! Check if node i is at the primal outflow boundary
          if (dnv .lt. -SYS_EPSREAL_DP) then
            
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
            Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)

            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at primal outflow boundary
            DvectorAtNode(inode) = -dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      

    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      print *, "Periodic boundary conditions are not implemented yet!"
      stop

    end select

  end subroutine transp_calcVecBdrConvP3d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcMatBdrConvD3d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local matrix entries
    ! DmatrixAtNode for the node $i$.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation in 3D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    !   DIMENSION(ncoeffs,nnodes)
    ! with ncoeffs the number of matrix coefficients at the node
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP) :: dtime,dnv
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcMatBdrConvP3d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(3), p_DvelocityZ)
    end if
    
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral (if any)

      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in node
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)

          ! Scale normal velocity by scaling parameter
          DmatrixAtNode(1,inode) = dscale*dnv
        end do
        
      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

      ! Loop over all noodes at the boundary
      do inode =  1, nnodes

        ! Impose Dirichlet boundary conditions via penalty method
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode))+&
            abs(DcoeffsAtNode(3,inode)) .gt. SYS_EPSREAL_DP) then
          DmatrixAtNode(1,inode) = dscale*BDRC_DIRICHLET_PENALTY
        else
          DmatrixAtNode(1,inode) = 0.0_DP
        end if
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form

      DmatrixAtNode = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcMatBdrConvP3d_sim')


    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions
      
      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes

          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)

          ! Check if node i is at the dual outflow boundary
          if (dnv .lt. -SYS_EPSREAL_DP) then
            DmatrixAtNode(1,inode) = dscale*dnv
          else
            DmatrixAtNode(1,inode) = 0.0_DP
          end if
        end do

      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if
      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcMatBdrConvD3d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_calcMatBdrConvD3d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcVecBdrConvD3d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DvectorAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local vector entries
    ! DvectorAtNode for the node $i$.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation in 3D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the vector for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(out) :: DvectorAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY,p_DvelocityZ
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime,dnv,dval
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcVecBdrConvD3d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(3), p_DvelocityZ)
    end if
        
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)
    
    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))
      
    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.

      DvectorAtNode = 0.0_DP
      
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcVecBdrConvP3d_sim')
      
      
    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Loop over all nodes at the boundary
      do inode =  1,  nnodes
        
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode))+&
            abs(DcoeffsAtNode(3,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
          Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Multiply by scaling coefficient
          DvectorAtNode(inode) = dscale*dval
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow ({\bf v}u)\cdot{\bf n}=({\bf v}g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Loop over all nodes at the boundary
      do inode = 1, nnodes
        
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode))+&
            abs(DcoeffsAtNode(3,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
          Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Impose Dirichlet value via penalty method
          DvectorAtNode(inode) = -dscale*dval*BDRC_DIRICHLET_PENALTY
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form
      !
      ! $$ ({\bf v}u-d\nabla u)\cdot{\bf n} = ({\bf v}g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      if (associated(p_rvelocity)) then
        
        ! Loop over all nodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)
          
          if (abs(dnv) .gt. SYS_EPSREAL_DP) then
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
            Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at Robin boundary
            DvectorAtNode(inode) = dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      
      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ ({\bf v}u-d\nabla u)\cdot{\bf n} = ({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      if (associated(p_rvelocity)) then
        
        ! Loop over all nodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)+&
                DcoeffsAtNode(3,inode)*p_DvelocityZ(i)

          ! Check if node i is at the dual outflow boundary
          if (dnv .gt. SYS_EPSREAL_DP) then
            
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM3D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM3D+2)
            Dvalue(3) = p_DdofCoords((i-1)*NDIM3D+3)

            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at dual outflow boundary
            DvectorAtNode(inode) = dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      

    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      print *, "Periodic boundary conditions are not implemented yet!"
      stop

    end select

  end subroutine transp_calcVecBdrConvD3d_sim

end module transport_callback3d
