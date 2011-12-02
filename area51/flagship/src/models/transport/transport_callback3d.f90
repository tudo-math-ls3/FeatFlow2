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
!#!# </purpose>
!##############################################################################

module transport_callback3d

  use collection
  use fsystem
  use genoutput
  use linearsystemscalar
  use linearsystemblock
  use storage

  implicit none

  private

  public :: transp_calcMatDiagConvP3d_sim
  public :: transp_calcMatGalConvP3d_sim
  public :: transp_calcMatUpwConvP3d_sim

  public :: transp_calcMatDiagConvD3d_sim
  public :: transp_calcMatGalConvD3d_sim
  public :: transp_calcMatUpwConvD3d_sim

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

end module transport_callback3d
