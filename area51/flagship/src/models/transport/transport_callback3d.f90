!##############################################################################
!# ****************************************************************************
!# <name> transport_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in 3D.
!#
!# The following routines are available:
!#
!# 1.) transp_hadaptCallback3d
!#      -> Performs application specific tasks in the adaptation algorithm in 3D
!#
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
  use flagship_callback
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemscalar
  use linearsystemblock
  use storage

  implicit none

  private

  public :: transp_hadaptCallback3d

  public :: transp_calcMatDiagConvP3d_sim
  public :: transp_calcMatGalConvP3d_sim
  public :: transp_calcMatUpwConvP3d_sim

  public :: transp_calcMatDiagConvD3d_sim
  public :: transp_calcMatGalConvD3d_sim
  public :: transp_calcMatUpwConvD3d_sim

contains
  
  !*****************************************************************************

!<subroutine>

  subroutine transp_hadaptCallback3d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   IquickAccess(1):     NEQ or ivt
    !   IquickAccess(2:5):   ivt1,...,ivt5
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution


    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion and set pointer
      rsolution => rcollection%p_rvectorQuickAccess1
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      p_Dsolution(rcollection%IquickAccess(1)) =&
          0.5_DP*(p_Dsolution(rcollection%IquickAccess(2))+&
                  p_Dsolution(rcollection%IquickAccess(3)))

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      p_Dsolution(rcollection%IquickAccess(1)) =&
          0.25_DP*(p_Dsolution(rcollection%IquickAccess(2))+&
                   p_Dsolution(rcollection%IquickAccess(3))+&
                   p_Dsolution(rcollection%IquickAccess(4))+&
                   p_Dsolution(rcollection%IquickAccess(5)))

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        p_Dsolution(rcollection%IquickAccess(1)) =&
            p_Dsolution(rcollection%IquickAccess(2))
      else
        p_Dsolution(rcollection%IquickAccess(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)
    end select

  end subroutine transp_hadaptCallback3d

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvP3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtNode
    
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
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
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
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode)&
          +p_DvelocityZ(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(3,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode)&
          +p_DvelocityZ(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(3,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvP3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvP3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

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
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
#endif

      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0
    end do

  end subroutine transp_calcMatGalConvP3d_sim
  
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvP3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

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
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

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
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
#endif

      ! Compute artificial diffusion coefficient 
      !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
      DcoefficientsAtEdge(1,iedge) =&
          max(-DcoefficientsAtEdge(2,iedge), 0.0_DP,&
              -DcoefficientsAtEdge(3,iedge))
    end do

  end subroutine transp_calcMatUpwConvP3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvD3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtNode

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
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
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
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode)&
          +p_DvelocityZ(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(3,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode)&
          +p_DvelocityZ(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(3,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvD3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

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
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
#endif

      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0
    end do
    
  end subroutine transp_calcMatGalConvD3d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

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
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

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
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_DvelocityZ(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_DvelocityZ(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
#endif

      ! Compute artificial diffusion coefficient
      !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
      DcoefficientsAtEdge(1,iedge) =&
          max(-DcoefficientsAtEdge(2,iedge), 0.0_DP,&
              -DcoefficientsAtEdge(3,iedge))
    end do

  end subroutine transp_calcMatUpwConvD3d_sim

end module transport_callback3d
