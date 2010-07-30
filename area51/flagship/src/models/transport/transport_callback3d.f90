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
!# 1.) transp_setVariable3d
!#     -> Sets global variables for external data, e.g., velocity fields in 3D
!#
!# 2.) transp_hadaptCallback3d
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
!4
!# 3.) transp_calcMatGalConvectionP3d
!#     -> Calculates the Galerkin transport coefficients
!#        for linear convection in 3D (primal formulation)
!#
!# 4.) transp_calcMatUpwConvP3d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 3D (primal formulation)
!#
!# 5.) transp_calcMatUpwConvectionP3d
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 3D (primal formulation)
!#
!# 6.) transp_calcMatDiagConvD3d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 3D (dual formulation)
!#
!# 7.) transp_calcMatGalConvD3d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 3D (dual formulation)
!4
!# 8.) transp_calcMatGalConvectionD3d
!#     -> Calculates the Galerkin transport coefficients
!#        for linear convection in 3D (dual formulation)
!#
!# 9.) transp_calcMatUpwConvD3d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 3D (dual formulation)
!#
!# 10.) transp_calcMatUpwConvectionD3d
!#      -> Calculates the Galerkin transport coefficients
!#         and applies scalar artificial diffusion (discrete upwinding)
!#         for linear convection in 3D (dual formulation)
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

  public :: transp_setVariable3d
  public :: transp_hadaptCallback3d

  public :: transp_calcMatDiagConvP3d_sim
  public :: transp_calcMatGalConvectionP3d
  public :: transp_calcMatGalConvP3d_sim
  public :: transp_calcMatUpwConvectionP3d
  public :: transp_calcMatUpwConvP3d_sim

  public :: transp_calcMatDiagConvD3d_sim
  public :: transp_calcMatGalConvectionD3d
  public :: transp_calcMatGalConvD3d_sim
  public :: transp_calcMatUpwConvectionD3d
  public :: transp_calcMatUpwConvD3d_sim

!<globals>

  !*****************************************************************
  ! Pointers to external data vectors.
  !
  ! Using global variables is not good programming style but it is the
  ! only way to allow for an efficient access to the velocity data
  ! from within the callback routines which are called repeatedly

  real(DP), dimension(:), pointer, save :: p_Dvariable1 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable2 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable3 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable4 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable5 => null()

!</globals>

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_setVariable3d(rvector, ivariable)

!<description>
    ! This subroutine sets one of the the global pointers to the given vector.
!</description>

!<input>
    ! scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! variable number
    integer, intent(in) :: ivariable
!</input>
!</subroutine>

    select case(ivariable)
    case (1)
      call lsyssc_getbase_double(rvector, p_Dvariable1)
    case (2)
      call lsyssc_getbase_double(rvector, p_Dvariable2)
    case (3)
      call lsyssc_getbase_double(rvector, p_Dvariable3)
    case (4)
      call lsyssc_getbase_double(rvector, p_Dvariable4)
    case (5)
      call lsyssc_getbase_double(rvector, p_Dvariable5)
    case DEFAULT
      call output_line('Invalid variable number!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_setVariable3d')
      call sys_halt()
    end select

  end subroutine transp_setVariable3d

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
    ! Collection
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


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)
    end select

  end subroutine transp_hadaptCallback3d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagConvP3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
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
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_Dvelocity
    integer :: inode

!!!    ! Set pointer to velocity vector
!!!    p_Dvelocity => collct_getvalue_vec(rcollection, 'velocity')
    
    do inode = 1, size(DcoefficientsAtNode,2)
      ! Compute convective coefficient  $-v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_Dvariable1(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_Dvariable2(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode)&
          +p_Dvariable3(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(3,inode))
    end do
    
  end subroutine transp_calcMatDiagConvP3d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalConvP3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
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
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_Dvelocity
    integer :: iedge

!!!    ! Set pointer to velocity vector
!!!    p_Dvelocity => collct_getvalue_vec(rcollection, 'velocity')
    
    do iedge = 1, size(DcoefficientsAtEdge,2)
      ! Compute convective coefficient  $-v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_Dvariable1(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_Dvariable2(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_Dvariable3(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient  $-v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_Dvariable1(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_Dvariable2(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_Dvariable3(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0
    end do

  end subroutine transp_calcMatGalConvP3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatGalConvectionP3d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D.
!</description>

!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)-p_Dvariable3(j)*C_ij(3)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)-p_Dvariable3(i)*C_ji(3)

    ! Set artificial diffusion to zero
    d_ij = 0.0_DP

  end subroutine transp_calcMatGalConvectionP3d
  
  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwConvP3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
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
    !</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_Dvelocity
    integer :: iedge

!!!    ! Set pointer to velocity vector
!!!    p_Dvelocity => collct_getvalue_vec(rcollection, 'velocity')
    
    do iedge = 1, size(DcoefficientsAtEdge,2)
      ! Compute convective coefficient  $-v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_Dvariable1(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_Dvariable2(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_Dvariable3(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient  $-v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_Dvariable1(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_Dvariable2(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_Dvariable3(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Compute artificial diffusion coefficient
      DcoefficientsAtEdge(1,iedge) =&
          max(-DcoefficientsAtEdge(2,iedge), 0.0_DP, -DcoefficientsAtEdge(3,iedge))
    end do

  end subroutine transp_calcMatUpwConvP3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatUpwConvectionP3d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)-p_Dvariable3(j)*C_ij(3)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)-p_Dvariable3(i)*C_ji(3)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)

  end subroutine transp_calcMatUpwConvectionP3d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagConvD3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
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
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_Dvelocity
    integer :: inode

!!!    ! Set pointer to velocity vector
!!!    p_Dvelocity => collct_getvalue_vec(rcollection, 'velocity')
    
    do inode = 1, size(DcoefficientsAtNode,2)
      ! Compute convective coefficient  $v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_Dvariable1(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_Dvariable2(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode)&
          +p_Dvariable3(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(3,inode))
    end do
    
  end subroutine transp_calcMatDiagConvD3d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalConvD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
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
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>
    
    ! local variable
    real(DP), dimension(:), pointer :: p_Dvelocity
    integer :: iedge

!!!    ! Set pointer to velocity vector
!!!    p_Dvelocity => collct_getvalue_vec(rcollection, 'velocity')
    
    do iedge = 1, size(DcoefficientsAtEdge,2)
      ! Compute convective coefficient  $v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_Dvariable1(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_Dvariable2(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_Dvariable3(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient  $v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_Dvariable1(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_Dvariable2(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_Dvariable3(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0
    end do
    
  end subroutine transp_calcMatGalConvD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatGalConvectionD3d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = p_Dvariable1(j)*C_ij(1)+p_Dvariable2(j)*C_ij(2)+p_Dvariable3(j)*C_ij(3)
    k_ji = p_Dvariable1(i)*C_ji(1)+p_Dvariable2(i)*C_ji(2)+p_Dvariable3(i)*C_ji(3)

    ! Set artificial diffusion to zero
    d_ij = 0.0_DP

  end subroutine transp_calcMatGalConvectionD3d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwConvD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
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
    !</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_Dvelocity
    integer :: iedge

!!!    ! Set pointer to velocity vector
!!!    p_Dvelocity => collct_getvalue_vec(rcollection, 'velocity')
    
    do iedge = 1, size(DcoefficientsAtEdge,2)
      ! Compute convective coefficient  $v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_Dvariable1(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_Dvariable2(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge)&
          +p_Dvariable3(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(3,1,iedge))
      ! Compute convective coefficient  $v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_Dvariable1(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_Dvariable2(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge)&
          +p_Dvariable3(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(3,2,iedge))
      ! Compute artificial diffusion coefficient
      DcoefficientsAtEdge(1,iedge) =&
          max(-DcoefficientsAtEdge(2,iedge), 0.0_DP, -DcoefficientsAtEdge(3,iedge))
    end do

  end subroutine transp_calcMatUpwConvD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatUpwConvectionD3d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = p_Dvariable1(j)*C_ij(1)+p_Dvariable2(j)*C_ij(2)+p_Dvariable3(j)*C_ij(3)
    k_ji = p_Dvariable1(i)*C_ji(1)+p_Dvariable2(i)*C_ji(2)+p_Dvariable3(i)*C_ji(3)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)

  end subroutine transp_calcMatUpwConvectionD3d

end module transport_callback3d
