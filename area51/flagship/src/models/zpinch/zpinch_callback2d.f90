!##############################################################################
!# ****************************************************************************
!# <name> hydro_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the simplified MHD equations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) zpinch_hadaptCallbackScalar2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 2.) zpinch_hadaptCallbackBlock2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in block format
!#
!# 3.) zpinch_calcMatDiagConvIntlP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#        for hydrodynamic systems stored in interleaved format
!#
!# 4.) zpinch_calcMatRusConvIntlP2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#        and applies scalar artificial viscosities of Rusanov-type
!#        for hydrodynamic systems stored in interleaved format
!#
!# 5.) zpinch_calcMatDiagConvIntlD2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#        for hydrodynamic systems stored in interleaved format
!#
!# 6.) zpinch_calcMatRusConvIntlD2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#        and applies scalar artificial viscosities of Rusanov-type
!#        for hydrodynamic systems stored in interleaved format
!#
!# 7.) zpinch_calcMatDiagConvBlockP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#        for hydrodynamic systems stored in block format
!#
!# 10.) zpinch_calcMatRusConvBlockP2d_sim
!#      -> Calculates the off-diagonal Galerkin transport coefficients
!#         for linear convection in 2D (primal formulation)
!#         and applies scalar artificial viscosities of Rusanov-type
!#         for hydrodynamic systems stored in block format
!#
!# 11.) zpinch_calcMatDiagConvBlockD2d_sim
!#      -> Calculates the diagonal Galerkin transport coefficients
!#         for linear convection in 2D (dual formulation)
!#         for hydrodynamic systems stored in block format
!#
!# 12.) zpinch_calcMatRusConvBlockD2d_sim
!#      -> Calculates the off-diagonal Galerkin transport coefficients
!#         for linear convection in 2D (dual formulation)
!#         and applies scalar artificial viscosities of Rusanov-type
!#         for hydrodynamic systems stored in block format
!#
!# </purpose>
!##############################################################################

module zpinch_callback2d

#include "hydro.h"

  use collection
  use hydro_basic
  use flagship_callback
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use storage

  implicit none

  private
  public :: zpinch_hadaptCallbackScalar2d
  public :: zpinch_hadaptCallbackBlock2d
  
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

  subroutine zpinch_hadaptCallbackScalar2d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in scalar interleave format.
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
    type(t_vectorBlock), pointer, save :: rsolutionHydro, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionHydro, p_DsolutionTransport
    integer :: ivar


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vectors from colletion and set pointer
      rsolutionHydro     => rcollection%p_rvectorQuickAccess1
      rsolutionTransport => rcollection%p_rvectorQuickAccess2

      ! Check if solution is stored in interleave format
      if (rsolutionHydro%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_hadaptCallbackScalar2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionHydro, p_DsolutionHydro)
      nullify(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .ne. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
            0.5_DP*(p_DsolutionHydro((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
                    p_DsolutionHydro((rcollection%IquickAccess(3)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.5_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                  p_DsolutionTransport(rcollection%IquickAccess(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)
      

    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
            0.25_DP*(p_DsolutionHydro((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
                     p_DsolutionHydro((rcollection%IquickAccess(3)-1)*NVAR2D+ivar)+&
                     p_DsolutionHydro((rcollection%IquickAccess(4)-1)*NVAR2D+ivar)+&
                     p_DsolutionHydro((rcollection%IquickAccess(5)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.25_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                   p_DsolutionTransport(rcollection%IquickAccess(3))+&
                   p_DsolutionTransport(rcollection%IquickAccess(4))+&
                   p_DsolutionTransport(rcollection%IquickAccess(5)))


      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the hydrodynamic model
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR2D
          p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
              p_DsolutionHydro((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)
        end do
      else
        do ivar = 1, NVAR2D
          p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (rcollection%IquickAccess(2) .ne. 0) then
        p_DsolutionTransport(rcollection%IquickAccess(1)) =&
            p_DsolutionTransport(rcollection%IquickAccess(2))
      else
        p_DsolutionTransport(rcollection%IquickAccess(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine zpinch_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_hadaptCallbackBlock2d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in block format.
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
    type(t_vectorBlock), pointer, save :: rsolutionHydro, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionHydro, p_DsolutionTransport
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vectors from colletion and set pointer
      rsolutionHydro     => rcollection%p_rvectorQuickAccess1
      rsolutionTransport => rcollection%p_rvectorQuickAccess2

      ! Check if solution is stored in interleave format
      if (rsolutionHydro%nblocks .ne. NVAR2D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_hadaptCallbackBlock2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionHydro, p_DsolutionHydro)
      nullify(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .ne. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro, NVAR2D&
            *rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      neq = rsolutionHydro%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5_DP*(p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.5_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                  p_DsolutionTransport(rcollection%IquickAccess(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro, NVAR2D&
            *rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      neq = rsolutionHydro%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) =&
            0.25_DP*(p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(2))+&
                     p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(3))+&
                     p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(4))+&
                     p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(5)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.25_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                   p_DsolutionTransport(rcollection%IquickAccess(3))+&
                   p_DsolutionTransport(rcollection%IquickAccess(4))+&
                   p_DsolutionTransport(rcollection%IquickAccess(5)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the hydrodynamic model
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolutionHydro%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolutionHydro%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (rcollection%IquickAccess(2) .ne. 0) then
        p_DsolutionTransport(rcollection%IquickAccess(1)) =&
            p_DsolutionTransport(rcollection%IquickAccess(2))
      else
        p_DsolutionTransport(rcollection%IquickAccess(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine zpinch_hadaptCallbackBlock2d

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvIntlP2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
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
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))

#else
      ! Compute convective coefficient  $k_{ii} = -v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvIntlP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvIntlP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient  $k_{ji} = v_i*Cx_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient  $k_{ij} = -v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient  $k_{ji} = -v_i*Cx_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Compute base indices
      idx = 4*(IverticesAtEdge(1,iedge)-1)
      jdx = 4*(IverticesAtEdge(2,iedge)-1)
      
      ! Compute auxiliary variables
      ui = p_DsolutionHydro(idx+2)/p_DsolutionHydro(idx+1)
      vi = p_DsolutionHydro(idx+3)/p_DsolutionHydro(idx+1)
      Ei = p_DsolutionHydro(idx+4)/p_DsolutionHydro(idx+1)
      uj = p_DsolutionHydro(jdx+2)/p_DsolutionHydro(jdx+1)
      vj = p_DsolutionHydro(jdx+3)/p_DsolutionHydro(jdx+1)
      Ej = p_DsolutionHydro(jdx+4)/p_DsolutionHydro(jdx+1)
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))**2)*ci )
#else
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(DmatrixCoeffsAtEdge(1,1,iedge)*uj +&
                   DmatrixCoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,iedge)*ui +&
                   DmatrixCoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvIntlP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvIntlD2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
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
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient  $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvIntlD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvIntlD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji] = -v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji] = v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

      ! Compute base indices
      idx = 4*(IverticesAtEdge(1,iedge)-1)
      jdx = 4*(IverticesAtEdge(2,iedge)-1)
      
      ! Compute auxiliary variables
      ui = p_DsolutionHydro(idx+2)/p_DsolutionHydro(idx+1)
      vi = p_DsolutionHydro(idx+3)/p_DsolutionHydro(idx+1)
      Ei = p_DsolutionHydro(idx+4)/p_DsolutionHydro(idx+1)
      uj = p_DsolutionHydro(jdx+2)/p_DsolutionHydro(jdx+1)
      vj = p_DsolutionHydro(jdx+3)/p_DsolutionHydro(jdx+1)
      Ej = p_DsolutionHydro(jdx+4)/p_DsolutionHydro(jdx+1)
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))**2)*ci )
#else 
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(DmatrixCoeffsAtEdge(1,1,iedge)*uj +&
                   DmatrixCoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,iedge)*ui +&
                   DmatrixCoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvIntlD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvBlockP2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
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
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvBlockP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvBlockP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Get number of equations for single variable
      neq = size(p_DvelocityX)

      ! Get equations i and j
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
      
      ! Compute auxiliary variables
      ui = p_DsolutionHydro(neq  +i)/p_DsolutionHydro(i)
      vi = p_DsolutionHydro(neq*2+i)/p_DsolutionHydro(i)
      Ei = p_DsolutionHydro(neq*3+i)/p_DsolutionHydro(i)
      uj = p_DsolutionHydro(neq  +j)/p_DsolutionHydro(j)
      vj = p_DsolutionHydro(neq*2+j)/p_DsolutionHydro(j)
      Ej = p_DsolutionHydro(neq*3+j)/p_DsolutionHydro(j)
            
! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))**2)*ci )
#else      
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(DmatrixCoeffsAtEdge(1,1,iedge)*uj +&
                   DmatrixCoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,iedge)*ui +&
                   DmatrixCoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvBlockP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatDiagConvBlockD2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
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
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine zpinch_calcMatDiagConvBlockD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvBlockD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)
    
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
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
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
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji] = -v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji] = v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Get number of equations for single variable
      neq = size(p_DvelocityX)
      
      ! Get equations i and j
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)

      ! Compute auxiliary variables
      ui = p_DsolutionHydro(neq  +i)/p_DsolutionHydro(i)
      vi = p_DsolutionHydro(neq*2+i)/p_DsolutionHydro(i)
      Ei = p_DsolutionHydro(neq*3+i)/p_DsolutionHydro(i)
      uj = p_DsolutionHydro(neq  +j)/p_DsolutionHydro(j)
      vj = p_DsolutionHydro(neq*2+j)/p_DsolutionHydro(j)
      Ej = p_DsolutionHydro(neq*3+j)/p_DsolutionHydro(j)
      
! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0_DP)*GAMMA*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))*uj+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))*vj)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,1,iedge)-&
                           DmatrixCoeffsAtEdge(1,2,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,1,iedge)-&
                           DmatrixCoeffsAtEdge(2,2,iedge))**2)*cj,&
               abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))*ui+&
                   0.5_DP*(DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))*vi)+&
              0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,2,iedge)-&
                           DmatrixCoeffsAtEdge(1,1,iedge))**2+&
                          (DmatrixCoeffsAtEdge(2,2,iedge)-&
                           DmatrixCoeffsAtEdge(2,1,iedge))**2)*ci )
#else
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max( abs(DmatrixCoeffsAtEdge(1,1,iedge)*uj +&
                   DmatrixCoeffsAtEdge(2,1,iedge)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,1,iedge)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,iedge)*ui +&
                   DmatrixCoeffsAtEdge(2,2,iedge)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,iedge)**2 +&
                        DmatrixCoeffsAtEdge(2,2,iedge)**2)*ci )
#endif
    end do
    
  end subroutine zpinch_calcMatRusConvBlockD2d_sim

end module zpinch_callback2d
