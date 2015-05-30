!##############################################################################
!# ****************************************************************************
!# <name> sse_callback_poisson </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the SSE problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!#
!#
!# 1.) sse_initMatVec_Poisson
!#     -> Initialises the matrices/vectors for the Poisson problem
!#
!# 2.) sse_postprocessing_Poisson
!#     -> Post-processes the solution (output and error computation)
!#
!# 3.) sse_outputTable_Poisson
!#     -> Output convergence table
!#
!# These callback routines are provided:
!#
!# 1.) coeff_Matrix_Poisson
!#     -> Returns analytic valyes for the system matrix.
!#
!# 2.) coeff_RHS_Poisson
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 3.) coeff_RHS_Bdr_Poisson
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 4.) getBoundaryValues_Poisson
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#
!# 5.) getReferenceFunction_Poisson
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# </purpose>
!##############################################################################

module sse_callback_poisson

  use basicgeometry
  use bilinearformevaluation
  use blockmatassemblybase
  use blockmatassemblystdop
  use boundary
  use collection
  use convergencetable
  use cubature
  use derivatives
  use discretebc
  use discretefbc
  use domainintegration
  use element
  use feevaluation2
  use fparser
  use fsystem
  use genoutput
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use pprocerror
  use pprocgradients
  use scalarpde
  use spatialdiscretisation
  use storage
  use triangulation
  use ucd

  use sse_base
  use sse_base_poisson

  implicit none

  private

  public :: sse_initMatVec_Poisson
  public :: sse_postprocessing_Poisson
  public :: sse_outputTable_Poisson

  public :: coeff_Matrix_Poisson
  public :: coeff_RHS_Poisson
  public :: coeff_RHS_Bdr_Poisson
  public :: getBoundaryValues_Poisson
  public :: getReferenceFunction_Poisson

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initMatVec_Poisson(rproblem)

!<description>
    ! Calculates the system matrix and RHS vector of the linear system
    ! by discretising the problem with the default discretisation structure
    ! in the problem structure.
    ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform

    select case(rproblem%cproblemType)
    case (POISSON_SCALAR)
      !---------------------------------------------------------------------------
      !
      ! Problem formulation in (test,trial)-notation
      !
      ! (grad_x,dpoisson*grad_x)*u + (grad_y,dpoisson*grad_y)*u = (func,f) + Neumann bc.s
      !

      ! Prepare quick-access arrays of collection structure
      rproblem%rcollection%IquickAccess(1) = rproblem%cproblemType

      do i=rproblem%ilvmin,rproblem%ilvmax

        ! Initialise the block matrix with default values based on
        ! the discretisation.
        call lsysbl_createMatrix(rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rmatrix)

        ! Save matrix to the collection.
        ! They may be used later, especially in nonlinear problems.
        call collct_setvalue_mat(rproblem%rcollection,"MATRIX",&
            rproblem%RlevelInfo(i)%rmatrix,.true.,i)

        ! Now as the discretisation is set up, we can start to generate
        ! the structure of the system matrix which is to solve.
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,1,LSYSSC_MATRIX9)

        ! (1,1)-block (grad_x w,grad_x u)+(grad_y w,grad_y w)
        rform%itermCount = 2
        rform%ballCoeffConstant = (cpoisson .eq. 0)
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Dcoefficients(1)  = dpoisson
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_Y
        rform%Dcoefficients(2)  = dpoisson

        ! Assemble matrix block(1,1)
        rproblem%rcollection%IquickAccess(2) = 11
        rproblem%rcollection%IquickAccess(3) = cpoisson
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_Matrix_Poisson,rproblem%rcollection)
      end do

      ! Next step: Create a RHS vector and a solution vector and a temporary
      ! vector. All are filled with zero.
      call lsysbl_createVector(&
          rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
          rproblem%rrhs,.true.)
      call lsysbl_createVector(&
          rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
          rproblem%rvector,.true.)

      ! Save the solution/RHS vector to the collection.
      call collct_setvalue_vec(rproblem%rcollection,"RHS",rproblem%rrhs,.true.)
      call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",rproblem%rvector,.true.)

      ! The vector structure is ready but the entries are missing.
      ! So the next thing is to calculate the content of that vector.
      !
      ! At first set up the corresponding linear form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC

      ! ... and then discretise the RHS to the first subvector of
      ! the block vector using the discretisation structure of the
      ! first block.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      rproblem%rcollection%IquickAccess(2) = 1
      rproblem%rcollection%IquickAccess(3) = crhs
      call linf_buildVectorScalar(&
          rlinform,.true.,rproblem%rrhs%RvectorBlock(1),&
          rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1),&
          coeff_RHS_Poisson,rproblem%rcollection)

      ! Loop over all boundary conditions
      do j=1,size(rproblem%RboundaryCondition)

        if (rproblem%RboundaryCondition(j)%cboundaryConditionType .eq. BDRCOND_NEUMANN) then

          ! Set position, equation number and boundary component and assemble linear form
          rproblem%rcollection%IquickAccess(2) = rproblem%RboundaryCondition(j)%iequation
          rproblem%rcollection%IquickAccess(3) = rproblem%RboundaryCondition(j)%iexpression
          rproblem%rcollection%IquickAccess(4) = rproblem%RboundaryCondition(j)%rboundaryRegion%iboundCompIdx

          ! Linear forms for Neumann boundary conditions        
          call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
              rproblem%rrhs%RvectorBlock(rproblem%RboundaryCondition(j)%iequation),&
              coeff_RHS_Bdr_Poisson,rproblem%RboundaryCondition(j)%rboundaryRegion,&
              rproblem%rcollection)
        end if
      end do

    case (POISSON_SYSTEM)
      !---------------------------------------------------------------------------
      !
      ! Problem formulation in (test,trial)-notation
      !
      !  /                                                           \   /       \   /           \
      ! |      0        (func,dpoisson*grad_x) (func,dpoisson*grad_y) | |    u    | | (func,-f)   |
      ! |                                                             | |         | |             |
      ! | (grad_x,func) (func,func)                 0                 |*| sigma_x |=| <func,g*nx> |
      ! |                                                             | |         | |             |
      ! | (grad_y,func)      0                 (func,func)            | | sigma_y | | <func,g*ny> |
      !  \                                                           /   \       /   \           /
      !
      ! where <.,.> denotes the surface integral over the Dirichlet(!) boundary parts

      ! Prepare quick-access arrays of collection structure
      rproblem%rcollection%IquickAccess(1) = rproblem%cproblemType

      do i=rproblem%ilvmin,rproblem%ilvmax

        ! Initialise the block matrix with default values based on
        ! the discretisation.
        call lsysbl_createMatrix(rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rmatrix)

        ! Save matrix to the collection.
        ! They may be used later, especially in nonlinear problems.
        call collct_setvalue_mat(rproblem%rcollection,"MATRIX",&
            rproblem%RlevelInfo(i)%rmatrix,.true.,i)

        ! Now as the discretisation is set up, we can start to generate
        ! the structure of the system matrix which is to solve.
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,2,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,3,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,1,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,2,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,1,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,3,LSYSSC_MATRIX9)

        ! (1,2)-block (w,-grad_x sigma_x)
        rform%itermCount = 1
        rform%ballCoeffConstant = (cpoisson .eq. 0)
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients(1)  = dpoisson

        ! Assemble matrix block(1,2)
        rproblem%rcollection%IquickAccess(2) = 12
        rproblem%rcollection%IquickAccess(3) = cpoisson
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_Matrix_Poisson,rproblem%rcollection)

        ! (1,3)-block (w,-grad_y sigma_y)
        rform%itermCount = 1
        rform%ballCoeffConstant = (cpoisson .eq. 0)
        rform%Idescriptors(1,1) = DER_DERIV_Y
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients(1)  = dpoisson

        ! Assemble matrix block(1,3)
        rproblem%rcollection%IquickAccess(2) = 13
        rproblem%rcollection%IquickAccess(3) = cpoisson
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,3),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_Matrix_Poisson,rproblem%rcollection)

        ! (2,1)-block (grad_x v_x,u)
        rform%itermCount = 1
        rform%ballCoeffConstant = (cpoisson .eq. 0)
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Dcoefficients(1)  = 1.0

        ! Assemble matrix block(2,1)
        rproblem%rcollection%IquickAccess(2) = 21
        rproblem%rcollection%IquickAccess(3) = 0
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(2),&
            coeff_Matrix_Poisson,rproblem%rcollection)

        ! (3,1)-block (grad_y v_y,u)
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_DERIV_Y
        rform%Dcoefficients(1)  = 1.0

        ! Assemble matrix block(3,1)
        rproblem%rcollection%IquickAccess(2) = 31
        rproblem%rcollection%IquickAccess(3) = 0
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(3),&
            coeff_Matrix_Poisson,rproblem%rcollection)

        ! (2,2)-block (v_x,sigma_x)
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients(1)  = 1.0

        ! Assemble matrix block(2,2)
        rproblem%rcollection%IquickAccess(2) = 22
        rproblem%rcollection%IquickAccess(3) = 0
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,2),&
            rproblem%RlevelInfo(i)%RcubatureInfo(2),&
            coeff_Matrix_Poisson,rproblem%rcollection)

        ! (3,3)-block (v_y,sigma_y)
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients(1)  = 1.0

        ! Assemble matrix block(3,3)
        rproblem%rcollection%IquickAccess(2) = 33
        rproblem%rcollection%IquickAccess(3) = 0
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,3),&
            rproblem%RlevelInfo(i)%RcubatureInfo(3),&
            coeff_Matrix_Poisson,rproblem%rcollection)
      end do

      ! Next step: Create a RHS vector and a solution vector and a temporary
      ! vector. All are filled with zero.
      call lsysbl_createVector(&
          rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
          rproblem%rrhs,.true.)
      call lsysbl_createVector(&
          rproblem%RlevelInfo(rproblem%ilvmax)%rdiscretisation,&
          rproblem%rvector,.true.)

      ! Save the solution/RHS vector to the collection.
      call collct_setvalue_vec(rproblem%rcollection,"RHS",rproblem%rrhs,.true.)
      call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",rproblem%rvector,.true.)

      ! The vector structure is ready but the entries are missing.
      ! So the next thing is to calculate the content of that vector.
      !
      ! At first set up the corresponding linear form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC

      ! ... and then discretise the RHS to the first subvector of
      ! the block vector using the discretisation structure of the
      ! first block.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      rproblem%rcollection%IquickAccess(2) = 1
      rproblem%rcollection%IquickAccess(3) = crhs
      call linf_buildVectorScalar(&
          rlinform,.true.,rproblem%rrhs%RvectorBlock(1),&
          rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1),&
          coeff_RHS_Poisson,rproblem%rcollection)

      ! Loop over all boundary conditions
      do j=1,size(rproblem%RboundaryCondition)

        if (rproblem%RboundaryCondition(j)%cboundaryConditionType .eq. BDRCOND_DIRICHLET) then

          ! Set position, equation number and boundary component and assemble linear form
          rproblem%rcollection%IquickAccess(2) = rproblem%RboundaryCondition(j)%iequation
          rproblem%rcollection%IquickAccess(3) = rproblem%RboundaryCondition(j)%iexpression
          rproblem%rcollection%IquickAccess(4) = rproblem%RboundaryCondition(j)%rboundaryRegion%iboundCompIdx

          ! Linear forms for Neumann boundary conditions        
          call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
              rproblem%rrhs%RvectorBlock(rproblem%RboundaryCondition(j)%iequation),&
              coeff_RHS_Bdr_Poisson,rproblem%RboundaryCondition(j)%rboundaryRegion,&
              rproblem%rcollection)
        end if
      end do

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initMatVec_Poisson")
      call sys_halt()
    end select

  end subroutine sse_initMatVec_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine sse_postprocessing_Poisson(rproblem,cgradType,cgradSubtype,&
      sucdDir,sucdFile,cucdType,rtable)

!<description>
    ! Exports the solution into UCD file. If the optional parameter
    ! rtable is given, then the error norms are also added to the
    ! convergence table.
!</description>

!<input>
    ! Type of gradient recovery
    integer, intent(in) :: cgradType

    ! Subtype of gradient recovery
    integer, intent(in) :: cgradSubtype

    ! Directory name for UCD output
    character(len=*), intent(in) :: sucdDir

    ! File name for UCD output
    character(len=*), intent(in) :: sucdFile

    ! Type of UCD output
    integer, intent(in) :: cucdType
!</input>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem

    ! OPTIONAL: Convergence table
    type(t_convergenceTable), intent(inout), optional :: rtable
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_ucdExport) :: rexport
    type(t_vectorBlock), target :: rvectorBlock,rvectorBlockX,rvectorBlockY
    type(t_vectorScalar), pointer :: p_rvectorDerivX,p_rvectorDerivY
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfoDerivX,p_rcubatureInfoDerivY
    real(DP) :: derror
    integer :: i

    select case(rproblem%cproblemType)
    case (POISSON_SCALAR,POISSON_SYSTEM)

      ! --- solution -------------------------------------------------------------

      if (present(rtable)) then
        ! Calculate the error to the reference function.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csol
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_x
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csol_y

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rproblem%rvector%RvectorBlock(1),&
            getReferenceFunction_Poisson,rproblem%rcollection,rcubatureInfo=&
            rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
        call ctab_addValue(rtable, "L2-error u", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rproblem%rvector%RvectorBlock(1),&
            getReferenceFunction_Poisson,rproblem%rcollection,rcubatureInfo=&
            rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
        call ctab_addValue(rtable, "H1-error u", derror)
      end if

      ! --- first derivative -----------------------------------------------------

      select case(rproblem%cproblemType)
      case (POISSON_SCALAR)
        ! Recover gradient by superconvergent patch recovery

        ! Create new block discretisation of size NDIM
        call spdiscr_initBlockDiscr(rblockDiscr,sse_getNDIM(rproblem%cproblemType),&
            rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
            rproblem%rvector%p_rblockDiscr%p_rboundary)

        ! Duplicate spatial discretisation from solution variable
        do i=1,sse_getNDIM(rproblem%cproblemType)
          call spdiscr_duplicateDiscrSc(&
              rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
              rblockDiscr%RspatialDiscr(i), .true.)
        end do

        ! Create block vector for the gradient
        call lsysbl_createVector(rblockDiscr, rvectorBlock, .false.)

        ! Recover the gradient vector
        call ppgrd_calcGradient(rproblem%rvector%RvectorBlock(1),&
            rvectorBlock, cgradType, cgradSubType)

        ! Set pointers to recovered gradient vector
        p_rvectorDerivX => rvectorBlock%RvectorBlock(1)
        p_rvectorDerivY => rvectorBlock%RvectorBlock(2)

        ! Set pointers to cubature structure from scalar field
        p_rcubatureInfoDerivX => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1)
        p_rcubatureInfoDerivY => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1)

      case (POISSON_SYSTEM)
        ! Set pointer to scalar solution components
        p_rvectorDerivX => rproblem%rvector%RvectorBlock(2)
        p_rvectorDerivY => rproblem%rvector%RvectorBlock(3)

        ! Set pointers to cubature structure from vector field
        p_rcubatureInfoDerivX => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(2)
        p_rcubatureInfoDerivY => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(3)
      end select

      if (present(rtable)) then
        ! Calculate the error to the reference x-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csol_x
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_xx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csol_xy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivX,&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivX)
        call ctab_addValue(rtable, "L2-error u_x", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivX,&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivX)
        call ctab_addValue(rtable, "H1-error u_x", derror)

        ! Calculate the error to the reference y-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csol_y
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_xy ! = csol_yx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csol_yy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivY,&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivY)
        call ctab_addValue(rtable, "L2-error u_y", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivY,&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivY)
        call ctab_addValue(rtable, "H1-error u_y", derror)
      end if

      ! --- second derivative ----------------------------------------------------

      select case(rproblem%cproblemType)
      case (POISSON_SYSTEM)
        ! Recover gradient by superconvergent patch recovery

        ! Create new block discretisation of size NDIM
        call spdiscr_initBlockDiscr(rblockDiscr,sse_getNDIM(rproblem%cproblemType),&
            rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
            rproblem%rvector%p_rblockDiscr%p_rboundary)

        ! Duplicate spatial discretisation from solution variable
        do i=1,sse_getNDIM(rproblem%cproblemType)
          call spdiscr_duplicateDiscrSc(&
              rproblem%rvector%RvectorBlock(1)%p_rspatialDiscr,&
              rblockDiscr%RspatialDiscr(i), .true.)
        end do
      end select

      ! Create block vectors for the second gradient
      call lsysbl_createVector(rblockDiscr, rvectorBlockX, .false.)
      call lsysbl_createVector(rblockDiscr, rvectorBlockY, .false.)

      ! Recover the gradient vector of u_x
      call ppgrd_calcGradient(p_rvectorDerivX, rvectorBlockX,&
          cgradType, cgradSubType)

      ! Recover the gradient vector of u_y
      call ppgrd_calcGradient(p_rvectorDerivY, rvectorBlockY,&
          cgradType, cgradSubType)

      if (present(rtable)) then
        ! Calculate the error to the reference xx-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csol_xx
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_xxx
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_xxy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX%RvectorBlock(1),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivX)
        call ctab_addValue(rtable, "L2-error u_xx", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX%RvectorBlock(1),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivX)
        call ctab_addValue(rtable, "H1-error u_xx", derror)

        ! Calculate the error to the reference xy-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csol_xy
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_xxy !! = csol_xyx
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_xyy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockX%RvectorBlock(2),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivX)
        call ctab_addValue(rtable, "L2-error u_xy", derror)

        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY%RvectorBlock(1),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivY)
        call ctab_addValue(rtable, "L2-error u_yx", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockX%RvectorBlock(2),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivX)
        call ctab_addValue(rtable, "H1-error u_xy", derror)

        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY%RvectorBlock(1),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivY)
        call ctab_addValue(rtable, "H1-error u_yx", derror)

        ! Calculate the error to the reference yy-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csol_yy
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_xyy !! = csol_yyx
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csol_yyy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockY%RvectorBlock(2),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivX)
        call ctab_addValue(rtable, "L2-error u_yy", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockY%RvectorBlock(2),&
            getReferenceFunction_Poisson,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivY)
        call ctab_addValue(rtable, "H1-error u_yy", derror)
      end if

      ! Start UCD export to file
      if ((trim(adjustl(sucdFile)) .ne. '') .and.&
          (cucdType .ne. UCD_FORMAT_NONE)) then

        select case(cucdType)
        case(UCD_FORMAT_GMV)
          call ucd_startGMV(rexport,UCD_FLAG_STANDARD,&
              rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
              trim(sucdDir)//"/"//trim(sucdFile)//"_"//&
              trim(sys_siL(rproblem%ilvmax,2))//".gmv")
        case(UCD_FORMAT_BGMV)
          call ucd_startBGMV(rexport,UCD_FLAG_STANDARD,&
              rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
              trim(sucdDir)//"/"//trim(sucdFile)//"_"//&
              trim(sys_siL(rproblem%ilvmax,2))//".gmv")
        case(UCD_FORMAT_AVS)
          call ucd_startAVS(rexport,UCD_FLAG_STANDARD,&
              rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
              trim(sucdDir)//"/"//trim(sucdFile)//"_"//&
              trim(sys_siL(rproblem%ilvmax,2))//".avs")
        case(UCD_FORMAT_VTK)
          call ucd_startVTK(rexport,UCD_FLAG_STANDARD,&
              rproblem%RlevelInfo(rproblem%ilvmax)%rtriangulation,&
              trim(sucdDir)//"/"//trim(sucdFile)//"_"//&
              trim(sys_siL(rproblem%ilvmax,2))//".vtk")
        case default
          call output_line("Invalid type of UCD output file.", &
              OU_CLASS_ERROR,OU_MODE_STD,"sse_postprocessing_Poisson")
          call sys_halt()
        end select

        ! Add the solution and its (recovered) gradient to the UCD exporter
        call ucd_addVectorByVertex(rexport, "u", UCD_VAR_STANDARD, &
            rproblem%rvector%RvectorBlock(1))
        call ucd_addVectorFieldByVertex(rexport, "grad u", UCD_VAR_STANDARD, &
            (/p_rvectorDerivX,p_rvectorDerivY/))
        call ucd_addVectorFieldByVertex(rexport, "grad u_x", UCD_VAR_STANDARD, &
            (/rvectorBlockX%RvectorBlock(1),rvectorBlockX%RvectorBlock(2)/))
        call ucd_addVectorFieldByVertex(rexport, "grad u_y", UCD_VAR_STANDARD, &
            (/rvectorBlockY%RvectorBlock(1),rvectorBlockY%RvectorBlock(2)/))

        ! Write the file to disc, that is it.
        call ucd_write(rexport)
        call ucd_release(rexport)
      end if

      ! Clean temporal structures
      call lsysbl_releaseVector(rvectorBlock)
      call lsysbl_releaseVector(rvectorBlockX)
      call lsysbl_releaseVector(rvectorBlockY)
      call spdiscr_releaseBlockDiscr(rblockDiscr)

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_postprocessing_Poisson")
      call sys_halt()
    end select

  end subroutine sse_postprocessing_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine sse_outputTable_Poisson(rproblem,rtable)

!<description>
    ! Exports the convergence table to file.
!</description>

!<input>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(in) :: rproblem
!</input>

!<inputoutput>
    ! A convergence table.
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    select case(rproblem%cproblemType)
    case (POISSON_SCALAR,POISSON_SYSTEM)

      ! Compute reduction rates
      call ctab_evalConvergenceRate(rtable,"L2-error u",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error u_x",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error u_y",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error u_xx",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error u_xy",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error u_yx",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error u_yy",CTAB_REDUCTION_RATE)
      
      call ctab_evalConvergenceRate(rtable,"H1-error u",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error u_x",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error u_y",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error u_xx",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error u_xy",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error u_yx",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error u_yy",CTAB_REDUCTION_RATE)

      ! Adjust format of convergence table
      call ctab_setPrecision(rtable,"L2-error u",3)
      call ctab_setPrecision(rtable,"L2-error u-convrate",3)
      call ctab_setPrecision(rtable,"L2-error u_x",3)
      call ctab_setPrecision(rtable,"L2-error u_y",3)
      call ctab_setPrecision(rtable,"L2-error u_x-convrate",3)
      call ctab_setPrecision(rtable,"L2-error u_y-convrate",3)
      call ctab_setPrecision(rtable,"L2-error u_xx",3)
      call ctab_setPrecision(rtable,"L2-error u_xy",3)
      call ctab_setPrecision(rtable,"L2-error u_yx",3)
      call ctab_setPrecision(rtable,"L2-error u_yy",3)
      call ctab_setPrecision(rtable,"L2-error u_xx-convrate",3)
      call ctab_setPrecision(rtable,"L2-error u_xy-convrate",3)
      call ctab_setPrecision(rtable,"L2-error u_yx-convrate",3)
      call ctab_setPrecision(rtable,"L2-error u_yy-convrate",3)

      call ctab_setPrecision(rtable,"H1-error u",3)
      call ctab_setPrecision(rtable,"H1-error u-convrate",3)
      call ctab_setPrecision(rtable,"H1-error u_x",3)
      call ctab_setPrecision(rtable,"H1-error u_y",3)
      call ctab_setPrecision(rtable,"H1-error u_x-convrate",3)
      call ctab_setPrecision(rtable,"H1-error u_y-convrate",3)
      call ctab_setPrecision(rtable,"H1-error u_xx",3)
      call ctab_setPrecision(rtable,"H1-error u_xy",3)
      call ctab_setPrecision(rtable,"H1-error u_yx",3)
      call ctab_setPrecision(rtable,"H1-error u_yy",3)
      call ctab_setPrecision(rtable,"H1-error u_xx-convrate",3)
      call ctab_setPrecision(rtable,"H1-error u_xy-convrate",3)
      call ctab_setPrecision(rtable,"H1-error u_yx-convrate",3)
      call ctab_setPrecision(rtable,"H1-error u_yy-convrate",3)

      ! Set scientific flag
      call ctab_setScientific(rtable,"L2-error u",.true.)
      call ctab_setScientific(rtable,"L2-error u_x",.true.)
      call ctab_setScientific(rtable,"L2-error u_y",.true.)
      call ctab_setScientific(rtable,"L2-error u_xx",.true.)
      call ctab_setScientific(rtable,"L2-error u_xy",.true.)
      call ctab_setScientific(rtable,"L2-error u_yx",.true.)
      call ctab_setScientific(rtable,"L2-error u_yy",.true.)

      call ctab_setScientific(rtable,"H1-error u",.true.)
      call ctab_setScientific(rtable,"H1-error u_x",.true.)
      call ctab_setScientific(rtable,"H1-error u_y",.true.)
      call ctab_setScientific(rtable,"H1-error u_xx",.true.)
      call ctab_setScientific(rtable,"H1-error u_xy",.true.)
      call ctab_setScientific(rtable,"H1-error u_yx",.true.)
      call ctab_setScientific(rtable,"H1-error u_yy",.true.)

      ! Set Tex captions
      call ctab_setTexCaption(rtable,"cells","\# cells")
      call ctab_setTexCaption(rtable,"dofs","\# dofs")

      call ctab_setTexCaption(rtable,"L2-error u","$L^2(u)$")
      call ctab_setTexCaption(rtable,"H1-error u","$H^1(u)$")
      select case(rproblem%cproblemType)
      case (POISSON_SCALAR)
        call ctab_setTexCaption(rtable,"L2-error u_x","$L^2(\partial_{x}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"L2-error u_y","$L^2(\partial_{y}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"L2-error u_xx","$L^2(\partial_{xx}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"L2-error u_xy","$L^2(\partial_{xy}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"L2-error u_yx","$L^2(\partial_{yx}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"L2-error u_yy","$L^2(\partial_{yy}^{\rm ZZ}u)$")

        call ctab_setTexCaption(rtable,"H1-error u_x","$H^1(\partial_{x}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"H1-error u_y","$H^1(\partial_{y}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"H1-error u_xx","$H^1(\partial_{xx}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"H1-error u_xy","$H^1(\partial_{xy}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"H1-error u_yx","$H^1(\partial_{yx}^{\rm ZZ}u)$")
        call ctab_setTexCaption(rtable,"H1-error u_yy","$H^1(\partial_{yy}^{\rm ZZ}u)$")

      case (POISSON_SYSTEM)
        call ctab_setTexCaption(rtable,"L2-error u_x","$L^2(\sigma_x)$")
        call ctab_setTexCaption(rtable,"L2-error u_y","$L^2(\sigma_y)$")
        call ctab_setTexCaption(rtable,"L2-error u_xx","$L^2(\partial_{x}^{\rm ZZ}\sigma_x)$")
        call ctab_setTexCaption(rtable,"L2-error u_xy","$L^2(\partial_{y}^{\rm ZZ}\sigma_x)$")
        call ctab_setTexCaption(rtable,"L2-error u_yx","$L^2(\partial_{x}^{\rm ZZ}\sigma_y)$")
        call ctab_setTexCaption(rtable,"L2-error u_yy","$L^2(\partial_{y}^{\rm ZZ}\sigma_y)$")

        call ctab_setTexCaption(rtable,"H1-error u_x","$H^1(\sigma_x)$")
        call ctab_setTexCaption(rtable,"H1-error u_y","$H^1(\sigma_y)$")
        call ctab_setTexCaption(rtable,"H1-error u_xx","$H^1(\partial_{x}^{\rm ZZ}\sigma_x)$")
        call ctab_setTexCaption(rtable,"H1-error u_xy","$H^1(\partial_{y}^{\rm ZZ}\sigma_x)$")
        call ctab_setTexCaption(rtable,"H1-error u_yx","$H^1(\partial_{x}^{\rm ZZ}\sigma_y)$")
        call ctab_setTexCaption(rtable,"H1-error u_yy","$H^1(\partial_{y}^{\rm ZZ}\sigma_y)$")
      end select

      call ctab_setTexCaption(rtable,"L2-error u-convrate","")
      call ctab_setTexCaption(rtable,"L2-error u_x-convrate","")
      call ctab_setTexCaption(rtable,"L2-error u_y-convrate","")
      call ctab_setTexCaption(rtable,"L2-error u_xx-convrate","")
      call ctab_setTexCaption(rtable,"L2-error u_xy-convrate","")
      call ctab_setTexCaption(rtable,"L2-error u_yx-convrate","")
      call ctab_setTexCaption(rtable,"L2-error u_yy-convrate","")

      call ctab_setTexCaption(rtable,"H1-error u-convrate","")
      call ctab_setTexCaption(rtable,"H1-error u_x-convrate","")
      call ctab_setTexCaption(rtable,"H1-error u_y-convrate","")
      call ctab_setTexCaption(rtable,"H1-error u_xx-convrate","")
      call ctab_setTexCaption(rtable,"H1-error u_xy-convrate","")
      call ctab_setTexCaption(rtable,"H1-error u_yx-convrate","")
      call ctab_setTexCaption(rtable,"H1-error u_yy-convrate","")

      ! Set Tex format
      call ctab_setTexFormat(rtable,"cells","r")
      call ctab_setTexFormat(rtable,"dofs","r")

      ! Hide all H1-columns
      call ctab_setHidden(rtable,"H1-error u",.true.)
      call ctab_setHidden(rtable,"H1-error u-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error u_x",.true.)
      call ctab_setHidden(rtable,"H1-error u_y",.true.)
      call ctab_setHidden(rtable,"H1-error u_x-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error u_y-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error u_xx",.true.)
      call ctab_setHidden(rtable,"H1-error u_xy",.true.)
      call ctab_setHidden(rtable,"H1-error u_yx",.true.)
      call ctab_setHidden(rtable,"H1-error u_yy",.true.)
      call ctab_setHidden(rtable,"H1-error u_xx-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error u_xy-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error u_yx-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error u_yy-convrate",.true.)

      select case(rproblem%cproblemType)
      case (POISSON_SCALAR)
        call ctab_setTexTableCaption(rtable,&
            "$L^2$-Convergence table: Poisson problem solved as scalar equation.")

      case (POISSON_SYSTEM)
        call ctab_setTexTableCaption(rtable,&
            "$L^2$-Convergence table: Poisson problem solved in mixed formulation.")
      end select
      call ctab_setTexTableLabel(rtable,"tab:l2_convergence_rate")

      ! Write convergence table to Tex file
      call ctab_outputTex(rtable,'./table_l2.tex')

      select case(rproblem%cproblemType)
      case (POISSON_SCALAR)
        call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table: Poisson problem solved as scalar equation.")
      case (POISSON_SYSTEM)
        call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table: Poisson problem solved in mixed formulation.")
      end select
      call ctab_setTexTableLabel(rtable,"tab:h1_convergence_rate")

      ! Unhide all H1-columns
      call ctab_setHidden(rtable,"H1-error u",.false.)
      call ctab_setHidden(rtable,"H1-error u-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error u_x",.false.)
      call ctab_setHidden(rtable,"H1-error u_y",.false.)
      call ctab_setHidden(rtable,"H1-error u_x-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error u_y-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error u_xx",.false.)
      call ctab_setHidden(rtable,"H1-error u_xy",.false.)
      call ctab_setHidden(rtable,"H1-error u_yx",.false.)
      call ctab_setHidden(rtable,"H1-error u_yy",.false.)
      call ctab_setHidden(rtable,"H1-error u_xx-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error u_xy-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error u_yx-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error u_yy-convrate",.false.)

      ! Hide all L2-columns
      call ctab_setHidden(rtable,"L2-error u",.true.)
      call ctab_setHidden(rtable,"L2-error u-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error u_x",.true.)
      call ctab_setHidden(rtable,"L2-error u_y",.true.)
      call ctab_setHidden(rtable,"L2-error u_x-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error u_y-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error u_xx",.true.)
      call ctab_setHidden(rtable,"L2-error u_xy",.true.)
      call ctab_setHidden(rtable,"L2-error u_yx",.true.)
      call ctab_setHidden(rtable,"L2-error u_yy",.true.)
      call ctab_setHidden(rtable,"L2-error u_xx-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error u_xy-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error u_yx-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error u_yy-convrate",.true.)

      ! Write convergence table to Tex file
      call ctab_outputTex(rtable,'./table_h1.tex')

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_outputTable_Poisson")
      call sys_halt()
    end select

  end subroutine sse_outputTable_Poisson
      
  ! ***************************************************************************

!<subroutine>

  subroutine coeff_Matrix_Poisson(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,IdofsTrial,IdofsTest,&
      rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection

!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    integer :: cproblemType,iexpression,imatrixpos,i,j
    real(DP) :: dvalue

    ! Get configuration from quick-access arrays
    cproblemType = rcollection%IquickAccess(1)
    imatrixpos   = rcollection%IquickAccess(2)
    iexpression  = rcollection%IquickAccess(3)

    select case(cproblemType)
    case(POISSON_SCALAR) !------------------------------------------------------

      select case(imatrixpos)
      case(11)
        do i=1,npointsPerElement
          do j=1,nelements
            call fparser_evalFunction(rfparser, iexpression,&
                Dpoints(:,i,j), dvalue)
            Dcoefficients(1:2,i,j) = dvalue
          end do
        end do

      case default
        call output_line("Invalid matrix position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_Matrix_Poisson")
        call sys_halt()
      end select

    case (POISSON_SYSTEM) !-----------------------------------------------------

      select case(imatrixpos)
      case(12,13)
        do i=1,npointsPerElement
          do j=1,nelements
            call fparser_evalFunction(rfparser, iexpression,&
                Dpoints(:,i,j), dvalue)
            Dcoefficients(1:2,i,j) = -dvalue
          end do
        end do

      case (21,31,22,33)
        ! But this case should be handled by constant coefficients
        Dcoefficients(1,:,:) = 1.0_DP

      case default
        call output_line("Invalid matrix position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_Matrix_Poisson")
        call sys_halt()
      end select

    case default !--------------------------------------------------------------
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"coeff_Matrix_Poisson")
      call sys_halt()
    end select

  end subroutine coeff_Matrix_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Poisson(rdiscretisation,rform, &
      nelements,npointsPerElement,Dpoints, &
      IdofsTest,rdomainIntSubset,&
      Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation

!<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection

!</input>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    integer :: cproblemType,iexpression,ivectorpos,i,j

    ! Get configuration from quick-access arrays
    cproblemType = rcollection%IquickAccess(1)
    ivectorpos   = rcollection%IquickAccess(2)
    iexpression  = rcollection%IquickAccess(3)

    select case(cproblemType)
    case(POISSON_SCALAR) !------------------------------------------------------

      select case(ivectorpos)
      case(1)
        do i=1,npointsPerElement
          do j=1,nelements
            call fparser_evalFunction(rfparser, iexpression,&
                Dpoints(:,i,j), Dcoefficients(1,i,j))
          end do
        end do

      case default
        call output_line("Invalid vector position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Poisson")
        call sys_halt()
      end select

    case(POISSON_SYSTEM) !------------------------------------------------------

      select case(ivectorpos)
      case(1)
        do i=1,npointsPerElement
          do j=1,nelements
            call fparser_evalFunction(rfparser, iexpression,&
                Dpoints(:,i,j), Dcoefficients(1,i,j))
            Dcoefficients(1,i,j) = -Dcoefficients(1,i,j)
          end do
        end do

      case default
        call output_line("Invalid vector position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Poisson")
        call sys_halt()
      end select

    case default !--------------------------------------------------------------
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Poisson")
      call sys_halt()
    end select

  end subroutine coeff_RHS_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Bdr_Poisson(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates. According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), allocatable :: Dnx,Dny
    real(DP) :: dvalue
    integer :: cproblemType,iboundCompIdx,iexpression,ivectorpos,i,j

    ! Get configuration from quick-access arrays
    cproblemType    = rcollection%IquickAccess(1)
    ivectorpos      = rcollection%IquickAccess(2)
    iexpression     = rcollection%IquickAccess(3)
    iboundCompIdx   = rcollection%IquickAccess(4)

    select case(cproblemType)
    case(POISSON_SCALAR) !------------------------------------------------------

      select case(ivectorpos)
      case(1)
        do i=1,npointsPerElement
          do j=1,nelements
            call fparser_evalFunction(rfparser, iexpression,&
                Dpoints(:,i,j), Dcoefficients(1,i,j))
          end do
        end do

      case default
        call output_line("Invalid vector position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
        call sys_halt()
      end select

    case (POISSON_SYSTEM) !-----------------------------------------------------

      ! Calculate normal vectors
      allocate(Dnx(npointsPerElement,nelements),Dny(npointsPerElement,nelements))
      call boundary_getNormalVec2D(rdiscretisation%p_rboundary, iboundCompIdx,&
          DpointPar, Dnx, Dny, cparType=BDR_PAR_LENGTH)

      select case(ivectorpos)
      case(2)
        do i=1,npointsPerElement
          do j=1,nelements
            call fparser_evalFunction(rfparser, iexpression,&
                Dpoints(:,i,j), dvalue)
            Dcoefficients(1,i,j) = Dnx(i,j)*dvalue
          end do
        end do

      case(3)
        do i=1,npointsPerElement
          do j=1,nelements
            call fparser_evalFunction(rfparser, iexpression,&
                Dpoints(:,i,j), dvalue)
            Dcoefficients(1,i,j) = Dny(i,j)*dvalue
          end do
        end do

      case default
        call output_line("Invalid vector position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
        call sys_halt()
      end select

    case default !--------------------------------------------------------------
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
      call sys_halt()
    end select

  end subroutine coeff_RHS_Bdr_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_Poisson(Icomponents,rdiscretisation,rboundaryRegion,&
      ielement,cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)

!<description>
    ! This subroutine is called during the discretisation of boundary
    ! conditions. It calculates a special quantity on the boundary, which is
    ! then used by the discretisation routines to generate a discrete
    ! "snapshot" of the (actually analytic) boundary conditions.
!</description>

!<input>
    ! Component specifier.
    ! For Dirichlet boundary:
    !   Icomponents(1) defines the number of the boundary component, the value
    !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry,
    !   2=2nd solution component, e.g. Y-velocity,...)
    integer, dimension(:), intent(in) :: Icomponents

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in) :: rboundaryRegion

    ! The element number on the boundary which is currently being processed
    integer, intent(in) :: ielement

    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in) :: cinfoNeeded

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDDERIV :
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN :
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    integer, intent(in) :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV :
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN :
    !   dwhere = 0 (not used)
    real(DP), intent(in) :: dwhere

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection

!</input>

!<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1).
    ! If multiple values are needed, they are collected here (e.g. for
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(NDIM2D) :: Dpoint
    integer :: cproblemType,iboundCompIdx,iexpression,ivectorpos

    ! Get configuration from quick-access arrays
    cproblemType    = rcollection%IquickAccess(1)
    ivectorpos      = rcollection%IquickAccess(2)
    iexpression     = rcollection%IquickAccess(3)
    iboundCompIdx   = rcollection%IquickAccess(4)

    ! Get xy-coordinates
    call boundary_getCoords(rdiscretisation%p_rboundary, iboundCompIdx,&
        dwhere, Dpoint(1), Dpoint(2), BDR_PAR_LENGTH)

    ! Get Dirichlet values
    call fparser_evalFunction(rfparser, iexpression,&
        Dpoint, Dvalues(1))

  end subroutine getBoundaryValues_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
      Dvalues,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration

!<description>
    ! This subroutine is called during the calculation of errors. It has to compute
    ! the (analytical) values of a function in a couple of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
!</description>

!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection

!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: iexpression,i,j

    ! Get configuration from quick-access arrays
    iexpression  = rcollection%IquickAccess(cderivative)

    if (iexpression .eq. 0) then

      ! No expression available
      Dvalues = 0.0_DP

    else

      ! Evaluate expression
      do i=1,npointsPerElement
        do j=1,nelements
          call fparser_evalFunction(rfparser, iexpression,&
              Dpoints(:,i,j), Dvalues(i,j))
        end do
      end do

    end if

  end subroutine getReferenceFunction_Poisson

end module sse_callback_poisson
