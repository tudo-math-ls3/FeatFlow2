
!##############################################################################
!# ****************************************************************************
!# <name> sse_callback_sse </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the SSE problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in "intf_xxxx.inc" files.
!#
!#
!# 1.) sse_initMatVec_SSE
!#     -> Initialises the matrices/vectors for the SSE problem
!#
!# 2.) sse_postprocessing_SSE
!#     -> Post-processes the solution of the SSE problem
!#
!# 3.) sse_outputTable_SSE
!#      -> Output the convergence table of the SSE problem
!#
!# These callback routines are provided:
!#
!# 1.) coeff_MatrixSc_SSE   /
!#     coeff_MatrixSys1_SSE /
!#     coeff_MatrixSys2_SSE
!#     -> Returns analytic valyes for the system matrix A.
!#
!# 2.) coeff_MatrixSc_Bdr_SSE   /
!#     coeff_MatrixSys1_Bdr_SSE /
!#     coeff_MatrixSys2_Bdr_SSE
!#     -> Returns analytic valyes for the system matrix A.
!#
!# 5.) coeff_RHSSc_SSE   /
!#     coeff_RHSSys1_SSE /
!#     coeff_RHSSys2_SSE
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 6.) coeff_RHSSc_Bdr_SSE   /
!#     coeff_RHSSys1_Bdr_SSE /
!#     coeff_RHSSys2_Bdr_SSE /
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 7.) getBoundaryValues_SSE
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#
!# 8.) getReferenceFunction_SSE
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 9.) getAnalyticValues_SSE
!#     -> Returns the values of the analytic function and its derivatives.
!#
!# 10.) getAnalyticVelocities_SSE
!#      -> Returns the values of the analytic velocities.
!#
!# 11.) 
!#
!# </purpose>
!##############################################################################

module sse_callback_sse

  use basicgeometry
  use bilinearformevaluation
  use boundary
  use collection
  use convergencetable
  use cubature
  use derivatives
  use discretebc
  use discretefbc
  use domainintegration
  use element
  use fparser
  use fsystem
  use genoutput
  use linearformevaluation
  use lineariser
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
  use sse_base_sse

  implicit none

  private

  public :: sse_initMatVec_SSE
  public :: sse_postprocessing_SSE
  public :: sse_outputTable_SSE
  public :: sse_calcVelocity_SSE

  public :: getReferenceFunction_SSE
  public :: getAnalyticValues_SSE
  public :: getAnalyticVelocities_SSE
  public :: getBoundaryValues_SSE

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initMatVec_SSE(rproblem)

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
    integer :: i,iboundarySeg

    ! A boundary segment
    type(t_boundaryRegion) :: rboundaryRegion

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform

    select case(rproblem%cproblemType)
    case (SSE_SCALAR)
      !---------------------------------------------------------------------------
      !
      ! Problem formulation in (test,trial)-notation
      !
      !  /                                                                         \
      ! | (grad,-Re(A)*grad)                   (grad,Im(A)*grad)-\omega*(func,func) |
      ! |                                                                           |
      ! | (grad,-Im(A)*grad)+\omega(func,func) (grad,-Re(A)*grad)                   |
      !  \                                                                         /
      !
      !   /     \   / \
      !  | Re(N) | | 0 |
      ! *|       |=|   |
      !  | Im(N) | | 0 |
      !   \     /   \ /

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
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,2,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,1,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,2,LSYSSC_MATRIX9)
        
        ! (1,1)-block (grad w,-Re(A)*grad SSE_Re)
        rform%itermCount = 4
        rform%ballCoeffConstant = .false.
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_X
        rform%Idescriptors(1,3) = DER_DERIV_X
        rform%Idescriptors(2,3) = DER_DERIV_Y
        rform%Idescriptors(1,4) = DER_DERIV_Y
        rform%Idescriptors(2,4) = DER_DERIV_Y

        ! Assemble matrix block(1,1)
        rproblem%rcollection%IquickAccess(2) = 11
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixSc_SSE,rproblem%rcollection)

        ! Assemble matrix block(2,2)
        rproblem%rcollection%IquickAccess(2) = 22
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,2),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixSc_SSE,rproblem%rcollection)
        
        ! Anisotropic diffusion matrix (imaginary part)
        ! plus consistent mass matrix:
        ! (grad,Im(A)*grad)-\omega*(func,func)
        rform%itermCount = 5
        rform%ballCoeffConstant = .false.
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_X
        rform%Idescriptors(1,3) = DER_DERIV_X
        rform%Idescriptors(2,3) = DER_DERIV_Y
        rform%Idescriptors(1,4) = DER_DERIV_Y
        rform%Idescriptors(2,4) = DER_DERIV_Y
        rform%Idescriptors(1,5) = DER_FUNC
        rform%Idescriptors(2,5) = DER_FUNC
        
        ! Assemble matrix block(1,2)
        rproblem%rcollection%IquickAccess(2) = 12
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixSc_SSE,rproblem%rcollection)
        
        ! Assemble matrix block(2,1)
        rproblem%rcollection%IquickAccess(2) = 21
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixSc_SSE,rproblem%rcollection)
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

    case (SSE_SYSTEM1)
      !---------------------------------------------------------------------------
      print *, "SSE system formulation 1 is not implemented yet!"
      stop

    case (SSE_SYSTEM2)
      !---------------------------------------------------------------------------
      !
      ! Problem formulation in matrix notation
      !
      !  /                                                \   /       \   /     \
      ! |     0     -omega*M  ReC_1  -ImC_1  ReC_2  -ImC_2 | |   ReN   | |   0   |
      ! | omega*M        0    ImC_1   ReC_1  ImC_2   ReC_2 | |   ImN   | |   0   |
      ! | B_x+D_1        0      A       0      0       0   | | ReTau_1 | | Reb_1 |
      ! |     0      B_x+D_1    0       A      0       0   |*| ImTau_1 |=| Imb_1 |
      ! | B_y+D_2        0      0       0      A       0   | | ReTau_2 | | Reb_2 |
      ! |     0      B_y+D_2    0       0      0       A   | | ImTau_2 | | Imb_2 |
      !  \                                                /   \       /   \     /
      !
      ! The definition of the matrices and vectors is given in the documentation.
      !
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
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,5,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,1,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,5,1,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,3,LSYSSC_MATRIX9)
        call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,5,5,LSYSSC_MATRIX9)

        ! (1,2)-block -\omega*(w,Re(N))
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients(1)  = -dtidalfreq

        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            rcollection=rproblem%rcollection)

        ! (2,1)-block \omega*(w,Im(N))
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        call lsyssc_scaleMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),-1.0_DP)

        ! (3,3)-block (v_x,Re(tau_x)
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients(1)  = 1.0

        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,3),&
            rproblem%RlevelInfo(i)%RcubatureInfo(2),&
            rcollection=rproblem%rcollection)

        ! (4,4)-block (v_x,Im(tau_x)
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,3),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(4,4),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        ! (5,5)-block (v_y,Re(tau_y)
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients(1)  = 1.0

        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,5),&
            rproblem%RlevelInfo(i)%RcubatureInfo(3),&
            rcollection=rproblem%rcollection)

        ! (6,6)-block (v_y,Im(tau_y)
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,5),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(6,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        ! (1,3)-block (w,grad_x Re(a_11)*Re(tau_x) + grad_y Re(a_21)*Re(tau_y))
        rform%itermCount = 2
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,3),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixC1_SSE2re,rproblem%rcollection)

        ! (2,4)-block (w,grad_x Re(a_11)*Re(tau_x) + grad_y Re(a_21)*Re(tau_y))
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,3),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,4),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        ! (1,4)-block -(w,grad_x Im(a_11)*Im(tau_x) - grad_y Im(a_21)*Im(tau_y))
        rform%itermCount = 2
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_FUNC
        rform%Dcoefficients     = -1.0

        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,3),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,4),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,4),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixC1_SSE2im,rproblem%rcollection)

        ! (2,3)-block (w,grad_x Im(a_11)*Im(tau_x) + grad_y Im(a_21)*Im(tau_y))
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,4),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,3),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        call lsyssc_scaleMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,3),-1.0_DP)

        ! (1,5)-block (w,grad_x Re(a_12)*Re(tau_x) + grad_y Re(a_22)*Re(tau_y))
        rform%itermCount = 2
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,5),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixC2_SSE2re,rproblem%rcollection)

        ! (2,6)-block (w,grad_x Re(a_12)*Re(tau_x) + grad_y Re(a_22)*Re(tau_y))
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,5),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        ! (1,6)-block -(w,grad_x Im(a_12)*Im(tau_x) - grad_y Im(a_22)*Im(tau_y))
        rform%itermCount = 2
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_FUNC
        rform%Dcoefficients     = -1.0

        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,5),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,6),&
            rproblem%RlevelInfo(i)%RcubatureInfo(1),&
            coeff_MatrixC2_SSE2im,rproblem%rcollection)

        ! (2,5)-block (w,grad_x Im(a_12)*Im(tau_x) + grad_y Im(a_22)*Im(tau_y))
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,6),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,5),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        call lsyssc_scaleMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,5),-1.0_DP)

        ! (3,1)-block (grad_x q_x,Re(N)) + boundary-term
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(2),&
            rcollection=rproblem%rcollection)

        ! (4,2)-block (grad_x q_x,Im(N)) + boundary-term
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(4,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        ! (5,1)-block (grad_y q_y,Re(N)) + boundary-term
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_DERIV_Y
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalar(rform,.true.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,1),&
            rproblem%RlevelInfo(i)%RcubatureInfo(3),&
            rcollection=rproblem%rcollection)

        ! (6,2)-block (grad_y q_y,Im(N)) + boundary-term
        call lsyssc_duplicateMatrix(&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,1),&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(6,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

#if defined(CASE_SSE_ALEX) || defined(CASE_SSE_WINANT)

        ! The contribution of the mass matrix evaluated along the
        ! Neumann boundary needs to be included in the bilinear form.

        do iboundarySeg = 1,3

          ! Create boundary region
          call boundary_createRegion(rproblem%rboundary,1,iboundarySeg,rboundaryRegion)
          select case(iboundarySeg)
          case(1)
            rboundaryRegion%iproperties = BDR_PROP_WITHEND
          case(2)
            rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
          case(3)
            rboundaryRegion%iproperties = BDR_PROP_WITHSTART
          end select

          ! (3,1)-block boundary term: (q_x,Re(N)*n_x)
          ! (4,2)-block boundary term: (q_x,Im(N)*n_x) shares all data
          rform%itermCount = 1
          rform%ballCoeffConstant = .true.
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          rform%Dcoefficients     = 1.0

          call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
              rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
              coeff_MatrixD1_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)

          ! (5,1)-block boundary term: (q_y,Re(N)*n_y)
          ! (6,2)-block boundary term: (q_y,Im(N)*n_y) shares all data
          rform%itermCount = 1
          rform%ballCoeffConstant = .true.
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          rform%Dcoefficients     = 1.0

          call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
              rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,1),&
              coeff_MatrixD2_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)
        end do

#elif defined(CASE_SSE_MARCHI)

        ! There are no Neumann boundary conditions that need to be
        ! included in the bilinear form.

#elif defined(CASE_SSE_WALTERS)

        ! Create boundary region - part 1
        call boundary_createRegion(rproblem%rboundary,1,1,rboundaryRegion)
        rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND

        ! (3,1)-block boundary term: (q_x,Re(N)*n_x)
        ! (4,2)-block boundary term: (q_x,Im(N)*n_x) shares all data
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
            coeff_MatrixD1_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)

        ! (5,1)-block boundary term: (q_y,Re(N)*n_y)
        ! (6,2)-block boundary term: (q_y,Im(N)*n_y) shares all data
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,1),&
            coeff_MatrixD2_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)

        ! Create boundary region - part 2
        call boundary_createRegion(rproblem%rboundary,1,2,rboundaryRegion)
        rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND

        ! (3,1)-block boundary term: (q_x,Re(N)*n_x)
        ! (4,2)-block boundary term: (q_x,Im(N)*n_x) shares all data
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
            coeff_MatrixD1_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)

        ! (5,1)-block boundary term: (q_y,Re(N)*n_y)
        ! (6,2)-block boundary term: (q_y,Im(N)*n_y) shares all data
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,1),&
            coeff_MatrixD2_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)

        ! Create boundary region - part 4
        call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
        rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND

        ! (3,1)-block boundary term: (q_x,Re(N)*n_x)
        ! (4,2)-block boundary term: (q_x,Im(N)*n_x) shares all data
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
            coeff_MatrixD1_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)

        ! (5,1)-block boundary term: (q_y,Re(N)*n_y)
        ! (6,2)-block boundary term: (q_y,Im(N)*n_y) shares all data
        rform%itermCount = 1
        rform%ballCoeffConstant = .true.
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Dcoefficients     = 1.0

        call bilf_buildMatrixScalarBdr2D(rform,CUB_G5_1D,.false.,&
            rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(5,1),&
            coeff_MatrixD2_Bdr_SSE2,rboundaryRegion,rproblem%rcollection)

#else
#error 'Test case is undefined.' 
#endif
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

#if defined(CASE_SSE_ALEX) || defined(CASE_SSE_WINANT)

      ! The contribution of the mass matrix evaluated along the
      ! Neumann boundary needs to be included in the bilinear form.

      ! Create boundary region
      call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = 0_I32!BDR_PROP_WITHSTART+BDR_PROP_WITHEND

      ! Initialise the linear form along the boundary
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC

      ! Assemble the linear forms
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(3),coeff_RHSb1_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(4),coeff_RHSb1_Bdr_SSEim,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(5),coeff_RHSb2_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(6),coeff_RHSb2_Bdr_SSEim,rboundaryRegion)

#elif defined(CASE_SSE_MARCHI)

      do iboundarySeg = 1,4

        ! Create boundary region
        call boundary_createRegion(rproblem%rboundary,1,iboundarySeg,rboundaryRegion)

        ! Initialise the linear form along the boundary
        rlinform%itermCount = 1
        rlinform%Idescriptors(1) = DER_FUNC

        ! Assemble the linear forms
        call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
            rproblem%rrhs%RvectorBlock(3),coeff_RHSb1_Bdr_SSEre,rboundaryRegion)
        call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
            rproblem%rrhs%RvectorBlock(4),coeff_RHSb1_Bdr_SSEim,rboundaryRegion)
        call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
            rproblem%rrhs%RvectorBlock(5),coeff_RHSb2_Bdr_SSEre,rboundaryRegion)
        call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
            rproblem%rrhs%RvectorBlock(6),coeff_RHSb2_Bdr_SSEim,rboundaryRegion)

      end do

#elif defined(CASE_SSE_WALTERS)

      ! Initialise the linear form along the boundary
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC

      ! Create boundary region - part 1
      call boundary_createRegion(rproblem%rboundary,1,1,rboundaryRegion)

      ! Assemble the linear forms
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(3),coeff_RHSb1_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(4),coeff_RHSb1_Bdr_SSEim,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(5),coeff_RHSb2_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(6),coeff_RHSb2_Bdr_SSEim,rboundaryRegion)

      ! Create boundary region - part 2
      call boundary_createRegion(rproblem%rboundary,1,2,rboundaryRegion)

      ! Assemble the linear forms
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(3),coeff_RHSb1_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(4),coeff_RHSb1_Bdr_SSEim,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(5),coeff_RHSb2_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(6),coeff_RHSb2_Bdr_SSEim,rboundaryRegion)

      ! Create boundary region - part 4
      call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)

      ! Assemble the linear forms
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(3),coeff_RHSb1_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(4),coeff_RHSb1_Bdr_SSEim,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(5),coeff_RHSb2_Bdr_SSEre,rboundaryRegion)
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(6),coeff_RHSb2_Bdr_SSEim,rboundaryRegion)

#else
#error 'Test case is undefined.' 
#endif
    end select

  end subroutine sse_initMatVec_SSE

  ! ***************************************************************************

!<subroutine>

  subroutine sse_postprocessing_SSE(rproblem,cgradType,cgradSubtype,&
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
    type(t_vectorBlock), target :: rvectorBlockReal,rvectorBlockImag
    type(t_vectorBlock), target :: rvectorBlockXReal,rvectorBlockXImag
    type(t_vectorBlock), target :: rvectorBlockyReal,rvectorBlockYImag
    type(t_vectorScalar), pointer :: p_rvectorDerivXReal,p_rvectorDerivXImag
    type(t_vectorScalar), pointer :: p_rvectorDerivYReal,p_rvectorDerivYImag
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfoDerivXReal,p_rcubatureInfoDerivXImag
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfoDerivYReal,p_rcubatureInfoDerivYImag
    
    real(DP) :: derror
    integer :: i

    select case(rproblem%cproblemType)
    case (SSE_SCALAR,SSE_SYSTEM1,SSE_SYSTEM2)

      ! --- solution -------------------------------------------------------------

      if (present(rtable)) then
        ! Calculate the error to the real part of the reference function.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolReal
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolReal_x
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolReal_y

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rproblem%rvector%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,rcubatureInfo=&
            rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
        call ctab_addValue(rtable, "L2-error Re(SSE)", derror)

        ! H1-errro
        call pperr_scalar(PPERR_H1ERROR,derror,rproblem%rvector%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,rcubatureInfo=&
            rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1))
        call ctab_addValue(rtable, "H1-error Re(SSE)", derror)

        ! Calculate the error to the imaginary part of the reference function.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolImag
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolImag_x
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolImag_y

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rproblem%rvector%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,rcubatureInfo=&
            rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(2))
        call ctab_addValue(rtable, "L2-error Im(SSE)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rproblem%rvector%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,rcubatureInfo=&
            rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(2))
        call ctab_addValue(rtable, "H1-error Im(SSE)", derror)
      end if

      ! --- first derivative -----------------------------------------------------

      select case(rproblem%cproblemType)
      case (SSE_SCALAR)
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

        ! Create block vectors for the gradient (real and imaginary parts)
        call lsysbl_createVector(rblockDiscr, rvectorBlockReal, .false.)
        call lsysbl_createVector(rblockDiscr, rvectorBlockImag, .false.)

        ! Recover the real part of the gradient vector
        call ppgrd_calcGradient(rproblem%rvector%RvectorBlock(1),&
            rvectorBlockReal, cgradType, cgradSubType)
        
        ! Recover the imaginary part of the gradient vector
        call ppgrd_calcGradient(rproblem%rvector%RvectorBlock(2),&
            rvectorBlockImag, cgradType, cgradSubType)

        ! Set pointers to recovered gradient vector
        p_rvectorDerivXReal => rvectorBlockReal%RvectorBlock(1)
        p_rvectorDerivYReal => rvectorBlockReal%RvectorBlock(2)
        p_rvectorDerivXImag => rvectorBlockImag%RvectorBlock(1)
        p_rvectorDerivYImag => rvectorBlockImag%RvectorBlock(2)

        ! Set pointers to cubature structure from scalar fields
        p_rcubatureInfoDerivXReal => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1)
        p_rcubatureInfoDerivYReal => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1)
        p_rcubatureInfoDerivXImag => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(2)
        p_rcubatureInfoDerivYImag => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(2)

      case (SSE_SYSTEM1,SSE_SYSTEM2)
        ! Set pointer to scalar solution components
        p_rvectorDerivXReal => rproblem%rvector%RvectorBlock(3)
        p_rvectorDerivXImag => rproblem%rvector%RvectorBlock(4)
        p_rvectorDerivYReal => rproblem%rvector%RvectorBlock(5)
        p_rvectorDerivYImag => rproblem%rvector%RvectorBlock(6)

        ! Set pointers to cubature structure from scalar fields
        p_rcubatureInfoDerivXReal => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(3)
        p_rcubatureInfoDerivYReal => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(4)
        p_rcubatureInfoDerivXImag => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(5)
        p_rcubatureInfoDerivYImag => rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(6)
      end select

      if (present(rtable)) then
        ! Calculate the error to the real part of the reference x-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolReal_x
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolReal_xx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolReal_xy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivXReal,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXReal)
        call ctab_addValue(rtable, "L2-error Re(SSE_x)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivXReal,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXReal)
        call ctab_addValue(rtable, "H1-error Re(SSE_x)", derror)

        ! Calculate the error to the imaginary part of the reference x-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolImag_x
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolImag_xx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolImag_xy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivXImag,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXImag)
        call ctab_addValue(rtable, "L2-error Im(SSE_x)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivXImag,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXImag)
        call ctab_addValue(rtable, "H1-error Im(SSE_x)", derror)

        ! Calculate the error to the real part of the reference y-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolReal_y
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolReal_xy  ! = csolReal_yx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolReal_yy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivYReal,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYReal)
        call ctab_addValue(rtable, "L2-error Re(SSE_y)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivYReal,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYReal)
        call ctab_addValue(rtable, "H1-error Re(SSE_y)", derror)

        ! Calculate the error to the imaginary part of the reference x-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolImag_y
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolImag_xy  ! = csolImag_yx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolImag_yy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,p_rvectorDerivYImag,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYImag)
        call ctab_addValue(rtable, "L2-error Im(SSE_y)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,p_rvectorDerivYImag,&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYImag)
        call ctab_addValue(rtable, "H1-error Im(SSE_y)", derror)
      end if

      ! --- second derivative ----------------------------------------------------

      select case(rproblem%cproblemType)
      case (SSE_SYSTEM1,SSE_SYSTEM2)
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
      call lsysbl_createVector(rblockDiscr, rvectorBlockXReal, .false.)
      call lsysbl_createVector(rblockDiscr, rvectorBlockXImag, .false.)
      call lsysbl_createVector(rblockDiscr, rvectorBlockYReal, .false.)
      call lsysbl_createVector(rblockDiscr, rvectorBlockYImag, .false.)

      ! Recover the gradient vector of Re(SSE_x)
      call ppgrd_calcGradient(p_rvectorDerivXReal, rvectorBlockXReal,&
          cgradType, cgradSubType)

      ! Recover the gradient vector of Im(SSE_x)
      call ppgrd_calcGradient(p_rvectorDerivXImag, rvectorBlockXImag,&
          cgradType, cgradSubType)

      ! Recover the gradient vector of Re(SSE_y)
      call ppgrd_calcGradient(p_rvectorDerivYReal, rvectorBlockYReal,&
          cgradType, cgradSubType)

      ! Recover the gradient vector of Im(SSE_y)
      call ppgrd_calcGradient(p_rvectorDerivYImag, rvectorBlockYImag,&
          cgradType, cgradSubType)

      if (present(rtable)) then
        ! Calculate the error to the real part of the reference xx-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolReal_xx
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolReal_xxx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolReal_xxy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockXReal%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXReal)
        call ctab_addValue(rtable, "L2-error Re(SSE_xx)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockXReal%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXReal)
        call ctab_addValue(rtable, "H1-error Re(SSE_xx)", derror)

        ! Calculate the error to the imaginary part of the reference xx-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolImag_xx
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolImag_xxx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolImag_xxy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockXImag%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXImag)
        call ctab_addValue(rtable, "L2-error Im(SSE_xx)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockXImag%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXImag)
        call ctab_addValue(rtable, "H1-error Im(SSE_xx)", derror)

        ! Calculate the error to the real part of the reference xy-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolReal_xy
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolReal_xxy  ! = csolReal_xyx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolReal_xyy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockXReal%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXReal)
        call ctab_addValue(rtable, "L2-error Re(SSE_xy)", derror)

        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockYReal%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYReal)
        call ctab_addValue(rtable, "L2-error Re(SSE_yx)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockXReal%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXReal)        
        call ctab_addValue(rtable, "H1-error Re(SSE_xy)", derror)

        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockYReal%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYReal)
        call ctab_addValue(rtable, "H1-error Re(SSE_yx)", derror)

        ! Calculate the error to the imaginary part of the reference xy-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolImag_xy
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolImag_xxy  ! = csolImag_xyx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolImag_xyy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockXImag%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXImag)
        call ctab_addValue(rtable, "L2-error Im(SSE_xy)", derror)

        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockYImag%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYImag)
        call ctab_addValue(rtable, "L2-error Im(SSE_yx)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockXImag%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivXImag)
        call ctab_addValue(rtable, "H1-error Im(SSE_xy)", derror)

        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockYImag%RvectorBlock(1),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYImag)
        call ctab_addValue(rtable, "H1-error Im(SSE_yx)", derror)

        ! Calculate the error to the real part of the reference yy-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolReal_yy
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolReal_xyy  ! = csolReal_yyx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolReal_yyy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockYReal%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYReal)
        call ctab_addValue(rtable, "L2-error Re(SSE_yy)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockYReal%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYReal)
        call ctab_addValue(rtable, "H1-error Re(SSE_yy)", derror)

        ! Calculate the error to the imaginary part of the reference yy-derivative.
        rproblem%rcollection%IquickAccess(DER_FUNC)    = csolImag_yy
        rproblem%rcollection%IquickAccess(DER_DERIV_X) = csolImag_xyy  ! = csolImag_yyx
        rproblem%rcollection%IquickAccess(DER_DERIV_Y) = csolImag_yyy

        ! L2-error
        call pperr_scalar(PPERR_L2ERROR,derror,rvectorBlockYImag%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYImag)
        call ctab_addValue(rtable, "L2-error Im(SSE_yy)", derror)

        ! H1-error
        call pperr_scalar(PPERR_H1ERROR,derror,rvectorBlockYImag%RvectorBlock(2),&
            getReferenceFunction_SSE,rproblem%rcollection,&
            rcubatureInfo=p_rcubatureInfoDerivYImag)
        call ctab_addValue(rtable, "H1-error Im(SSE_yy)", derror)
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
        call ucd_addVectorByVertex(rexport, "Re(SSE)", UCD_VAR_STANDARD, &
            rproblem%rvector%RvectorBlock(1))
        call ucd_addVectorByVertex(rexport, "Im(SSE)", UCD_VAR_STANDARD, &
            rproblem%rvector%RvectorBlock(2))
        call ucd_addVectorFieldByVertex(rexport, "Re(grad SSE)", UCD_VAR_STANDARD, &
            (/p_rvectorDerivXReal,p_rvectorDerivYReal/))
        call ucd_addVectorFieldByVertex(rexport, "Im(grad SSE)", UCD_VAR_STANDARD, &
            (/p_rvectorDerivXImag,p_rvectorDerivYImag/))
        call ucd_addVectorFieldByVertex(rexport, "Re(grad SSE_x)", UCD_VAR_STANDARD, &
            (/rvectorBlockXReal%RvectorBlock(1),rvectorBlockXReal%RvectorBlock(2)/))
        call ucd_addVectorFieldByVertex(rexport, "Im(grad SSE_x)", UCD_VAR_STANDARD, &
            (/rvectorBlockXImag%RvectorBlock(1),rvectorBlockXImag%RvectorBlock(2)/))
        call ucd_addVectorFieldByVertex(rexport, "Re(grad SSE_y)", UCD_VAR_STANDARD, &
            (/rvectorBlockYReal%RvectorBlock(1),rvectorBlockYReal%RvectorBlock(2)/))
        call ucd_addVectorFieldByVertex(rexport, "Im(grad SSE_y)", UCD_VAR_STANDARD, &
            (/rvectorBlockYImag%RvectorBlock(1),rvectorBlockYImag%RvectorBlock(2)/))

        ! Write the file to disc, that is it.
        call ucd_write(rexport)
        call ucd_release(rexport)
      end if

      ! Clean temporal structures
      call lsysbl_releaseVector(rvectorBlockReal)
      call lsysbl_releaseVector(rvectorBlockImag)
      call lsysbl_releaseVector(rvectorBlockXReal)
      call lsysbl_releaseVector(rvectorBlockXImag)
      call lsysbl_releaseVector(rvectorBlockYReal)
      call lsysbl_releaseVector(rvectorBlockYImag)
      call spdiscr_releaseBlockDiscr(rblockDiscr)

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_postprocessing_sse")
      call sys_halt()
    end select

  end subroutine sse_postprocessing_SSE

  ! ***************************************************************************

!<subroutine>

  subroutine sse_outputTable_SSE(rproblem,rtable)

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
      case (SSE_SCALAR,SSE_SYSTEM1,SSE_SYSTEM2)

      ! Compute reduction rates
      call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_x)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_y)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_xx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_xy)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_yx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Re(SSE_yy)",CTAB_REDUCTION_RATE)

      call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_x)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_y)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_xx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_xy)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_yx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Re(SSE_yy)",CTAB_REDUCTION_RATE)

      call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_x)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_y)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_xx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_xy)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_yx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"L2-error Im(SSE_yy)",CTAB_REDUCTION_RATE)

      call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_x)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_y)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_xx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_xy)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_yx)",CTAB_REDUCTION_RATE)
      call ctab_evalConvergenceRate(rtable,"H1-error Im(SSE_yy)",CTAB_REDUCTION_RATE)

      ! Adjust format of convergence table
      call ctab_setPrecision(rtable,"L2-error Re(SSE)",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_x)",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_y)",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_x)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_y)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_xx)",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_xy)",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_yx)",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_yy)",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_xx)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_xy)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_yx)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Re(SSE_yy)-convrate",3)

      call ctab_setPrecision(rtable,"L2-error Im(SSE)",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_x)",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_y)",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_x)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_y)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_xx)",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_xy)",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_yx)",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_yy)",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_xx)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_xy)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_yx)-convrate",3)
      call ctab_setPrecision(rtable,"L2-error Im(SSE_yy)-convrate",3)

      call ctab_setPrecision(rtable,"H1-error Re(SSE)",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_x)",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_y)",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_x)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_y)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_xx)",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_xy)",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_yx)",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_yy)",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_xx)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_xy)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_yx)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Re(SSE_yy)-convrate",3)

      call ctab_setPrecision(rtable,"H1-error Im(SSE)",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_x)",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_y)",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_x)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_y)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_xx)",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_xy)",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_yx)",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_yy)",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_xx)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_xy)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_yx)-convrate",3)
      call ctab_setPrecision(rtable,"H1-error Im(SSE_yy)-convrate",3)

      ! Set scientific flag
      call ctab_setScientific(rtable,"L2-error Re(SSE)",.true.)
      call ctab_setScientific(rtable,"L2-error Re(SSE_x)",.true.)
      call ctab_setScientific(rtable,"L2-error Re(SSE_y)",.true.)
      call ctab_setScientific(rtable,"L2-error Re(SSE_xx)",.true.)
      call ctab_setScientific(rtable,"L2-error Re(SSE_xy)",.true.)
      call ctab_setScientific(rtable,"L2-error Re(SSE_yx)",.true.)
      call ctab_setScientific(rtable,"L2-error Re(SSE_yy)",.true.)

      call ctab_setScientific(rtable,"L2-error Im(SSE)",.true.)
      call ctab_setScientific(rtable,"L2-error Im(SSE_x)",.true.)
      call ctab_setScientific(rtable,"L2-error Im(SSE_y)",.true.)
      call ctab_setScientific(rtable,"L2-error Im(SSE_xx)",.true.)
      call ctab_setScientific(rtable,"L2-error Im(SSE_xy)",.true.)
      call ctab_setScientific(rtable,"L2-error Im(SSE_yx)",.true.)
      call ctab_setScientific(rtable,"L2-error Im(SSE_yy)",.true.)

      call ctab_setScientific(rtable,"H1-error Re(SSE)",.true.)
      call ctab_setScientific(rtable,"H1-error Re(SSE_x)",.true.)
      call ctab_setScientific(rtable,"H1-error Re(SSE_y)",.true.)
      call ctab_setScientific(rtable,"H1-error Re(SSE_xx)",.true.)
      call ctab_setScientific(rtable,"H1-error Re(SSE_xy)",.true.)
      call ctab_setScientific(rtable,"H1-error Re(SSE_yx)",.true.)
      call ctab_setScientific(rtable,"H1-error Re(SSE_yy)",.true.)

      call ctab_setScientific(rtable,"H1-error Im(SSE)",.true.)
      call ctab_setScientific(rtable,"H1-error Im(SSE_x)",.true.)
      call ctab_setScientific(rtable,"H1-error Im(SSE_y)",.true.)
      call ctab_setScientific(rtable,"H1-error Im(SSE_xx)",.true.)
      call ctab_setScientific(rtable,"H1-error Im(SSE_xy)",.true.)
      call ctab_setScientific(rtable,"H1-error Im(SSE_yx)",.true.)
      call ctab_setScientific(rtable,"H1-error Im(SSE_yy)",.true.)

      ! Set Tex captions
      call ctab_setTexCaption(rtable,"cells","\# cells")
      call ctab_setTexCaption(rtable,"dofs","\# dofs")

      call ctab_setTexCaption(rtable,"L2-error Re(SSE)","$L^2(Re(N))$")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE)","$L^2(Im(N))$")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE)","$H^1(Re(N))$")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE)","$H^1(Im(N))$")
      select case(rproblem%cproblemType)
      case (SSE_SCALAR)
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_x)","$L^2(\partial_{x}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_y)","$L^2(\partial_{y}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_xx)","$L^2(\partial_{xx}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_xy)","$L^2(\partial_{xy}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_yx)","$L^2(\partial_{yx}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_yy)","$L^2(\partial_{yy}^{\rm ZZ}Re(N))$")

        call ctab_setTexCaption(rtable,"L2-error Im(SSE_x)","$L^2(\partial_{x}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_y)","$L^2(\partial_{y}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_xx)","$L^2(\partial_{xx}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_xy)","$L^2(\partial_{xy}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_yx)","$L^2(\partial_{yx}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_yy)","$L^2(\partial_{yy}^{\rm ZZ}Im(N))$")

        call ctab_setTexCaption(rtable,"H1-error Re(SSE_x)","$H^1(\partial_{x}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_y)","$H^1(\partial_{y}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_xx)","$H^1(\partial_{xx}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_xy)","$H^1(\partial_{xy}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_yx)","$H^1(\partial_{yx}^{\rm ZZ}Re(N))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_yy)","$H^1(\partial_{yy}^{\rm ZZ}Re(N))$")

        call ctab_setTexCaption(rtable,"H1-error Im(SSE_x)","$H^1(\partial_{x}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_y)","$H^1(\partial_{y}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_xx)","$H^1(\partial_{xx}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_xy)","$H^1(\partial_{xy}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_yx)","$H^1(\partial_{yx}^{\rm ZZ}Im(N))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_yy)","$H^1(\partial_{yy}^{\rm ZZ}Im(N))$")

      case (SSE_SYSTEM1,SSE_SYSTEM2)
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_x)","$L^2(Re(\sigma_x))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_y)","$L^2(Re(\sigma_y))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_xx)","$L^2(\partial_{x}^{\rm ZZ}Re(\sigma_x))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_xy)","$L^2(\partial_{y}^{\rm ZZ}Re(\sigma_x))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_yx)","$L^2(\partial_{x}^{\rm ZZ}Re(\sigma_y))$")
        call ctab_setTexCaption(rtable,"L2-error Re(SSE_yy)","$L^2(\partial_{y}^{\rm ZZ}Re(\sigma_y))$")

        call ctab_setTexCaption(rtable,"L2-error Im(SSE_x)","$L^2(Im(\sigma_x))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_y)","$L^2(Im(\sigma_y))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_xx)","$L^2(\partial_{x}^{\rm ZZ}Im(\sigma_x))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_xy)","$L^2(\partial_{y}^{\rm ZZ}Im(\sigma_x))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_yx)","$L^2(\partial_{x}^{\rm ZZ}Im(\sigma_y))$")
        call ctab_setTexCaption(rtable,"L2-error Im(SSE_yy)","$L^2(\partial_{y}^{\rm ZZ}Im(\sigma_y))$")

        call ctab_setTexCaption(rtable,"H1-error Re(SSE_x)","$H^1(Re(\sigma_x))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_y)","$H^1(Re(\sigma_y))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_xx)","$H^1(\partial_{x}^{\rm ZZ}Re(\sigma_x))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_xy)","$H^1(\partial_{y}^{\rm ZZ}Re(\sigma_x))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_yx)","$H^1(\partial_{x}^{\rm ZZ}Re(\sigma_y))$")
        call ctab_setTexCaption(rtable,"H1-error Re(SSE_yy)","$H^1(\partial_{y}^{\rm ZZ}Re(\sigma_y))$")

        call ctab_setTexCaption(rtable,"H1-error Im(SSE_x)","$H^1(Im(\sigma_x))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_y)","$H^1(Im(\sigma_y))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_xx)","$H^1(\partial_{x}^{\rm ZZ}Im(\sigma_x))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_xy)","$H^1(\partial_{y}^{\rm ZZ}Im(\sigma_x))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_yx)","$H^1(\partial_{x}^{\rm ZZ}Im(\sigma_y))$")
        call ctab_setTexCaption(rtable,"H1-error Im(SSE_yy)","$H^1(\partial_{y}^{\rm ZZ}Im(\sigma_y))$")
      end select

      call ctab_setTexCaption(rtable,"L2-error Re(SSE)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_x)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_y)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_xx)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_xy)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_yx)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Re(SSE_yy)-convrate","")

      call ctab_setTexCaption(rtable,"L2-error Im(SSE)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_x)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_y)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_xx)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_xy)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_yx)-convrate","")
      call ctab_setTexCaption(rtable,"L2-error Im(SSE_yy)-convrate","")

      call ctab_setTexCaption(rtable,"H1-error Re(SSE)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_x)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_y)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_xx)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_xy)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_yx)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Re(SSE_yy)-convrate","")

      call ctab_setTexCaption(rtable,"H1-error Im(SSE)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_x)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_y)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_xx)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_xy)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_yx)-convrate","")
      call ctab_setTexCaption(rtable,"H1-error Im(SSE_yy)-convrate","")

      ! Set Tex format
      call ctab_setTexFormat(rtable,"cells","r")
      call ctab_setTexFormat(rtable,"dofs","r")

      ! Hide all H1-columns
      call ctab_setHidden(rtable,"H1-error Re(SSE)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_x)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_y)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_x)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_y)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xx)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xy)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yx)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yy)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xx)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xy)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yx)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yy)-convrate",.true.)

      call ctab_setHidden(rtable,"H1-error Im(SSE)",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_x)",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_y)",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_x)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_y)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xx)",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xy)",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yx)",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yy)",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xx)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xy)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yx)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yy)-convrate",.true.)

      ! Hide imaginary L2-columns
      call ctab_setHidden(rtable,"L2-error Im(SSE)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_x)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_y)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_x)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_y)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xx)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xy)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yx)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yy)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xy)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yy)-convrate",.true.)

      call ctab_setTexTableCaption(rtable,"$L^2$-Convergence table, real part")
      call ctab_setTexTableLabel(rtable,"tab:l2_convergence_rate_real")

      ! Write convergence table to Tex file
      call ctab_outputTex(rtable,'./table_l2.tex')

      ! Unhide imaginary parts of L2-columns
      call ctab_setHidden(rtable,"L2-error Im(SSE)",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE)-convrate",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_x)",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_y)",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_x)-convrate",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_y)-convrate",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xx)",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xy)",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yx)",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yy)",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xx)-convrate",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xy)-convrate",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yx)-convrate",.false.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yy)-convrate",.false.)

      ! Hide real parts of L2-columns
      call ctab_setHidden(rtable,"L2-error Re(SSE)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_x)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_y)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_x)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_y)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xx)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xy)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yx)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yy)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xy)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yy)-convrate",.true.)

      call ctab_setTexTableCaption(rtable,"$L^2$-Convergence table, imaginary part")
      call ctab_setTexTableLabel(rtable,"tab:l2_convergence_rate_aimag")

      ! Write convergence table to Tex file
      call ctab_outputTex(rtable,'./table_l2.tex',SYS_APPEND)

      ! Hide all L2-columns
      call ctab_setHidden(rtable,"L2-error Re(SSE)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_x)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_y)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_x)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_y)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xx)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xy)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yx)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yy)",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_xy)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Re(SSE_yy)-convrate",.true.)

      call ctab_setHidden(rtable,"L2-error Im(SSE)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_x)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_y)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_x)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_y)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xx)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xy)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yx)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yy)",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_xy)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yx)-convrate",.true.)
      call ctab_setHidden(rtable,"L2-error Im(SSE_yy)-convrate",.true.)

      ! Unhide real parts of H1-columns
      call ctab_setHidden(rtable,"H1-error Re(SSE)",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_x)",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_y)",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_x)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_y)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xx)",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xy)",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yx)",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yy)",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xx)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xy)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yx)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yy)-convrate",.false.)

      call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table, real part")
      call ctab_setTexTableLabel(rtable,"tab:h1_convergence_rate_real")

      ! Write convergence table to Tex file
      call ctab_outputTex(rtable,'./table_h1.tex')

      ! Hide real parts of H1-columns
      call ctab_setHidden(rtable,"H1-error Re(SSE)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_x)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_y)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_x)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_y)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xx)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xy)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yx)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yy)",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xx)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_xy)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yx)-convrate",.true.)
      call ctab_setHidden(rtable,"H1-error Re(SSE_yy)-convrate",.true.)

      ! Unhide imaginary parts of H1-columns
      call ctab_setHidden(rtable,"H1-error Im(SSE)",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_x)",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_y)",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_x)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_y)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xx)",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xy)",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yx)",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yy)",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xx)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_xy)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yx)-convrate",.false.)
      call ctab_setHidden(rtable,"H1-error Im(SSE_yy)-convrate",.false.)

      call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table, imaginary part")
      call ctab_setTexTableLabel(rtable,"tab:h1_convergence_rate_aimag")

      ! Write convergence table to Tex file
      call ctab_outputTex(rtable,'./table_h1.tex', SYS_APPEND)

    case default
      call output_line("Invalid type of problem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_outputTable_SSE")
      call sys_halt()
    end select
    
  end subroutine sse_outputTable_SSE
    
  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixSc_SSE(rdiscretisationTrial,rdiscretisationTest,&
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
    integer :: cproblemType,imatrixpos

    ! Get configuration from quick-access arrays
    cproblemType = rcollection%IquickAccess(1)
    imatrixpos   = rcollection%IquickAccess(2)

    select case(cproblemType)
    case(SSE_SCALAR)

      select case(imatrixpos)
        
      case(11,22)

      case(12)

      case(21)

      case default
        call output_line("Invalid matrix position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_MatrixSc_SSE")
        call sys_halt()
      end select
      
    case default
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"coeff_MatrixSc_SSE")
      call sys_halt()
    end select
    
  end subroutine coeff_MatrixSc_SSE
    
  ! ***************************************************************************
  
!<subroutine>

  subroutine coeff_MatrixA_SSEre(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) = -0.5_DP * real(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = -0.5_DP * real( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(3,ipoint,iel) = -0.5_DP * real(-cimg*(cCalpha1-cCalpha2)) ! C3
        Dcoefficients(4,ipoint,iel) = -0.5_DP * real(       cCalpha1+cCalpha2 ) ! C4
      end do
    end do

  end subroutine coeff_MatrixA_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_SSEim(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(3,ipoint,iel) = 0.5_DP * aimag(-cimg*(cCalpha1-cCalpha2)) ! C3
        Dcoefficients(4,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C4
        Dcoefficients(5,ipoint,iel) =        - dtidalfreq
      end do
    end do

  end subroutine coeff_MatrixA_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC1_SSE2re(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * real(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * real(-cimg*(cCalpha1-cCalpha2)) ! C3
      end do
    end do

  end subroutine coeff_MatrixC1_SSE2re

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC1_SSE2im(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag(-cimg*(cCalpha1-cCalpha2)) ! C3
      end do
    end do

  end subroutine coeff_MatrixC1_SSE2im

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC2_SSE2re(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * real( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(2,ipoint,iel) = 0.5_DP * real(       cCalpha1+cCalpha2 ) ! C4
      end do
    end do

  end subroutine coeff_MatrixC2_SSE2re

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC2_SSE2im(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C4
      end do
    end do

  end subroutine coeff_MatrixC2_SSE2im

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_Bdr_SSEre(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,ibct,DpointPar,IdofsTrial,&
      IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use triangulation
    use spatialdiscretisation, only: t_spatialDiscretisation

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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) =&
            0.5_DP * real(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,1) ! C1
        Dcoefficients(2,ipoint,iel) =&
            0.5_DP * real( cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,2) ! C2
        Dcoefficients(3,ipoint,iel) =&
            0.5_DP * real(-cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,1) ! C3
        Dcoefficients(4,ipoint,iel) =&
            0.5_DP * real(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,2) ! C4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_MatrixA_Bdr_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_Bdr_SSEim(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,ibct,DpointPar,IdofsTrial,&
      IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use triangulation
    use spatialdiscretisation, only: t_spatialDiscretisation

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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
            ds*sinh(-calpha1*dh)-&
            calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
            ds*sinh(-calpha2*dh)-&
            calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) =&
            -0.5_DP * aimag(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,1) ! C1
        Dcoefficients(2,ipoint,iel) =&
            -0.5_DP * aimag( cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,2) ! C2
        Dcoefficients(3,ipoint,iel) =&
            -0.5_DP * aimag(-cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,1) ! C3
        Dcoefficients(4,ipoint,iel) =&
            -0.5_DP * aimag(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,2) ! C4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_MatrixA_Bdr_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD1_Bdr_SSE2(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,ibct,DpointPar,IdofsTrial,&
      IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use triangulation
    use spatialdiscretisation, only: t_spatialDiscretisation

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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Store first normal component
        Dcoefficients(1,ipoint,iel) = Dnormal(ipoint,iel,1)
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_MatrixD1_Bdr_SSE2

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD2_Bdr_SSE2(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,ibct,DpointPar,IdofsTrial,&
      IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use triangulation
    use spatialdiscretisation, only: t_spatialDiscretisation

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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Store second normal component
        Dcoefficients(1,ipoint,iel) = Dnormal(ipoint,iel,2)
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_MatrixD2_Bdr_SSE2

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_SSEre(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
            ds*sinh(calpha1*dz)+&
            calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
            ds*sinh(calpha2*dz)+&
            calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) = -0.5_DP * real(       cDalpha1+cDalpha2 ) ! D1
        Dcoefficients(2,ipoint,iel) = -0.5_DP * real( cimg*(cDalpha1-cDalpha2)) ! D2
        Dcoefficients(3,ipoint,iel) = -0.5_DP * real(-cimg*(cDalpha1-cDalpha2)) ! D3
        Dcoefficients(4,ipoint,iel) = -0.5_DP * real(       cDalpha1+cDalpha2 ) ! D4
      end do
    end do

  end subroutine coeff_MatrixD_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_SSEim(rdiscretisationTrial,rdiscretisationTest,&
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
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
            ds*sinh(calpha1*dz)+&
            calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
            ds*sinh(calpha2*dz)+&
            calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag(       cDalpha1+cDalpha2 ) ! D1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag( cimg*(cDalpha1-cDalpha2)) ! D2
        Dcoefficients(3,ipoint,iel) = 0.5_DP * aimag(-cimg*(cDalpha1-cDalpha2)) ! D3
        Dcoefficients(4,ipoint,iel) = 0.5_DP * aimag(       cDalpha1+cDalpha2 ) ! D4
        Dcoefficients(5,ipoint,iel) =        - dtidalfreq
      end do
    end do

  end subroutine coeff_MatrixD_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_Bdr_SSEre(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,ibct,DpointPar,IdofsTrial,&
      IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use triangulation
    use spatialdiscretisation, only: t_spatialDiscretisation

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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
            ds*sinh(calpha1*dz)+&
            calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
            ds*sinh(calpha2*dz)+&
            calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) =&
            0.5_DP * real(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,1) ! D1
        Dcoefficients(2,ipoint,iel) =&
            0.5_DP * real( cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,2) ! D2
        Dcoefficients(3,ipoint,iel) =&
            0.5_DP * real(-cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,1) ! D3
        Dcoefficients(4,ipoint,iel) =&
            0.5_DP * real(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,2) ! D4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_MatrixD_Bdr_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_Bdr_SSEim(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,ibct,DpointPar,IdofsTrial,&
      IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use triangulation
    use spatialdiscretisation, only: t_spatialDiscretisation

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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
            ds*sinh(calpha1*dz)+&
            calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
            ds*sinh(calpha2*dz)+&
            calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) =&
            -0.5_DP * aimag(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,1) ! D1
        Dcoefficients(2,ipoint,iel) =&
            -0.5_DP * aimag( cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,2) ! D2
        Dcoefficients(3,ipoint,iel) =&
            -0.5_DP * aimag(-cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,1) ! D3
        Dcoefficients(4,ipoint,iel) =&
            -0.5_DP * aimag(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,2) ! D4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_MatrixD_Bdr_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_SSEre(rdiscretisation,rform, &
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine coeff_RHS_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_SSEim(rdiscretisation,rform, &
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine coeff_RHS_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Bdr_SSEre(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = dforcing

  end subroutine coeff_RHS_Bdr_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Bdr_SSEim(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine coeff_RHS_Bdr_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb1_Bdr_SSEre(rdiscretisation, rform,&
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
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        Dcoefficients (1,ipoint,iel) = dforcing * Dnormal(ipoint,iel,1)
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_RHSb1_Bdr_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb1_Bdr_SSEim(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine coeff_RHSb1_Bdr_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb2_Bdr_SSEre(rdiscretisation, rform,&
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
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        Dcoefficients (1,ipoint,iel) = dforcing * Dnormal(ipoint,iel,2)
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine coeff_RHSb2_Bdr_SSEre

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb2_Bdr_SSEim(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine coeff_RHSb2_Bdr_SSEim

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_SSE(Icomponents,rdiscretisation,rboundaryRegion,&
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
    
  end subroutine getBoundaryValues_SSE

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues3_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)

    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)

    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
        ds*sinh(-calpha1*dh)-&
        calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
        ds*sinh(-calpha2*dh)-&
        calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4

    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!real(cc2*dnx+cc4*dny)

  end subroutine getBoundaryValues3_SSE2

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues4_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)

    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)

    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
        ds*sinh(-calpha1*dh)-&
        calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
        ds*sinh(-calpha2*dh)-&
        calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4

    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!aimag(cc2*dnx+cc4*dny)

  end subroutine getBoundaryValues4_SSE2

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues5_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)

    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)

    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
        ds*sinh(-calpha1*dh)-&
        calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
        ds*sinh(-calpha2*dh)-&
        calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4

    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!real(-cc1*dnx-cc3*dny)

  end subroutine getBoundaryValues5_SSE2

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues6_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)

    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)

    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
        ds*sinh(-calpha1*dh)-&
        calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
        ds*sinh(-calpha2*dh)-&
        calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4

    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!aimag(-cc1*dnx-cc3*dny)

  end subroutine getBoundaryValues6_SSE2

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_SSE(cderivative,rdiscretisation, &
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
    
  end subroutine getReferenceFunction_SSE

  ! ***************************************************************************

!<subroutine>

  elemental subroutine getAnalyticValues_SSE(dx,dy,&
      cN,cN_x,cN_y,cN_xx,cN_xy,cN_yy)

!<description>
    ! This function computes the values of the analytic function and
    ! its derivatives.
!</description>

!<input>
    ! Coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<output>
    ! Solution value
    complex(DP), intent(out) :: cN

    ! Values of the first derivatives
    complex(DP), intent(out) :: cN_x,cN_y

    ! Values of the second derivatives
    complex(DP), intent(out) :: cN_xx,cN_xy,cN_yy
!</output>
!</subroutine>

#if defined(CASE_SSE_ALEX)
    ! local variables
    complex(DP) :: cC,calpha,cr1,cr2
    real(DP) :: dAv,dh,ds

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)

    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)

    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients
    calpha = sqrt(-cimg*dtidalfreq/dAv)
    cC     = dgravaccel/(dAv*(calpha**3))*&
        ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
    cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
    cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

    ! Solution values    
    cN = (cr1*exp(cr1*dlength+cr2*dx)-&
          cr2*exp(cr1*dx+cr2*dlength))/&
         (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength))

    ! Values of first derivatives
    cN_x = cr1*cr2*(cr2*exp(cr1*dlength+cr2*dx)-&
           cr1*exp(cr1*dx+cr2*dlength))/&
          (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength))
    cN_y = cmplx(0.0_DP,0.0_DP)

    ! Values of second derivatives
    cN_xx = (cr1**2)*(cr2**2)*(exp(cr1*dlength+cr2*dx)-&
        exp(cr1*dx+cr2*dlength))/&
        (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength))
    cN_xy = cmplx(0.0_DP,0.0_DP)
    cN_yy = cmplx(0.0_DP,0.0_DP)
#else
#error 'Test case is undefined.'
#endif

  end subroutine getAnalyticValues_SSE

  ! ***************************************************************************

!<subroutine>

  elemental subroutine getAnalyticVelocities_SSE(dx,dy,dz,cU,cV,cW)

!<description>
    ! This function computes the values of the analytic function and
    ! its derivatives.
!</description>

!<input>
    ! Coordinates
    real(DP), intent(in) :: dx,dy,dz
!</input>

!<output>
    ! Velocity values
    complex(DP), intent(out) :: cU,cV,cW
!</output>
!</subroutine>

#if defined(CASE_SSE_ALEX)
    ! local variables
    complex(DP) :: calpha,cbeta,calpha1,calpha2,cr1,cr2
    complex(DP) :: cN,cN_x,cN_y,cN_xx,cN_xy,cN_yy
    real(DP) :: dAv,dh,ds

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)

    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)

    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute sea-surface elevation and its derivatives
    call getAnalyticValues_SSE(dx,dy,cN,cN_x,cN_y,cN_xx,cN_xy,cN_yy)

    ! Compute coefficients
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel/dAv))
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel/dAv))

    cr1 = dgravaccel/(dAv*calpha1**2)*((ds*cosh(calpha1*dz))/&
        (calpha1*dAv*sinh(calpha1*dh) + ds*cosh(calpha1*dh))-1.0_DP)*(cN_x+cimg*cN_y)
    cr2 = dgravaccel/(dAv*calpha2**2)*((ds*cosh(calpha2*dz))/&
        (calpha2*dAv*sinh(calpha2*dh) + ds*cosh(calpha2*dh))-1.0_DP)*(cN_x-cimg*cN_y)

    ! Horizontal velocities
    cU = 0.5_DP*(cr1+cr2)
    cV = 0.5_DP*(cr1-cr2)/cimg

    ! Vertical velocity
    cbeta  = sqrt(cimg*dtidalfreq/dAv);
    calpha = ds/(dAv*cbeta*sinh(cbeta*dheight)+ds*cosh(cbeta*dheight))
    cW     = dgravaccel/(cimg*dtidalfreq)*(cN_x/dlengthB - cN_xx)*&
        (dz-calpha/cbeta*sinh(cbeta*dz)) + cimg*dtidalfreq*cN

#else
#error 'Test case is undefined.'
#endif

  end subroutine getAnalyticVelocities_SSE

  ! ***************************************************************************

!<subroutine>

  subroutine sse_calcVelocity_SSE(rvelocity,dz,rvector,rvectorGrad_SSEre,&
      rvectorGrad_SSEim,rvectorHessX_SSEre,rvectorHessX_SSEim,&
      rvectorHessY_SSEre,rvectorHessY_SSEim)

!<description>
    ! Calculates the vertical and horizontal velocities (u,v,w) from the
    ! first and (if present) second derivatives of the sea surface
    ! elevation N provided via rvector. The discretisation structure is
    ! provided via the problem structure rproblem.
!</description>

!<input>
    ! Sea surface elevation
    type(t_vectorBlock), intent(in) :: rvector

    ! Z-coordinate, where the velocity should be calculated
    real(DP), intent(in) :: dz

    ! Gradients of the sea surface elevation
    type(t_vectorBlock), intent(in) :: rvectorGrad_SSEre
    type(t_vectorBlock), intent(in) :: rvectorGrad_SSEim

    ! OPTIONAL: second derivatives of the sea surface elevation
    type(t_vectorBlock), intent(in), optional :: rvectorHessX_SSEre
    type(t_vectorBlock), intent(in), optional :: rvectorHessX_SSEim
    type(t_vectorBlock), intent(in), optional :: rvectorHessY_SSEre
    type(t_vectorBlock), intent(in), optional :: rvectorHessY_SSEim
!</input>

!<output>
    ! Velocity vector
    type(t_vectorBlock), intent(out) :: rvelocity
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_vectorBlock) :: rcoordsDOF
    real(DP), dimension(:), pointer :: p_DcoordsX,p_DcoordsY
    real(DP), dimension(:), pointer :: p_DsseX_SSEre,p_DsseX_SSEim
    real(DP), dimension(:), pointer :: p_DsseY_SSEre,p_DsseY_SSEim
    real(DP), dimension(:), pointer :: p_DvelU_SSEre,p_DvelU_SSEim
    real(DP), dimension(:), pointer :: p_DvelV_SSEre,p_DvelV_SSEim
    complex(DP) :: calpha1,calpha2,cr1,cr2,cSSEx,cSSEY,cvelU,cvelV
    real(DP) :: dAv,dh,ds
    integer, dimension(NDIM2D) :: Isize
    integer :: ipoint,i

    ! Compute coordinates of DOFs
    Isize = rvector%RvectorBlock(1)%NEQ
    call lsysbl_createVector(rcoordsDOF, Isize, .false.)
    call lin_calcDofCoords(rvector%RvectorBlock(1)%p_rspatialDiscr, rcoordsDOF)
    call lsyssc_getbase_double(rcoordsDOF%RvectorBlock(1), p_DcoordsX)
    call lsyssc_getbase_double(rcoordsDOF%RvectorBlock(2), p_DcoordsY)

    ! Create block discretisation structure
    call spdiscr_initBlockDiscr(rblockDiscr,6,&
        rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
        rvector%p_rblockDiscr%p_rboundary)
    do i=1,6
      call spdiscr_duplicateDiscrSc(rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(i), .true.)
    end do

    ! Set pointer to horizontal velocity values
    call lsysbl_createVector(rblockDiscr, rvelocity, .true.)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(1), p_DvelU_SSEre)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(2), p_DvelU_SSEim)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(3), p_DvelV_SSEre)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(4), p_DvelV_SSEim)

    ! Set pointers to gradient values
    call lsyssc_getbase_double(rvectorGrad_SSEre%RvectorBlock(1), p_DsseX_SSEre)
    call lsyssc_getbase_double(rvectorGrad_SSEre%RvectorBlock(2), p_DsseY_SSEre)
    call lsyssc_getbase_double(rvectorGrad_SSEim%RvectorBlock(1), p_DsseX_SSEim)
    call lsyssc_getbase_double(rvectorGrad_SSEim%RvectorBlock(2), p_DsseY_SSEim)

    ! Compute horizontal velocities (U,V) from analytical expressions
    do ipoint = 1, size(p_DcoordsX)

      ! Compute bottom profile
      dh = sse_bottomProfile(p_DcoordsX(ipoint),p_DcoordsY(ipoint))

      ! Compute bottom stress
      ds = sse_bottomStress(p_DcoordsX(ipoint),p_DcoordsY(ipoint))

      ! Compute vertical eddy viscosity
      dAv = sse_eddyViscosity(p_DcoordsX(ipoint),p_DcoordsY(ipoint))

      ! Compute coefficients calpha1 and calpha2
      calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
      calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

      ! Compute complex sea surface elevation
      cSSEx = cmplx(p_DsseX_SSEre(ipoint),p_DsseX_SSEim(ipoint))
      cSSEy = cmplx(p_DsseY_SSEre(ipoint),p_DsseY_SSEim(ipoint))

      ! Compute coefficient cr1
      cr1 = dgravaccel/(dAv*(calpha1**2))*&
          ((ds*cosh(calpha1*dz))/(calpha1*dAv*sinh(calpha1*dh)+&
          ds*cosh(calpha1*dh))-1.0_DP)*&
          (cSSEx+cimg*cSSEy)

      ! Compute coefficient cr2
      cr2 = dgravaccel/(dAv*(calpha2**2))*&
          ((ds*cosh(calpha2*dz))/(calpha2*dAv*sinh(calpha2*dh)+&
          ds*cosh(calpha2*dh))-1.0_DP)*&
          (cSSEx-cimg*cSSEy)

      ! Compute complex velocities
      cvelU = 0.5_DP*(cr1+cr2)
      cvelV = 0.5_DP*(cr1-cr2)/cimg

      ! And separate them into real and imaginary parts
      p_DvelU_SSEre(ipoint) = real(cvelU)
      p_DvelV_SSEre(ipoint) = real(cvelV)
      p_DvelU_SSEim(ipoint) = aimag(cvelU)
      p_DvelV_SSEim(ipoint) = aimag(cvelV)
    end do

    ! Release DOF coordinates
    call lsysbl_releaseVector(rcoordsDOF)

  end subroutine sse_calcVelocity_SSE
  
  ! ***************************************************************************

!<function>

  elemental function sse_bottomProfile(dx,dy) result(dh)

!<description>
    ! This function computes the bathymetry as a function of the
    ! Cartisian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Height of the bottom profile
    real(DP) :: dh
!</result>

!</function>

    select case(ibathymetryType)
    case(SSE_BOTTOMPROFILE_LINEAR)
      
      ! linear bottom profile
      dh = dheight0 + (dheight-dheight0)*(1-dx/dlength)
      
    case(SSE_BOTTOMPROFILE_CONSTANT)
      
      ! constant bottom profile
      dh = dheight
      
    case (SSE_BOTTOMPROFILE_PARABOLIC)
      
      ! parabolic bottom profile
      dh = dheight*(dheightRatio+(1.0_DP-dheightRatio)*(1.0-(dy/dwidth)**2))
      
    case default
      dh = 0.0_DP
      
    end select

  end function sse_bottomProfile

  ! ***************************************************************************

!<function>

  elemental function sse_bottomStress(dx,dy) result(ds)

!<description>
    ! This function computes the bottom stress as a function of the
    ! Cartesian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Bottom stress
    real(DP) :: ds
!</result>

!</function>

    ! local parameters
    real(DP), parameter :: dx1 = -15000.0_DP
    real(DP), parameter :: dx2 = -10000.0_DP
    real(DP), parameter :: ds1 =  0.1_DP
    real(DP), parameter :: ds2 =  0.00049_DP

    real(DP), parameter :: da = (ds1-ds2)/(dx1-dx2)
    real(DP), parameter :: db = (ds2*dx1-ds1*dx2)/(dx1-dx2)

    select case(istressType)

    case(SSE_STRESS_VARIABLE)
      ! variable stress
      if (dx .le. dx1) then
        ds = ds1
      elseif (dx .ge. dx2) then
        ds = ds2
      elseif ((dx .gt. dx1) .and. (dx .lt. dx2)) then
        ds = da*dx+db
      else
        ds = 0.0_DP
      end if

    case(SSE_STRESS_CONSTANT)
      ! constant stress
      ds = dstress

    case(SSE_STRESS_PROPORTIONAL)
      ! stress proportional to bathymetry
      ds = dstress * sse_bottomProfile(dx,dy) / dheight

    case default
      ds = 0.0_DP
    end select

  end function sse_bottomStress

  ! ***************************************************************************

!<function>

  elemental function sse_eddyViscosity(dx,dy) result(dAv)

!<description>
    ! This function computes the vertical eddy viscosity as a function
    ! of the Cartesian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Vertical eddy viscosity
    real(DP) :: dAv
!</result>

!</function>

    ! local parameters
    real(DP), parameter :: dx1 = -15000.0_DP
    real(DP), parameter :: dx2 = -10000.0_DP
    real(DP), parameter :: dAv1 = 1.0_DP
    real(DP), parameter :: dAv2 = 0.012_DP

    real(DP), parameter :: da = (dAv1-dAv2)/(dx1-dx2)
    real(DP), parameter :: db = (dAv2*dx1-dAv1*dx2)/(dx1-dx2)

    select case(iviscosityType)
    case(SSE_VISCOSITY_VARIABLE)
      ! variable eddy viscosity
      if (dx .le. dx1) then
        dAv = dAv1
      elseif (dx .ge. dx2) then
        dAv = dAv2
      elseif ((dx .gt. dx1) .and. (dx .lt. dx2)) then
        dAv = da*dx+db
      else
        dAv = 0.0_DP
      end if

    case(SSE_VISCOSITY_CONSTANT)
      ! constant eddy viscosity
      dAv = dviscosity

    case(SSE_VISCOSITY_PROPORTIONAL)
      ! eddy viscosity proportional to bathymetry
      dAv = dviscosity * sse_bottomProfile(dx,dy) / dheight

    case default
      dAv = 0.0_DP
    end select

  end function sse_eddyViscosity
    
end module sse_callback_sse
