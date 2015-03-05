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
!# 2.) sse_initDiscreteBC_Poisson
!#     -> Initialises the discrete version of the boundary conditions
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
!# 6.) getReferenceDerivX_Poisson
!#     -> Returns the values of the analytic derivatives in x-direction.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the recovered FE gradient in comparison to the
!#        analytic derivative
!#
!# 7.) getReferenceDerivY_Poisson
!#      -> Returns the values of the analytic derivatives in x-direction.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#         $H_1$-error of the recovered FE gradient in comparison to the
!#         analytic derivative
!#
!# 8.) getReferenceDerivXX_Poisson
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 9.) getReferenceDerivXY_Poisson
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#       $H_1$-error of the FE function in comparison to the analytic
!#       function
!#
!# 10.) getReferenceDerivYX_Poisson
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 11.) getReferenceDerivYY_Poisson
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 12.) getAnalyticValues_Poisson
!#      -> Returns the values of the analytic function and its derivatives.
!#
!# </purpose>
!##############################################################################

module sse_callback_poisson

  use bcassembly
  use bilinearformevaluation
  use blockmatassemblybase
  use blockmatassemblystdop
  use boundary
  use collection
  use cubature
  use derivatives
  use discretebc
  use discretefbc
  use domainintegration
  use element
  use feevaluation2
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
  
  use sse_base
  use sse_base_poisson

  implicit none

  private

  public :: sse_initMatVec_Poisson
  public :: sse_initDiscreteBC_Poisson
  
  public :: coeff_Matrix_Poisson
  public :: coeff_RHS_Poisson
  public :: coeff_RHS_Bdr_Poisson
  public :: getBoundaryValues_Poisson
  public :: getReferenceFunction_Poisson
  public :: getReferenceDerivX_Poisson
  public :: getReferenceDerivY_Poisson
  public :: getReferenceDerivXX_Poisson
  public :: getReferenceDerivXY_Poisson
  public :: getReferenceDerivYX_Poisson
  public :: getReferenceDerivYY_Poisson
  public :: getAnalyticValues_Poisson

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
  integer :: i,j,k

  ! A boundary segment
  type(t_boundaryRegion) :: rboundaryRegion

  ! A bilinear and linear form describing the analytic problem to solve
  type(t_bilinearForm) :: rform
  type(t_linearForm) :: rlinform

  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR)
    !---------------------------------------------------------------------------
    !
    ! Problem formulation in (test,trial)-notation
    !
    ! (grad_x,grad_x)*u + (grad_y,grad_y)*u = (func,f) + Neumann bc.s
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
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,1,1,LSYSSC_MATRIX9)

      ! Prepare quick-access arrays of collection structure
      rproblem%rcollection%IquickAccess(1) = rproblem%cproblemtype
      rproblem%rcollection%IquickAccess(2) = rproblem%cproblemsubtype
      
      ! (1,1)-block (grad_x w,grad_x u)+(grad_y w,grad_y w)
      rform%itermCount = 2
      rform%ballCoeffConstant = .false.
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Dcoefficients(1)  = 1.0
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_Y
      rform%Dcoefficients(2)  = 1.0

      ! Assemble matrix entry (1,1)
      rproblem%rcollection%IquickAccess(3) = 11
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
    rproblem%rcollection%IquickAccess(3) = 1
    call linf_buildVectorScalar(&
        rlinform,.true.,rproblem%rrhs%RvectorBlock(1),&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1),&
        coeff_RHS_Poisson,rproblem%rcollection)

    select case(rproblem%cproblemsubtype)
    case (POISSON_DIRICHLET)
      ! No natural boundary conditions in the vector
      
    case (POISSON_DIRICHLET_NEUMANN)
      ! Natural boundary conditions in the vector at regions 2
      call boundary_createRegion(rproblem%rboundary,1,2,rboundaryRegion)
      rproblem%rcollection%IquickAccess(4) = 1
      rproblem%rcollection%IquickAccess(5) = 2
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(1),coeff_RHS_Bdr_Poisson,&
          rboundaryRegion,rproblem%rcollection)
      
      ! Natural boundary conditions in the vector at regions 4
      call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
      rproblem%rcollection%IquickAccess(4) = 1
      rproblem%rcollection%IquickAccess(5) = 4
      call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
          rproblem%rrhs%RvectorBlock(1),coeff_RHS_Bdr_Poisson,&
          rboundaryRegion,rproblem%rcollection)
      
    case (POISSON_NEUMANN)
      ! Natural boundary conditions in the vector at regions 1-4
      do j=1,4
        call boundary_createRegion(rproblem%rboundary,1,j,rboundaryRegion)
        rproblem%rcollection%IquickAccess(4) = 1
        rproblem%rcollection%IquickAccess(5) = j
        call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
            rproblem%rrhs%RvectorBlock(1),coeff_RHS_Bdr_Poisson,&
            rboundaryRegion,rproblem%rcollection)
      end do
      
    case default
      call output_line("Invalid type of subproblem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initMatVec_Poisson")
      call sys_halt()
    end select

  case (POISSON_SYSTEM)
    !---------------------------------------------------------------------------
    !
    ! Problem formulation in (test,trial)-notation
    !
    !  /                                         \   /       \   /         \
    ! |      0        (func,grad_x) (func,grad_y) | |    u    | | (func,-f) |
    ! |                                           | |         | |           |
    ! | (grad_x,func)  (func,func)       0        |*| sigma_x |=|     0     |
    ! |                                           | |         | |           |
    ! | (grad_y,func)       0        (func,func)  | | sigma_y | |     0     |
    !  \                                         /   \       /   \         /
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
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,1,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,2,2,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,1,LSYSSC_MATRIX9)
      call bilf_createMatrixStructure(rproblem%RlevelInfo(i)%rmatrix,3,3,LSYSSC_MATRIX9)
      
      ! Prepare quick-access arrays of collection structure
      rproblem%rcollection%IquickAccess(1) = rproblem%cproblemtype
      rproblem%rcollection%IquickAccess(2) = rproblem%cproblemsubtype
      
      ! (1,2)-block (w,grad_x sigma_x)
      rform%itermCount = 1
      rform%ballCoeffConstant = .false.
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0

      rproblem%rcollection%IquickAccess(3) = 12
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          coeff_Matrix_Poisson,rproblem%rcollection)

      ! (1,3)-block (w,grad_y sigma_y)
      rform%itermCount = 1
      rform%ballCoeffConstant = .false.
      rform%Idescriptors(1,1) = DER_DERIV_Y
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0

      rproblem%rcollection%IquickAccess(3) = 13
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,3),&
          rproblem%RlevelInfo(i)%RcubatureInfo(1),&
          coeff_Matrix_Poisson,rproblem%rcollection)

      ! (2,1)-block (grad_x v_x,u)
      rform%itermCount = 1
      rform%ballCoeffConstant = .false.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Dcoefficients(1)  = 1.0

      rproblem%rcollection%IquickAccess(3) = 21
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(2),&
          coeff_Matrix_Poisson,rproblem%rcollection)

      ! (3,1)-block (grad_y v_y,u)
      rform%itermCount = 1
      rform%ballCoeffConstant = .false.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y
      rform%Dcoefficients(1)  = 1.0

      rproblem%rcollection%IquickAccess(3) = 31
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(3,1),&
          rproblem%RlevelInfo(i)%RcubatureInfo(3),&
          coeff_Matrix_Poisson,rproblem%rcollection)

      ! (2,2)-block (v_x,sigma_x)
      rform%itermCount = 1
      rform%ballCoeffConstant = .false.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0

      rproblem%rcollection%IquickAccess(3) = 22
      call bilf_buildMatrixScalar(rform,.true.,&
          rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(2,2),&
          rproblem%RlevelInfo(i)%RcubatureInfo(2),&
          coeff_Matrix_Poisson,rproblem%rcollection)

      ! (3,3)-block (v_y,sigma_y)
      rform%itermCount = 1
      rform%ballCoeffConstant = .false.
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      rform%Dcoefficients(1)  = 1.0
      
      rproblem%rcollection%IquickAccess(3) = 22
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
    rproblem%rcollection%IquickAccess(3) = 1
    call linf_buildVectorScalar(&
        rlinform,.true.,rproblem%rrhs%RvectorBlock(1),&
        rproblem%RlevelInfo(rproblem%ilvmax)%RcubatureInfo(1),&
        coeff_RHS_Poisson,rproblem%rcollection)

    select case(rproblem%cproblemsubtype)
    case (POISSON_DIRICHLET)
      ! Essential boundary conditions in the vector at regions 1-4
      do j=1,4
        call boundary_createRegion(rproblem%rboundary,1,j,rboundaryRegion)
        do k=2,3   ! second and third equation
          rproblem%rcollection%IquickAccess(3) = k
          rproblem%rcollection%IquickAccess(4) = 1
          rproblem%rcollection%IquickAccess(5) = j
          call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
              rproblem%rrhs%RvectorBlock(k),coeff_RHS_Bdr_Poisson,&
              rboundaryRegion,rproblem%rcollection)
        end do
      end do
      
    case (POISSON_DIRICHLET_NEUMANN)
      ! Essential boundary conditions in the vector at regions 1
      call boundary_createRegion(rproblem%rboundary,1,1,rboundaryRegion)
      do k=2,3   ! second and third equation
        rproblem%rcollection%IquickAccess(3) = k
        rproblem%rcollection%IquickAccess(4) = 1
        rproblem%rcollection%IquickAccess(5) = 1
        call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
            rproblem%rrhs%RvectorBlock(k),coeff_RHS_Bdr_Poisson,&
            rboundaryRegion,rproblem%rcollection)
      end do
      
      ! Essential boundary conditions in the vector at regions 3
      call boundary_createRegion(rproblem%rboundary,1,3,rboundaryRegion)
      do k=2,3   ! second and third equation
        rproblem%rcollection%IquickAccess(3) = k
        rproblem%rcollection%IquickAccess(4) = 1
        rproblem%rcollection%IquickAccess(5) = 3
        call linf_buildVectorScalarBdr2d(rlinform,CUB_G5_1D,.false.,&
            rproblem%rrhs%RvectorBlock(k),coeff_RHS_Bdr_Poisson,&
            rboundaryRegion,rproblem%rcollection)
      end do
      
    case (POISSON_NEUMANN)
      ! No essential boundary conditions in the vector
      
    case default
      call output_line("Invalid type of subproblem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initMatVec_Poisson")
      call sys_halt()
    end select

  case default
    call output_line("Invalid type of problem.", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_initMatVec_Poisson")
    call sys_halt()
  end select
    
  end subroutine sse_initMatVec_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initDiscreteBC_Poisson(rproblem)

!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_boundaryRegion) :: rboundaryRegion
  integer :: i,j,k

  ! Prepare quick-access arrays of collection structure
  rproblem%rcollection%IquickAccess(1) = rproblem%cproblemtype
  rproblem%rcollection%IquickAccess(2) = rproblem%cproblemsubtype

  ! Create a t_discreteBC structure where we store all discretised
  ! boundary conditions.
  do i=rproblem%ilvmin,rproblem%ilvmax
    call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
  end do
  
  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR) !-------------------------------------------------------
    
    ! Essential boundary conditions go into the first vector
    rproblem%rcollection%IquickAccess(3) = 1
    
    select case(rproblem%cproblemsubtype)
    case (POISSON_DIRICHLET)
      ! Essential boundary conditions at regions 1-4
      do j=1,4
        rproblem%rcollection%IquickAccess(4) = 1
        rproblem%rcollection%IquickAccess(5) = j
        ! Create boundary region and discretize boundary conditions
        call boundary_createRegion(rproblem%rboundary,1,j,rboundaryRegion)
        do i=rproblem%ilvmin,rproblem%ilvmax
          call bcasm_newDirichletBConRealBD(&
              rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
              rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
              getBoundaryValues_Poisson,rproblem%rcollection)
        end do
      end do
      
    case (POISSON_DIRICHLET_NEUMANN)
      ! Essential boundary conditions at regions 1
      rproblem%rcollection%IquickAccess(4) = 1
      rproblem%rcollection%IquickAccess(5) = 1
      ! Create boundary region and discretize boundary conditions
      call boundary_createRegion(rproblem%rboundary,1,1,rboundaryRegion)
      do i=rproblem%ilvmin,rproblem%ilvmax
        call bcasm_newDirichletBConRealBD(&
            rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
            rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
            getBoundaryValues_Poisson,rproblem%rcollection)
      end do

      ! Essential boundary conditions at regions 3
      rproblem%rcollection%IquickAccess(4) = 1
      rproblem%rcollection%IquickAccess(5) = 3
      ! Create boundary region and discretize boundary conditions
      call boundary_createRegion(rproblem%rboundary,1,3,rboundaryRegion)
      do i=rproblem%ilvmin,rproblem%ilvmax
        call bcasm_newDirichletBConRealBD(&
            rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,1,&
            rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
            getBoundaryValues_Poisson,rproblem%rcollection)
      end do
      
    case (POISSON_NEUMANN)
      ! No essential boundary conditions
      
    case default
      call output_line("Invalid type of subproblem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscreteBC_Poisson")
      call sys_halt()
    end select
    
  case (POISSON_SYSTEM) !-------------------------------------------------------

    select case(rproblem%cproblemsubtype)
    case (POISSON_DIRICHLET)
      ! No essential boundary conditions
      
    case (POISSON_DIRICHLET_NEUMANN)

      ! Essential boundary conditions at regions 2
      rproblem%rcollection%IquickAccess(4) = 1
      rproblem%rcollection%IquickAccess(5) = 2
      ! Create boundary region and discretize boundary conditions
      call boundary_createRegion(rproblem%rboundary,1,2,rboundaryRegion)
      do i=rproblem%ilvmin,rproblem%ilvmax
        do k=2,3   ! second and third equation
          rproblem%rcollection%IquickAccess(3) = k
          call bcasm_newDirichletBConRealBD(&
              rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,k,&
              rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
              getBoundaryValues_Poisson,rproblem%rcollection)
        end do
      end do

      ! Essential boundary conditions at regions 4
      rproblem%rcollection%IquickAccess(4) = 1
      rproblem%rcollection%IquickAccess(5) = 4
      ! Create boundary region and discretize boundary conditions
      do i=rproblem%ilvmin,rproblem%ilvmax
        call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
        do k=2,3   ! second and third equation
          rproblem%rcollection%IquickAccess(3) = k
          call bcasm_newDirichletBConRealBD(&
              rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,k,&
              rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
              getBoundaryValues_Poisson,rproblem%rcollection)
        end do
      end do
      
    case (POISSON_NEUMANN)
      ! Essential boundary conditions at regions 1-4
      do j=1,4
        rproblem%rcollection%IquickAccess(4) = 1
        rproblem%rcollection%IquickAccess(5) = j
        ! Create boundary region and discretize boundary conditions
        call boundary_createRegion(rproblem%rboundary,1,j,rboundaryRegion)
        do i=rproblem%ilvmin,rproblem%ilvmax
          do k=2,3   ! second and third equation
            rproblem%rcollection%IquickAccess(3) = k
            call bcasm_newDirichletBConRealBD(&
                rproblem%RlevelInfo(i)%rmatrix%p_rblockDiscrTest,k,&
                rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
                getBoundaryValues_Poisson,rproblem%rcollection)
          end do
        end do
      end do
      
    case default
      call output_line("Invalid type of subproblem.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscreteBC_Poisson")
      call sys_halt()
    end select
    
  case default
    call output_line("Invalid type of problem.", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_initDiscreteBC_Poisson")
    call sys_halt()
  end select 
  
  end subroutine sse_initDiscreteBC_Poisson
  
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
    integer :: cproblemtype,cproblemsubtype,imatrixpos

    ! Get configuration from quick-access arrays
    cproblemtype    = rcollection%IquickAccess(1)
    cproblemsubtype = rcollection%IquickAccess(2)
    imatrixpos      = rcollection%IquickAccess(3)

    select case(cproblemtype)
    case(POISSON_SCALAR) !------------------------------------------------------
      ! The problem subtype has no influence on the matrix
      
      select case(imatrixpos)
      case(11)
        Dcoefficients(1:2,:,:) = dpoisson
        
      case default
        call output_line("Invalid matrix position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_Matrix_Poisson")
        call sys_halt()
      end select

    case (POISSON_SYSTEM) !-----------------------------------------------------
      ! The problem subtype has no influence on the matrix
      
      select case(imatrixpos)
      case(12,13)
        Dcoefficients(1,:,:) = dpoisson

      case (21,31,22,33)
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
    integer :: cproblemtype,cproblemsubtype,ivectorpos
    
    ! Get configuration from quick-access arrays
    cproblemtype    = rcollection%IquickAccess(1)
    cproblemsubtype = rcollection%IquickAccess(2)
    ivectorpos      = rcollection%IquickAccess(3)

    select case(cproblemtype)
    case(POISSON_SCALAR) !------------------------------------------------------
      ! The problem subtype has no influence on the matrix
      
      select case(ivectorpos)
      case(1)
        !    u(x,y) = SIN(PI * x) * SIN(PI * y)
        ! => f(x,y) = 2 * PI^2 * SIN(PI * x) * SIN(PI * y)
        Dcoefficients (1,:,:) = 2.0_DP * SYS_PI * SYS_PI &
                              * sin(SYS_PI * Dpoints(1,:,:)) &
                              * sin(SYS_PI * Dpoints(2,:,:))
        
      case default
        call output_line("Invalid vector position.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Poisson")
        call sys_halt()
      end select
      
    case (POISSON_SYSTEM) !-----------------------------------------------------
      ! The problem subtype has no influence on the matrix
      
      select case(ivectorpos)
      case(1)
        !     u(x,y) = SIN(PI * x) * SIN(PI * y)
        ! => -f(x,y) = -2 * PI^2 * SIN(PI * x) * SIN(PI * y)
        Dcoefficients (1,:,:) = - 2.0_DP * SYS_PI * SYS_PI &
                              * sin(SYS_PI * Dpoints(1,:,:)) &
                              * sin(SYS_PI * Dpoints(2,:,:))

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
    integer :: cproblemtype,cproblemsubtype,ivectorpos,icomp,isegment
    
    ! Get configuration from quick-access arrays
    cproblemtype    = rcollection%IquickAccess(1)
    cproblemsubtype = rcollection%IquickAccess(2)
    ivectorpos      = rcollection%IquickAccess(3)
    icomp           = rcollection%IquickAccess(4)
    isegment        = rcollection%IquickAccess(5)

    select case(cproblemtype)
    case(POISSON_SCALAR) !------------------------------------------------------
      
      select case(cproblemsubtype)
      case(POISSON_DIRICHLET)
        call output_line("Invalid boundary segment.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
        call sys_halt()
        
      case(POISSON_DIRICHLET_NEUMANN)
        if (ivectorpos .eq. 1) then
          if ((icomp .eq. 1)      .and.&
              (isegment .eq. 2 .or. isegment .eq. 4)) then
            Dcoefficients(1,:,:) = dneumann
          else
            call output_line("Invalid boundary segment.", &
                OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
            call sys_halt()
          end if
        else
          call output_line("Invalid vector position.", &
              OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
          call sys_halt()
        end if
        
      case(POISSON_NEUMANN)
        if (ivectorpos .eq. 1) then
          Dcoefficients(1,:,:) = dneumann
        else
          call output_line("Invalid vector position.", &
              OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
          call sys_halt()
        end if
        
      case default
        call output_line("Invalid problem subtype.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
        call sys_halt()
      end select
      
    case (POISSON_SYSTEM) !-----------------------------------------------------
      
      select case(cproblemsubtype)
      case(POISSON_DIRICHLET)
        ! Calculate normal vectors
        allocate(Dnx(npointsPerElement,nelements),Dny(npointsPerElement,nelements))
        call boundary_getNormalVec2D(rdiscretisation%p_rboundary, icomp,&
            DpointPar, Dnx, Dny, cparType=BDR_PAR_LENGTH)
        
        if    (ivectorpos .eq. 2) then
          Dcoefficients(1,:,:) = dneumann*Dnx(:,:)
        elseif(ivectorpos .eq. 3) then
          Dcoefficients(1,:,:) = dneumann*Dny(:,:)
        else
          call output_line("Invalid vector position.", &
              OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
          call sys_halt()
        end if
        deallocate(Dnx,Dny)
        
      case(POISSON_DIRICHLET_NEUMANN)
        if ((icomp .eq. 1) .and.&
            (isegment .eq. 1 .or. isegment .eq. 3)) then

          ! Calculate normal vectors
          allocate(Dnx(npointsPerElement,nelements),Dny(npointsPerElement,nelements))
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary, icomp,&
              DpointPar, Dnx, Dny, cparType=BDR_PAR_LENGTH)

          if    (ivectorpos .eq. 2) then
            Dcoefficients(1,:,:) = dneumann*Dnx(:,:)
          elseif(ivectorpos .eq. 3) then
            Dcoefficients(1,:,:) = dneumann*Dny(:,:)
          else
            call output_line("Invalid vector position.", &
                OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
            call sys_halt()
          end if
          deallocate(Dnx,Dny)
          
        else
          call output_line("Invalid boundary segment.", &
              OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
          call sys_halt()
        end if
        
      case(POISSON_NEUMANN)
        call output_line("Invalid boundary segment.", &
            OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
        call sys_halt()
        
      case default
        call output_line("Invalid problem subtype.", &
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
  integer :: cproblemtype,cproblemsubtype,ivectorpos,icomp,isegment
  
  ! Get configuration from quick-access arrays
  cproblemtype    = rcollection%IquickAccess(1)
  cproblemsubtype = rcollection%IquickAccess(2)
  ivectorpos      = rcollection%IquickAccess(3)
  icomp           = rcollection%IquickAccess(4)
  isegment        = rcollection%IquickAccess(5)

  
  
  
  Dvalues(1) = ddirichlet
  return
  
#if defined(CASE_POISSON_DIRICHLET)

    select case(rcollection%IquickAccess(1))
    case(1)
      
      ! Return Dirichlet boundary values for all situations.
      Dvalues(1) = ddirichlet
      
    case default
      call output_line("There is no boundary integral in this benchmark", &
          OU_CLASS_ERROR,OU_MODE_STD,"getBoundaryValues_Poisson")
      call sys_halt()
    end select

#elif defined(CASE_POISSON_DIRICHLET_NEUMANN)

    ! local variables
    real(DP) :: dnx,dny

    select case(rcollection%IquickAccess(1))
    case(1)

      ! Return Dirichlet boundary values for all situations.
      Dvalues(1) = ddirichlet

    case(2)

      ! Compute the normal vector in the point on the boundary
      call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

      ! Return Neumann boundary values for all situations.
      Dvalues(1) = dneumann*dnx

    case(3)

      ! Compute the normal vector in the point on the boundary
      call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

      ! Return Neumann boundary values for all situations.
      Dvalues(1) = dneumann*dny

    case default
      call output_line("There is no boundary integral in this benchmark", &
          OU_CLASS_ERROR,OU_MODE_STD,"getBoundaryValues_Poisson")
      call sys_halt()
    end select

#else
#error 'Test case is undefined.'
#endif
  
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

  ! local variable
  integer :: cproblemtype,cproblemsubtype,cderivativebase

  select case (cderivative)
  case (DER_FUNC)
    ! u(x,y) = SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case default
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select

  return
  
  ! Get configuration from quick-access arrays
  cproblemtype    = rcollection%IquickAccess(1)
  cproblemsubtype = rcollection%IquickAccess(2)
  cderivativebase = rcollection%IquickAccess(3)

  ! u(x,y) = SIN(PI * x) * SIN(PI * y)
  select case(cderivativebase)
  case (DER_FUNC)
    ! We start from the solution
    ! u(x,y) = SIN(PI * x) * SIN(PI * y)
    ! ----------------------------------

    select case (cderivative)
    case (DER_FUNC)
      ! u(x,y) = SIN(PI * x) * SIN(PI * y)
      Dvalues (:,:) = sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
      
    case (DER_DERIV_X)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
      Dvalues (:,:) = SYS_PI * &
          cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
      
    case (DER_DERIV_Y)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
      Dvalues (:,:) = SYS_PI * &
          sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))
      
    case default
      ! Unknown. Set the result to 0.0.
      Dvalues = 0.0_DP
    end select
    
  case (DER_DERIV_X)
    ! We start from the x-derivative
    ! u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
    ! -----------------------------------------
    
    select case (cderivative)
    case (DER_FUNC)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
      Dvalues (:,:) = SYS_PI * &
          cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
      
    case (DER_DERIV_X)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_xx(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
      Dvalues (:,:) = -SYS_PI*SYS_PI * &
          sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
      
    case (DER_DERIV_Y)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_xy(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
      Dvalues (:,:) = SYS_PI*SYS_PI * &
          cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))
      
    case DEFAULT
      ! Unknown. Set the result to 0.0.
      Dvalues = 0.0_DP
    end select
    
  case (DER_DERIV_Y)
    ! We start from the y-derivative
    ! u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
    ! -----------------------------------------
    
    select case (cderivative)
    case (DER_FUNC)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
      Dvalues (:,:) = SYS_PI * &
          sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))
      
    case (DER_DERIV_Y)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_yy(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
      Dvalues (:,:) = -SYS_PI*SYS_PI * &
          sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
      
    case (DER_DERIV_X)
      !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
      ! => u_yx(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
      Dvalues (:,:) = SYS_PI*SYS_PI * &
          cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))
      
    case DEFAULT
      ! Unknown. Set the result to 0.0.
      Dvalues = 0.0_DP
    end select
    
  end select
  
  end subroutine getReferenceFunction_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivX_Poisson(cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xx(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xy(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine getReferenceDerivX_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivY_Poisson(cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yy(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yx(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine getReferenceDerivY_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXX_Poisson(cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xx(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xxx(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xxy(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine getReferenceDerivXX_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXY_Poisson(cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xy(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xyx(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xyy(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine getReferenceDerivXY_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYX_Poisson(cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yx(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yxx(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yxy(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine getReferenceDerivYX_Poisson

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYY_Poisson(cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yy(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yyx(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yyy(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine getReferenceDerivYY_Poisson

  ! ***************************************************************************

!<subroutine>

  elemental subroutine getAnalyticValues_Poisson(dx,dy,&
      du,du_x,du_y,du_xx,du_xy,du_yy)

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
    real(DP), intent(out) :: du

    ! Values of the first derivatives
    real(DP), intent(out) :: du_x,du_y

    ! Values of the second derivatives
    real(DP), intent(out) :: du_xx,du_xy,du_yy
!</output>
!</subroutine>

    ! Solution values
    du = sin(SYS_PI*dx) * sin(SYS_PI*dy)

    ! Values of first derivatives
    du_x = SYS_PI * cos(SYS_PI*dx) * sin(SYS_PI*dy)
    du_y = SYS_PI * sin(SYS_PI*dx) * cos(SYS_PI*dy)

    ! Values of second derivatives
    du_xx = -SYS_PI*SYS_PI * sin(SYS_PI*dx) * sin(SYS_PI*dy)
    du_xy =  SYS_PI*SYS_PI * cos(SYS_PI*dx) * cos(SYS_PI*dy)
    du_yy = -SYS_PI*SYS_PI * sin(SYS_PI*dx) * sin(SYS_PI*dy)

  end subroutine getAnalyticValues_Poisson

end module sse_callback_poisson
