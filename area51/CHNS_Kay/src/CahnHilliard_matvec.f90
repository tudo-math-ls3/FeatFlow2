!##############################################################################
!# ****************************************************************************
!# <name> CahnHilliard_matvec </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains matrix-vector assembly routines for the heat conduction
!# problem. The following routines can be found here:
!# 0.)
!# 1.) CH_allocMatVec
!#     -> Allocate memory for matrices/vectors, generate 'static' matrices that
!#        don't change in time.
!#
!# 2') CH_generateStaticMatrices
!#     -> generate static matrix for CH problem including Laplacian and Mass matrix
!#
!# 2.) CH_calcRHS
!#     -> Calculate RHS vector. (doesn't implement BC's.)
!#
!# 3.) CH_doneMatVec
!#     -> Release memory allocated in CH_initMatVec.
!#
!# 4.) CH_nonlinearMatMul
!#     -> Calculate multiplication of nonlinear matrix times a vector
!#
!# 5.) CH_assembleMatrix
!#     -> Assemble linear and (nonliear) Matrix.
!#
!# </purpose>
!##############################################################################

module CahnHilliard_matvec

  use fsystem
  use storage
  use basicgeometry
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use sortstrategy
  use coarsegridcorrection
  use ucd
  use timestepping
  use genoutput

  use spdiscprojection
  use nonlinearsolver
  use scalarpde
  use stdoperators

 use trilinearformevaluation
  use collection
  use paramlist
    
!~~~~~~~~NS modules~~~~~~~~~~~~~~~~~~~
   use ccbasic
!   use cccallback
! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   use CahnHilliard_callback
   use CahnHilliard_basic
!
  IMPLICIT NONE

!</constantblock>

!<constantblock description="Matrix type ID's specifying the general matrix class to set up.">

  ! Standard matrix.
  integer, parameter :: CHMASM_MTP_AUTOMATIC         = 0
  
  ! Standard matrix with decoupled velocity blocks
  integer, parameter :: CHMASM_MTP_DECOUPLED         = 1
  
  ! Extended 'full-tensor' matrix with submatrices A11, A12, A21, A22, all independent from
  ! each other.
  integer, parameter :: CHMASM_MTP_FULLTENSOR        = 2

!<types>

!<typeblock>

  ! This routine describes the nonlinear system matrix. The system matrix
  ! does actually not exist in memory -- since it's nonlinear! Therefore,
  ! this structure contains all parameters and settings which are necessary
  ! do apply(!) the matrix to a vector or to evaluate it.
  ! ('Evaluate a nonlinear matrix' means: Using a given FE-function $y$,
  ! assemble the linear matrix D(y).)
!~~~~~~~~~~~~~~~~~~~~~~~Nonlinear Matrix for CH problem~~~~~~~~~~~~~~~~
!  (  A     B  )  [c] = [f]
!  (  D     C  )  [w] = [g]
!
!  with A = \alpha M + \gamma {\bf u} grad, where M <--> Identity operator in c
!       B = \eta - div(grad w), where \eta including time stepping, b(c) and \rho(c)
!       D = \tau N(c) + \theta div(grad c)
!       C = \beta M   where M <--> Identity operator in w
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  type t_nonlinearCHMatrix
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    real(DP) :: dalpha = 0.0_DP
    
	! GAMMA-parameter that controls the weight in front of the
    ! Convetive term...
    real(DP) :: dgamma = 0.0_DP

    ! THETA-parameter that controls the weight of the A matrix
    ! in the core equation. =1.0 for stationary simulations.
    real(DP) :: dtheta = 0.0_DP

    ! ETA-parameter that switch the B-term on/off.
    real(DP) :: deta = 0.0_DP
    
    ! TAU-parameter that switch the D-term on/off
    real(DP) :: dtau = 0.0_DP

    ! BETA-parameter that switch the B-term on/off.
    real(DP) :: dbeta = 0.0_DP

    ! Weight for the Newton matrix N*(u).
    ! = 0.0 deactivates the Newton part.
    real(DP) :: dnewton = 0.0_DP

    ! STABILISATION: so far we don't need it. For stabilization, refer to Feng's paper
    
    ! STABILISATION: Viscosity parameter. Used for stabilisation schemes when
    ! a nonlinearity is set up.
    real(DP) :: dnu = 0.0_DP
    
    ! STABILISATION: Stabilisation parameter for streamline diffusion, upwind and
    ! edge-oriented stabilisation. If iupwind=CCMASM_STAB_STREAMLINEDIFF, a value of
    ! 0.0 deactivates any stabilisation.
    real(DP) :: dupsam = 0.0_DP
    
    ! MATRIX RESTRICTION: Parameter to activate matrix restriction.
    ! Can be used to generate parts of the matrices on coarse grids where the
    ! aspect ratio of the cells is large. Only applicable for $\tilde Q_1$
    ! discretisations.
    ! Standard = 0 = deactivate matrix restriction
    integer :: iadaptiveMatrices = 0
    
    ! MATRIX RESTRICTION: Threshold parameter for adaptive matrix generation
    ! of coarse grid matrices (here: maximum aspect ratio).
    ! Only applicable if iadaptiveMatrices <> 0.
    ! Standard = 20.0
    real(DP) :: dadmatthreshold = 20.0_DP
    
    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...).
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! pointer to a template FEM matrix that defines the structure of
    ! Laplace/... matrices.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM => null()

    ! pointer to Laplace matrix.
    type(t_matrixScalar), pointer :: p_rmatrixLaplace => null()
    ! pointer to a Mass matrix.
    ! May point to NULL() during matrix creation.
    type(t_matrixScalar), pointer :: p_rmatrixMass => null()
    ! pointer to a Conv matrix.
    ! May point to NULL() during matrix creation.
    type(t_matrixScalar), pointer :: p_rmatrixConv => null()

    ! pointer to a A-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixA => null()

    ! pointer to a B-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB => null()

    ! pointer to a C-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixC => null()

    ! pointer to a D-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixD => null()

  end type
  
CONTAINS

  ! ***************************************************************************
!~~~~~~~We may first remark this part, Shall we assemble matrix for all levels?
! MCai, assembling matrix part need to be rewritten, if we want to use multigrid
! code to solve the system.  How to do that? for all levels?
! rewrite this part according to cc2d

!<subroutine>

  subroutine CH_assembleMatrix (rCHproblem, rmatrix,&
                      rnonlinearCHMatrix,rCHvector,rCHparams,&
                                    rNSproblem, rNSvector)
  !<description>
  ! Calculates the matrices of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system, allocates memory
  ! for a RHS vector.
  !</description>

  !<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem

  ! A parameter list with informations from the DAT file.
  type(t_parlist), intent(IN) :: rCHparams

  ! A t_nonlinearCHMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, p_rmatrixTemplateXXXX
  ! must be initialised as well as p_rdiscretisation!
  ! The new matrix is created based p_rmatrixTemplateXXXX as well as
  ! p_rdiscretisation. Memory is allocated automatically if it's missing.
  type(t_nonlinearCHMatrix), intent(INOUT) :: rnonlinearCHMatrix

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity.
  type(t_vectorBlock), intent(IN), optional :: rCHvector

  ! The destination matrix which should be set up.
  ! If not initialised, a new matrix is created (as if CCMASM_ALLOCxxxx
  ! was specified).
  ! If initialised, the existing matrix is updated or recreated, depending on
  ! coperation.
  type(t_matrixBlock), intent(INOUT) :: rmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~NS  problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform,rformmass
    type(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), pointer :: p_rmatrixLaplace,p_rmatrixMass,p_rmatrixConv, p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM

	! A temp matrix for convection term
    type(t_matrixScalar), target :: rmatrixConvTmp
   

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Arrays for the Cuthill McKee renumbering strategy
    !integer, dimension(1) :: H_Iresort
    !integer(PREC_VECIDX), dimension(:), pointer :: p_Iresort
 
    ! Parameters from the DAT file
    real(DP) :: Pe, eps, gamma
    real(DP) :: dgamma_tmp
    CHARACTER(LEN=10) :: Sstr
 
     ! get the coefficients from the parameter list
 !    call parlst_getvalue_string (rCHparams, 'EQUATION', 'Pe', Sstr, '1.0')
 !    READ(Sstr,*) Pe
 !    call parlst_getvalue_string (rCHparams, 'EQUATION', 'Epsilon', Sstr, '1.0')
 !    READ(Sstr,*) eps
    Pe=1000.0_DP
    eps=0.02_DP
!	gamma=0.1_DP


!    print *, rmatrix%RmatrixBlock(1,1)%h_Da

 ! For (1,1) block <------> matrix A

      ! Release the matrix if present.
    call lsysbl_releaseMatrix (rmatrix)
    p_rmatrixTemplateFEM => rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixTemplateFEM

!      p_rmatrixTemplateFEM => rnonlinearCHMatrix%p_rmatrixTemplateFEM

    if (associated(p_rdiscretisation)) then
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
    else
      ! No discretisation structure; create the matrix directly as 3x3 matrix.
      call lsysbl_createEmptyMatrix (rmatrix,NDIM2D)
    end if
        
    ! the entries.,here we need to use 'empty' rather than 'remove'
    call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    call lsysbl_updateMatStrucInfo (rmatrix)

    ! (1,1) Block include both variable density based Mass + Conv
!	if (rnonlinearCHMatrix%dalpha .ne. 0.0_DP) then
!	  call CH_generateVarMass (rnonlinearCHmatrix,rCHvector, rmatrix%RmatrixBlock(1,1))
!	  call lsyssc_scaleMatrix(rmatrix%RmatrixBlock(1,1), &
!	                             rnonlinearCHMatrix%dalpha)
!	end if

!   ! (1,1) Block, mass matrix does not depend on c.
    if (rnonlinearCHMatrix%dalpha .ne. 0.0_DP) then
      call lsyssc_matrixLinearComb (&
        rnonlinearCHMatrix%p_rmatrixMass,rnonlinearCHMatrix%dalpha,&
        rmatrix%RmatrixBlock(1,1),0.0_DP,&
        rmatrix%RmatrixBlock(1,1),&
        .false.,.false.,.true.,.true.)
    end if

! For (1,1) block, we also have the convective term.
! For {\bf u}(t_n) \cdot grad c, we need to first create a scalar matrix for this
! How can we incorporate in the velocity in previous step for assembling convetion
! term?
! Do we need upwind?
! call the upwind method to calculate the nonlinear matrix.
    if (rnonlinearCHMatrix%dgamma .ne. 0.0_DP) then

! Q: hot incorporate dgamma into the matrix?
!      call conv_upwind2d (rvector, rvector, &
!                         dvecWeight, 0.0_DP,&
!                         rupwind, CONV_MODMATRIX, &
!                         rConvmatrixTmp)

!      ! the convection depending on other outer vector
      call CH_generateConvMatrix(rCHproblem,rCHvector, rmatrixConvTmp, &
	                             rNSproblem,rNSvector)
      ! Add the convective matrix to (1,1) block
      call lsyssc_matrixLinearComb (&
        rmatrixConvTmp,rnonlinearCHMatrix%dgamma,&
        rmatrix%RmatrixBlock(1,1),1.0_DP,&
        rmatrix%RmatrixBlock(1,1),&
        .false.,.false.,.true.,.true.)

  ! MCai
  ! We use rmatrixConvTmp to evaluate p_rMatrixConv,
  ! Question: how to pass a matrix to ...p_rMatrixConv,
  ! rnonlinearCHMatrix is just input parameter.
  !         rnonlinearCHMatrix%p_rMatrixConv=>rmatrixConvTmp

    end if
!    dgamma_tmp=rnonlinearCHMatrix%dgamma*..
!    call lsyssc_matrixLinearComb (&
!      rConvmatrixTmp,dgamma_tmp,&
!      rmatrix%RmatrixBlock(1,1),1.0_DP,&
!      rmatrix%RmatrixBlock(1,1),&
!      .false.,.false.,.true.,.true.)

! For (1,2) block <----------> B, we need c(t_n)
    if (rnonlinearCHMatrix%deta .ne. 0.0_DP) then
      ! Concentration dependent Mobility
	  call CH_generateVarLaplace (rnonlinearCHMatrix,rCHvector, rmatrix%RmatrixBlock(1,2))
	  call lsyssc_scaleMatrix(rmatrix%RmatrixBlock(1,2), &
	                            rnonlinearCHMatrix%deta/Pe)

      ! Constant Mobility
!      call lsyssc_matrixLinearComb (&
!            rnonlinearCHMatrix%p_rmatrixLaplace,gamma*rnonlinearCHMatrix%deta/Pe,&
!            rmatrix%RmatrixBlock(1,2),0.0_DP,&
!            rmatrix%RmatrixBlock(1,2),&
!            .false.,.false.,.true.,.true.)
	end if

! For (2,1) block <----------> D, this one is nonlinear part, how to do this?
! Nonlinear part : 1/eps*2 N(.)
    if (rnonlinearCHMatrix%dtau .ne. 0.0_DP) then
! MCai, we first remark this part, because we are only constructing preconditioner

! For N'(.)/eps: here we use Jacobian N', and eps is in callback.
      call CH_generateNonlinearMat (rCHproblem, rCHvector, &
	                             rmatrix%RmatrixBlock(2,1))
	end if

! Add Lapace term to nonlinear matrix.
! for Laplacian term in (2,1) block: Laplace term: - eps* Laplace
    if (rnonlinearCHMatrix%dtheta .ne. 0.0_DP) then
      call lsyssc_matrixLinearComb (rnonlinearCHMatrix%p_rmatrixLaplace,&
	        eps*rnonlinearCHMatrix%dtheta,&
            rmatrix%RmatrixBlock(2,1),1.0_DP,&
            rmatrix%RmatrixBlock(2,1),&
            .false.,.false.,.true.,.true.)
    end if

!-------------------------------------------------------------------------------
! For (2,2) block <----------> C, this one corresponds to Mass matrix
    if (rnonlinearCHMatrix%dbeta .ne. 0.0_DP) then
      call lsyssc_matrixLinearComb (&
	        rnonlinearCHMatrix%p_rmatrixMass,rnonlinearCHMatrix%dbeta,&
	    	rmatrix%RmatrixBlock(2,2), 0.0_DP,&
	    	rmatrix%RmatrixBlock(2,2),&
            .false.,.false.,.true.,.true.)
	end if

    ! That's it, all submatrices are basically set up.
    !
    ! Update the structural information of the block matrix, as we manually
    ! changed the submatrices:
    call lsysbl_updateMatStrucInfo(rMatrix)


    ! release tmp matrix
    call lsyssc_releaseMatrix(rmatrixConvTmp)
! Perhaps, we need to release some tmp matrix....
!    call lsyssc_releaseMatrix (rConvmatrixTmp)
! We need to write this ....   MCai

  end subroutine
  ! ***************************************************************************

!<subroutine>

  subroutine CH_nonlinearMatMul (rnonlinearCHMatrix,rx,rd,dcx,dcd,ry, &
                                 rNSproblem, rNSvector)

!<description>
  ! This routine performs a matrix-vector multiplication with a nonlinear
  ! matrix:
  !      rd := dcx A(ry) rx + dcd rd
  ! with the system matrix A(.) defined by the configuration in rnonlinearCHMatrix.
  ! The caller must initialise the rnonlinearCHMatrix according to how the
  ! matrix should look like.
  !
  ! The parameter ry is optional. If specified, this parameter defines where to
  ! evaluate the nonlinearity (if the system matrix $A$ contains a nonlinearity).
  ! If not specified, ry=rx is assumed.
  !
  ! The routine will not include any boundary conditions in the defect.
!</description>

  ! A t_nonlinearCHMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  type(t_nonlinearCHMatrix), intent(INOUT) :: rnonlinearCHMatrix

  ! This vector specifies the 'x' that is multiplied to the matrix.
  type(t_vectorBlock), intent(IN), target :: rx

  ! Multiplication factor in front of the term 'A(ry) rx'.
  real(DP), intent(IN) :: dcx

  ! Multiplication factor in front of the term 'rd'.
  real(DP), intent(IN) :: dcd

  ! OPTIONAL: Point where to evaluate the nonlinearity. If not specified,
  ! ry=rx is assumed.
  type(t_vectorBlock), intent(IN), target, optional :: ry

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  type(t_vectorBlock), intent(INOUT) :: rd
  
!~~~~~~~~~~~~~~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(IN), optional :: rNSproblem
  type(t_vectorBlock), intent(IN), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_ry
    type(t_matrixBlock) :: rmatrix
    
    ! DEBUG!!!
    !real(dp), dimension(:), pointer :: p_DdataX,p_DdataD
    
    !call lsysbl_getbase_double (rx,p_DdataX)
    !call lsysbl_getbase_double (rd,p_DdataD)
    
    p_ry => rx
    if (present(ry)) p_ry => ry

    ! Probably weight the input vector.
    if (dcd .ne. 1.0_DP) then
      call lsysbl_scaleVector (rd,dcd)
    end if
    
    ! The system matrix looks like:
    !
    !    ( A   B  )
    !    ( D   C  )
    ! Create a temporary matrix that covers this structure.
    call lsysbl_createEmptyMatrix (rmatrix,NDIM2D)

    ! Put references to the Laplace- and B-matrices to Aij. assembleDefect
    ! needs this template matrix to provide the structure for the stabilisation
    ! routines! The B-matrices are needed later.
!     call lsyssc_duplicateMatrix (rnonlinearCHMatrix%p_rmatrixA,&
!         rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!     call lsyssc_duplicateMatrix (rnonlinearCHMatrix%p_rmatrixC,&
!         rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!     call lsyssc_duplicateMatrix (rnonlinearCHMatrix%p_rmatrixB,&
!         rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
! ! MCai, we may also duplicate content(for nonlinear matrix..)
! !    call lsyssc_duplicateMatrix (rnonlinearCHMatrix%p_rmatrixB,&
! !        rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!     call lsyssc_duplicateMatrix (rnonlinearCHMatrix%p_rmatrixD,&
!         rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    ! Update the structural information of the block matrix, as we manually
    ! changed the submatrices:
    call lsysbl_updateMatStrucInfo (rmatrix)
    
    ! First, we assemble the defect that arises in the velocity components.
    ! assembleDefect handles exactly these submatrices.
    call assembleDefect (rnonlinearCHMatrix,rmatrix,rx,rd,p_ry, -dcx, &
	                     rNSproblem, rNSvector)

! MCai, ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MCai, Note that here it is different from cc2d

    ! Release the temporary matrix, we don't need it anymore.
    call lsysbl_releaseMatrix (rmatrix)
    
  contains

    subroutine assembleDefect (rnonlinearCHMatrix,&
        rmatrix,rvector,rdefect,ry,dvectorWeight, rNSproblem, rNSvector)

    ! Calculate defect:
    !       rdefect = rdefect - dvectorWeight * rmatrix(ry) * rvector
    ! Assembles the velocity defect in the block matrix rmatrix at position
    ! itop..itop+1 in the velocity vector. rdefect must have been initialised
    ! with the right hand side vector.
    !
    ! With a matrix of the theoretical form:
    !       [ A    B ]
    !       [ D    C ]
    !
    !       A := dalpha*M + dgamma Conv
    !       B := deta Lap
    !       C := dbeta Mass
    !       D := dtau N(.) + dtheta Lap
    !
    ! and c=dvectorWeight, the routine will construct
    !
    ! A t_nonlinearCHMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCHMatrix), intent(INOUT) :: rnonlinearCHMatrix

    ! Reference to the system matrix. Only the structure of the matrix
    ! is used to reconstruct the structure of the discretisation.
    ! The content of the matrix is not changed or used.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Solution vector.
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Weight for the velocity vector; usually = 1.0
    real(DP), intent(IN) :: dvectorWeight
    
    ! the vector for nonlinearity.
    type(t_vectorBlock), intent(IN) :: ry

!~~~~~~~~~~~~~~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    type(t_problem), intent(IN), optional :: rNSproblem
    type(t_vectorBlock), intent(IN), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! local variables
      logical :: bshared

  	  real(DP) :: Pe=1000.0_DP
	  real(DP) :: eps=0.02_DP
!	  real(DP) :: gamma=0.1_DP

      ! A tmp Mass matrix
	  type(t_matrixScalar) :: rmatrixMassTmp
	  ! A tmp Laplace matrix
	  type(t_matrixScalar) :: rmatrixLaplaceTmp
      ! A tmp conv matrix
	  type(t_matrixScalar) :: rmatrixConvTmp

    ! DEBUG!!!
     real(dp), dimension(:), pointer :: p_DdataX,p_DdataD
      ! DEBUG
      call lsysbl_getbase_double (rvector,p_DdataX)
      call lsysbl_getbase_double (rdefect,p_DdataD)

      !MCai~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------------------Start to cal dcx * rmatrix * rx---------------------
	  ! (1,1) block
      ! Subtract the mass matrix stuff?
      if (rnonlinearCHMatrix%dalpha .ne. 0.0_DP) then
        ! MCai, for CH model, there is no concentration dependent mass matrix
!	    ! Concentration dependent Mass matrix
!		call CH_generateVarMass (rnonlinearCHMatrix,rvector, rmatrixMassTmp)
 
!	    call lsyssc_scalarMatVec (rmatrixMassTmp, &
!        rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
!         -rnonlinearCHMatrix%dalpha, 1.0_DP)
!		call lsyssc_releaseMatrix(rmatrixMassTmp)

        ! Constant Mass matrix
    	call lsyssc_scalarMatVec (rnonlinearCHMatrix%p_rmatrixMass, &
           rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
           -rnonlinearCHMatrix%dalpha, 1.0_DP)
      end if

      if (rnonlinearCHMatrix%dgamma .ne. 0.0_DP) then
         ! MCai, we need to use convective matrix, \gamma ConvMat
!MCai, we need to evaluate p_rmatrixConv first, we can not do this so far.
!         call lsyssc_scalarMatVec (rnonlinearCHMatrix%p_rmatrixConv, &
!           rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
!           -rnonlinearCHMatrix%gamma, 1.0_DP)
     
         ! Alternative way
    	call CH_generateConvMat(rnonlinearCHMatrix,rvector,rmatrixConvTmp,&
	                               rNSproblem,rNSvector)
        call lsyssc_scalarMatVec (rmatrixConvTmp, &
           rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rnonlinearCHMatrix%dgamma, 1.0_DP)
        ! release temporary convection matrix
        call lsyssc_releaseMatrix(rmatrixConvTmp)

      end if
	      
      ! ---------------------------------------------------
	  ! (1,2) block
      ! Subtract the Laplace matrix stuff from chemical potential
      ! Take care, here the input of the following subroutine is different from above
      if (rnonlinearCHMatrix%deta .ne. 0.0_DP) then
	    ! Concentration dependent mobility
		call CH_generateVarLaplace (rnonlinearCHMatrix,rvector,rmatrixLaplaceTmp)
		call lsyssc_scalarMatVec (rmatrixLaplaceTmp, &
               rvector%RvectorBlock(2), rdefect%RvectorBlock(1), &
               -rnonlinearCHMatrix%deta/Pe, 1.0_DP)
	    call lsyssc_releaseMatrix(rmatrixLaplaceTmp)

!        call lsyssc_scalarMatVec (rmatrix%RmatrixBlock(1,2), &
!               rvector%RvectorBlock(2), rdefect%RvectorBlock(1), &
!               -rnonlinearCHMatrix%deta/Pe, 1.0_DP)

        ! constant mobility
!        call lsyssc_scalarMatVec (rnonlinearCHMatrix%p_rmatrixLaplace, &
!               rvector%RvectorBlock(2), rdefect%RvectorBlock(1), &
!               -gamma*rnonlinearCHMatrix%deta/Pe, 1.0_DP)
      end if

!-----------------Then the second equation--------------------------
      ! (2,1) block

	  ! then the Laplace term:  -eps * Lap
      if (rnonlinearCHMatrix%dtheta .ne. 0.0_DP) then
	    call lsyssc_scalarMatVec (rnonlinearCHMatrix%p_rmatrixLaplace, &
               rvector%RvectorBlock(1), rdefect%RvectorBlock(2), &
               -eps*rnonlinearCHMatrix%dtheta, 1.0_DP)
	  end if
      ! That was not easy -- the adventure begins now... The nonlinearity!
      ! It is different from preconditioner
      ! Since this subroutine is for setting up right hand side, where we need it again?
      ! first the nonlinear term: 1/eps * f(c)
      if (rnonlinearCHMatrix%dtau .ne. 0.0_DP) then
!	    print *, 'so far, we have no nonlinear matrix, use ry to assemble nonlinear part'
!	    call lsyssc_scalarMatVec (rnonlinearCHMatrix%p_rmatrixD, &
!               ry%RvectorBlock(1), rdefect%RvectorBlock(2), &
!               -rnonlinearCHMatrix%dtau/eps, 1.0_DP)

!MCai, we first remark this to make sure nonlinear solver is correct.
      ! nonlinear term: rd=rd - (f(c), \ki)/eps
         call nonlinear_calcDefect (rnonlinearCHMatrix, rvector%RvectorBlock(1), &
		      rdefect%RvectorBlock(2), -rnonlinearCHMatrix%dtau/(eps), 1.0_DP)
	  end if
      
	  ! (2,2) block
      if (rnonlinearCHMatrix%dbeta .ne. 0.0_DP) then
	    call lsyssc_scalarMatVec (rnonlinearCHMatrix%p_rmatrixMass, &
               rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
               -rnonlinearCHMatrix%dbeta, 1.0_DP)
      end if

    end subroutine

  end subroutine

  ! ***************************************************************************
!<subroutine>
!MCai
  subroutine nonlinear_calcDefect (rnonlinearMatrix, rvectorScalar, &
		         rdefect, dcx, dcd)
!<description>
  ! Calculates:
  !      rd=dcd * rd + dcx * N(rx)
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_nonlinearCHMatrix), intent(INOUT) :: rnonlinearMatrix

  ! A pointer to the RHS vector.
  type(t_vectorScalar), intent(IN)  :: rvectorScalar
  type(t_vectorScalar), intent(INOUT) :: rdefect
   
  real(DP), intent(in) :: dcx
  real(DP), intent(in) :: dcd
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!</inputoutput>

!</subroutine>

  ! local variables
    INTEGER :: i
  
    ! A linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! temp rproblem
    type(t_CHproblem) :: rproblem

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

	! A tmp vector created from rvectorScalar
    type(t_vectorBlock), target :: rvectorBlocktmp
    type(t_vectorScalar) :: rdefecttmp
    ! A tmp collection
    type(t_collection) :: rcollectionTmp

!~~~Mcai,
!   Notice that there is a nonlinear term f(\phi)
    !
    ! At first set up the corresponding linear form (f,Phi_j):

    ! Create two temp vectors.
     call lsysbl_createVecFromScalar (rvectorScalar,rvectorBlocktmp)
     call lsyssc_duplicateVector (rdefect,rdefecttmp,&
	                       LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC

    ! get the discretisation
    p_rdiscretisation => rnonlinearMatrix%p_rdiscretisation

    call collct_init(rcollectionTmp)
  
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is

    rcollectionTmp%p_rvectorQuickAccess2 => rvectorBlocktmp

    ! We use rcollection to pass rvectorBlocktmp into nonlinear coeff.
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
              rdefecttmp,coeff_nonlinear,&
              rcollectionTmp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! rd=dcd * rd + dcx * rdefecttmp
    call lsyssc_vectorLinearComb (rdefecttmp, &
	          rdefect,dcx,dcd)

	! Release tem vectors.
    call lsysbl_releaseVector(rvectorBlocktmp)
    call lsyssc_releaseVector(rdefecttmp)

    call collct_done(rcollectionTmp)
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine CH_calcRHS (rCHproblem,rCHrhs, &
	            rNSproblem, rNSvector, rNSrhs)
  
!<description>
  ! Calculates the RHS vector at the current point in time.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem

  ! A pointer to the RHS vector.
  type(t_vectorBlock), intent(INOUT) :: rCHrhs
!~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !optional
  type(t_problem), intent(IN), optional :: rNSproblem
  type(t_vectorBlock), intent(IN), target, optional :: rNSvector
  type(t_vectorBlock), intent(IN), target, optional :: rNSrhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
  
    ! A linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    p_rdiscretisation => rCHrhs%p_rblockDiscr

    ! Put the current simulation time as parameter "TIME" into the collection.
    ! Also set Dquickaccess (1) to the simulation time for faster access by the
    ! callback routine.
    rCHproblem%rcollection%Dquickaccess (1) = rCHproblem%rtimedependence%dtime

    call collct_setvalue_real(rCHproblem%rcollection,'TIME',&
         rCHproblem%rtimedependence%dtime,.TRUE.)
    ! The vector structure is done but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC

    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the first block

    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
	
    ! for the first component \phi
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.TRUE.,&
              rCHrhs%rvectorBlock(1),coeff_RHS_phi,&
              rCHproblem%rcollection)

    ! for the first component: chemical potential w
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.TRUE.,&
              rCHrhs%rvectorBlock(2),coeff_RHS_w,&
              rCHproblem%rcollection)
!~~~~~~~~~NS problem need to be used to generate RHS~~~~~~~~~~~~~~~~~~~
! Time dependent cases: coupled with Navier-Stokes equations

! In this model, we treat convective term explicitly, so we do not need to worry
! about right hand side, (no matter whether we have rNSvector)

    if(present(rNSvector)) then


  ! Then, determine we use explicit convection term, if explicit, we need to
  ! treat the convective term as source term.
  ! Debug
  !  write(*,*)'check whether rNSvector is sorted?'
  !  mcai=lsysbl_isVectorSorted (rACvector)
  !  print *, mcai
  !  mcai=lsysbl_isVectorSorted (rNSvector)
  !  print *, mcai

  !      Implicit_Convective=0

 !	  if (Implicit_Convective .eq. 0) then

! we need to calculate ({\bf u} \cdot \nabla \phi, \psi)
!        rlinform%itermCount = 1
!        rlinform%Idescriptors(1) = DER_FUNC

   !Implicit_Convective=1, we don't need rNSvector in the RHS of AC problem

! Take care of rcollection%...., do not abuse it.
!        rACproblem%rcollection%p_rvectorQuickAccess2 => rNSvector

!        call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
!              rACrhs%RvectorBlock(1),coeff_RHS,&
!              rACproblem%rcollection)

    	! First, deal with the source term from AC problem
	    !MCai, here we set bclear=.false. and call coeff_RHS0
        ! in the collection.
        ! if (rACproblem%rtimedependence%dtime .eq. 0.0_DP) then
!	! at the initial step, the right hand side is 0.
!	! Debug,
!	  write(*,*)'Make sure that this step is called at the very beginning'
!      call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
!              rACrhs%RvectorBlock(1),coeff_RHS0)
!	else

! MCai, we do not need to treat the nonlinear term and mass conservation term(AllenCah problem)
! in the right hand side, but we keep it here if needed.
!      ! first calculate (f(phi(t_n)), psi)
!        call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
!              rACrhs%RvectorBlock(1),coeff_RHS1,rACproblem%rcollection)

! We need the following subroutine to guarantee mass conservation, i.e. add constraint
! We need innerVector=rACvector.
!
!        call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
!              rACrhs%RvectorBlock(1),coeff_RHS3,&
!              rACproblem%rcollection)
!	  end if

    end if

!*******************************************************************************

    ! Remove the "TIME"-parameter from the collection again.
    call collct_deletevalue (rCHproblem%rcollection,'TIME')
    
  end subroutine

  ! ***************************************************************************
!MCai, March 19, 2009, we use the following code in assemble_Defect, rather than
! CH_assembleMatrix
  subroutine CH_generateConvMat(rNonlinearCHMatrix,rvector,rmatrixConvTmp, &
	                               rNSproblem,rNSvector)


 !<description>
  ! Calculates entries of Convetive matrix, it is done by using trilinear form eval
  !
  ! the matrix should first be generated by
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_nonlinearCHMatrix), intent(IN) :: rNonlinearCHMatrix
  
  type(t_vectorBlock), intent(IN) :: rvector

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_matrixScalar), intent(INOUT) :: rmatrixConvTmp

!~~~~~~~~~~~~~~~from NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(IN), optional :: rNSproblem
  type(t_vectorBlock), intent(IN), target, optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!</inputoutput>

!</inputoutput>

!</subroutine>

    ! local variables
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
	type(t_bilinearform) :: rform
    ! The trilinear form specifying the underlying PDE of the discretisation.
!    type(t_trilinearForm) :: rtriform
 
    type(t_collection) :: rcollection
	
    ! Ask the problem structure to give us the discretisation structure
!    p_rdiscretisation => rvector%p_rblockDiscr
    p_rdiscretisation => rnonlinearCHMatrix%p_rdiscretisation


    ! Create matrix by discretisation
    call bilf_createMatrixStructure (p_rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrixConvTmp)

    ! Prepare a collection structure to be passed to the callback
    ! routine. We attach the vector T in the quick-access variables
    ! so the callback routine can access it.
    call collct_init(rcollection)
    rcollection%p_rvectorQuickAccess1 => rNSvector

	rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X   !(1,i) for trial function
    rform%Idescriptors(2,1) = DER_FUNC      !(2,i) for test function
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_FUNC

    ! In this case, we have nonconstant coefficients.
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff(:) = .false.
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    ! The collection is passed as additional parameter. That's the way
    ! how we get the vector to the callback routine.
    call bilf_buildMatrixScalar (rform,.true.,rmatrixConvTmp,&
                             coeff_Conv,rcollection)
    

!   ! Alternative way, use trilinearform
!	! We use p_rvectorQuickAccess2 to get rNSvector
!    rcollection%p_rvectorQuickAccess2 => rNSvector
!    rtriform%itermcount=1
!    rtriform%BconstantCoeff= .FALSE.
!    rtriform%ballcoeffconstant = .FALSE.
!    rtriform%Dcoefficients(1)= 1.0_DP
!    rtriform%Idescriptors(1,1)=DER_FUNC
!    rtriform%Idescriptors(2,1)=DER_DERIV_X
!    rtriform%Idescriptors(3,1)=DER_FUNC

!    ! Now we can build the matrix entries. We use coeff_Conv1
!    CALL trilf_buildMatrixScalar (rtriform,.FALSE.,rmatrixConvTmp,&
!             rvector%rvectorBlock(1),coeff_Conv1,rcollection)

!    ! For second term of velocity, we use coeff_Conv2
!    rtriform%Idescriptors(2,1)=DER_DERIV_Y
 
!    !MCai, question, shall we use rvector%rvectorBlock(2)?
!    ! No, there is only 1 block is needed.
!    CALL trilf_buildMatrixScalar (rtriform,.FALSE.,rmatrixConvTmp,&
!             rvector%rvectorBlock(1),coeff_Conv2,rcollection)
    

!MCai, or we should use the following code
!    !MCai, we need the rVeclvector block to evaluate the coefficient.
!	! because velocity act as convection coefficient.
!    rproblem%rcollection%p_rvectorQuickAccess1=>rNSvector

!    ! use trilinear form evaluation.
!	! use coeff_Conv (in call_back) to pass the coefficient.
!    call trilf_buildMatrixScalar (rform,.TRUE.,rmatrix,rvector,&
!                                  coeff_Conv,rcollection)


    call collct_done(rcollection)

  end subroutine


  ! ***************************************************************************

!MCai, March 18, 2009,
  subroutine CH_generateConvMatrix(rproblem,rvector, rmatrixConvTmp, &
	                               rNSproblem,rNSvector)


 !<description>
  ! Calculates entries of Convetive matrix, it is done by using trilinear form eval
  !
  ! the matrix should first be generated by
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT) :: rproblem

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_matrixScalar), intent(INOUT) :: rmatrixConvTmp

  type(t_vectorBlock), intent(IN) :: rvector

!~~~~~~~~~~~~~~~from NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(IN), optional :: rNSproblem
  type(t_vectorBlock), intent(IN), target, optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!</inputoutput>

!</subroutine>

    ! local variables
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
	type(t_bilinearform) :: rform
    ! The trilinear form specifying the underlying PDE of the discretisation.
!    type(t_trilinearForm) :: rtriform
 
    type(t_collection) :: rcollection
	
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rproblem%rlevelInfo(rproblem%NLMAX)%rdiscretisation

    ! Create matrix by discretisation
    call bilf_createMatrixStructure (p_rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrixConvTmp)

    ! Prepare a collection structure to be passed to the callback
    ! routine. We attach the vector T in the quick-access variables
    ! so the callback routine can access it.
    call collct_init(rcollection)

    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X   !(1,i) for trial function
    rform%Idescriptors(2,1) = DER_FUNC      !(2,i) for test function
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_FUNC

    ! In this case, we have nonconstant coefficients.
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff(:) = .false.
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0

    ! Prepare a collection structure to be passed to the callback
    ! routine. We attach the vector T in the quick-access variables
    ! so the callback routine can access it.
    call collct_init(rcollection)
    rcollection%p_rvectorQuickAccess1 => rNSvector

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    ! The collection is passed as additional parameter. That's the way
    ! how we get the vector to the callback routine.
    call bilf_buildMatrixScalar (rform,.true.,rmatrixConvTmp,&
                             coeff_Conv,rcollection)
    

!   ! Alternative way, use trilinearform
!	! We use p_rvectorQuickAccess2 to get rNSvector
!    rcollection%p_rvectorQuickAccess2 => rNSvector
!    rtriform%itermcount=1
!    rtriform%BconstantCoeff= .FALSE.
!    rtriform%ballcoeffconstant = .FALSE.
!    rtriform%Dcoefficients(1)= 1.0_DP
!    rtriform%Idescriptors(1,1)=DER_FUNC
!    rtriform%Idescriptors(2,1)=DER_DERIV_X
!    rtriform%Idescriptors(3,1)=DER_FUNC

!    ! Now we can build the matrix entries. We use coeff_Conv1
!    CALL trilf_buildMatrixScalar (rtriform,.FALSE.,rmatrixConvTmp,&
!             rvector%rvectorBlock(1),coeff_Conv1,rcollection)

!    ! For second term of velocity, we use coeff_Conv2
!    rtriform%Idescriptors(2,1)=DER_DERIV_Y
 
!    !MCai, question, shall we use rvector%rvectorBlock(2)?
!    ! No, there is only 1 block is needed.
!    CALL trilf_buildMatrixScalar (rtriform,.FALSE.,rmatrixConvTmp,&
!             rvector%rvectorBlock(1),coeff_Conv2,rcollection)
    

!MCai, or we should use the following code
!    !MCai, we need the rVeclvector block to evaluate the coefficient.
!	! because velocity act as convection coefficient.
!    rproblem%rcollection%p_rvectorQuickAccess1=>rNSvector

!    ! use trilinear form evaluation.
!	! use coeff_Conv (in call_back) to pass the coefficient.
!    call trilf_buildMatrixScalar (rform,.TRUE.,rmatrix,rvector,&
!                                  coeff_Conv,rcollection)


    call collct_done(rcollection)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine CH_generateNonlinearMat (rproblem,rvector, rmatrix)
  
!<description>
  ! Calculates entries of Nonlinear mass matrix, it is of polynomial form
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT) :: rproblem

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_matrixScalar), intent(INOUT) :: rmatrix

  type(t_vectorBlock), intent(IN), target :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: j

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Structure for the bilinear form for assembling Laplacian, Mass
    type(t_bilinearForm) :: rform
 
    type(t_collection) :: rcollectionTmp
    ! -----------------------------------------------------------------------
    ! Basic CH problem
    ! -----------------------------------------------------------------------
    
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rproblem%rlevelInfo(rproblem%NLMAX)%rdiscretisation
    
	! Create matrix by discretisation
    call bilf_createMatrixStructure (p_rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrix)

    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    
    ! In this case, we have nonconstant coefficients.
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff(:) = .false.
    rform%Dcoefficients(1)  = 1.0
!    rform%Dcoefficients(2)  = 1.0

    ! Prepare a collection structure to be passed to the callback
    ! routine. We attach the vector T in the quick-access variables
    ! so the callback routine can access it.
    call collct_init(rcollectionTmp)

!MCai, we need the first block to evaluate the coefficient.
    rcollectionTmp%p_rvectorQuickAccess1=>rvector

    ! coeff_NonlinearPrec is in callback subroutine.
    call bilf_buildMatrixScalar (rform,.TRUE.,&
                                 rmatrix,coeff_NonlinearPrec,&
                                 rcollectionTmp)
    call collct_done(rcollectionTmp)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine CH_generateVarMass (rnonlinearCHmatrix,rvector, rmatrix)
  
!<description>
  ! Calculates entries of Laplace matrix with variable coefficient, the coefficient
  ! depending on rvector. This subroutine is similar to variable coeff Laplacian
  ! problem.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_nonlinearCHmatrix), intent(INOUT) :: rnonlinearCHmatrix

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_matrixScalar), intent(INOUT) :: rmatrix

  type(t_vectorBlock), intent(IN), target :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: j

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Structure for the bilinear form for assembling Laplacian, Mass
    type(t_bilinearForm) :: rform
 
    ! A temporary collection
    type(t_collection) :: rcollectionTmp
    ! -----------------------------------------------------------------------
    ! Basic CH problem
    ! -----------------------------------------------------------------------
    
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rnonlinearCHmatrix%p_rdiscretisation
    
	! Create matrix by discretisation
    call bilf_createMatrixStructure (p_rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrix)

    ! For Mass matrix itermCount =1.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
   

    ! In this case, we have nonconstant coefficients.
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff(:) = .false.
    rform%Dcoefficients(1)  = 1.0_DP

    ! Prepare a collection structure to be passed to the callback
    ! routine. We attach the vector T in the quick-access variables
    ! so the callback routine can access it.
    call collct_init(rcollectionTmp)

!MCai, we need the first block to evaluate the coefficient.
    rcollectionTmp%p_rvectorQuickAccess1=>rvector

    ! coeff_VarMass is in callback subroutine.
    call bilf_buildMatrixScalar (rform,.TRUE.,&
                                 rmatrix,coeff_VarMass,&
                                 rcollectionTmp)

    call collct_done(rcollectionTmp)

  end subroutine


  ! ***************************************************************************


!<subroutine>

  subroutine CH_generateVarLaplace (rnonlinearCHMatrix,rvector, rmatrix)
  
!<description>
  ! Calculates entries of Laplace matrix with variable coefficient, the coefficient
  ! depending on rvector. This subroutine is similar to variable coeff Laplacian
  ! problem.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_nonlinearCHMatrix), intent(INOUT) :: rnonlinearCHMatrix

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_matrixScalar), intent(INOUT) :: rmatrix

  type(t_vectorBlock), intent(IN), target :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: j

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Structure for the bilinear form for assembling Laplacian, Mass
    type(t_bilinearForm) :: rform
    
	! A tmp collcetion
    type(t_collection) :: rcollectionTmp
    ! -----------------------------------------------------------------------
    ! Basic CH problem
    ! -----------------------------------------------------------------------
    
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rnonlinearCHMatrix%p_rdiscretisation
    
	! Create matrix by discretisation
    call bilf_createMatrixStructure (p_rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrix)

    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y

    ! In this case, we have nonconstant coefficients.
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff(:) = .false.
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0

    ! Prepare a collection structure to be passed to the callback
    ! routine. We attach the vector T in the quick-access variables
    ! so the callback routine can access it.
    call collct_init(rcollectionTmp)

!MCai, we need the first block to evaluate the coefficient.
    rcollectionTmp%p_rvectorQuickAccess1=>rvector

    ! coeff_VarLaplace is in callback subroutine.
    call bilf_buildMatrixScalar (rform,.TRUE.,&
                                 rmatrix,coeff_VarLaplace,&
                                 rcollectionTmp)

    call collct_done(rcollectionTmp)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine CH_allocMatVec (rCHproblem, rCHparams)
!<description>
  ! Allocate memory for rCHproblem system matrix.
  ! Allocates memory for all matrices and vectors of the problem on the heap
  ! by evaluating the parameters in the problem structure.
  ! Matrices/vectors of global importance are added to the collection
  ! structure of the problem, given in rCHproblem. Matrix/vector entries
  ! are not calculated, they are left uninitialised.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), target :: rCHproblem

  type(t_parlist), intent(IN) :: rCHparams

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,cmatBuildtype,ielementtype,icubtemp
    integer(I32) :: icubA
    character(LEN=SYS_NAMELEN) :: sstr
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM
   
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_vectorBlock), pointer :: p_rrhs
  

    call parlst_getvalue_int (rCHparams,'CH-DISCRETISATION',&
                              'iElementtype',ielementtype,0)

    call parlst_getvalue_string (rCHparams,'CH-DISCRETISATION',&
                                 'scubA',sstr,'')
    if (sstr .eq. '') then
    	icubtemp = CUB_G2X2
      call parlst_getvalue_int (rCHparams,'CH-DISCRETISATION',&
                                'icubA',icubtemp,icubtemp)
      icubA = icubtemp
    else
      icubA = cub_igetID(sstr)
    end if

    ! Initialise all levels...
    do i=rCHproblem%NLMIN,rCHproblem%NLMAX

      ! -----------------------------------------------------------------------
      ! Basic Cahn Hilliard problem
      ! -----------------------------------------------------------------------

      ! Ask the problem structure to give us the discretisation structure
      ! (init_discretisation should be called before this subroutine)
      p_rdiscretisation => rCHproblem%RlevelInfo(i)%rdiscretisation

      ! The global system looks as follows:
      !
      !    ( A    B )
      !    ( D    C )
      !
      ! Get a pointer to the template FEM matrix. This is used for the
      ! Laplace/Stokes matrix and probably for the mass matrix.

      p_rmatrixTemplateFEM => rCHproblem%RlevelInfo(i)%rmatrixTemplateFEM
!MCai, whether rmatrixTemplateFEM need to be evaluated first?

      ! Create the matrix structure: p_rmatrixTemplateFEM is the output

      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixTemplateFEM,cconstrtype=BILF_MATC_ELEMENTBASED)

      ! We use Ciarlet-Raviart mixed FEM: RspatialDiscr(1)=RspatialDiscr(2)
     ! call bilf_createMatrixStructure (&
     !           p_rdiscretisation%RspatialDiscr(2),LSYSSC_MATRIX9,&
     !           p_rmatrixTemplateFEM,cconstrtype=BILF_MATC_ELEMENTBASED)

!       call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                    rCHproblem%RlevelInfo(i)%rmatrixA,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!       call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                    rCHproblem%RlevelInfo(i)%rmatrixB,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!       call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                    rCHproblem%RlevelInfo(i)%rmatrixC,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!       call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                    rCHproblem%RlevelInfo(i)%rmatrixD,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

!~~~~~It seems that one the following codes is not needed~~~~~~~~~~~~~~~~~~~~
      ! Allocate memory for the entries; don't initialise the memory.
!       call lsyssc_allocEmptyMatrix ( rCHproblem%RlevelInfo(i)%rmatrixA,LSYSSC_SETM_UNDEFINED)
!       call lsyssc_allocEmptyMatrix ( rCHproblem%RlevelInfo(i)%rmatrixB,LSYSSC_SETM_UNDEFINED)
!       call lsyssc_allocEmptyMatrix ( rCHproblem%RlevelInfo(i)%rmatrixC,LSYSSC_SETM_UNDEFINED)
!       call lsyssc_allocEmptyMatrix ( rCHproblem%RlevelInfo(i)%rmatrixD,LSYSSC_SETM_UNDEFINED)

!If for other type of problems, we have Dirichlet BC, we need also alloc the golbal matrix
! by the following way
      ! Initialise the block matrix with default values based on
      ! the discretisation.
!      if (associated(p_rdiscretisation)) then
!        call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
!      else
!        ! No discretisation structure; create the matrix directly as 3x3 matrix.
!        call lsysbl_createEmptyMatrix (rmatrix,NDIM2D)
!      end if
!      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!   rCHproblem%RlevelInfo(i)%rmatrix%RmatirxBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!   rCHproblem%RlevelInfo(i)%rmatrix%RmatirxBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!   rCHproblem%RlevelInfo(i)%rmatrix%RmatirxBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!   rCHproblem%RlevelInfo(i)%rmatrix%RmatirxBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
 
     ! -----------------------------------------------------------------------
	  ! MCai, we may need temporary vector
      ! Temporary vectors
      !
      ! Now on all levels except for the maximum one, create a temporary
      ! vector on that level, based on the block discretisation structure.
      ! It's used for building the matrices on lower levels.

        call lsysbl_createVecBlockByDiscr (&
            rCHproblem%RlevelInfo(i)%rdiscretisation,&
            rCHproblem%RlevelInfo(i)%rtempVector,.true.)

    end do
    
    ! (Only) on the finest level, we need to have to allocate a RHS vector
    ! and a solution vector.
    !
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template.
    ! Initialise the vectors with 0.
    call lsysbl_createVecBlockByDiscr (rCHproblem%RlevelInfo(rCHproblem%NLMAX)%&
	              rdiscretisation,rCHproblem%rrhs,.true.)

! Do we need the following code?
    p_rrhs    => rCHproblem%rrhs
! Save the solution/RHS vector to the collection. Might be used
! later (e.g. in nonlinear problems)
    call collct_setvalue_vec(rCHproblem%rcollection,'RHS',p_rrhs,.TRUE.)
 
  end subroutine


  ! ***************************************************************************
!<subroutine>

  subroutine CH_generateStaticMatrices (rCHproblem,rlevelInfo)

  
!<description>
  ! Calculates entries of all static matrices: Laplace and Mass matrix
  ! in the specified problem structure, i.e. the entries of all matrices
  ! that do not change during the computation or which serve as a template for
  ! generating other matrices.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT) :: rCHproblem

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_CHproblem_lvl), intent(INOUT),target :: rlevelInfo
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: j

    ! A pointer to the Laplace and mass matrix. Remark: Usually, p_rmatrixConv
	! is not a static matrix, because it will depend velocity, changing with time.
    type(t_matrixScalar), pointer :: p_rmatrixLaplace, p_rmatrixMass, p_rmatrixConv

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

!MCai
!	! We may need bilinear form if we use it
!	type(t_bilinearform) :: rform

    ! The trilinear form specifying the underlying PDE of the discretisation.
    type(t_trilinearForm) :: rtriform
  
    !
    ! -----------------------------------------------------------------------
    ! Basic CH problem
    ! -----------------------------------------------------------------------
    
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%rdiscretisation
    
    ! Get a pointer to the (scalar) Laplace matrix and Mass matrix:
    p_rmatrixLaplace => rlevelInfo%rmatrixLaplace
    p_rmatrixMass => rlevelInfo%rmatrixMass
!    p_rmatrixConv => rlevelInfo%rmatrixConv
    
    ! The global system looks as follows:
    !
    !    ( A    B )
    !    ( D    C )
    !
    ! with A = ? B=? C=? D=?
    
    ! If there is an existing Laplace matrix, release it.
    call lsyssc_releaseMatrix (rlevelInfo%rmatrixLaplace)

    call lsyssc_duplicateMatrix (rlevelInfo%rmatrixTemplateFEM,&
                rlevelInfo%rmatrixLaplace,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

! We assemble Laplace matrix, we can use two ways.
!     rform%itermCount = 2
!     rform%Idescriptors(1,1) = DER_DERIV_X
!     rform%Idescriptors(2,1) = DER_DERIV_X
!     rform%Idescriptors(1,2) = DER_DERIV_Y
!     rform%Idescriptors(2,2) = DER_DERIV_Y
!
!     ! In the standard case, we have constant coefficients:
!     rform%ballCoeffConstant = .TRUE.
!     rform%BconstantCoeff = .TRUE.
!     rform%Dcoefficients(1)  = 1.0_DP
!     rform%Dcoefficients(2)  = 1.0_DP
!
!     ! Now we can build the matrix entries.
!     ! We specify the callback function coeff_Stokes for the coefficients.
!     ! As long as we use constant coefficients, this routine is not used.
!     ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!     ! the framework will call the callback routine to get analytical data.
!     !
!     ! We pass our collection structure as well to this routine,
!     ! so the callback routine has access to everything what is
!     ! in the collection.
!
!     call bilf_buildMatrixScalar (rform,.TRUE.,&
!                                  p_rmatrixLaplace,coeff_CahnHilliard,&
!                                  rCHproblem%rcollection)
!
! ! Alternatively,
     call stdop_assembleLaplaceMatrix (p_rmatrixLaplace,.true.,1.0_DP)

    ! -----------------------------------------------------------------------
    ! Mass matrices. They are used in so many cases, it's better we always
    ! have them available.
    ! -----------------------------------------------------------------------

    ! If there is an existing mass matrix, release it.
    call lsyssc_releaseMatrix (rlevelInfo%rmatrixMass)

    ! Generate mass matrix. The matrix has basically the same structure as
    ! our template FEM matrix, so we can take that.
    call lsyssc_duplicateMatrix (rlevelInfo%rmatrixTemplateFEM,&
                rlevelInfo%rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    ! Change the discretisation structure of the mass matrix to the
    ! correct one; at the moment it points to the discretisation structure
    ! of the Laplace matrix...
    call lsyssc_assignDiscrDirectMat (rlevelInfo%rmatrixMass,&
        rlevelInfo%rdiscretisationMass)
!MCai
!or
!    call lsyssc_assignDiscretDirectMat (rlevelInfo%rmatrixMass,&
!        rlevelInfo%rdiscretisation(1))
!or
!    call lsyssc_assignDiscretDirectMat (rlevelInfo%rmatrixMass,&
!        rlevelInfo%rdiscretisation)

    ! call the standard matrix setup routine to build the matrix.
    call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixMass,DER_FUNC,DER_FUNC)

! For convective matrix, we need to assemble it each timestep
! !    ! -----------------------------------------------------------------------
! !    ! Conv matrices. if the velocity field does not change with time
! !    ! -----------------------------------------------------------------------
!
! !    ! If there is an existing Conv matrix, release it.
!     call lsyssc_releaseMatrix (rlevelInfo%rmatrixConv)
!
! !    ! Generate mass matrix. The matrix has basically the same structure as
! !    ! our template FEM matrix, so we can take that.
!     call lsyssc_duplicateMatrix (rlevelInfo%rmatrixTemplateFEM,&
!                 rlevelInfo%rmatrixConv,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
! !    ! Change the discretisation structure of the mass matrix to the
! !    ! correct one; at the moment it points to the discretisation structure
! !    ! of the Laplace matrix...
!     call lsyssc_assignDiscretDirectMat (rlevelInfo%rmatrixConv,&
!         rlevelInfo%rdiscretisationMass)

    ! We assemble Convective matrix, we can use two ways.
! MCai
! where to get rVelvector????????, we do not have it in generate_static
!     rCHproblem%rcollection%p_rvectorQuickAccess1 => rVelvector
!     rtriform%itermcount=1
!     rtriform%BconstantCoeff= .FALSE.
!     rtriform%ballcoeffconstant = .FALSE.
!     rtriform%Dcoefficients(1)= 1.0_DP
!     rtriform%Idescriptors(1,1)=DER_FUNC
!     rtriform%Idescriptors(2,1)=DER_DERIV_X
!     rtriform%Idescriptors(3,1)=DER_FUNC
!
!     ! Now we can build the matrix entries. We use coeff_Conv1
!     CALL trilf_buildMatrixScalar (rtriform,.FALSE.,p_rmatrixConv,&
!              rvector%rvectorBlock(1),coeff_Conv1,rCHproblem%rcollection)
!
!     ! For second term of velocity, we use coeff_Conv2
!     rtriform%Idescriptors(2,1)=DER_DERIV_Y
!
!     CALL trilf_buildMatrixScalar (rtriform,.FALSE.,p_rmatrixConv,&
!              rvector%rvectorBlock(2),coeff_Conv2,rCHproblem%rcollection)
   

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine CH_doneMatVec (rCHproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem
!</inputoutput>

!</subroutine>

    integer :: ihandle,i

    ! Release matrices and vectors on all levels
    do i=rCHproblem%NLMAX,rCHproblem%NLMIN,-1

      ! Delete the variables from the collection.
!      call collct_deletevalue (rCHproblem%rcollection,'SYSTEM',i)
!      call collct_deletevalue (rCHproblem%rcollection,'LAPLACE',i)
!      call collct_deletevalue (rCHproblem%rcollection,'MASS',i)

      ! Delete the system matrix
!       call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixA)
!       call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixB)
!       call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixC)
!       call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixD)

! Delete the static matrices
      call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixLaplace)
      call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixMass)
!      call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixConv)
      call lsyssc_releaseMatrix (rCHproblem%RlevelInfo(i)%rmatrixTemplateFEM)
    
      call lsysbl_releasevector (rCHproblem%RlevelInfo(i)%rtempVector)
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rCHproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue (rCHproblem%rcollection,'RHS')
    call collct_deletevalue (rCHproblem%rcollection,'TIME')

  end subroutine


end module
