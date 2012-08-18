!##############################################################################
!# ****************************************************************************
!# <name> AllenCahn_matvec </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains matrix-vector assembly routines for the Allen_Cahn
!# problem. The following routines can be found here:
!#
!# 1.) AC_initMatVec
!#     -> Allocate memory for matrices/vectors, generate 'static' matrices that
!#        don't change in time.
!#
!# 2.) AC_calcRHS
!#     -> Calculate RHS vector. (Doesn't implement BC's.)
!#
!# 3.) AC_doneMatVec
!#     -> Release memory allocated in AC_initMatVec.
!# </purpose>
!##############################################################################

module AllenCahn_matvec

  use fsystem
  use storage
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

  use collection
  use paramlist

!~~~~~~~~NS modules~~~~~~~~~~~~~~~~~~~
  use ccbasic
  use cccallback
!  use ccnonstationary   ! do not use this, otherwise, wrong.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  use AllenCahn_callback
  use AllenCahn_basic
  
  IMPLICIT NONE
  
CONTAINS

  !****************************************************************************
  subroutine AC_assembleMatPoly(rACproblem, rACvector)
! for assembling the matrix from nonlinear term: we do it this way:
!  ( (\phi(t_n)^2 -1) \phi(t_n+1), \psi)
  type(t_ACproblem), intent(INOUT) :: rACproblem
  type(t_vectorBlock), intent(INOUT) :: rACvector

	! local variables
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Structure for the bilinear form for assembling Laplacian, Mass
    type(t_bilinearForm) :: rform
    ! A tmp collection
    type(t_collection) :: rcollectionTmp
    ! A tmp vector
    type(t_vectorblock), target :: rACvectorTemp
    type(t_matrixBlock), POINTER :: p_rmatrixPoly
    integer :: i


    call lsysbl_copyvector (rACvector,rACvectorTemp)
   
    do i = rACproblem%NLMAX, rACproblem%NLMIN,-1
      p_rdiscretisation => rACproblem%RlevelInfo(i)%p_rdiscretisation
      p_rmatrixPoly => rACproblem%RlevelInfo(i)%rmatrixPoly

      ! One and only term for the mass matrix
      rform%itermCount=1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixPoly)

      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixPoly%RmatrixBlock(1,1))

      ! Update the structural information
	  call lsysbl_updateMatStrucInfo (p_rmatrixPoly)

      ! In this case, we have nonconstant coefficients.
      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff(:) = .false.
      rform%Dcoefficients(1)  = 1.0

      ! Prepare a collection structure to be passed to the callback
      ! routine. We attach the vector T in the quick-access variables
      ! so the callback routine can access it.
      call collct_init(rcollectionTmp)

      !MCai, we need the first block to evaluate the coefficient.
      rcollectionTmp%p_rvectorQuickAccess1=>rACvectorTemp

      ! Build matrix, similar to variable coeff Mass matrix
      call bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrixPoly%RmatrixBlock(1,1),coeff_Poly,&
                                   rcollectionTmp)

      if(i.GT.rACproblem%NLMIN) then
        call prjF2C(rACvectorTemp,rACproblem%Rlevelinfo(i-1)%p_rdiscretisation)
     end if

      call collct_done(rcollectionTmp)
      
    end do

    ! release temp vector block
    call lsysbl_releasevector(rACvectorTemp)
    
  end subroutine

  !****************************************************************************
  subroutine AC_assembleMatConv(rACproblem, rACvector, rNSproblem, rNSvector)
! for assembling the matrix from nonlinear term: we do it this way:
!  ( (\phi(t_n)^2 -1) \phi(t_n+1), \psi)
  type(t_ACproblem), intent(INOUT) :: rACproblem
  type(t_vectorBlock), intent(INOUT) :: rACvector

  ! information from NS problem, we need the discretisation from NS problem
  type(t_problem), intent(INOUT) :: rNSproblem
  type(t_vectorBlock), intent(INOUT) :: rNSvector

	! local variables
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Structure for the bilinear form for assembling Laplacian, Mass
    type(t_bilinearForm) :: rform
    ! A tmp collection
    type(t_collection) :: rcollectionTmp
    ! A tmp vector
    type(t_vectorblock), target :: rNSvectorTemp
    type(t_matrixBlock), POINTER :: p_rmatrixConv
    integer :: i

! For stabilization, if necessary
!    type(t_jumpStabilisation) :: rjumpstabil



    call lsysbl_copyvector (rNSvector,rNSvectorTemp)

    do i = rACproblem%NLMAX, rACproblem%NLMIN,-1
      p_rdiscretisation => rACproblem%RlevelInfo(i)%p_rdiscretisation
      p_rmatrixConv => rACproblem%RlevelInfo(i)%rmatrixConv

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixConv)
	  
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixConv%RmatrixBlock(1,1),cconstrtype=BILF_MATC_ELEMENTBASED)

      ! Update the structural information
       call lsysbl_updateMatStrucInfo (p_rmatrixConv)

      ! In this case, we have nonconstant coefficients.
	  ! Bilinear form for convective matrix
	  ! The first coefficient is from the term of x derivative
      rform%itermCount=1
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_FUNC
      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff(:) = .false.
      rform%Dcoefficients(1)  = 1.0_DP

      ! Prepare a collection structure to be passed to the callback
      ! routine. We attach the vector T in the quick-access variables
      ! so the callback routine can access it.
      call collct_init(rcollectionTmp)

      !MCai, we need the first block to evaluate the coefficient.
      rcollectionTmp%p_rvectorQuickAccess1=>rNSvectorTemp

      ! Build matrix, similar to variable coeff Mass matrix
      call bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrixConv%RmatrixBlock(1,1),coeff_Conv1,&
                                   rcollectionTmp)

      ! The first coefficient is from the term of y derivative
      ! take care, our itermCound is only 1.
      rform%Idescriptors(1,1) = DER_DERIV_Y
      call bilf_buildMatrixScalar (rform,.false.,&
                                  p_rmatrixConv%RmatrixBlock(1,1), coeff_Conv2,&
                                  rcollectionTmp)

! For jumpstabilization
!     rjumpstabil%dnu = 0.0_DP
!     rjumpstabil%dgamma = 0.01_DP
!     rjumpstabil%dgammastar = 0.01_DP
!     rjumpstabil%dtheta = rtimestepping%dweightMatrixLHS
    
!MCai, think about which discretisation we should use? I do this similar to Cui
!     call conv_JumpStabilisation2d ( &
!       rjumpstabil, conv_modmatrix, p_rmatrixConv%rmatrixBlock(1,1),&
!         p_rdiscretisation)

      IF(i.GT.rACproblem%NLMIN) &
        call prjF2C(rNSvectorTemp,rNSproblem%rlevelinfo(i-1)%rdiscretisation)

      call collct_done(rcollectionTmp)

    end do
    
    ! release temp vector block
	call lsysbl_releasevector(rNSvectorTemp)
  end subroutine


!**********************************************************************************
 subroutine prjF2C(rvector,rdiscretisation)
! This subroutine is wrote by Zejun Cui, the purpose is to use fine grid vector to get
! coarse grid vector. Projection from fine to coarse.
 
    type(t_vectorBlock), INTENT(INOUT) :: rvector
    !type(t_problem_lvl), INTENT(IN),TARGET :: rlevelInfo
    type(t_blockDiscretisation),intent(in):: rdiscretisation
    
    !local
    type(t_interlevelProjectionBlock) :: rprojection
    type(t_vectorBlock):: rtmpvector
    type(t_vectorScalar) :: rvectorTemp
    INTEGER:: NEQ

    call lsysbl_createVecBlockByDiscr (&
            rdiscretisation,rtmpvector,.false.)
  
    call mlprj_initProjectionVec (rprojection,rtmpVector)
    NEQ = mlprj_getTempMemoryVec (rprojection,rtmpVector,rvector)
   
    IF (NEQ .NE. 0) call lsyssc_createVector (rvectorTemp,NEQ,.FALSE.)
   
    call mlprj_performInterpolation (rprojection,rtmpvector,rvector, &
                                       rvectorTemp)
       
    IF (NEQ .NE. 0) call lsyssc_releaseVector (rvectorTemp)
    call lsysbl_swapVectors (rvector,rtmpvector)
    call lsysbl_releaseVector (rtmpvector)
    call mlprj_doneProjection (rprojection)
  
  end subroutine

  !****************************************************************************
!~~~~~~~~~~This is main module, we need to modify many stuff~~~~~~~~~~~~~~~~~~~

!<subroutine>

  subroutine AC_initMatVec (rACproblem,rparams, &
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
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem

  ! A parameter list with informations from the DAT file.
  type(t_parlist), INTENT(IN) :: rparams

!~~~~~~NS problem, optional, it depends on whether we treat the convective term
!~~~~~~~~~implicitly or explicitly. If implicitly, we need them~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(in), optional :: rNSproblem
  type(t_vectorBlock), intent(in), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform,rformmass
    type(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), POINTER :: p_rmatrixStatic,p_rmatrixMass,p_rmatrix
    type(t_vectorBlock), POINTER :: p_rrhs

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! Arrays for the Cuthill McKee renumbering strategy
!    INTEGER, DIMENSION(1) :: H_Iresort
!    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iresort

    ! Parameters from the DAT file
    REAL(DP) :: alpha11,alpha12,alpha21,alpha22,beta1,beta2,gamma
    CHARACTER(LEN=10) :: Sstr
  
  
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.
    
    rform%itermCount = 7
    
    ! alpha * Laplace(u)
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_X
    
    rform%Idescriptors(1,3) = DER_DERIV_X
    rform%Idescriptors(2,3) = DER_DERIV_Y
    
    rform%Idescriptors(1,4) = DER_DERIV_Y
    rform%Idescriptors(2,4) = DER_DERIV_Y
    
    ! (beta1, beta2)^T * grad(u)
    rform%Idescriptors(1,5) = DER_DERIV_X
    rform%Idescriptors(2,5) = DER_FUNC
    
    rform%Idescriptors(1,6) = DER_DERIV_Y
    rform%Idescriptors(2,6) = DER_FUNC
    
    ! gamma * u
    rform%Idescriptors(1,7) = DER_FUNC
    rform%Idescriptors(2,7) = DER_FUNC

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    
    ! get the coefficients from the parameter list
    call parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA11', Sstr, '1.0')
    READ(Sstr,*) alpha11
    call parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA12', Sstr, '0.0')
    READ(Sstr,*) alpha12
    call parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA21', Sstr, '0.0')
    READ(Sstr,*) alpha21
    call parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA22', Sstr, '1.0')
    READ(Sstr,*) alpha22
    call parlst_getvalue_string (rparams, 'EQUATION', 'BETA1', Sstr, '0.0')
    READ(Sstr,*) beta1
    call parlst_getvalue_string (rparams, 'EQUATION', 'BETA2', Sstr, '0.0')
    READ(Sstr,*) beta2
    call parlst_getvalue_string (rparams, 'EQUATION', 'GAMMA', Sstr, '0.0')
    READ(Sstr,*) gamma
    
    rform%Dcoefficients(1)  = alpha11
    rform%Dcoefficients(2)  = alpha12
    rform%Dcoefficients(3)  = alpha21
    rform%Dcoefficients(4)  = alpha22
    rform%Dcoefficients(5)  = beta1
    rform%Dcoefficients(6)  = beta2
    rform%Dcoefficients(7)  = gamma

    ! For the time dependent problem, we also need a mass matrix. We set up another
    ! bilinear form for that.

    rformmass%itermCount = 1

    ! One and only term for the mass matrix
    rformmass%Idescriptors(1,1) = DER_FUNC
    rformmass%Idescriptors(2,1) = DER_FUNC

    ! The coefficient in front of the term of the mass matrix
    rformmass%Dcoefficients(1)  = 1.0_DP

    do i = rACproblem%NLMIN, rACproblem%NLMAX
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rACproblem%RlevelInfo(i)%p_rdiscretisation
    
      ! -------------------------------------------------------------
      ! Laplace matrix
    
      p_rmatrixStatic => rACproblem%RlevelInfo(i)%rmatrixStatic
    
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixStatic)

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rACproblem%rcollection,'LAPLACE',p_rmatrixStatic,.TRUE.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixStatic%RmatrixBlock(1,1))
    
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      call lsysbl_updateMatStrucInfo (p_rmatrixStatic)
    
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_AllenCahn for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.

!~~~~~build the matrix entries, we may do other thing stuff here.
      call bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrixStatic%RmatrixBlock(1,1),coeff_AllenCahn,&
                                   rACproblem%rcollection)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! -------------------------------------------------------------
      ! Mass matrix, similar to Laplace matrix
      
      p_rmatrixMass => rACproblem%RlevelInfo(i)%rmatrixMass
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixMass)

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rACproblem%rcollection,'MASS',p_rmatrixMass,.TRUE.,i)
      
      ! The structure of the mass matrix is the same as the system matrix.
      ! Initialise the structure as "shared" with the system matrix.
      ! Reserve memory for the entries.
      call lsyssc_duplicateMatrix(p_rmatrixStatic%RmatrixBlock(1,1),&
           p_rmatrixMass%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      call lsysbl_updateMatStrucInfo (p_rmatrixMass)

      ! Now we can build the matrix entries of the mass matrix.
      call bilf_buildMatrixScalar (rformmass,.TRUE.,&
                                   p_rmatrixMass%RmatrixBlock(1,1))

      ! -------------------------------------------------------------
      ! System matrix.

      p_rmatrix => rACproblem%RlevelInfo(i)%rmatrix

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rACproblem%rcollection,'SYSTEM',p_rmatrix,.TRUE.,i)
      
      ! The structure of the mass matrix is the same as the system matrix.
      ! Initialise the structure as "shared" with the static matrix.
      ! Reserve memory for the entries.
      call lsyssc_duplicateMatrix(p_rmatrixStatic%RmatrixBlock(1,1),&
           p_rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      call lsysbl_updateMatStrucInfo (p_rmatrix)

    end do

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rACproblem%rrhs
    p_rmatrixStatic => rACproblem%RlevelInfo(rACproblem%NLMAX)%rmatrixStatic

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    call collct_setvalue_vec(rACproblem%rcollection,'RHS',p_rrhs,.TRUE.)

    ! Reserve memory for all the vectors on the finest level.
    call lsysbl_createVecBlockIndMat (p_rmatrixStatic,p_rrhs, .FALSE.)

    ! -------------------------------------------------------------
    ! Matrix resorting
    
    ! Finally, sort the matrices on all levels. We dfo this after the
    ! creation of the vectors, so the vectors stay unsorted!
!    do i = rACproblem%NLMIN, rACproblem%NLMAX
    
!      p_rmatrixStatic => rACproblem%RlevelInfo(i)%rmatrixStatic
!      p_rmatrixMass => rACproblem%RlevelInfo(i)%rmatrixMass
!      p_rmatrix => rACproblem%RlevelInfo(i)%rmatrix
      
      ! Allocate an array for holding the resorting strategy.
!      call storage_new ('AC_initMatVec', 'Iresort', &
!            p_rmatrixStatic%RmatrixBlock(1,1)%NEQ*2, ST_INT, h_Iresort(1), &
!            ST_NEWBLOCK_ZERO)
!      call storage_getbase_int(h_Iresort(1),p_Iresort)
      
      ! Calculate the resorting strategy.
!      call sstrat_calcCuthillMcKee (p_rmatrixStatic%RmatrixBlock(1,1),p_Iresort)
      
      ! Save the handle of the resorting strategy to the collection.
!      call collct_setvalue_int(rACproblem%rcollection,'LAPLACE-CM',h_Iresort(1),.TRUE.,i)
      
      ! Resort the matrices according to the sorting strategy.
      ! Note that as we share the structure between all matrices, we first
      ! have to sort the 'child' matrices...
!      call lsyssc_sortMatrix (p_rmatrixMass%RmatrixBlock(1,1),.TRUE.,&
!                              SSTRAT_CM,h_Iresort(1))
!      call lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.TRUE.,&
!                              SSTRAT_CM,h_Iresort(1))

      ! ...before sorting the 'parent' matrix!
!      call lsyssc_sortMatrix (p_rmatrixStatic%RmatrixBlock(1,1),.TRUE.,&
!                              SSTRAT_CM,h_Iresort(1))
                              
!    end do

  end subroutine

 ! ***************************************************************************

!<subroutine>
!MCai
  subroutine AC_calcRHS (rACproblem, rACvector, rACrhs,&
                          rNSproblem, rNSvector, rNSrhs)
  
!<description>
  ! Calculates the RHS vector at the current point in time.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem

  ! A pointer to the RHS vector.
  type(t_vectorBlock), INTENT(INOUT), target :: rACvector

  ! A pointer to the RHS vector.
  type(t_vectorBlock), INTENT(INOUT) :: rACrhs


!~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !optional
  type(t_problem), intent(IN), optional :: rNSproblem
  
  ! OPTIONAL:
  type(t_vectorBlock), intent(IN), target, optional :: rNSvector
  type(t_vectorBlock), intent(IN), target, optional :: rNSrhs

! Debug
!  logical :: mcai

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!</inputoutput>

!</subroutine>

  ! local variables
    INTEGER :: i
  
    ! A linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation

!~~~Mcai,
 	! A parameter to determine, whether we use implicit scheme or explicit scheme for
	! Allen-Cahn problem, if it is 0, we treat convective term explictly. Otherwise,implicitly
	integer :: Implicit_Convective=1

    ! Put the current simulation time as parameter "TIME" into the collection.
    ! Also set Dquickaccess (1) to the simulation time for faster access by the
    ! callback routine.

    rACproblem%rcollection%Dquickaccess (1) = rACproblem%rtimedependence%dtime
!~~~Mcai, do we need the following code in cc2d(determine time dependent)?
	call collct_setvalue_real(rACproblem%rcollection,'TIME',&
         rACproblem%rtimedependence%dtime,.TRUE.)

!~~~Mcai,
! Notice that there is a source term f(\phi)
 
    ! The vector structure is done but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
 
    p_rdiscretisation => rACproblem%RlevelInfo(rACproblem%NLMAX)%p_rdiscretisation

    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is

    rACproblem%rcollection%p_rvectorQuickAccess1 => rACvector

!MCai, here we set bclear=.false. and call coeff_RHS0
    ! in the collection.
!	if (rACproblem%rtimedependence%dtime .eq. 0.0_DP) then
!	! at the initial step, the right hand side is 0.
!	! Debug,
!	  write(*,*)'Make sure that this step is called at the very beginning'
!      call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
!              rACrhs%RvectorBlock(1),coeff_RHS0)
!	else
    
	! The following is for initialization, there is no time evolution.

      call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
              rACrhs%RvectorBlock(1),coeff_RHS0)

!	end if

    ! first calculate (f(phi(t_n)), psi)
        call linf_buildVectorScalar (&
               p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
               rACrhs%RvectorBlock(1),coeff_RHS1,rACproblem%rcollection)

!~~~~~~~~~NS problem need to be used to generate RHS~~~~~~~~~~~~~~~~~~~
! Time dependent cases: coupled with Navier-Stokes equations

    if(present(rNSvector)) then


  ! Then, determine we use explicit convection term, if explicit, we need to
  ! treat the convective term as source term.
  ! Debug

      Implicit_Convective=0

	  if (Implicit_Convective .eq. 0) then

  ! we need to calculate ({\bf u} \cdot \nabla \phi, \psi)
!        rlinform%itermCount = 1
!        rlinform%Idescriptors(1) = DER_FUNC

  !Implicit_Convective=1, we don't need rNSvector in the RHS of AC problem

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

!****************************************************************************************************

	  end if

    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! Remove the "TIME"-parameter from the collection again.
    call collct_deletevalue (rACproblem%rcollection,'TIME')
    
  end subroutine

!  ! *************Mcai, the following one is for conservative RHS******************************************
!<subroutine>
  subroutine AC_conservativeRHS(rACproblem, rACvector, rACrhs)

  
!<description>
  ! Calculates the RHS vector at the current point in time.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem

  ! A pointer to the RHS vector.
  type(t_vectorBlock), INTENT(INOUT), target :: rACvector

  ! A pointer to the RHS vector.
  type(t_vectorBlock), INTENT(INOUT) :: rACrhs
!</inputoutput>

!</subroutine>

  ! local variables
    INTEGER :: i
  
    ! A linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation

!~~~Mcai,
! Notice that there is a source term f(\phi)
 
    ! The vector structure is done but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
 
    p_rdiscretisation => rACproblem%RlevelInfo(rACproblem%NLMAX)%p_rdiscretisation

    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is

    rACproblem%rcollection%p_rvectorQuickAccess1 => rACvector

! MCai,
! We need the following subroutine to guarantee mass conservation, i.e. add constraint
! We need innerVector=rACvector.


           call linf_buildVectorScalar (&
                 p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
                 rACrhs%RvectorBlock(1),coeff_RHS3,&
                 rACproblem%rcollection)

    ! Remove the "TIME"-parameter from the collection again.
    call collct_deletevalue (rACproblem%rcollection,'TIME')
    
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine AC_doneMatVec (rACproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
!</inputoutput>

!</subroutine>

    INTEGER :: ihandle,i

    ! Release matrices and vectors on all levels
    do i=rACproblem%NLMAX,rACproblem%NLMIN,-1

      ! Delete the variables from the collection.
      call collct_deletevalue (rACproblem%rcollection,'SYSTEM',i)
      call collct_deletevalue (rACproblem%rcollection,'LAPLACE',i)
      call collct_deletevalue (rACproblem%rcollection,'MASS',i)

      ! Delete the system matrix
      call lsysbl_releaseMatrix (rACproblem%RlevelInfo(i)%rmatrix)
      
      ! Delete the mass matrix
      call lsysbl_releaseMatrix (rACproblem%RlevelInfo(i)%rmatrixMass)

      ! Delete the matrix
      call lsysbl_releaseMatrix (rACproblem%RlevelInfo(i)%rmatrixStatic)

      ! Release the permutation for sorting matrix/vectors
!      ihandle = collct_getvalue_int (rACproblem%rcollection,'LAPLACE-CM',i)
!      call storage_free (ihandle)
!      call collct_deletevalue (rACproblem%rcollection,'LAPLACE-CM',i)
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rACproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue (rACproblem%rcollection,'RHS')

  end subroutine

end module
