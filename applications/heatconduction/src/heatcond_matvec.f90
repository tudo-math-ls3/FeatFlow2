!##############################################################################
!# ****************************************************************************
!# <name> heatcond_matvec </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains matrix-vector assembly routines for the heat conduction
!# problem. The following routines can be found here:
!#
!# 1.) hc5_initMatVec
!#     -> Allocate memory for matrices/vectors, generate "static" matrices that
!#        do not change in time.
!#
!# 2.) hc5_calcRHS
!#     -> Calculate RHS vector. (Does not implement BC`s.)
!#
!# 3.) hc5_doneMatVec
!#     -> Release memory allocated in hc5_initMatVec.
!# </purpose>
!##############################################################################

module heatcond_matvec

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use sortstrategy
  use coarsegridcorrection
  use scalarpde
  use ucd
  use timestepping
  use genoutput
  
  use collection
  use paramlist
    
  use heatcond_callback
  
  use heatcond_basic
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_initMatVec (rproblem,rparams)
  
!<description>
  ! Calculates the matrices of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system, allocates memory
  ! for a RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! A parameter list with informations from the DAT file.
  type(t_parlist), intent(in) :: rparams
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform,rformmass
    type(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), pointer :: p_rmatrixStatic,p_rmatrixMass,p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Parameters from the DAT file
    real(DP) :: alpha11,alpha12,alpha21,alpha22,beta1,beta2,gamma
    character(LEN=10) :: Sstr
  
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
  
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
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    
    ! get the coefficients from the parameter list
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA11", Sstr, "1.0")
    read(Sstr,*) alpha11
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA12", Sstr, "0.0")
    read(Sstr,*) alpha12
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA21", Sstr, "0.0")
    read(Sstr,*) alpha21
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA22", Sstr, "1.0")
    read(Sstr,*) alpha22
    call parlst_getvalue_string (rparams, "EQUATION", "BETA1", Sstr, "0.0")
    read(Sstr,*) beta1
    call parlst_getvalue_string (rparams, "EQUATION", "BETA2", Sstr, "0.0")
    read(Sstr,*) beta2
    call parlst_getvalue_string (rparams, "EQUATION", "GAMMA", Sstr, "0.0")
    read(Sstr,*) gamma
    
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

    do i = rproblem%ilvmin, rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
    
      ! Get the cubature information structure
      p_rcubatureInfo => rproblem%RlevelInfo(i)%rcubatureInfo
    
      ! -------------------------------------------------------------
      ! Laplace matrix
    
      p_rmatrixStatic => rproblem%RlevelInfo(i)%rmatrixStatic
    
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixStatic)

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"LAPLACE",p_rmatrixStatic,.true.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixStatic%RmatrixBlock(1,1))
    
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_heatcond for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      call bilf_buildMatrixScalar (rform,.true.,&
          p_rmatrixStatic%RmatrixBlock(1,1),p_rcubatureInfo,&
          coeff_heatcond,rproblem%rcollection)

      ! -------------------------------------------------------------
      ! Mass matrix
      
      p_rmatrixMass => rproblem%RlevelInfo(i)%rmatrixMass
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixMass)

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"MASS",p_rmatrixMass,.true.,i)
      
      ! The structure of the mass matrix is the same as the system matrix.
      ! Initialise the structure as "shared" with the system matrix.
      ! Reserve memory for the entries.
      call lsyssc_duplicateMatrix(p_rmatrixStatic%RmatrixBlock(1,1),&
           p_rmatrixMass%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! Now we can build the matrix entries of the mass matrix.
      call bilf_buildMatrixScalar (rformmass,.true.,&
          p_rmatrixMass%RmatrixBlock(1,1),p_rcubatureInfo)

      ! -------------------------------------------------------------
      ! System matrix.

      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"SYSTEM",p_rmatrix,.true.,i)
      
      ! The structure of the mass matrix is the same as the system matrix.
      ! Initialise the structure as "shared" with the static matrix.
      ! Reserve memory for the entries.
      call lsyssc_duplicateMatrix(p_rmatrixStatic%RmatrixBlock(1,1),&
           p_rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
    end do

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs
    p_rmatrixStatic => rproblem%RlevelInfo(rproblem%ilvmax)%rmatrixStatic

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    call collct_setvalue_vec(rproblem%rcollection,"RHS",p_rrhs,.true.)

    ! Reserve memory for all the vectors on the finest level.
    p_rdiscretisation => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscretisation
    call lsysbl_createVectorBlock(p_rdiscretisation,p_rrhs,.true.)

    ! -------------------------------------------------------------
    ! Matrix resorting
    
    ! Finally, sort the matrices on all levels. We dfo this after the
    ! creation of the vectors, so the vectors stay unsorted!
    do i = rproblem%ilvmin, rproblem%ilvmax
    
      p_rmatrixStatic => rproblem%RlevelInfo(i)%rmatrixStatic
      p_rmatrixMass => rproblem%RlevelInfo(i)%rmatrixMass
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Create a sort strategy structure for our discretisation
      call sstrat_initBlockSorting (rproblem%RlevelInfo(i)%rsortStrategy,p_rdiscretisation)

      ! Calculate the resorting strategy.
      call sstrat_initCuthillMcKee (rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
          p_rmatrix%RmatrixBlock(1,1))
      
      ! Attach the sorting strategy to the matrices. The matrix is not yet sorted.
      call lsysbl_setSortStrategy (p_rmatrixStatic,&
          rproblem%RlevelInfo(i)%rsortStrategy,&
          rproblem%RlevelInfo(i)%rsortStrategy)

      call lsysbl_setSortStrategy (p_rmatrixMass,&
          rproblem%RlevelInfo(i)%rsortStrategy,&
          rproblem%RlevelInfo(i)%rsortStrategy)
      
      call lsysbl_setSortStrategy (p_rmatrix,&
          rproblem%RlevelInfo(i)%rsortStrategy,&
          rproblem%RlevelInfo(i)%rsortStrategy)

      ! Resort the system- and the mass matrices
      ! For the system matrix, we only resort the structure.
      ! The entries are not sorted, they are undefined anyway.
      call lsysbl_sortMatrix (p_rmatrixMass,.true.)
      call lsysbl_sortMatrix (p_rmatrix,.true.,.false.)

      ! Finally, sort the Laplace matrix.
      ! Be aware that this matrix is the "owner" of the structure, i.e.,
      ! it has to be sorted after the other matrices. Otherwise,
      ! The other matrices would try to resort an already resorted structure!
      call lsysbl_sortMatrix (p_rmatrixStatic,.true.)
      
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_calcRHS (rproblem,rrhs)
  
!<description>
  ! Calculates the RHS vector at the current point in time.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! A pointer to the RHS vector.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
    ! A linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Put the current simulation time as parameter "TIME" into the collection.
    ! Also set Dquickaccess (1) to the simulation time for faster access by the
    ! callback routine.
    rproblem%rcollection%Dquickaccess (1) = rproblem%rtimedependence%dtime
    call collct_setvalue_real(rproblem%rcollection,"TIME",&
         rproblem%rtimedependence%dtime,.true.)

    ! The vector structure is done but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! Get discretisation and cubature information structure
    p_rdiscretisation => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscretisation
    p_rcubatureInfo => rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo

    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    call linf_buildVectorScalar (&
        rlinform,.true.,rrhs%RvectorBlock(1),p_rcubatureInfo,&
        coeff_RHS,rproblem%rcollection)
    
    ! Remove the "TIME"-parameter from the collection again.
    call collct_deletevalue (rproblem%rcollection,"TIME")
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    integer :: ihandle,i

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue (rproblem%rcollection,"RHS")

    ! Release matrices and vectors on all levels
    do i=rproblem%ilvmax,rproblem%ilvmin,-1

      ! Delete the variables from the collection.
      call collct_deletevalue (rproblem%rcollection,"SYSTEM",i)
      call collct_deletevalue (rproblem%rcollection,"LAPLACE",i)
      call collct_deletevalue (rproblem%rcollection,"MASS",i)

      ! Delete the system matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)
      
      ! Delete the mass matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixMass)

      ! Delete the matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixStatic)

      ! Release the sorting strategy
      call sstrat_doneBlockSorting (rproblem%RlevelInfo(i)%rsortStrategy)
    end do

  end subroutine

end module
