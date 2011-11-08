!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2discretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic spatial discretisation related routines for
!# CC2D. Here, matrix and RHS creation routines can be found as well as
!# routines to initialise/clean up discretisation structures.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_initDiscretisation
!#     -> Initialise the discretisation structure inside of the problem
!#        structure using the parameters from the INI/DAT files.
!#
!# 2.) c2d2_initMatVec
!#     -> Basic initialisation of matrices/vectors
!#
!# 3.) c2d2_doneMatVec
!#     -> Cleanup of matrices/vectors
!#
!# 4.) c2d2_doneDiscretisation
!#     -> Cleanup of the underlying discretisation structures
!#
!# </purpose>
!##############################################################################

module cc2dminim2discretisation

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use scalarpde
  use derivatives
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  
  use collection
  use convection
    
  use cc2dminim2basic
  use cc2dmini_callback
  
  implicit none
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: I,ielementType,icubtemp
  integer(I32) :: icubA,icubB,icubF
  
  ! Number of equations in our problem. velocity+velocity+pressure = 3
  integer, parameter :: nequations = 3
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Which discretisation is to use?
    ! Which cubature formula should be used?
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iElementType',ielementType,3)
		
		icubtemp = CUB_G2X2
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'icubLaplace',icubtemp,icubtemp)
		icubA = icubtemp
		
		icubtemp = CUB_G2X2
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'icubB',icubtemp,icubtemp)
		icubB = icubtemp
                              
		icubtemp = CUB_G2X2
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'icubF',icubtemp,icubtemp)
		icubF = icubtemp

    ! Now set up discrezisation structures on all levels:

    do i=rproblem%NLMIN,rproblem%NLMAX
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector. In this simple problem, we only have one block.
      allocate(p_rdiscretisation)
      call spdiscr_initBlockDiscr (p_rdiscretisation,nequations,&
                                   p_rtriangulation, p_rboundary)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      select case (ielementType)
      case (3)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the
        ! velocity...
        call spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_EM30,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory,
        ! as both structures will share the same dynamic information afterwards.
        call spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        call spdiscr_deriveSimpleDiscrSc (p_rdiscretisation%RspatialDiscr(1),&
            EL_Q0,icubB,p_rdiscretisation%RspatialDiscr(3))
                    
      case DEFAULT
        print *,'Unknown discretisation: iElementType = ',ielementType
        stop
      end select

    end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Rremove the block discretisation structure and all the substructures.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! Remove the discretisation from the heap.
      deallocate(rproblem%RlevelInfo(i)%p_rdiscretisation)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_initMatVec (rproblem)
  
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
  integer :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector,p_rtempVector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    do i=rproblem%NLMIN,rproblem%NLMAX
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! The global system looks as follows:
      !
      !    ( A         B1 )
      !    (      A    B2 )
      !    ( B1^T B2^T    )
      !
      ! with A = L + nonlinear Convection. We compute in advance
      ! a standard Laplace matrix L which can be added later to the
      ! convection matrix, resulting in the nonlinear system matrix.
      !
      ! Get a pointer to the (scalar) Laplace matrix:
      p_rmatrixLaplace => rproblem%RlevelInfo(i)%rmatrixLaplace
      
      ! and save it to the collection for later use.
      call collct_setvalue_matsca(rproblem%rcollection,'LAPLACE',&
                                  p_rmatrixLaplace,.true.,i)
      
      ! Create the matrix structure of the Laplace matrix:
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixLaplace)
      
      ! And now to the entries of the matrix. For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = rproblem%dnu
      rform%Dcoefficients(2)  = rproblem%dnu

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Stokes for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      call bilf_buildMatrixScalar (rform,.true.,&
                                   p_rmatrixLaplace,coeff_Stokes,&
                                   rproblem%rcollection)
      
      ! In the global system, there are two coupling matrices B1 and B2.
      ! Both have the same structure.
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
                rproblem%RlevelInfo(i)%rmatrixB1,&
                p_rdiscretisation%RspatialDiscr(1))
                
      ! Duplicate the B1 matrix structure to the B2 matrix, so use
      ! lsyssc_duplicateMatrix to create B2. Share the matrix
      ! structure between B1 and B2 (B1 is the parent and B2 the child).
      ! Don't create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
                  rproblem%RlevelInfo(i)%rmatrixB2,LSYSSC_DUP_COPY,LSYSSC_DUP_REMOVE)

      ! Build the first pressure matrix B1.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,&
                                  rproblem%RlevelInfo(i)%rmatrixB1,coeff_Pressure,&
                                  rproblem%rcollection)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,&
                                  rproblem%RlevelInfo(i)%rmatrixB2,coeff_Pressure,&
                                  rproblem%rcollection)
                                  
      ! Now let's come to the main system matrix, which is a block matrix.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)
      
      ! Save the system matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,'SYSTEMMAT',p_rmatrix,.true.,i)

      ! Inform the matrix that we build a saddle-point problem.
      ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
      ! but probably some solvers can use the special structure later.
      p_rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Let's consider the global system in detail:
      !
      !    ( A         B1 ) = ( A11  A12  A13 )
      !    (      A    B2 )   ( A21  A22  A23 )
      !    ( B1^T B2^T    )   ( A31  A32  A33 )
      !
      ! The matrices A11 and A22 of the global system matrix have exactly
      ! the same structure as the original Laplace matrix from above!
      ! Initialise them with the same structure, i.e. A11, A22 and the
      ! Laplace matrix L share(!) the same structure.
      !
      ! For this purpose, use the "duplicate matric" routine.
      ! The structure of the matrix is shared with the Laplace matrix.
      ! For the content, a new empty array is allocated which will later receive
      ! the entries.
      call lsyssc_duplicateMatrix (p_rmatrixLaplace,&
                  p_rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                  
      ! The matrix A22 is identical to A11! So mirror A11 to A22 sharing the
      ! structure and the content.
      call lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                  p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      ! Manually change the discretisation structure of the Y-velocity
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscrDirectMat (p_rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))
                                  
      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The reason is
      ! that these matrices are modified for bondary conditions later.
      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(1,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(2,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(3,1),&
                                   LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(3,2),&
                                   LSYSSC_TR_VIRTUAL)

      ! Now on all levels except for the maximum one, create a temporary
      ! vector on that level, based on the matrix template.
      ! It's used for building the matrices on lower levels.
      if (i .lt. rproblem%NLMAX) then
        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
        call lsysbl_createVecBlockIndMat (p_rmatrix,p_rtempVector,.false.)
        
        ! Add the temp vector to the collection on level i
        ! for use in the callback routine
        call collct_setvalue_vec(rproblem%rcollection,'RTEMPVEC',p_rtempVector,&
                                .true.,i)
      end if

    end do

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    call lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .false.)
    call lsysbl_createVecBlockIndMat (p_rmatrix,p_rvector, .false.)

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    call collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.true.)
    call collct_setvalue_vec(rproblem%rcollection,'SOLUTION',p_rvector,.true.)
    
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
    !
    ! Note that the vector is unsorted after calling this routine!
    !
    ! Discretise the X-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
              p_rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
              p_rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(p_rrhs%RvectorBlock(3))
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    integer :: i

    ! Release matrices and vectors on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Delete the matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)

      ! Delete the variables from the collection.
      call collct_deletevalue (rproblem%rcollection,'SYSTEMMAT',i)
      call collct_deletevalue (rproblem%rcollection,'LAPLACE',i)
      
      ! Release Laplace, B1 and B2 matrix
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixLaplace)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      if (i .lt. rproblem%NLMAX) then
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
        call collct_deletevalue(rproblem%rcollection,'RTEMPVEC',i)
      end if
      
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rproblem%rvector)
    call lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue (rproblem%rcollection,'RHS')
    call collct_deletevalue (rproblem%rcollection,'SOLUTION')

  end subroutine

end module
