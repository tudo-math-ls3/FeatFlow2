!##############################################################################
!# ****************************************************************************
!# <name> vanca </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the implementations of the VANCA preconditioner.
!# These are more or less auxiliary routines called by the VANCA 
!# preconditioner in the linearsolver.f90 solver library.
!#
!# The following routines can be found here:
!#
!# 1.) vanca_initGeneralVanca
!#     -> Initialise the general Vanca solver
!#
!# 2.) vanca_general
!#     -> Perform one step of the full VANCA solver for general block systems.
!#
!# 3.) vanca_doneGeneralVanca
!#     -> Clean up the general VANCA solver
!#
!# 4.) vanca_init2DSPQ1TQ0simple
!#     -> Initialise specialised VANCA solver for 2D saddle point problems
!#        with $\tilde Q_1/Q_0$ discretisation
!#
!# 5.) vanca_2DSPQ1TQ0simple
!#     -> Perform one step of the specialised VANCA solver for 2D saddle point 
!#        problems with $\tilde Q_1/Q_0$ discretisation
!#
!# </purpose>
!##############################################################################

MODULE vanca

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock

  IMPLICIT NONE

!<types>
  
!<typeblock>
  
  ! A structure that accepts a pointer to the column/row/data arrays
  ! of a structure-7/structure 9 matrix. This is usually passed to
  ! the VANCA preconditioner(s) to specify the matrices to handle.
  
  TYPE t_matrixPointer79Vanca
    ! Is set to FALSE if the matrix does not exist/is empty.
    ! In this case, the pointers below are undefined!
    LOGICAL :: bexists
    
    ! TRUE if the matrix is saved transposed
    LOGICAL :: btransposed
    
    ! The scaling factor of the matrix; from the matrix structure.
    REAL(DP) :: dscaleFactor
  
    ! Pointer to the data - currently only double precision 
    REAL(DP), DIMENSION(:), POINTER :: p_DA
    
    ! Pointer to the column structure
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    
    ! Pointer to the row structure
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A structure for saving precalculated information for the general VANCA.
  ! This is initialised by vanca_initGeneralVanca and released by
  ! vanca_doneGeneralVanca.
  
  TYPE t_vancaGeneral
  
    ! Number of blocks in the global matrix
    INTEGER                                              :: nblocks
    
    ! Pointer to the block matrix
    TYPE(t_matrixBlock), POINTER                         :: p_rmatrix
  
    ! Pointers to t_matrixPointer79Vanca structures specifying
    ! the submatrices and their properties.
    TYPE(t_matrixPointer79Vanca), DIMENSION(:,:),POINTER :: p_Rmatrices
    
    ! Maximum number of local DOF's.
    INTEGER(PREC_DOFIDX)                              :: nmaxLocalDOFs
    
    ! Total number of local DOF's
    INTEGER                                           :: ndofsPerElement

    ! Number of local DOF's in the element distributions of all blocks.
    ! Note that this VANCA supports only uniform discretisations, so
    ! each entry corresponds to one block in the solution vector.
    INTEGER(PREC_DOFIDX), DIMENSION(LSYSBL_MAXBLOCKS) :: InDofsLocal
    
    ! Pointers to the spatial discretisation of all the blocks
    TYPE(t_spatialDiscretisation), DIMENSION(LSYSBL_MAXBLOCKS) :: p_Rdiscretisation
    
    ! Offset indices of the blocks in the solution vector. IblockOffset(i)
    ! points to the beginning of the i'th block of the solution vector.
    INTEGER(PREC_DOFIDX), DIMENSION(LSYSBL_MAXBLOCKS+1) :: IblockOffset
    
    ! Temporary array that saves the DOF's that are in processing when
    ! looping over an element set.
    ! DIMENSION(nmaxLocalDOFs,VANCA_NELEMSIM,nblocks)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), POINTER     :: IelementDOFs

  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D-Saddle-Point
  ! VANCA method with Q1~/Q0 discretisation support.
  
  TYPE t_vancaPointer2DSPQ1TQ0
    ! Pointer to the column structure of the velocity matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    
    ! Pointer to the row structure of the velocity matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA
    
    ! Pointer to diagonal entries in the velocity matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA

    ! Pointer to the matrix entries of the velocity matrix.
    REAL(DP), DIMENSION(:), POINTER             :: p_DA

    ! Pointer to the column structure of the B/D-matrices.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    
    ! Pointer to the row structure of the B/D-matrices
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    
    ! Pointer to the entries of the B1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1

    ! Pointer to the entries of the B2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    
    ! Pointer to the entries of the D1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1

    ! Pointer to the entries of the D2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
  END TYPE
  
!</typeblock>


!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D-Saddle-Point
  ! VANCA method with Q2/QP1 discretisation support.
  
  TYPE t_vancaPointer2DSPQ2QP1
    ! Pointer to the column structure of the velocity matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    
    ! Pointer to the row structure of the velocity matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA
    
    ! Pointer to diagonal entries in the velocity matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalA

    ! Pointer to the matrix entries of the velocity matrix.
    REAL(DP), DIMENSION(:), POINTER             :: p_DA

    ! Pointer to the column structure of the B/D-matrices.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    
    ! Pointer to the row structure of the B/D-matrices
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    
    ! Pointer to the entries of the B1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1

    ! Pointer to the entries of the B2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    
    ! Pointer to the entries of the D1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1

    ! Pointer to the entries of the D2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
  END TYPE
  
!</typeblock>


!</types>

!<constants>
!<constantblock description="Constants defining the blocking of element sets in VANCA">

  ! Number of elements to handle simultaneously inb general VANCA
  INTEGER :: VANCA_NELEMSIM   = 1000
  
!</constantblock>
!</constants>

CONTAINS

  ! ***************************************************************************
  ! GENERAL VANCA
  ! Supports (non-transposed) matrices of any kind.
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_initGeneralVanca (rmatrix,rvancaGeneral)

!<description>
  ! This routine initialises the general VANCA solver and allocates
  ! necessary memory for the iteration.
!</description>

!<input>
  ! The system matrix that is used during the VANCA iteration.
  ! Remark: VANCA saves a pointer to this matrix, so the structure must exist
  !  until the system is solved! (Usually this points to the system matrix in
  !  the corresponding solver structure of the underlying linear solver...)
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: rmatrix
!</input>
  
!<output>
  ! VANCA spiecific structure. Contains internal data and allocated memory.
  TYPE(t_vancaGeneral), INTENT(OUT)       :: rvancaGeneral
!</output>

!</subroutine>

    ! local variables
    LOGICAL :: bfirst
    INTEGER :: nblocks,i,j,nmaxLocalDOFs,ndofsPerElement
    TYPE(t_spatialDiscretisation), POINTER            :: p_rdiscretisation
    
    nblocks = rmatrix%ndiagBlocks
    nmaxLocalDOFs = 0
    ndofsPerElement = 0
    
    ! Allocate memory for the matrix structures
    rvancaGeneral%nblocks = nblocks
    ALLOCATE(rvancaGeneral%p_Rmatrices(nblocks,nblocks))
    
    ! This type of VANCA only supports a uniform discretisation
    ! and matrix format 7 or 9. Transposed matrices are not allowed.
    !
    ! Check all matrices that this is the case.
    ! Build the Rmatrices structure array. We manually create a
    ! (nblock,nblock) array inside of a 1-dimensional array and cast a rank
    ! change of the array on call to the actual VANCA subroutine later.
    !
    ! Loop through the columns of the block matrix.
    !
    ! Offset position of the first block is = 0.
    rvancaGeneral%IblockOffset(1) = 0
    
    DO i=1,nblocks
    
      ! Note this block as 'not processed'
      bfirst = .TRUE.
      
      ! Loop through the rows of the current matrix column.
      DO j=1,nblocks
        IF (rmatrix%RmatrixBlock(j,i)%NCOLS .NE. 0) THEN
          ! Get a/the discretisation structure of the current block/matrix column
          p_rdiscretisation => rmatrix%RmatrixBlock(j,i)%p_rspatialDiscretisation
          
          IF (p_rdiscretisation%ccomplexity .NE. &
              SPDISC_UNIFORM) THEN
            PRINT *,'General VANCA supports only uniform triangulation!'
            STOP
          END IF
          
          IF ((rmatrix%RmatrixBlock(j,i)%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND. &
              (rmatrix%RmatrixBlock(j,i)%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
            PRINT *,'General VANCA supports only matrix structure 7 and 9!'
            STOP
          END IF
          
          rvancaGeneral%p_Rmatrices(j,i)%bexists = .TRUE.
          rvancaGeneral%p_Rmatrices(j,i)%dscaleFactor = &
            rmatrix%RmatrixBlock(j,i)%dscaleFactor
          rvancaGeneral%p_Rmatrices(j,i)%btransposed = &
            IAND(rmatrix%RmatrixBlock(j,i)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0
            
          CALL storage_getbase_double (rmatrix%RmatrixBlock(j,i)%h_Da,&
                                       rvancaGeneral%p_Rmatrices(j,i)%p_DA)
          CALL storage_getbase_int (rmatrix%RmatrixBlock(j,i)%h_Kcol,&
                                    rvancaGeneral%p_Rmatrices(j,i)%p_Kcol)
          CALL storage_getbase_int (rmatrix%RmatrixBlock(j,i)%h_Kld,&
                                    rvancaGeneral%p_Rmatrices(j,i)%p_Kld)
                                    
          IF (bfirst) THEN
            ! This block has not yet been processed.
            !
            ! Get the NEQ of the current block and save it as offset position
            ! in the global vector for the next block.
            rvancaGeneral%IblockOffset(i+1) = &
              rvancaGeneral%IblockOffset(i) + rmatrix%RmatrixBlock(j,i)%NCOLS

            ! We need some information for calculating DOF's later.
            ! Get the number of local DOF's in the current block.
            ! Note that we restrict to uniform discretisations!
            rvancaGeneral%InDofsLocal(i) = elem_igetNDofLoc(p_rdiscretisation% &
                                                RelementDistribution(1)%itrialElement)
            
            ! Calculate the maximum number of local DOF's
            nmaxLocalDOFs = MAX(nmaxLocalDOFs,rvancaGeneral%InDofsLocal(i))
            
            ! Calculate the total number of local DOF's
            ndofsPerElement = ndofsPerElement + rvancaGeneral%InDofsLocal(i)
            
            bfirst = .FALSE.
          
          END IF
          
        ELSE
          rvancaGeneral%p_Rmatrices(j,i)%bexists = .FALSE.
        END IF
        
      END DO
    END DO

    ! Save the max. and total number of local DOF's
    rvancaGeneral%nmaxLocalDOFs = nmaxLocalDOFs
    rvancaGeneral%ndofsPerElement = ndofsPerElement
    
    ! We know the maximum number of DOF's now. For the later over the elements,
    ! allocate memory for storing the DOF's of an element set.
    ALLOCATE(rvancaGeneral%IelementDOFs(nmaxLocalDOFs,VANCA_NELEMSIM,nblocks))
    
    ! Remember the matrix
    rvancaGeneral%p_rmatrix => rmatrix
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_doneGeneralVanca (rvancaGeneral)

!<description>
  ! This routine cleans up a general VANCA solver. All memory allocated in
  ! rvancaGeneral is released.
!</description>
  
!<inputoutput>
  ! The general-VANCA structure to be cleaned up.
  TYPE(t_vancaGeneral), INTENT(INOUT)       :: rvancaGeneral
!</inputoutput>

!</subroutine>

    ! Release memory allocated in the init-routine
    
    IF (ASSOCIATED(rvancaGeneral%IelementDOFs)) &
      DEALLOCATE(rvancaGeneral%IelementDOFs)
    IF (ASSOCIATED(rvancaGeneral%p_Rmatrices)) &
      DEALLOCATE(rvancaGeneral%p_Rmatrices)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_general (rvancaGeneral, rvector, rrhs, domega)

!<description>
  ! This routine applies the general-VANCA algorithm to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvancaGeneral structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix A.
!</description>

!<input>
  
  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The general-VANCA structure. Must have been initialised with 
  ! vanca_initGeneralVanca before.
  TYPE(t_vancaGeneral), INTENT(INOUT)     :: rvancaGeneral

  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,j
  INTEGER(PREC_ELEMENTIDX) :: IELmax, IELset, iel
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementList
  REAL(DP), DIMENSION(:), POINTER                 :: p_Drhs,p_Dvector
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER     :: p_Ipermutation
    
  ! Saved matrix and the vector(s) must be compatible!
  CALL lsysbl_isMatrixCompatible(rvector,rvancaGeneral%p_rmatrix)
  CALL lsysbl_isVectorCompatible(rvector,rrhs)
    
  ! Get the data arrays of the vector/rhs
  CALL lsysbl_getbase_double (rvector,p_Dvector)
  CALL lsysbl_getbase_double (rrhs,p_Drhs)
  
  ! p_IelementList must point to our set of elements in the discretisation
  ! with that combination of trial/test functions.
  CALL storage_getbase_int (rvector%RvectorBlock(1)%p_rspatialDiscretisation% &
                            RelementDistribution(1)%h_IelementList, &
                            p_IelementList)
    
  ! Loop over the elements - blockwise.
  DO IELset = 1, rvector%RvectorBlock(1)%p_rspatialDiscretisation%&
                         p_rtriangulation%NEL, VANCA_NELEMSIM
  
    ! We always handle LINF_NELEMSIM elements simultaneously.
    ! How many elements have we actually here?
    ! Get the maximum element number, such that we handle at most VANCA_NELEMSIM
    ! elements simultaneously.
    IELmax = MIN(rvector%RvectorBlock(1)%p_rspatialDiscretisation%&
                        p_rtriangulation%NEL,IELset-1+VANCA_NELEMSIM)
  
    ! Loop over the blocks in the block vector to get the DOF's everywhere.
    
    DO i=1,rvector%nblocks

      ! Calculate the global DOF's of all blocks.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our VANCA_NELEMSIM elements simultaneously.
      CALL dof_locGlobMapping_mult(rvector%RvectorBlock(i)%p_rspatialDiscretisation,&
                                   p_IelementList(IELset:IELmax), &
                                   .FALSE.,rvancaGeneral%IelementDOFs(:,:,i))

      ! If the vector is sorted, push the DOF's through the permutation to get
      ! the actual DOF's.
      IF (rvector%RvectorBlock(i)%isortStrategy .GT. 0) THEN
      
        CALL storage_getbase_int(rvector%RvectorBlock(i)%h_IsortPermutation,&
                                 p_Ipermutation)

        DO iel=1,IELmax-IELset+1
          DO j=1,rvancaGeneral%InDofsLocal(i)
            ! We are not resorting the vector but determining the 'sorted'
            ! DOF's - this needs the 2nd half of the permutation.
            rvancaGeneral%IelementDOFs(j,iel,i) = &
              p_Ipermutation(rvancaGeneral%IelementDOFs(j,iel,i)+rvector%RvectorBlock(i)%NEQ)
          END DO
        END DO
      END IF
    
    END DO  
  
    ! Now, IdofsTotal contains all DOF's on each element, over all discretisations.
    !
    ! Call the actual VANCA to process the DOF's on each element.
    CALL vanca_general_double_mat79 (p_Dvector, p_Drhs, domega, &
         rvancaGeneral%p_Rmatrices,IELmax-IELset+1_PREC_ELEMENTIDX,&
         rvancaGeneral%IblockOffset,rvancaGeneral%nblocks,&
         rvancaGeneral%InDofsLocal,rvancaGeneral%ndofsPerElement,&
         rvancaGeneral%IelementDOFs)
                 
  END DO
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE vanca_general_double_mat79 (Dvector, Drhs, domega, Rmatrices,&
             nelements,IblockOffset,nblocks,InDofsLocal,ndofsPerElement,&
             IelementDofs)

!<description>
  ! This routine applies one step of the VANCA solver (local block 
  ! gauss-seidel) to a given set of solution/RHS vector. 
  ! All given matrices are assumed to be format 7/9,
  ! double precision. The given vector is assumed to be double precision.
!</description>

!<input>
  ! The (block) RHS vector, given as one large array.
  REAL(DP), DIMENSION(:), INTENT(IN)             :: Drhs
  
  ! A relaxation parameter. Standard = 1.0_DP.
  REAL(DP), INTENT(IN)                           :: domega
  
  ! A list of matrices to handle; directly specified by pointers
  ! to the substructures (data/columns/rows).
  TYPE(t_matrixPointer79Vanca), DIMENSION(:,:),&
                                INTENT(IN)       :: Rmatrices

  ! Number of elements that should be processed in this sweep.
  INTEGER(PREC_ELEMENTIDX), INTENT(IN)           :: nelements
  
  ! Number of blocks in the vectors
  INTEGER, INTENT(IN)                            :: nblocks
  
  ! Offset position of the blocks in the vector.
  ! Block i starts at position IblockOffset(i)+1 in Dvector / Drhs.
  ! IblockOffset(nblocks+1) gives the number of equations in Dvector/Drhs.
  INTEGER(PREC_DOFIDX), DIMENSION(nblocks+1), INTENT(IN) :: IblockOffset
  
  ! Number of local DOF's in each block.
  INTEGER(PREC_DOFIDX), DIMENSION(nblocks), INTENT(IN)   :: InDofsLocal
  
  ! Total number of local DOF's per element
  INTEGER, INTENT(IN)                                    :: ndofsPerElement
  
  ! List of DOF's on every element for every block.
  ! DIMENSION(nmaxDOFs,nelements,nblocks)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), INTENT(IN)     :: IelementDOFs
  
!</input>

!<inputoutput>
  ! The initial (block) solution vector. Is overwritten by the new (block) 
  ! solution vector.
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dvector
!</inputoutput>

!</subroutine>

    ! One iteration of vanka smother (block gauss-seidel) on the system
    !
    !    A11 A12 ... A1n | U1    F1
    !    A21 A22 ... A2n | U2  = F2
    !     :   :  ...  :  |  :     :
    !    An1 An2 ... Ann | Un  = Fn
    !
    ! Let's first describe the method used here.
    !
    ! The above block system can be written as a general block system of 
    ! the form
    !
    !             A x = f
    !
    ! Consider a general defect correction approach for this system:
    !
    !     x_{n+1}  =  x_n  +  \omega C^{-1} (b - A x_n)
    !
    ! With some damping parameter \omega and some preconditioner C^{-1}.
    ! Normally, this iteration is used globally - and the usual linear
    ! solver that calls this routine usually realises this approach.
    ! VANCA now applies this defect correction loop *locally*, i.e.
    ! not for the full system at once.
    !
    ! This local approach is based on a geometric point of view.
    ! In general, one could imagine a set of connected cells in the 
    ! global domain \Omega where to apply the algorithm to (the LMPSC 
    ! approach), but for simplicity, consider only one cell. Again for
    ! simplicity, imagine that our FE-spaces is Q1~/Q0. 
    !
    ! We loop over each cell in the domain, one after the other, and
    ! change the DOF's in the solution vector there:
    !
    ! We fetch all the data (e.g. velocity, pressure) on that cell. On the
    ! first cell, we have only "old" velocity entries. These values
    ! are updated and then the calculation proceeds with the 2nd cell.
    !
    !        old                      new
    !     |---X---|                |---X---|
    !     |       |                |       |
    ! old X   1   X old   -->  new X   1   X new
    !     |       |                |       |
    !     |---X---|                |---X---|
    !        old                      new
    !
    ! From the second cell on, there might be "old" data and "new" 
    ! data on that cell - the old data that has not been updated and
    ! perhaps some already updated velocity data from a neighbor cell.
    ! 
    !        new     old                   new     new
    !     |---X---|---X---|             |---X---|---X---|
    !     |       |       |             |       |       |
    ! new X   1   X   2   X old --> new X   1   X   2   X new
    !     |       |new    |             |       |newer  |
    !     |---X---|---X---|             |---X---|---X---|
    !        new     old                   new     new
    !
    ! These values are updated and then the calculation proceeds
    ! with the next cell.
    ! As can be seen in the above picture, the "new" node in the
    ! middle is even going to be a "newer" node when handled again
    ! for the 2nd cell. This is meant by "Gauss-Seldel" character:
    ! Information is updated subsequently by using "old" data and
    ! "new" data from a previous calculation.
    !
    ! So how to do this in detail?
    !
    ! Again consider the problem:
    !
    !             A x = f
    !
    ! We assume, that all components in the vector x are      
    ! given - except for the unknowns current element; these unknowns
    ! are located anywhere in the vector x. The idea is to    
    ! shift "everything known" to the right hand side to obtain   
    ! a system for only these unknowns!  
    !                         
    ! We extract all the lines of the system that correspond to
    ! our 'unknown' DOF's on our element; this results in a rectangular      
    ! system of the form                                            
    !                                                               
    !    [ === A~ === ] x  = (f~)
    !
    ! So #rows(A~)=#rows(f~)=#DOF's on the element! Furthermore we set
    ! Now we make a defect-correction approach for this system:
    !
    !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
    !                                     -----------
    !                                        =d~
    !
    ! Here the 'projection' operator simply converts the small
    ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
    ! of the same size as x - what is easy using the number of
    ! the DOF's on the element.
    !
    ! The only question now will be: What is C^{-1}?
    !
    ! Well, here we have different choices. Everything depends on the
    ! matrix A~, which is unfortunately rectangular: It's a
    ! (#DOF's on the element, #DOF's in the space) matrix - so
    ! in case of a Q1~/Q0 discretisation, it's a (5,NEQ) matrix!
    !
    ! For full linear systems, one would choose C=A, which ist the
    ! theoretically best preconditioner. What we do here is simply
    ! extracting all columns of A~ that correspond to the DOF's
    ! on the current element: 
    !
    !   C:=delete columns of A~ that don't belong to DOF's on the element
    !
    ! This then leads to a square preconditioner C^{-1} - and that's the
    ! full method, because C^{-1} can be applied directly using Lapack e.g.!
    !
    !
    ! So all in all we write our algorithm in short:
    !
    ! loop over all elements in the given element set
    !   extract local matrix aa
    !   extract the local rhs ff
    !   compute local residuum: ff := ff - aa * x
    !   solve: xx := aa^(-1) ff
    !   update global solution vector: x := x + omega*xx
    ! end loop

    ! local variables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER, DIMENSION(nblocks+1) :: IlocalIndex
    INTEGER(PREC_DOFIDX), DIMENSION(ndofsPerElement) :: IlocalDOF,IglobalDOF
    INTEGER :: i,j,k,iidx,iminiDOF
    INTEGER(PREC_VECIDX) :: irow,idof
    INTEGER(PREC_MATIDX) :: icol
    
    ! Memory for our local system; let's hope it's not too big :)
    REAL(DP), DIMENSION(ndofsPerElement,ndofsPerElement) :: Daa
    REAL(DP), DIMENSION(ndofsPerElement)                 :: Dff
    REAL(DP) :: dscale
    INTEGER(I32), DIMENSION(ndofsPerElement) :: Ipiv
    INTEGER(I32) :: iinfo
    
    ! Quickly check the matrices if one of them is saved transposed.
    ! VANCA does not support transposed matrices; would kill computational
    ! time, as we cannot extract columns from a structure-9 matrix with
    ! reasonable effort!
    DO i=1,SIZE(Rmatrices,1)
      DO j=1,SIZE(Rmatrices,2)
        IF (Rmatrices(i,j)%bexists .AND. Rmatrices(i,j)%btransposed) THEN
          PRINT *,'General VANCA does not support transposed matrices!'
          STOP
        END IF 
      END DO ! j
    END DO ! i
        
    ! Build an index pointer for accessing local DOF's
    IlocalIndex(1) = 0
    DO i=2,nblocks+1
      IlocalIndex(i) = IlocalIndex(i-1)+InDOFsLocal(i-1)
    END DO
        
    ! Ok, let's start with the loop over the elements in the given
    ! element list.
    DO iel = 1,nelements
    
      ! IelementDOFs (.,iel,.) gives now for every block in the system
      ! the DOF's on this element.
      !
      ! First copy the RHS entries of f to f~.
      
      iidx = 1
      DO i=1,nblocks
        DO j=1,InDOFsLocal(i)
          iidx = IlocalIndex(i)+j
          
          ! Get the DOF on the element:
          idof = IelementDOFs(j,iel,i)
          
          ! Number of DOF relative to this block
          IlocalDOF(iidx) = idof
          
          ! Calculate corresponding global DOF:
          IglobalDOF(iidx) = idof+IblockOffset(i)
          
          Dff(iidx) = Drhs(IglobalDOF(iidx))
        END DO
      END DO
      
      ! Compute  ff := ff - A x  for the local unknowns
      ! to build the local residuum Dff = f~-A~x.
      !
      ! Simultaneously extract the local matrix into one array Daa(:,:).
      ! But first clear the matrix - maybe that there are unused subblocks!
      Daa = 0
      
      ! Loop through the rows of the block matrix
      DO i=1,SIZE(Rmatrices,1)
        
        ! Loop through the columns of the block matrix
        DO j=1,SIZE(Rmatrices,2)
        
          ! Is there a matrix saved at this position?
          IF (Rmatrices(i,j)%bexists) THEN

            dscale = Rmatrices(i,j)%dscaleFactor
        
            ! Block column j in the global matrix corresponds to
            ! block column j in the small local matrix we have to fill with data.
            !
            !   ....      :           
            !   ....      :           a11 a12         a15
            !                         a21 a22         a25
            !        .... :   ==>             a33 a34 a35
            !        .... :                   a43 a44 a45
            !   .... ....             a51 a52 a53 a54
            !   ==== ==== =           ======= ======= ===
            !    1    2   3 corr. to     1       2     3
            !
            !     n x n-mat.          5x5-mat or so
            !
            ! From IlocalIndex, get the starting address of the j'th block
            ! in the local solution vector. In the above example, for j=2 e.g.
            ! this gives the starting address of the two global DOF's
            ! that correspond to the columns 3 and 4 in the local matrix.
            iidx = IlocalIndex(j)

            ! Loop through the DOF's that correspond to this block:
            DO irow = IlocalIndex(i)+1,IlocalIndex(i+1)
              
              ! Get the actual DOF, relative to this block.
              ! This is the row of the matrix that must be multiplied by x
              ! and be subtracted from the local RHS.
              idof = IlocalDOF(irow)
              
              ! Loop through the row of the matrix to its contribution to "b-Ax".
              DO k = Rmatrices(i,j)%p_Kld(idof) , Rmatrices(i,j)%p_Kld(idof+1)-1

                ! Get the column number in the global matrix. This gives
                ! the global DOF = number of the element in x that must be multiplied
                ! with that matrix entry.             
                icol = Rmatrices(i,j)%p_Kcol(k)+IblockOffset(j)

                ! Build the defect
                Dff(irow) = Dff(irow) &
                          - dscale * Rmatrices(i,j)%p_Da(k) * Dvector(icol)

                ! icol is the number of a DOF.
                ! Check if this DOF belongs to the DOF's we have to
                ! extract to our local system.
                ! Loop through the DOF's corresponding to column j of the block system.
                ! In the above example, this checks only the two DOF's corresponding
                ! to a?3 and a?4.
                DO iminiDOF = 1,InDOFsLocal(j)
                  IF (icol .EQ. IglobalDOF(iidx+iminiDOF)) THEN
                    ! Yes. Get the matrix entry, write it to the local matrix
                    Daa(irow,iidx+iminiDOF) = Rmatrices(i,j)%p_Da(k)
                    EXIT
                  END IF
                END DO ! idof
                
              END DO ! k
            
            END DO ! irow
                        
          END IF ! exists
        
        END DO ! j
        
      END DO ! i
      
      ! Ok, we now have our local matrix and our local defect vector d~.
      ! Apply LAPACK to the local system to solve C^{-1} d~.
      
      CALL DGETRF( ndofsPerElement, ndofsPerElement, Daa, ndofsPerElement, &
                   Ipiv, iinfo )
                  
      ! Note: It may happen that the matrix is singular!
      !  That is the case if all DOF's are Dirichlet-DOF's - for example
      !  if the element is completely inside of a rigid fictitious boundary
      !  object.
      ! What to do in this case? Nothing! Ignore the system!
      ! Why? 
      !  - The values for the 'velocity' DOF's are fixed, so it's no use
      !    to try to calculate them.
      !  - If there's a zero-block involved in a saddle-point problem,
      !    the pressure (that corresponds to the zero-block) is not
      !    connected to the velocity - it's completely free!
      !    By ignoring this system, we let the pressure as it is.
      !
      ! One can theoretically also solve a least-squares system with
      ! LAPACK. This calculate an y with ||f~ - C y|| -> min and |y|->min.
      ! But this is (because of the 0-block in A and therefore also in C
      ! and because of the unit-vectors in the other rows of C)
      ! exactly the case if the 'velocity' y matches f~ and the 'pressure'
      ! components are zero - which means nothing else than that there's
      ! no contribution in the preconditioned defect to correct the
      ! pressure entries of x.
                   
      IF (iinfo .EQ. 0) THEN
        CALL DGETRS('N', ndofsPerElement, 1, Daa, ndofsPerElement, &
                    Ipiv, Dff, ndofsPerElement, iinfo )
        IF (iinfo .EQ. 0) THEN
          
          ! Dff is in-situ replaced by the solution - the preconditioned
          ! defect. Add this to our original solution vector to complete
          ! the 'local' defect-correction.
          
          DO i=1,ndofsPerElement
            j = IglobalDOF(i)
            Dvector(j) = Dvector(j) + domega*Dff(i)
          END DO
        
        END IF
      END IF
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D Saddle-Point VANCA, simple diagonal version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_init2DSPQ1TQ0simple (rmatrix,rvanca)
  
!<description>
  ! Checks if the "2D-Saddle-Point-Q1T-Q0" VANCA variant can be applied to
  ! the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>

!<output>
  ! t_vancaPointer2DSPQ1TQ0 structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DSPQ1TQ0), INTENT(OUT) :: rvanca
!</output>

!</subroutine>

    INTEGER :: i,j

    ! Matrix must be 3x3.
    IF (rmatrix%ndiagBlocks .NE. 3) THEN
      PRINT *,'vanca_check2DSPQ1TQ0: System matrix is not 3x3.'
      STOP
    END IF
    
    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(3,:) must be (virtually) transposed. All matrices must be double precision.
    DO i=1,3
      DO j=1,3
      
        IF (rmatrix%RmatrixBlock(i,j)%NEQ .NE. 0) THEN
        
          IF (i .LE. 2) THEN
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .NE. 0) THEN
              PRINT *,'vanca_check2DSPQ1TQ0: Transposed submatrices not supported.'
              STOP
            END IF
          ELSE
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .EQ. 0) THEN
              PRINT *,'vanca_check2DSPQ1TQ0: B1/B2 submatrices must be virtually'
              PRINT *,'transposed (LSYSSC_MSPEC_TRANSPOSED)!'
              STOP
            END IF
          END IF
          
          IF ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX7) .AND. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
            PRINT *,'vanca_check2DSPQ1TQ0: Only format 7 and 9 matrices supported.'
            STOP
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            PRINT *,'vanca_check2DSPQ1TQ0: Only double precision matrices supported.'
            STOP
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) THEN
            PRINT *,'vanca_check2DSPQ1TQ0: Scaled matrices supported.'
            STOP
          END IF
          
        END IF ! neq != 0
      END DO
    END DO
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    IF ((rmatrix%RmatrixBlock(1,3)%NA .NE. rmatrix%RmatrixBlock(3,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .NE. rmatrix%RmatrixBlock(3,1)%NCOLS)) THEN
      PRINT *,'vanca_check2DSPQ1TQ0: Structure of B1 and B1^T different!'
      STOP
    END IF

    IF ((rmatrix%RmatrixBlock(2,3)%NA .NE. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .NE. rmatrix%RmatrixBlock(3,2)%NCOLS)) THEN
      PRINT *,'vanca_check2DSPQ1TQ0: Structure of B2 and B2^T different!'
      STOP
    END IF
    
    ! Fill the output structure with data of the matrices.
    CALL storage_getbase_double(rmatrix%RmatrixBlock(1,1)%h_Da,rvanca%p_DA )
    CALL storage_getbase_double(rmatrix%RmatrixBlock(1,3)%h_Da,rvanca%p_DB1)
    CALL storage_getbase_double(rmatrix%RmatrixBlock(2,3)%h_Da,rvanca%p_DB2)
    CALL storage_getbase_double(rmatrix%RmatrixBlock(3,1)%h_Da,rvanca%p_DD1)
    CALL storage_getbase_double(rmatrix%RmatrixBlock(3,2)%h_Da,rvanca%p_DD2)
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,3)%h_Kcol,rvanca%p_KcolB)
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,3)%h_Kld, rvanca%p_KldB )
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,1)%h_Kcol,rvanca%p_KcolA)
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,1)%h_Kld, rvanca%p_KldA )
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL storage_getbase_int(rmatrix%RmatrixBlock(1,1)%h_Kdiagonal, &
                               rvanca%p_KdiagonalA)
    ELSE
      rvanca%p_KdiagonalA => rvanca%p_KldA
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ1TQ0simple (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Saddle-Point problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vancaPointer2DSPQ1TQ0 structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DSPQ1TQ0), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_POINTIDX)   :: NVT
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8
    REAL(DP) :: daux
    
    ! Local arrays for informations about one element
    REAL(DP), DIMENSION(4) :: AA,BB1,BB2,DD1,DD2
    REAL(DP), DIMENSION(9) :: FF,UU
    INTEGER(PREC_VECIDX), DIMENSION(4) :: idofGlobal
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NVT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscretisation% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    DO iel=1,NEL
    
      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF's    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      DO inode=1,4
      
        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)-NVT

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof
        
        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the four velocity unknowns and the         
        ! pressure unknown on the current element; these five unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity      
        ! components and B being an (2*4) x 1 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |        |             |        |        |                
        !   --|---X----|--           |--------|--------|                
        !                                                               
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to  
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have       
        ! 2 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        DO ib = ib1,ib2
          IF (p_KcolB(ib) .EQ. IEL) THEN
          
            J = p_KcolB(ib)
          
            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            
            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            
            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)
          
            ! Quit the loop - the other possible entry belongs to another 
            ! element, not to the current one
            EXIT
          END IF
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.
      
      CALL vanca_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)
    
      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,4
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ************************************************************************

!<subroutine>

  PURE SUBROUTINE vanca_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)
  
!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two 
  ! velocity vectors and one pressure vector. It's used as auxiliary 
  ! routine in the simple VANCA solver to calculate an update vector
  ! for velocity/pressure.
!</description>

!<input>
  ! Diagonal elements of the local system matrix.
  REAL(DP), DIMENSION(4), INTENT(IN) :: AA
  
  ! Entries in the submatrix B1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB1

  ! Entries in the submatrix B2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: BB2
  
  ! Entries in the submatrix D1.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD1

  ! Entries in the submatrix D2.
  REAL(DP), DIMENSION(4), INTENT(IN) :: DD2

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(IN) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  REAL(DP), DIMENSION(9), INTENT(OUT) :: UU
!</output>

!</subroutine>

    ! local variables

    INTEGER :: inode
    REAL(DP) :: PP,dpres
    REAL(DP), DIMENSION(9) :: AI,dff
    
    INTEGER, PARAMETER :: lofsv = 4
    INTEGER, PARAMETER :: lofsp = 8

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA(1)                                                   BB1(1) )
    !     (        AA(2)                                            BB1(2) )
    !     (               AA(3)                                     BB1(3) )
    !     (                      AA(4)                              BB1(4) )
    !     (                             AA(1)                       BB2(1) )
    !     (                                    AA(2)                BB2(2) )
    !     (                                           AA(3)         BB2(3) )
    !     (                                                  AA(4)  BB2(4) )
    !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the 
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial! 
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity 
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    DO inode=1,4
    
      ! Quick check if everything is ok - we don't want to divide by 0.
      IF (AA(inode)*AA(inode) .LT. 1E-20_DP) THEN
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        RETURN
      END IF

      ! AI(.) saves the diagonal matrix S^-1:
    
      AI(inode)=1E0_DP/AA(inode)
        
    END DO

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 4x1 matrix 
    ! b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = 0.0_DP
    dff = FF
      
    DO inode = 1,4
      dpres        = dpres &
                   - AI(inode)*(DD1(inode)*BB1(inode)+DD2(inode)*BB2(inode))
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI(inode)*(DD1(inode)*dff(inode)+DD2(inode)*dff(inode+lofsv))
    END DO

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF's are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    IF (dpres*dpres .LT. 1E-20_DP)  THEN
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      RETURN
    ENDIF
      
    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:
      
    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP
      
    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:
      
    DO inode=1,4
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  ! 2D Saddle-Point VANCA, simple diagonal and full version.
  ! Supports Q2/QP1 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_init2DSPQ2QP1 (rmatrix,rvanca)
  
!<description>
  ! Checks if the "2D-Saddle-Point-Q2-QP1" VANCA variant can be applied to
  ! the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>

!<output>
  ! t_vancaPointer2DSPQ1TQ0 structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DSPQ2QP1), INTENT(OUT) :: rvanca
!</output>

!</subroutine>

    INTEGER :: i,j

    ! Matrix must be 3x3.
    IF (rmatrix%ndiagBlocks .NE. 3) THEN
      PRINT *,'vanca_check2DSPQ1TQ0: System matrix is not 3x3.'
      STOP
    END IF
    
    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(3,:) must be (virtually) transposed. All matrices must be double precision.
    DO i=1,3
      DO j=1,3
      
        IF (rmatrix%RmatrixBlock(i,j)%NEQ .NE. 0) THEN
        
          IF (i .LE. 2) THEN
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .NE. 0) THEN
              PRINT *,'vanca_init2DSPQ2QP1simple: Transposed submatrices not supported.'
              STOP
            END IF
          ELSE
            IF (IAND(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .EQ. 0) THEN
              PRINT *,'vanca_init2DSPQ2QP1simple: B1/B2 submatrices must be virtually'
              PRINT *,'transposed (LSYSSC_MSPEC_TRANSPOSED)!'
              STOP
            END IF
          END IF
          
          IF ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX7) .AND. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
            PRINT *,'vanca_init2DSPQ2QP1simple: Only format 7 and 9 matrices supported.'
            STOP
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            PRINT *,'vanca_init2DSPQ2QP1simple: Only double precision matrices supported.'
            STOP
          END IF

          IF (rmatrix%RmatrixBlock(i,j)%dscaleFactor .NE. 1.0_DP) THEN
            PRINT *,'vanca_init2DSPQ2QP1simple: Scaled matrices supported.'
            STOP
          END IF
          
        END IF ! neq != 0
      END DO
    END DO
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    IF ((rmatrix%RmatrixBlock(1,3)%NA .NE. rmatrix%RmatrixBlock(3,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .NE. rmatrix%RmatrixBlock(3,1)%NCOLS)) THEN
      PRINT *,'vanca_init2DSPQ2QP1simple: Structure of B1 and B1^T different!'
      STOP
    END IF

    IF ((rmatrix%RmatrixBlock(2,3)%NA .NE. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .NE. rmatrix%RmatrixBlock(3,2)%NCOLS)) THEN
      PRINT *,'vanca_init2DSPQ2QP1simple: Structure of B2 and B2^T different!'
      STOP
    END IF
    
    ! Fill the output structure with data of the matrices.
    CALL storage_getbase_double(rmatrix%RmatrixBlock(1,1)%h_Da,rvanca%p_DA )
    CALL storage_getbase_double(rmatrix%RmatrixBlock(1,3)%h_Da,rvanca%p_DB1)
    CALL storage_getbase_double(rmatrix%RmatrixBlock(2,3)%h_Da,rvanca%p_DB2)
    CALL storage_getbase_double(rmatrix%RmatrixBlock(3,1)%h_Da,rvanca%p_DD1)
    CALL storage_getbase_double(rmatrix%RmatrixBlock(3,2)%h_Da,rvanca%p_DD2)
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,3)%h_Kcol,rvanca%p_KcolB)
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,3)%h_Kld, rvanca%p_KldB )
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,1)%h_Kcol,rvanca%p_KcolA)
    CALL storage_getbase_int(rmatrix%RmatrixBlock(1,1)%h_Kld, rvanca%p_KldA )
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL storage_getbase_int(rmatrix%RmatrixBlock(1,1)%h_Kdiagonal, &
                               rvanca%p_KdiagonalA)
    ELSE
      rvanca%p_KdiagonalA => rvanca%p_KldA
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ2QP1simple (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the specialised diagonal VANCA algorithm for
  ! 2D Saddle-Point problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vancaPointer2DSPQ2QP1 structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DSPQ2QP1), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_POINTIDX)   :: NVT
    INTEGER(PREC_EDGEIDX)    :: NMT
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 9      ! Q2 = 9 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 3 ! QP1 = 3 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,j,isubdof
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscretisation% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscretisation% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    DO iel=1,NEL
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          8    9    6      
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Put on AA(.) the diagonal entry of matrix A -- the 1st and the
        ! 2nd block
        AA(inode,inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+nnvel,inode+nnvel) = p_DA(p_KdiagonalA(idof))
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity      
        ! components and B~ being an (2*9) x 12 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most  
        ! 4*3 pressure elements on the neighbour cell, so we have       
        ! 12 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        ! or
        !
        !              IEL   
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to  
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF's on that cell, so we have       
        ! 3 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of the original B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF's), 6 (if the velocity DOF is an edge with 
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on 
        ! the boundary and there is no neighbour, or if it's the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.
        
        DO ib = ib1,ib2
        
          IF (p_KcolB(ib) .EQ. IEL) THEN
            isubdof = 1
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL) THEN
            isubdof = 2
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL*2) THEN
            isubdof = 3
          ELSE
            ! Cycle the loop - the entry belongs to another 
            ! element, not to the current one
            CYCLE
          END IF
          
          J = p_KcolB(ib)
          
          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)                                                   :::::: )
      !     (          ..                                               :AA :: )
      !     (               ..                                          :(B1): )
      !     (                     AA(9,9)                               :::::: )
      !     (                            AA(10,10)                      :::::: )
      !     (                                      ..                   :AA :: )
      !     (                                           ..              :(B2): )
      !     (                                                AA(18,18)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      IF (ilapackInfo .NE. 0) PRINT *,'ERROR: LAPACK(DGESV) solver'
      
      ! Ok, we got the update vector in FF. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,nnvel
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                domega * FF(1+lofsp)
      p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                    domega * FF(2+lofsp)
      p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                      domega * FF(3+lofsp)
    
    END DO ! iel

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanca_2DSPQ2QP1full (rvanca, rvector, rrhs, domega)
  
!<description>
  ! This routine applies the specialised full local system VANCA algorithm for
  ! 2D Saddle-Point problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanca structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vancaPointer2DSPQ2QP1 structure that saves algorithm-specific parameters.
  TYPE(t_vancaPointer2DSPQ2QP1), INTENT(IN) :: rvanca

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
!</inputoutput>

!</subroutine>

    ! local vairables
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: inode,idof
    
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KdiagonalA
    REAL(DP), DIMENSION(:), POINTER             :: p_DA
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldB
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2
    
    ! Triangulation information
    INTEGER(PREC_ELEMENTIDX) :: NEL
    INTEGER(PREC_POINTIDX)   :: NVT
    INTEGER(PREC_EDGEIDX)    :: NMT
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs,p_Dvector
    
    ! Local arrays for informations about one element
    INTEGER, PARAMETER :: nnvel = 9      ! Q2 = 9 DOF's per velocity
    INTEGER, PARAMETER :: nnpressure = 3 ! QP1 = 3 DOF's per pressure
    INTEGER, PARAMETER :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF's per element
    INTEGER(PREC_VECIDX), DIMENSION(nnvel) :: IdofGlobal
    REAL(DP), DIMENSION(nnld,nnld) :: AA
    REAL(DP), DIMENSION(nnld) :: FF
    
    ! LAPACK temporary space
    INTEGER :: Ipiv(nnld),ilapackInfo
    
    ! offset information in arrays
    INTEGER(PREC_VECIDX)     :: ioffsetv,ioffsetp,j
    INTEGER :: ia1,ia2,ib1,ib2,ia,ib,isubdof,k
    INTEGER, PARAMETER :: lofsv = nnvel
    INTEGER, PARAMETER :: lofsp = 2*nnvel
    REAL(DP) :: daux
    
    ! Get pointers to the system matrix, so we don't have to write
    ! so much - and it's probably faster.
    
    p_KcolA => rvanca%p_KcolA
    p_KldA => rvanca%p_KldA
    p_KdiagonalA => rvanca%p_KdiagonalA
    p_DA => rvanca%p_DA
    p_KcolB => rvanca%p_KcolB
    p_KldB => rvanca%p_KldB
    p_DB1 => rvanca%p_DB1
    p_DB2 => rvanca%p_DB2
    p_DD1 => rvanca%p_DD1
    p_DD2 => rvanca%p_DD2
    
    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation%NEL
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscretisation% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    CALL storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscretisation% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    CALL lsysbl_getbase_double (rvector,p_Dvector)
    CALL lsysbl_getbase_double (rrhs,p_Drhs)
    
    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    
    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform 
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new" 
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !    
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    DO iel=1,NEL
    
      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP
      
      ! We now have the element
      !                                               
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF's          8    9    6      
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !                                               
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF's coincide with the definition
      ! in dofmapping.f90!
    
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)
      
      ! Get the velocity DOF's on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That's the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      DO inode=1,nnvel
      
        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)
        
        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        
        ! What do we have at this point?                           
        ! FF     : "local" RHS vector belonging to the DOF's on the
        !          current element                                 
        ! AA     : Diagonal entries of A belonging to these DOF's  
        !                                                          
        ! And at the moment:                                       
        ! idof      : number of current DOF on element IEL            
        ! inode     : "local" number of DOF on element IEL, i.e.      
        !              number of the edge         
        !                     
        ! Now comes the crucial point with the "update": How to         
        ! subsequently update the vertex values, such that the whole    
        ! thing still converges to the solution, even if a node         
        ! is updated more than once? Here, we use a typical             
        ! matrix-decomposition approach:                                
        !                                                               
        ! Again consider the problem:                                   
        !                                                               
        !    [ A   B ] (u) = (f)                                        
        !    [ B^t 0 ] (p)   (g)                                        
        !                                                               
        ! We assume, that all components in the vector (u,p) are        
        ! given - except for the velocity and pressure unknowns 
        ! on the current element; these 21 unknowns  
        ! are located anywhere in the (u,p) vector. The idea is to      
        ! shift "everything known" to the right hand side to obtain     
        ! a system for only these unknowns!                             
        !                                                               
        ! Extracting all the lines of the system that correspond to     
        ! DOF's on our single element IEL results in a rectangular      
        ! system of the form                                            
        !                                                               
        !    [ === A^ === B~ ] (|) = (f1)                                
        !    [ B~^t       0  ] (u)   (f2)                                
        !                      (|)   (g )                                   
        !                      (p)                                      
        !                                                               
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity      
        ! components and B~ being an (2*9) x 12 matrix that couples the   
        ! velocities to the pressure on our current element.            
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most  
        ! 4*3 pressure elements on the adjacent cells, so we have       
        ! 12 columns in the B-matrix.                                    
        !                                                               
        !        IEL                              IEL                   
        !     |--------|             |--------|--------|                
        !     |        |             |        |        |                
        !     |   P    |      or     |   Q    X   P    |                
        !     |   X    |             |        |        |                
        !   --|--------|--           |--------|--------|                
        !
        ! or
        !
        !              IEL   
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!  
        !                                                               
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)                               
        !   (d2)    (f2)   [ B~^t       0  ] (u2)                             
        !   (dp)    (g )                     (p)
        !                                                                 
        !
        ! That way, A^ is reduced to a square matrix with two square    
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to  
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF's on that cell, so we have       
        ! 3 columns in the B-matrix).                                   
        !
        ! At first build: fi = fi-Aui
        
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        DO ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          
          ! Whereever we find a DOF that couples to another DOF on the 
          ! same element, we put that to both A-blocks of our local matrix.
          DO k=1,nnvel
            IF (j .EQ. IdofGlobal(k)) THEN
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              EXIT
            END IF
          END DO          
        END DO
        
        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        DO ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        END DO
        
        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF's
        ! and so to element numbers. 
        !
        ! Btw: Each row of B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF's), 6 (if the velocity DOF is an edge with 
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on 
        ! the boundary and there is no neighbour, or if it's the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.
        
        DO ib = ib1,ib2
        
          IF (p_KcolB(ib) .EQ. IEL) THEN
            isubdof = 1
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL) THEN
            isubdof = 2
          ELSE IF (p_KcolB(ib) .EQ. IEL+NEL*2) THEN
            isubdof = 3
          ELSE
            ! Cycle the loop - the entry belongs to another 
            ! element, not to the current one
            CYCLE
          END IF
          
          J = p_KcolB(ib)
          
          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        END DO ! ib
        
      END DO ! inode
    
      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF's on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices. 
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF's on the current element. We already set up the preconditioner 
      ! in the above variables. It has the form:
      ! 
      ! C = ( AA(1,1)  ..............                                   :::::: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(9,9)                               :::::: )
      !     (                            AA(10,10) ..............       :::::: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA(18,18)  :::::: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.
      
      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'
      
      CALL DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      IF (ilapackInfo .NE. 0) PRINT *,'ERROR: LAPACK(DGESV) solver'
      
      ! Ok, we got the update vector in FF. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!
      
      DO inode=1,nnvel
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
      END DO
      
      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                domega * FF(1+lofsp)
      p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                    domega * FF(2+lofsp)
      p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                      domega * FF(3+lofsp)
    
    END DO ! iel

  END SUBROUTINE

END MODULE 
