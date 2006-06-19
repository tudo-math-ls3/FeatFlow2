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
!# Currently, there is only the implementation of the  
!#
!# The following routines can be found here:
!#
!# 1.) vanca_initGeneralVanca
!#     -> Initialise the general Vanca solver
!#
!# 2.) vanca_general
!#     -> Perform one step of the full VANCA solver for general block systems.
!#
!# 3.) vanca_dsoneGeneralVanca
!#     -> Clean up the general Vanca solver
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

!</types>

!<constants>
!<constantblock description="Constants defining the complexity of the discretisation">

  ! Number of elements to handle simultaneously when building vectors
  INTEGER :: VANCA_NELEMSIM   = 1000
  
!</constantblock>
!</constants>

CONTAINS

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
    ! Loop through the columns of the block matrix:
    DO i=1,nblocks
    
      ! Note this block as 'not processed'
      bfirst = .TRUE.
      
      ! Loop through the rows of the current matrix column.
      DO j=1,nblocks
        IF (rmatrix%RmatrixBlock(j,i)%NEQ .NE. 0) THEN
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
            rvancaGeneral%IblockOffset(i+1) = rmatrix%RmatrixBlock(j,i)%NEQ

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

    ! Offset position of the first block is = 0.
    rvancaGeneral%IblockOffset(1) = 0
    
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
  ! vector. THe rvancaGeneral structure has to be initialised before calling
  ! this routine, as this gets a reference to the system matrix A.
!</description>

!<input>
  ! The initial solution vector
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
  
  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
!</input>

!<inputoutput>
  ! The general-VANCA structure. Must have been initialised with 
  ! vanca_initGeneralVanca before.
  TYPE(t_vancaGeneral), INTENT(INOUT)     :: rvancaGeneral
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
      CALL dof_locGlobMapping_mult(rvector%RvectorBlock(1)%p_rspatialDiscretisation,&
                                   p_IelementList(IELset:IELmax), &
                                   .TRUE.,rvancaGeneral%IelementDOFs(:,:,i))

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
         rvancaGeneral%p_Rmatrices,IELmax-IELset+1,&
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
  INTEGER, DIMENSION(nblocks), INTENT(IN)                :: InDofsLocal
  
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
    INTEGER(PREC_DOFIDX), DIMENSION(ndofsPerElement+1) :: IlocalDOF,IglobalDOF
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
        
            ! Loop through the DOF's that correspond to this block:
            iidx = IlocalIndex(i)
            DO irow = IlocalIndex(i)+1,IlocalIndex(i+1)
              
              ! Get the actual DOF, relative to this block
              idof = IlocalDOF(irow)
              
              ! Loop through the row of the matrix to its contribution to "b-Ax".
              DO k = Rmatrices(i,j)%p_Kld(idof) , Rmatrices(i,j)%p_Kld(idof+1)-1

                ! Get the column number in the global matrix:              
                icol = Rmatrices(i,j)%p_Kcol(k)+IblockOffset(i)

                ! Build the defect
                Dff(irow) = Dff(irow) &
                          - dscale * Rmatrices(i,j)%p_Da(k) * Dvector(icol)

                ! icol is the number of a DOF.
                ! Check if this DOF belongs to the DOF's we have to
                ! extract to our local system.
                ! Loop through the DOF's corresponding to column j of the block system
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

END MODULE 

