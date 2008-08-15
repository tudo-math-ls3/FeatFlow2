!##############################################################################
!# ****************************************************************************
!# <name> multileveloperators </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the creation of extended multi-level
!# operators, e.g. Matrix-Creation for L2-Projection-operators.
!#
!# This module contains the following routines:
!#
!#  1.) mlop_create2LvlMatrixStruct
!#      -> Creates a scalar matrix structure in a specified matrix format
!#         for multi-level-based operations between two discretisations.
!#
!#  2.) mlop_build2LvlMassMatrix
!#      -> Assembles the entries of a 2-level-mass matrix that is needed for
!#         L2-Prolongation and restriction operators.
!# </purpose>
!##############################################################################

MODULE multileveloperators

  USE fsystem
  USE genoutput
  USE linearsystemscalar
  USE spatialdiscretisation
  USE scalarpde
  USE derivatives
  USE cubature
  USE collection
  USE domainintegration
  USE element
  USE elementpreprocessing
  
  IMPLICIT NONE

!<types>
  
!<typeblock>
  
  ! A private type block. Used as memory block during the creation of matrices.
  TYPE t_matrixmem
    ! A pointer to a 1D memory block of integers - receives the
    ! column identifiers of a matrix
    INTEGER(I32), DIMENSION(:), POINTER :: p_Icol
    
    ! A pointer to a 1D memory block of integers - receives 
    ! indices of next entries in the list of each line in the matrix
    INTEGER(I32), DIMENSION(:), POINTER :: p_Iindx
  END TYPE
  
!</typeblock>

  PRIVATE :: t_matrixmem

!</types>

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building matrices
  INTEGER, PARAMETER :: MLOP_NELEMSIM   = 100
  
!</constantblock>
!</constants>

CONTAINS


  !****************************************************************************

!<subroutine>

  SUBROUTINE mlop_create2LvlMatrixStruct (rdiscretisationCoarse,&
                      rdiscretisationFine, iformat, rmatrixScalar, imemguess)
  
!<description>
  ! This routine allows to calculate the structure of a finite-element matrix
  ! on the heap. The size of the matrix is determined dynamically.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisationFine

  ! Format of the matrix structure to be created. One of the LSYSSC_xxxx
  ! constants.
  INTEGER, INTENT(IN) :: iformat

  ! OPTIONAL: An initial guess about how much memory the matrix needs. If set 
  ! to 0 or not given, an initial guess of 16*NEQ (but at least 10000 matrix 
  ! entries) is assumed.
  INTEGER(I32), INTENT(IN), OPTIONAL :: imemGuess
  
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation
  ! combination. Memory fo the structure is allocated dynamically on the heap.
  TYPE(t_matrixScalar), INTENT(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  INTEGER(I32) :: imem
  
    imem = 0
    IF (PRESENT(imemguess)) THEN
      imem = MAX(0,imemguess)
    END IF
    
    ! Let's make sure that the discretisations are uniform - we don't support
    ! anything else right now.
    IF ((rdiscretisationCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
        (rdiscretisationFine%ccomplexity .NE. SPDISC_UNIFORM)) THEN
        
      CALL output_line ('Discretisations must be uniform!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      CALL sys_halt()
      
    END IF
    
    ! Which matrix structure do we have to create?
    SELECT CASE (iformat) 
    
    CASE (LSYSSC_MATRIX9)
    
      ! Call the creation routine for structure 9:
      CALL mlop_create2LvlMatStruct9_uni (rdiscretisationCoarse,&
                           rdiscretisationFine,rmatrixScalar,imem)
     
    CASE (LSYSSC_MATRIX7)
    
      ! Call the creation routine for structure 9:
      CALL mlop_create2LvlMatStruct9_uni (rdiscretisationCoarse,&
                           rdiscretisationFine,rmatrixScalar,imem)

      ! Translate to matrix structure 7:
      CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
      
    CASE DEFAULT
      CALL output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      CALL sys_halt()
      
    END SELECT

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE mlop_build2LvlMassMatrix (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
  
!<description>
  ! This routine calculates the entries of a 2-Level mass matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  TYPE(t_matrixScalar) :: rmatrixBackup
  
    ! The matrix must be unsorted, otherwise we can't set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there's a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    IF (rmatrixScalar%isortStrategy .GT. 0) THEN
      CALL output_line ('Matrix-structure must be unsorted!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      CALL sys_halt()
    END IF

    ! Let's make sure that the discretisations are uniform - we don't support
    ! anything else right now.
    IF ((rdiscretisationCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
        (rdiscretisationFine%ccomplexity .NE. SPDISC_UNIFORM)) THEN
        
      CALL output_line ('Discretisations must be uniform!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      CALL sys_halt()
      
    END IF

    ! Which matrix structure do we have?
    SELECT CASE (rmatrixScalar%cmatrixFormat) 
    CASE (LSYSSC_MATRIX9)
      CALL mlop_build2LvlMass9_uni (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
    
    CASE (LSYSSC_MATRIX7)
      ! Convert structure 7 to structure 9.For that purpose, make a backup of
      ! the original matrix...
      CALL lsyssc_duplicateMatrix (rmatrixScalar,rmatrixBackup,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      ! Convert the matrix 
      CALL lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX9)
      
      ! Create the matrix in structure 9
      CALL mlop_build2LvlMass9_uni (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
                                     
      ! Convert back to structure 7
      CALL lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX7)
      
      ! Copy the entries back to the original matrix and release memory.
      CALL lsyssc_duplicateMatrix (rmatrixBackup,rmatrixScalar,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
          
      CALL lsyssc_releaseMatrix (rmatrixBackup)
                                     
    CASE DEFAULT
      CALL output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      CALL sys_halt()
      
    END SELECT


  END SUBROUTINE
    
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE mlop_create2LvlMatStruct9_uni (rdiscrCoarse,rdiscrFine,&
                                            rmatrixScalar,imemGuess)
  
!<description>
  ! This routine creates according to a given discretisation the matrix 
  ! structure of a structure-9 matrix. The discretisation is assumed to be
  ! conformal, i.e. the DOF's of different FE spaces in the trial space
  ! fit together. The function space for trial and test functions 
  ! may be different.
!</description>

!<input>
  
  ! The underlying discretisation structures which are to be used to
  ! create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscrCoarse
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscrFine
  
  ! An initial guess about how much memory the matrix needs. If set to 0,
  ! an initial guess of 16*NEQ (but at least 10000 matrix entries) is assumed.
  INTEGER(I32), INTENT(IN) :: imemGuess
  
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation.
  ! Memory fo rthe structure is allocated dynamically on the heap.
  TYPE(t_matrixScalar), INTENT(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  INTEGER(PREC_DOFIDX) :: NEQ, IEQ, IROW, JCOL, IPOS, istartIdx, NA, nmaxCol
  INTEGER :: IDOFE, JDOFE, i, IHELP
  INTEGER(PREC_ELEMENTIDX) :: IELC,IELF,IELIDX,NELREF
  LOGICAL :: BSORT
  
  ! An allocateable list of handles for memory blocks. Size is dynamically 
  ! increased if there are too many columns in the matrix.
  INTEGER, DIMENSION(:), POINTER :: p_Ihcol, p_Ihindx, p_IhTmp
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_Isize, p_ISizeTmp
  
  ! An allocateable list of pointers to the memory blocks - corresponds
  ! to the handles in p_Ihcol/p_Ihindx
  TYPE(t_matrixmem), DIMENSION(:), POINTER :: Rmemblock, RmemblockTmp
  
  ! Pointer to current KCOL memory block,
  ! pointer to current index memory block
  INTEGER(I32), DIMENSION(:), POINTER :: p_Icol, p_Iindx
  
  ! Number of currently allocated pointers in Ihmemblock
  INTEGER :: iblocks
  
  ! Currently active memory block
  INTEGER :: icurrentblock
  
  ! Number of elements in the current coarse/fine mesh element distribution
  INTEGER(PREC_ELEMENTIDX) :: NELC, NELF
  
  ! Size of memory blocks
  INTEGER(PREC_DOFIDX) :: imemblkSize
  
  ! Blocksize in terms of NEQ for guessing memory.
  ! The initial guess for memory is iblkSize*iblkSize*NEQ and every time
  ! memory is needed, another iblkSize*NEQ elements are added.
  INTEGER, PARAMETER :: iblkSize = 4
  
  ! Number of memory blocks to allocate
  INTEGER, PARAMETER :: NmemBlkCount = 5

  ! Pointer to KLD, KCOL
  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL
  
  ! Size of memory currently allocated
  INTEGER(PREC_DOFIDX) :: iallocated
  
  ! An allocateable array accepting the DOF's of a set of elements.
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,MLOP_NELEMSIM), TARGET :: &
  !                  IdofsTest, IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial, p_IdofsTest
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList,p_IelementRef
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_relementDistribution

  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  INTEGER :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  INTEGER :: nelementsTrial, nelementsTest
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  INTEGER(I32), DIMENSION(:), POINTER :: p_IrefPatchIdx, p_IrefPatch

    ! The algorithm is: Test every DOF on one element against each other
    ! DOF on the same element and save the combination into a matrix
    ! in structure 9!
    !
    ! At first, initialise the structure-9 matrix:
    
    !rmatrixScalar%p_rspatialDiscretisation => rdiscretisation
    rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9
    
    ! Get the #DOF's of the test space - as #DOF's of the test space is
    ! the number of equations in our matrix. The #DOF's in the trial space
    ! gives the number of columns of our matrix.
    rmatrixScalar%NCOLS         = dof_igetNDofGlob(rdiscrCoarse)
    rmatrixScalar%NEQ           = dof_igetNDofGlob(rdiscrFine)
    
    ! and get a pointer to the triangulation.
    p_rtriaCoarse => rdiscrCoarse%p_rtriangulation
    p_rtriaFine => rdiscrFine%p_rtriangulation
    
    ! Get the refinement patch arrays from the fine triangulation
    CALL storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    CALL storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Get NEQ - we need it for guessing memory...
    NEQ = rmatrixScalar%NEQ
    
    IF (NEQ .EQ. 0) THEN
      CALL output_line ('Empty matrix!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatStruct9_uni')
      CALL sys_halt()
    END IF
    
    ! Allocate KLD...
    CALL storage_new1D ('mlop_create2LvlMatStruct9_uni', 'KLD', &
                        NEQ+1_I32, ST_INT, rmatrixScalar%h_KLD, &
                        ST_NEWBLOCK_NOINIT)
    ! This must be a storage_getbase, no lsyssc_getbase, since this is the
    ! matrix construction routine!
    CALL storage_getbase_int(rmatrixScalar%h_Kld,p_KLD)
    
    ! Allocate a list of handles and a list of pointers corresponding to it.
    ! Initially allocate NmemBlkCount pointers
    ALLOCATE(p_Ihcol(NmemBlkCount))
    ALLOCATE(p_Ihindx(NmemBlkCount))
    ALLOCATE(p_Isize(NmemBlkCount))
    ALLOCATE(Rmemblock(NmemBlkCount))
    
    ! Allocate the first memory block that receives a part of the
    ! temporary matrix structure.
    ! We make an initial guess of iblkSize*iblkSize*NEQ elements in the matrix,
    ! if imemguess is not given.
    IF (imemguess .NE. 0) THEN 
      ! at least one element per line!
      iallocated = MAX(NEQ,imemguess) 
    ELSE  
      iallocated = MAX(10000,iblkSize*iblkSize*NEQ)
    END IF
    iblocks = 1
    imemblkSize = iblkSize*NEQ
    p_Isize(1) = iallocated
    
    ! imemblkSize = iallocated is necessary at the moment to simplify
    ! whether we leave a block or not.
    CALL storage_new1D ('mlop_create2LvlMatStruct9_uni', 'Ihicol', &
                        p_Isize(1), ST_INT, p_Ihcol(1), ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int (p_Ihcol(1),p_Icol)

    ! The new index array must be filled with 0 - otherwise
    ! the search routine below won't work!
    CALL storage_new1D ('mlop_create2LvlMatStruct9_uni', 'p_Ihindx', &
                        p_Isize(1), ST_INT, p_Ihindx(1), ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int (p_Ihindx(1),p_Iindx)
    
    Rmemblock(1)%p_Icol => p_Icol
    Rmemblock(1)%p_Iindx => p_Iindx
    
    ! Initialise Iindx and Icol.
    ! Initially, we have only diagonal elements in our matrix.
    !
    ! The basic idea behind the building of the matrix is a linked
    ! list of column numbers in each row!
    ! We collect all upcoming columns in the whole matrix in the
    ! array Icol, i.e. each new entry is attached to that
    ! (-> the array is resorted to KCOL at the end).
    ! Iindx points for every entry in the matrix to the position
    ! inside of Icol of the next entry in the line.
    !
    ! At the beginning, we no entries in the matrix. 
    ! We initialise the "head" of this collection of linked
    ! lists with 0 to indicate this. 
    ! Iindx(IEQ) is set to 0 to indicate that there is no following
    ! element in each line, i.e. we only have diagonal elements.
    ! Later, we fill Icol(1..NEQ) with the first column number in
    ! each row. When more entries appear in a row, they are appended
    ! to KCOL at position NEQ+1..*. Iindx keeps track of the entries
    ! in each row by storing (starting from the head in Icol(IEQ))
    ! the positions of the corresponding next entry inside of KCOL1
    ! in each row - so the column numbers in each row are to be
    ! found in KCOL1 at positions IEQ,Iindx(IEQ),Iindx(Iindx(IEQ)),...
    !
    ! Example: We want to add entry (1,3) to the matrix. Then
    ! we enlarge the lists as follows:
    ! - "2" (the column number" is added to the end of Icol, i.e.
    !   Icol(NEQ+1) = 3
    ! - We add a "follower" for the diagonal element by setting
    !   Iindx(1) = NEQ+1. Furthermore we set KIND(NEQ+1)=0
    ! So we have a linked list of matrix entries for the rows:
    !   Icol:     1   2   3   ...   NEQ     3
    !   Iindx:  NEQ+1 0   0   ...    0      0
    !             |                        /:\
    !             +-------------------------|
    ! i.e. row 1 can be computed as:
    !   Icol(1)               (=1),
    !   Icol(Iindx(1))        (=3),
    !   Icol(Iindx(Iindx(1))  -> not defined, as Iindx(Iindx(1))=0, 
    !                            line finished
    
    DO IEQ=1,NEQ
      p_Iindx(IEQ) = 0
      p_Icol(IEQ)  = 0
    END DO
    
    ! The first NEQ elements are reserved. The matrix is assumed
    ! to have at least one element per line.
    NA = NEQ
    
    ! Activate the current coarse mesh element distribution
    p_relementDistribution => rdiscrCoarse%RelementDistr(1)

    ! Cancel if this element distribution is empty.
    IF (p_relementDistribution%NEL .EQ. 0) THEN
      CALL output_line('Element space is empty!')
      CALL sys_halt()
    END IF

    ! Get the number of local DOF's for trial and test functions
    indofTrial = elem_igetNDofLoc(rdiscrCoarse%RelementDistr(1)%celement)
    indofTest = elem_igetNDofLoc(rdiscrFine%RelementDistr(1)%celement)
    
    ! Calculate the number of coarse mesh elements we want to process
    ! in one run.
    nelementsTrial = MIN(MLOP_NELEMSIM,p_relementDistribution%NEL) 
    
    ! Now calculate the number of fine mesh elements we want to process
    ! in one run. 
    SELECT CASE(p_rtriaFine%ndim)
    CASE (1)
      nelementsTest = 2*nelementsTrial + 10
    CASE (2)
      nelementsTest = 4*nelementsTrial + 20
    CASE (3)
      nelementsTest = 8*nelementsTrial + 40
    END SELECT
    
    ! Allocate an array saving a couple of DOF's for trial and test functions
    ALLOCATE(p_IdofsTrial(indofTrial,nelementsTrial))
    ALLOCATE(p_IdofsTest(indofTest,nelementsTest))
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that the trial functions
    CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)
    
    ! And allocate the refinemed element list for the test functions
    ALLOCATE(p_IelementRef(nelementsTest))
    
    ! Get the number of coarse mesh elements there.
    NELC = p_relementDistribution%NEL

    ! Set the pointers/indices to the initial position. During the
    ! search for new DOF's, these might be changed if there's not enough
    ! memory in the first block.    
    icurrentblock = 1
    istartidx = 0
    p_Icol => Rmemblock(1)%p_Icol
    p_Iindx => Rmemblock(1)%p_Iindx
    
    ! Loop over the elements. 
    nelementsDone = 0
    DO WHILE(nelementsDone .LT. NELC)
    
      ! We always try to handle nelementsTrial elements simultaneously.
      ! Of course, we will not handle more elements than the coarse
      ! mesh discretisation has.
      nelementsToDo = MIN(NELC-nelementsDone, nelementsTrial)
      
      ! Now comes the interesting part - we have to ensure that the DOF-mapping
      ! of the fine mesh discretisation fits into our DOF-array.
      ! If, for example, a coarse mesh quad was refined into more than 4 fine
      ! mesh quads, then it might happen that we cannot handle nelementsToDo
      ! coarse mesh elements at once, but we need to decrease nelementsToDo.
      
      NELF = 0
      DO IELC = 1, nelementsToDo
      
        ! Get the index of the coarse mesh element
        IELIDX = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that have been refined from the
        ! currently processed coarse mesh element.
        NELREF = p_IrefPatchIdx(IELIDX+1) - p_IrefPatchIdx(IELIDX)
        
        ! Now if (NELF+NELREF) is greater than nelementsTest, then we need
        ! to decrease nelementsToDo and exit the loop...
        ! This case should never happen if the coarse mesh was refined using
        ! the 2-Level-Ordering algorithm, but might happen for more freaky
        ! refinement techniques...
        IF((NELF+NELREF) .GT. nelementsTest) THEN
          nelementsToDo = IELC-1
          EXIT
        END IF
        
        ! Copy the indices of the elements into the element list for the
        ! fine mesh discretisation
        DO IELF = 1, NELREF
          p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IELIDX)+IELF-1)
        END DO
        
        ! Add the number of refined elements to the counter
        NELF = NELF + NELREF
      
      END DO
      
      ! If nelementsToDo is 0, then we have a serious problem...
      IF (nelementsToDo .LE. 0) THEN
        PRINT *, "ERROR: mlop_create2LvlMatStruct9_uni"
        PRINT *, "nelementsToDo = 0 !!!"
        CALL sys_halt()
      END IF
      
      ! Call the DOF-mapping routine for the coarse and fine mesh
      CALL dof_locGlobMapping_mult(rdiscrCoarse, &
          p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
          p_IdofsTrial)
      CALL dof_locGlobMapping_mult(rdiscrFine, p_IelementRef(1:NELF), &
          p_IdofsTest)
      
      ! Reset the counter
      NELF = 0
      
      ! Loop through all the elements of the coarse mesh set
      DO IELC = 1, nelementsToDo
      
        ! Get the index of the currently processed coarse mesh element
        IELIDX = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that are refined from the
        ! current coarse mesh element
        NELREF = p_IrefPatchIdx(IELIDX+1) - p_IrefPatchIdx(IELIDX)

        ! And loop through all elements of the current refinement patch
        DO IELF = 1, NELREF
        
          ! For building the local matrices, we have first to
          ! loop through the test functions (the "O"'s), as these
          ! define the rows in the matrix.
          DO IDOFE=1,indofTest

            ! The DOF IDOFE is now our "O".
            ! This global DOF gives us the row we have to build.
            IROW = p_IdofsTest(IDOFE,NELF+IELF)
            
            ! Now we loop through the other DOF's on the current element
            ! (the "X"'s).
            ! All these have common support with our current basis function
            ! and will therefore give an additive value to the global
            ! matrix.

            DO JDOFE=1,indofTrial
              
              ! Get the global DOF - our "X". This gives the column number
              ! in the matrix where an entry occurs in row IROW (the line of 
              ! the current global DOF "O").
              JCOL = p_IdofsTrial(JDOFE,IELC)

              ! This JCOL has to be inserted into line IROW.
              ! But first check, whether the element is already in that line,
              ! i.e. whether element (IROW,JCOL) is already in the matrix.
              ! This may happen because of an earlier loop in a neighbour
              ! element (imagine Q1, where there are two vertices on an edge
              ! and the edge is shared between two elements)...
              !
              ! We start walking through the linked list of row IROW to 
              ! look for column JCOL. IPOS is the position in the KCOL1 
              ! array of the column in row IROW we want to test. 
              ! KINDX(IPOS) is =0 if we reach the end of the linked list.
              !
              ! Start searching at the "head" of the list, which is the
              ! diagonal element at Icol(IROW).
              ! This is always found in the first memory block.
        
              ! Remark: This IF command gains 8% performance due to slow
              ! pointer handling!
              IF (icurrentblock .NE. 1) THEN
                icurrentblock = 1
                istartidx = 0
                p_Icol => Rmemblock(1)%p_Icol
                p_Iindx => Rmemblock(1)%p_Iindx
              END IF
              
              ! Is the list empty?

              IF (p_Icol(IROW).EQ.0) THEN

                ! Yes, row IROW is empty at the moment. Add the column as
                ! head of the list of row IROW.

                p_Icol(IROW) = JCOL
                
              ELSE

                ! No, the list is not empty, we have a "head".
                !
                ! We start walking through the linked list of row IROW to 
                ! look for column JCOL. IPOS is the position in the KCOL1 
                ! array of the column in row IROW we want to test. 
                ! KINDX(IPOS) is =0 if we reach the end of the linked list.
                !
                ! Start searching at the "head" of the list, which is the
                ! diagonal element at p_Icol(IROW).

                IPOS=IROW
              
                ! Loop through the elements in the list until we find
                ! column JCOL - or the end of the list!
                ! IPOS must be corrected by istartidx, which is the number
                ! of elements in all blocks before the current one.
                
                searchloop: DO WHILE ( (p_Icol(IPOS-istartIdx)) .NE. JCOL)
              
                  ! Did we reach the end of the list? Then we have to insert
                  ! a new element...
                  IF (p_Iindx(IPOS-istartIdx) .EQ. 0) THEN
                  
                    ! Increase NA, which is the actual length of the matrix -
                    ! and at the same time tells us how much memory we actually
                    ! use.
                    NA = NA+1
                    
                    ! Let p_Iindx of the last element of the row point to our
                    ! new element. The new element is now the last in the row.
                  
                    p_Iindx(IPOS-istartIdx) = NA
                    
                    ! Before really appending JCOL, first test
                    ! NA is now larger than the marimum amount
                    ! of storage we allocated!
                    
                    IF (NA .GT. iallocated) THEN
                     
                      ! Hmmm, we have to allocate more memory.
                      ! Do we have enough pointers left or do we have
                      ! to enlarge our list?
                      
                      IF (iblocks .GE. SIZE(p_Ihcol)) THEN 
                      
                        ! Not enough blocks, we have to reallocate the pointer lists!
                        ALLOCATE (p_IhTmp(iblocks+NmemBlkCount))
                        p_IhTmp(1:iblocks) = p_Ihcol(1:iblocks)
                        DEALLOCATE(p_Ihcol)
                        p_Ihcol => p_IhTmp

                        ALLOCATE (p_IhTmp(iblocks+NmemBlkCount))
                        p_IhTmp(1:iblocks) = p_Ihindx(1:iblocks)
                        DEALLOCATE(p_Ihindx)
                        p_Ihindx => p_IhTmp
                      
                        ALLOCATE (p_IsizeTmp(iblocks+NmemBlkCount))
                        p_IsizeTmp(1:iblocks) = p_Isize(1:iblocks)
                        DEALLOCATE(p_Isize)
                        p_Isize => p_IsizeTmp

                        ALLOCATE (RmemblockTmp(iblocks+NmemBlkCount))
                        RmemblockTmp(1:iblocks) = Rmemblock(1:iblocks)
                        DEALLOCATE(Rmemblock)
                        Rmemblock => RmemblockTmp
                        
                        ! Now we have enough blocks again.
                      END IF

                      ! Add a new block

                      iblocks = iblocks + 1
                      p_Isize (iblocks) = imemblkSize
                      
                      ! Move the start index behind the last completely
                      ! occupied block
                               
                      istartIdx = iallocated
                      icurrentblock = iblocks
                    
                      ! Allocate a new memory block of size imemblkSize
                      !
                      ! Use p_Icol and p_Iindx - they are not used anymore.
                      ! Allocate another imemblkSize elements for column numbers and
                      ! list pointers.

                      CALL storage_new1D ('mlop_create2LvlMatStruct9_uni', 'Ihicol', &
                                          p_Isize (iblocks), ST_INT, p_Ihcol(iblocks), &
                                          ST_NEWBLOCK_NOINIT)
                      CALL storage_getbase_int (p_Ihcol(iblocks),p_Icol)

                      ! The new index array must be filled with 0 - otherwise
                      ! the search routine below won't work!
                      CALL storage_new1D ('mlop_create2LvlMatStruct9_uni', 'p_Ihindx', &
                                          p_Isize (iblocks), ST_INT, p_Ihindx(iblocks), &
                                          ST_NEWBLOCK_ZERO)
                      CALL storage_getbase_int (p_Ihindx(iblocks),p_Iindx)
                      
                      Rmemblock(iblocks)%p_Icol => p_Icol
                      Rmemblock(iblocks)%p_Iindx => p_Iindx

                      iallocated = iallocated + p_Isize (iblocks)
                      
                    ELSE
                      ! Be careful when leaving the current memory block
                      ! for insertion of an element at position NA!
                      !
                      ! If the new position is not in the current block,
                      ! it's in the last block... and so set the pointer
                      ! and indices appropriately!
                      
                      IF ( NA .GT. (istartidx+p_Isize (icurrentblock))) THEN
                        istartidx = iallocated-p_Isize(iblocks)
                        icurrentblock = iblocks
                        p_Icol => Rmemblock(iblocks)%p_Icol
                        p_Iindx => Rmemblock(iblocks)%p_Iindx
                      END IF
                      
                    END IF
                  
                    ! Append JCOL to p_Icol
                    p_Icol(NA-istartIdx) = JCOL
                    
                    ! We have to make sure that p_Indx(NA)=0 to indicate the end of
                    ! the list. Ok, this is trivial because we allocated it with 
                    ! storage_new, configured to fill the memory with 0, so it is 0.
                    !
                    ! The searchloop ends here, continue with next JDOFE
                    
                    EXIT 
                  
                  ELSE
                  
                    ! No, the list does not end here.
                    ! Take the next element in the list
                    IPOS = p_Iindx(IPOS-istartidx)
                    
                    ! Be careful when leaving the current memory block
                    DO WHILE ( IPOS .GT. (istartidx+p_Isize (icurrentblock)) )
                    
                      ! go to the next memory block and search there
                      istartidx = istartidx+p_Isize(icurrentblock)
                      icurrentblock = icurrentblock+1
                      p_Icol => Rmemblock(icurrentblock)%p_Icol
                      p_Iindx => Rmemblock(icurrentblock)%p_Iindx
                      
                    END DO ! IPOS .GT. (istartidx+p_Isize (iblocks))
                  
                  END IF ! p_Iindx(IPOS) = 0
                  
                END DO searchloop
              
              END IF ! p_Icol(IROW) = 0
                  
            END DO ! JDOFE
          
          END DO ! IDOFE
        
        END DO ! IELF

        ! Add the number of refined elements to the counter
        NELF = NELF + NELREF
      
      END DO ! IELC
      
      ! Add the number of processed coarse mesh elements
      nelementsDone = nelementsDone + nelementsToDo
    
    END DO ! WHILE(nelementsDone .LT. NEL)
    
    ! Release the fine mesh element list
    DEALLOCATE(p_IelementRef)

    ! Clean up the DOF's arrays    
    DEALLOCATE(p_IdofsTest)
    DEALLOCATE(p_IdofsTrial)
    
    ! Ok, p_Icol is built. The hardest part is done!
    ! Now build KCOL by collecting the entries in the linear lists of 
    ! each row.
    !
    ! At first, as we now NA, we can allocate the real KCOL now!
    
    CALL storage_new1D ('mlop_create2LvlMatStruct9_uni', 'KCOL', &
                        NA, ST_INT, rmatrixScalar%h_KCOL, &
                        ST_NEWBLOCK_NOINIT)
    ! This must be a storage_getbase, no lsyssc_getbase, since this is the
    ! matrix construction routine!
    CALL storage_getbase_int (rmatrixScalar%h_Kcol,p_KCOL)
    
    ! Save NA in the matrix structure
    rmatrixScalar%NA = NA
    
    ! Set back NA to 0 at first.
    NA=0
        
    ! Loop through all of the NEQ linear lists:
    DO IEQ=1,NEQ
    
      ! We are at the head of the list, now we have to walk
      ! through it to append the entries to KCOL.
      ! We always start in the first memory block.
      IF (icurrentblock .NE. 1) THEN
        icurrentblock = 1
        istartidx = 0
        p_Icol => Rmemblock(1)%p_Icol
        p_Iindx => Rmemblock(1)%p_Iindx
      END IF
      
      ! Add the head of the list to KCOL:
      NA=NA+1
      p_KCOL(NA)=p_Icol(IEQ)

      ! Set KLD appropriately:
      p_KLD(IEQ) = NA
      
      IPOS = IEQ
      
      DO WHILE (p_Iindx(IPOS-istartidx).NE.0)
      
        ! Get the position of the next entry in p_Icol:
        IPOS = p_Iindx(IPOS-istartidx)
        
        ! Be careful when leaving the current memory block
        DO WHILE ( IPOS .GT. (istartidx+p_Isize (icurrentblock)) )
        
          ! go to the next memory block and search there
          istartidx = istartidx+p_Isize(icurrentblock)
          icurrentblock = icurrentblock+1
          p_Icol => Rmemblock(icurrentblock)%p_Icol
          p_Iindx => Rmemblock(icurrentblock)%p_Iindx
          
        END DO ! IPOS .GT. (istartidx+p_Isize (iblocks))
        
        ! Add the column number to the row in KCOL:
        NA=NA+1
        p_KCOL(NA)=p_Icol(IPOS-istartidx)
      
      END DO ! KINDX(IPOS) <> 0

    END DO ! IEQ
    
    ! Append the final entry to KLD:
    p_KLD(NEQ+1)=NA+1
    
    ! Sort entries on KCOL separately for each row.
    ! This is a small bubble-sort...
    !
    ! Loop through all rows:
    nmaxCol = 0

    DO IEQ=1,NEQ

      ! Repeat until everything is sorted.
      
      BSORT=.FALSE.
      DO WHILE (.NOT. BSORT)
      
        BSORT=.TRUE.

        !  Loop through the line 

        DO JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-2
        
          ! If the next element is larger...
        
          IF (p_KCOL(JCOL) .GT. p_KCOL(JCOL+1)) THEN
          
            ! Change position of the current and next element
          
            IHELP=p_KCOL(JCOL)
            p_KCOL(JCOL)=p_KCOL(JCOL+1)
            p_KCOL(JCOL+1)=IHELP
            
            ! And repeat the sorting of that line
            
            BSORT=.FALSE.
            
          END IF
          
        END DO ! JCOL
        
      END DO ! (not BSORT)      

      ! Grab the largest column number. As the current line is sorted,
      ! we can find this using the end of the line.
      nmaxCol = MAX(nmaxCol,p_Kcol(p_Kld(IEQ+1)-1))

    END DO ! IEQ
    
    ! HOORAY, THAT'S IT!
    ! Deallocate all temporary memory...
    
    DO i=iblocks,1,-1
      CALL storage_free(p_Ihcol(i))
      CALL storage_free(p_Ihindx(i))
    END DO
    
    DEALLOCATE(Rmemblock)
    DEALLOCATE(p_Isize)
    DEALLOCATE(p_Ihindx)
    DEALLOCATE(p_Ihcol)
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE mlop_build2LvlMass9_uni (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
  
!<description>
  ! This routine calculates the entries of a 2-Level mass matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>


  ! local variables
  INTEGER :: i,k,JDFG, ICUBP, NELC,NELF
  INTEGER(I32) :: IELC,IELF, IDXC, NELREF, IDOFE, JDOFE
  INTEGER(PREC_DOFIDX) :: JCOL0,JCOL
  REAL(DP) :: OM, DB
  
  ! Array to tell the element which derivatives to calculate
  LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp,ncubpc
  
  ! Pointer to KLD, KCOL, DA
  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsCoarse, IdofsFine
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasCoarse,DbasFine
  
  ! Number of entries in the matrix - for quicker access
  INTEGER(PREC_DOFIDX) :: NA
  INTEGER(I32) :: NEQ
  
  ! Type of transformation from the reference to the real element 
  INTEGER :: ctrafoCoarse, ctrafoFine
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  INTEGER(I32) :: cevalTagCoarse, cevalTagFine
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofCoarse, indofFine
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList, p_IelementRef
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dentry

  ! An array that takes coordinates of the cubature formula on the reference element
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_DcubPtsRefFine, p_DcubPtsRefCoarse
  
  ! Pointer to the jacobian determinants
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj

  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_relemDistCoarse, p_relemDistFine
  
  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  INTEGER :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  INTEGER :: nelementsCoarse, nelementsFine
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  INTEGER(I32), DIMENSION(:), POINTER :: p_IrefPatchIdx, p_IrefPatch
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines as well as element evaluation routines.
  TYPE(t_domainIntSubset) :: rintSubsetCoarse, rintSubsetFine

    ! We only need the function values as we want to assemble a mass matrix.
    Bder = .FALSE.
    Bder(DER_FUNC) = .TRUE.
    
    ! Get information about the matrix:
    NA = rmatrixScalar%NA
    NEQ = rmatrixScalar%NEQ
    
    ! We need KCOL/KLD of our matrix
    IF ((rmatrixScalar%h_KCOL .EQ. ST_NOHANDLE) .OR. &
        (rmatrixScalar%h_KLD .EQ. ST_NOHANDLE)) THEN
      CALL output_line ('No discretisation structure! Cannot assemble matrix!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMass9_uni')
      CALL sys_halt()
    END IF
    
    CALL lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
    CALL lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    
    ! Check if the matrix entries exist. If not, allocate the matrix.
    IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN

      ! Clear the entries in the matrix - we need to start with zero
      ! when assembling a new matrix!
      CALL storage_new1D ('mlop_build2LvlMass9_uni', 'DA', &
                          NA, ST_DOUBLE, rmatrixScalar%h_DA, &
                          ST_NEWBLOCK_ZERO)
      CALL lsyssc_getbase_double (rmatrixScalar,p_DA)

    ELSE
    
      CALL lsyssc_getbase_double (rmatrixScalar,p_DA)

      ! If desired, clear the matrix before assembling.
      IF (bclear) THEN
        CALL lalg_clearVectorDble (p_DA)
      END IF
      
    END IF
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriaCoarse => rdiscretisationCoarse%p_rtriangulation
    p_rtriaFine => rdiscretisationFine%p_rtriangulation

    ! Get the refinement patch arrays from the fine triangulation
    CALL storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    CALL storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Activate the current element distributions
    p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(1)
    p_relemDistFine => rdiscretisationFine%RelementDistr(1)
  
    ! Get the number of local DOF's for trial and test functions
    indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
    indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
      
    ! Calculate the number of coarse mesh elements we want to process
    ! in one run.
    nelementsCoarse = MIN(MLOP_NELEMSIM,p_relemDistCoarse%NEL) 
    
    ! Now calculate the number of fine mesh elements we want to process
    ! in one run. 
    SELECT CASE(p_rtriaFine%ndim)
    CASE (1)
      nelementsFine = 2*nelementsCoarse
    CASE (2)
      nelementsFine = 4*nelementsCoarse
    CASE (3)
      nelementsFine = 8*nelementsCoarse
    END SELECT
    
    ! Allocate an array saving a couple of DOF's for trial and test functions
    ALLOCATE(IdofsCoarse(indofCoarse,nelementsCoarse))
    ALLOCATE(IdofsFine(indofFine,nelementsFine))
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that the trial functions
    CALL storage_getbase_int (p_relemDistCoarse%h_IelementList, p_IelementList)
    
    ! And allocate the refinemed element list for the test functions
    ALLOCATE(p_IelementRef(nelementsFine))
    
    ! Get the number of coarse mesh elements there.
    NELC = p_relemDistCoarse%NEL
      
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
    ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
    
    ! Initialise the cubature formula, get cubature weights and point
    ! coordinates on the reference element of the fine mesh
    CALL cub_getCubPoints(p_relemDistFine%ccubTypeBilForm, ncubp, Dxi, Domega)
    
    ! Allocate some memory to hold the cubature points on the fine mesh
    ALLOCATE(p_DcubPtsRefFine(trafo_igetReferenceDimension(ctrafoFine),ncubp))

    ! Reformat the cubature points; they are in the wrong shape!
    DO i=1,ncubp
      DO k=1,UBOUND(p_DcubPtsRefFine,1)
        p_DcubPtsRefFine(k,i) = Dxi(i,k)
      END DO
    END DO
    
    ! Now we need to transform the points from the fine mesh into the coarse mesh
    ! Please note that the following trick does only work for 2-level ordered
    ! meshes!
    SELECT CASE(p_rtriaFine%ndim)
    CASE (1)
      ncubpc = 2*ncubp
      ALLOCATE(p_DcubPtsRefCoarse(trafo_igetReferenceDimension(ctrafoCoarse),ncubpc))
      CALL trafo_mapCubPtsRef2LvlEdge1D(ncubp,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
      
    CASE (2)
      ncubpc = 4*ncubp
      ALLOCATE(p_DcubPtsRefCoarse(trafo_igetReferenceDimension(ctrafoCoarse),ncubpc))
      CALL trafo_mapCubPtsRef2LvlQuad2D(ncubp,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

    CASE (3)
      ncubpc = 8*ncubp
      ALLOCATE(p_DcubPtsRefCoarse(trafo_igetReferenceDimension(ctrafoCoarse),ncubpc))
      CALL trafo_mapCubPtsRef2LvlHexa3D(ncubp,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
      
    END SELECT
    
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    ALLOCATE(DbasFine(indofFine,&
             elem_getMaxDerivative(p_relemDistFine%celement),&
             ncubp,nelementsFine))
    ALLOCATE(DbasCoarse(indofCoarse,&
             elem_getMaxDerivative(p_relemDistCoarse%celement), &
             ncubpc,nelementsCoarse))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    ALLOCATE(Kentry(indofCoarse,indofFine,nelementsFine))
    ALLOCATE(Dentry(indofCoarse,indofFine,nelementsFine))
    
    ! Loop over the elements - blockwise.
    nelementsDone = 0
    DO WHILE(nelementsDone .LT. NELC)
    
      ! We always try to handle nelementsTrial elements simultaneously.
      ! Of course, we will not handle more elements than the coarse
      ! mesh discretisation has.
      nelementsToDo = MIN(NELC-nelementsDone, nelementsCoarse)
      
      ! Now comes the interesting part - we have to ensure that the DOF-mapping
      ! of the fine mesh discretisation fits into our DOF-array.
      ! If, for example, a coarse mesh quad was refined into more than 4 fine
      ! mesh quads, then it might happen that we cannot handle nelementsToDo
      ! coarse mesh elements at once, but we need to decrease nelementsToDo.
      NELF = 0
      DO IELC = 1, nelementsToDo
      
        ! Get the index of the coarse mesh element
        IDXC = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that have been refined from the
        ! currently processed coarse mesh element.
        NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)
        
        ! Now if (NELF+NELREF) is greater than nelementsFine, then we need
        ! to decrease nelementsToDo and exit the loop...
        ! This case should never happen if the coarse mesh was refined using
        ! the 2-Level-Ordering algorithm, but might happen for more freaky
        ! refinement techniques...
        IF((NELF+NELREF) .GT. nelementsFine) THEN
          nelementsToDo = IELC-1
          EXIT
        END IF
        
        ! Copy the indices of the elements into the element list for the
        ! fine mesh discretisation
        DO IELF = 1, NELREF
          p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IDXC)+IELF-1)
        END DO
        
        ! Add the number of refined elements to the counter
        NELF = NELF + NELREF
      
      END DO
      
      ! If nelementsToDo is 0, then we have a serious problem...
      IF (nelementsToDo .LE. 0) THEN
        PRINT *, "ERROR: mlop_build2LvlMass9_uni"
        PRINT *, "nelementsToDo = 0 !!!"
        CALL sys_halt()
      END IF
      
      ! Call the DOF-mapping routine for the coarse and fine mesh
      CALL dof_locGlobMapping_mult(rdiscretisationCoarse, &
          p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
          IdofsCoarse)
      CALL dof_locGlobMapping_mult(rdiscretisationFine, p_IelementRef(1:NELF), &
          IdofsFine)
      
      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
      NELF = 0
      DO IELC = 1, nelementsToDo
      
        ! Get the index of the currently processed coarse mesh element
        IDXC = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that are refined from the
        ! current coarse mesh element
        NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

        ! And loop through all elements of the current refinement patch
        DO IELF = 1, NELREF
        
          ! Get the index of the currently processed fine mesh element
          !IDXF = p_IelementRef(NELF+IELF)
      
          ! For building the local matrices, we have first to
          ! loop through the test functions (the "O"'s), as these
          ! define the rows in the matrix.
          DO IDOFE=1,indofFine
          
            ! Row IDOFE of the local matrix corresponds 
            ! to row=global DOF KDFG(IDOFE) in the global matrix.
            ! This is one of the the "O"'s in the above picture.
            ! Get the starting position of the corresponding row
            ! to JCOL0:
            JCOL0 = p_KLD(IdofsFine(IDOFE,NELF+IELF))
            
            ! Now we loop through the other DOF's on the current element
            ! (the "O"'s).
            ! All these have common support with our current basis function
            ! and will therefore give an additive value to the global
            ! matrix.
            DO JDOFE=1,indofCoarse
              
              ! Get the global DOF of the "X" which interacts with 
              ! our "O".
              JDFG = IdofsCoarse(JDOFE,IELC)
              
              ! Starting in JCOL0 (which points to the beginning of
              ! the line initially), loop through the elements in
              ! the row to find the position of column IDFG.
              ! Jump out of the DO loop if we find the column.
              DO JCOL=JCOL0,NA
                IF (p_KCOL(JCOL) .EQ. JDFG) EXIT
              END DO

              ! Because columns in the global matrix are sorted 
              ! ascendingly (except for the diagonal element),
              ! the next search can start after the column we just found.
              
              ! JCOL0=JCOL+1
              
              ! Save the position of the matrix entry into the local
              ! matrix.
              ! Note that a column in Kentry corresponds to a row in
              ! the real matrix. We aligned Kentry/DENTRY this way to get
              ! higher speed of the assembly routine, since this leads
              ! to better data locality.
              
              Kentry(JDOFE,IDOFE,NELF+IELF)=JCOL
              
            END DO ! IDOFE
            
          END DO ! JDOFE
        
        END DO ! IELF
        
        NELF = NELF + NELREF
        
      END DO ! IELC
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevalTagCoarse = elem_getEvaluationTag(p_relemDistCoarse%celement)
      cevalTagFine = elem_getEvaluationTag(p_relemDistFine%celement)
                      
      cevalTagCoarse = IOR(cevalTagCoarse,EL_EVLTAG_REFPOINTS)
      cevalTagFine   = IOR(cevalTagFine,EL_EVLTAG_REFPOINTS)

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      CALL elprep_prepareSetForEvaluation (rintSubsetCoarse%revalElementSet,&
          cevalTagCoarse, p_rtriaCoarse, &
          p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
          ctrafoCoarse, p_DcubPtsRefCoarse(:,1:ncubpc))

      CALL elprep_prepareSetForEvaluation (rintSubsetFine%revalElementSet,&
          cevalTagFine, p_rtriaFine, p_IelementRef(1:NELF), &
          ctrafoFine, p_DcubPtsRefFine(:,1:ncubp))
      p_Ddetj => rintSubsetFine%revalElementSet%p_Ddetj
      
      ! Calculate the values of the basis functions.
      CALL elem_generic_sim2 (p_relemDistCoarse%celement, &
          rintSubsetCoarse%revalElementSet, Bder, DbasCoarse)
      CALL elem_generic_sim2 (p_relemDistFine%celement, &
          rintSubsetFine%revalElementSet, Bder, DbasFine)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      
      ! Clear the local matrix
      Dentry = 0.0_DP

      ! Loop over the elements in the current set.
      NELF = 0
      DO IELC=1,nelementsToDo
      
        ! Get the index of the currently processed coarse mesh element
        IDXC = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that are refined from the
        ! current coarse mesh element
        NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

        ! And loop through all elements of the current refinement patch
        DO IELF = 1, NELREF
          
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,NELF+IELF))

            ! Now loop through all possible combinations of DOF's
            ! in the current cubature point. The outer loop
            ! loops through the "O"'s in the above picture,
            ! the test functions:
            DO IDOFE=1,indofFine
            
              ! Get the value of the (test) basis function 
              ! phi_i (our "O") in the cubature point:
              DB = DbasFine(IDOFE,DER_FUNC,ICUBP,NELF+IELF)*OM
              
              ! Perform an inner loop through the other DOF's
              ! (the "X"). 
              DO JDOFE=1,indofCoarse
              
                Dentry(JDOFE,IDOFE,NELF+IELF) = Dentry(JDOFE,IDOFE,NELF+IELF) + &
                         DB*DbasCoarse(JDOFE,DER_FUNC,&
                         ICUBP + (IELF-1)*ncubp,IELC)
              
              END DO ! JDOFE
            
            END DO ! IDOFE

          END DO ! ICUBP 
          
        END DO ! IELF
        
        NELF = NELF + NELREF

      END DO ! IELC

      ! Incorporate the local matrix into the global one.
      ! Kentry gives the position of the additive contributions in Dentry.
      DO IELF = 1, NELF
        DO IDOFE=1,indofFine
          DO JDOFE=1,indofCoarse
            p_DA(Kentry(JDOFE,IDOFE,IELF)) = &
              p_DA(Kentry(JDOFE,IDOFE,IELF)) + Dentry(JDOFE,IDOFE,IELF)
          END DO
        END DO
      END DO

      ! Release the element sets here
      CALL elprep_releaseElementSet(rintSubsetFine%revalElementSet)
      CALL elprep_releaseElementSet(rintSubsetCoarse%revalElementSet)
      
      ! Increase the number of done elements
      nelementsDone = nelementsDone + nelementsToDo
      
    END DO ! WHILE(nelementsDone .LE. NELC)
    !%OMP END DO
    
    ! Release memory
    DEALLOCATE(p_DcubPtsRefCoarse)
    DEALLOCATE(p_DcubPtsRefFine)
    DEALLOCATE(IdofsCoarse)
    DEALLOCATE(IdofsFine)
    DEALLOCATE(DbasCoarse)
    DEALLOCATE(DbasFine)
    DEALLOCATE(Kentry)
    DEALLOCATE(Dentry)
    DEALLOCATE(p_IelementRef)

    !%OMP END PARALLEL

  ! That's it
  
  END SUBROUTINE

END MODULE
