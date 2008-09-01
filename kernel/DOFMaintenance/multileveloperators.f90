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

module multileveloperators

  use fsystem
  use genoutput
  use linearsystemscalar
  use spatialdiscretisation
  use scalarpde
  use derivatives
  use cubature
  use collection
  use domainintegration
  use element
  use elementpreprocessing
  
  implicit none

!<types>
  
!<typeblock>
  
  ! A private type block. Used as memory block during the creation of matrices.
  type t_matrixmem
    ! A pointer to a 1D memory block of integers - receives the
    ! column identifiers of a matrix
    integer(I32), dimension(:), pointer :: p_Icol
    
    ! A pointer to a 1D memory block of integers - receives 
    ! indices of next entries in the list of each line in the matrix
    integer(I32), dimension(:), pointer :: p_Iindx
  end type
  
!</typeblock>

  private :: t_matrixmem

!</types>

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building matrices
  integer, parameter :: MLOP_NELEMSIM   = 100
  
!</constantblock>
!</constants>

contains


  !****************************************************************************

!<subroutine>

  subroutine mlop_create2LvlMatrixStruct (rdiscretisationCoarse,&
                      rdiscretisationFine, iformat, rmatrixScalar, imemguess)
  
!<description>
  ! This routine allows to calculate the structure of a finite-element matrix
  ! on the heap. The size of the matrix is determined dynamically.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationFine

  ! Format of the matrix structure to be created. One of the LSYSSC_xxxx
  ! constants.
  integer, intent(IN) :: iformat

  ! OPTIONAL: An initial guess about how much memory the matrix needs. If set 
  ! to 0 or not given, an initial guess of 16*NEQ (but at least 10000 matrix 
  ! entries) is assumed.
  integer(I32), intent(IN), optional :: imemGuess
  
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation
  ! combination. Memory fo the structure is allocated dynamically on the heap.
  type(t_matrixScalar), intent(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  integer(I32) :: imem
  
    imem = 0
    if (present(imemguess)) then
      imem = max(0,imemguess)
    end if
    
    ! Let's make sure that the discretisations are uniform - we don't support
    ! anything else right now.
    if ((rdiscretisationCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
        (rdiscretisationFine%ccomplexity .ne. SPDISC_UNIFORM)) then
        
      call output_line ('Discretisations must be uniform!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      call sys_halt()
      
    end if
    
    ! Which matrix structure do we have to create?
    select case (iformat) 
    
    case (LSYSSC_MATRIX9)
    
      ! Call the creation routine for structure 9:
      call mlop_create2LvlMatStruct9_uni (rdiscretisationCoarse,&
                           rdiscretisationFine,rmatrixScalar,imem)
     
    case (LSYSSC_MATRIX7)
    
      ! Call the creation routine for structure 9:
      call mlop_create2LvlMatStruct9_uni (rdiscretisationCoarse,&
                           rdiscretisationFine,rmatrixScalar,imem)

      ! Translate to matrix structure 7:
      call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
      
    case DEFAULT
      call output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      call sys_halt()
      
    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine mlop_build2LvlMassMatrix (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
  
!<description>
  ! This routine calculates the entries of a 2-Level mass matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(IN) :: bclear
  
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_matrixScalar) :: rmatrixBackup
  
    ! The matrix must be unsorted, otherwise we can't set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there's a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrixScalar%isortStrategy .gt. 0) then
      call output_line ('Matrix-structure must be unsorted!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      call sys_halt()
    end if

    ! Let's make sure that the discretisations are uniform - we don't support
    ! anything else right now.
    if ((rdiscretisationCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
        (rdiscretisationFine%ccomplexity .ne. SPDISC_UNIFORM)) then
        
      call output_line ('Discretisations must be uniform!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      call sys_halt()
      
    end if

    ! Which matrix structure do we have?
    select case (rmatrixScalar%cmatrixFormat) 
    case (LSYSSC_MATRIX9)
      call mlop_build2LvlMass9_uni (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
    
    case (LSYSSC_MATRIX7)
      ! Convert structure 7 to structure 9.For that purpose, make a backup of
      ! the original matrix...
      call lsyssc_duplicateMatrix (rmatrixScalar,rmatrixBackup,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      ! Convert the matrix 
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX9)
      
      ! Create the matrix in structure 9
      call mlop_build2LvlMass9_uni (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
                                     
      ! Convert back to structure 7
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX7)
      
      ! Copy the entries back to the original matrix and release memory.
      call lsyssc_duplicateMatrix (rmatrixBackup,rmatrixScalar,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
          
      call lsyssc_releaseMatrix (rmatrixBackup)
                                     
    case DEFAULT
      call output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      call sys_halt()
      
    end select


  end subroutine
    
  !****************************************************************************
  
!<subroutine>
  
  subroutine mlop_create2LvlMatStruct9_uni (rdiscrCoarse,rdiscrFine,&
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
  type(t_spatialDiscretisation), intent(IN), target :: rdiscrCoarse
  type(t_spatialDiscretisation), intent(IN), target :: rdiscrFine
  
  ! An initial guess about how much memory the matrix needs. If set to 0,
  ! an initial guess of 16*NEQ (but at least 10000 matrix entries) is assumed.
  integer(I32), intent(IN) :: imemGuess
  
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation.
  ! Memory fo rthe structure is allocated dynamically on the heap.
  type(t_matrixScalar), intent(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  integer(PREC_DOFIDX) :: NEQ, IEQ, IROW, JCOL, IPOS, istartIdx, NA, nmaxCol
  integer :: IDOFE, JDOFE, i, IHELP
  integer(PREC_ELEMENTIDX) :: IELC,IELF,IELIDX,NELREF
  logical :: BSORT
  
  ! An allocateable list of handles for memory blocks. Size is dynamically 
  ! increased if there are too many columns in the matrix.
  integer, dimension(:), pointer :: p_Ihcol, p_Ihindx, p_IhTmp
  integer(PREC_DOFIDX), dimension(:), pointer :: p_Isize, p_ISizeTmp
  
  ! An allocateable list of pointers to the memory blocks - corresponds
  ! to the handles in p_Ihcol/p_Ihindx
  type(t_matrixmem), dimension(:), pointer :: Rmemblock, RmemblockTmp
  
  ! Pointer to current KCOL memory block,
  ! pointer to current index memory block
  integer(I32), dimension(:), pointer :: p_Icol, p_Iindx
  
  ! Number of currently allocated pointers in Ihmemblock
  integer :: iblocks
  
  ! Currently active memory block
  integer :: icurrentblock
  
  ! Number of elements in the current coarse/fine mesh element distribution
  integer(PREC_ELEMENTIDX) :: NELC, NELF
  
  ! Size of memory blocks
  integer(PREC_DOFIDX) :: imemblkSize
  
  ! Blocksize in terms of NEQ for guessing memory.
  ! The initial guess for memory is iblkSize*iblkSize*NEQ and every time
  ! memory is needed, another iblkSize*NEQ elements are added.
  integer, parameter :: iblkSize = 4
  
  ! Number of memory blocks to allocate
  integer, parameter :: NmemBlkCount = 5

  ! Pointer to KLD, KCOL
  integer(I32), dimension(:), pointer :: p_KLD, p_KCOL
  
  ! Size of memory currently allocated
  integer(PREC_DOFIDX) :: iallocated
  
  ! An allocateable array accepting the DOF's of a set of elements.
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,MLOP_NELEMSIM), TARGET :: &
  !                  IdofsTest, IdofsTrial
  integer(PREC_DOFIDX), dimension(:,:), pointer :: p_IdofsTrial, p_IdofsTest
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  integer(I32), dimension(:), pointer :: p_IelementList,p_IelementRef
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  integer :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  integer :: nelementsTrial, nelementsTest
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  integer(I32), dimension(:), pointer :: p_IrefPatchIdx, p_IrefPatch

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
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Get NEQ - we need it for guessing memory...
    NEQ = rmatrixScalar%NEQ
    
    if (NEQ .eq. 0) then
      call output_line ('Empty matrix!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatStruct9_uni')
      call sys_halt()
    end if
    
    ! Allocate KLD...
    call storage_new1D ('mlop_create2LvlMatStruct9_uni', 'KLD', &
                        NEQ+1_I32, ST_INT, rmatrixScalar%h_KLD, &
                        ST_NEWBLOCK_NOINIT)
    ! This must be a storage_getbase, no lsyssc_getbase, since this is the
    ! matrix construction routine!
    call storage_getbase_int(rmatrixScalar%h_Kld,p_KLD)
    
    ! Allocate a list of handles and a list of pointers corresponding to it.
    ! Initially allocate NmemBlkCount pointers
    allocate(p_Ihcol(NmemBlkCount))
    allocate(p_Ihindx(NmemBlkCount))
    allocate(p_Isize(NmemBlkCount))
    allocate(Rmemblock(NmemBlkCount))
    
    ! Allocate the first memory block that receives a part of the
    ! temporary matrix structure.
    ! We make an initial guess of iblkSize*iblkSize*NEQ elements in the matrix,
    ! if imemguess is not given.
    if (imemguess .ne. 0) then 
      ! at least one element per line!
      iallocated = max(NEQ,imemguess) 
    else  
      iallocated = max(10000,iblkSize*iblkSize*NEQ)
    end if
    iblocks = 1
    imemblkSize = iblkSize*NEQ
    p_Isize(1) = iallocated
    
    ! imemblkSize = iallocated is necessary at the moment to simplify
    ! whether we leave a block or not.
    call storage_new1D ('mlop_create2LvlMatStruct9_uni', 'Ihicol', &
                        p_Isize(1), ST_INT, p_Ihcol(1), ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (p_Ihcol(1),p_Icol)

    ! The new index array must be filled with 0 - otherwise
    ! the search routine below won't work!
    call storage_new1D ('mlop_create2LvlMatStruct9_uni', 'p_Ihindx', &
                        p_Isize(1), ST_INT, p_Ihindx(1), ST_NEWBLOCK_ZERO)
    call storage_getbase_int (p_Ihindx(1),p_Iindx)
    
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
    
    do IEQ=1,NEQ
      p_Iindx(IEQ) = 0
      p_Icol(IEQ)  = 0
    end do
    
    ! The first NEQ elements are reserved. The matrix is assumed
    ! to have at least one element per line.
    NA = NEQ
    
    ! Activate the current coarse mesh element distribution
    p_relementDistribution => rdiscrCoarse%RelementDistr(1)

    ! Cancel if this element distribution is empty.
    if (p_relementDistribution%NEL .eq. 0) then
      call output_line('Element space is empty!')
      call sys_halt()
    end if

    ! Get the number of local DOF's for trial and test functions
    indofTrial = elem_igetNDofLoc(rdiscrCoarse%RelementDistr(1)%celement)
    indofTest = elem_igetNDofLoc(rdiscrFine%RelementDistr(1)%celement)
    
    ! Calculate the number of coarse mesh elements we want to process
    ! in one run.
    nelementsTrial = min(MLOP_NELEMSIM,p_relementDistribution%NEL) 
    
    ! Now calculate the number of fine mesh elements we want to process
    ! in one run. 
    select case(p_rtriaFine%ndim)
    case (1)
      nelementsTest = 2*nelementsTrial + 10
    case (2)
      nelementsTest = 4*nelementsTrial + 20
    case (3)
      nelementsTest = 8*nelementsTrial + 40
    end select
    
    ! Allocate an array saving a couple of DOF's for trial and test functions
    allocate(p_IdofsTrial(indofTrial,nelementsTrial))
    allocate(p_IdofsTest(indofTest,nelementsTest))
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that the trial functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)
    
    ! And allocate the refinemed element list for the test functions
    allocate(p_IelementRef(nelementsTest))
    
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
    do while(nelementsDone .lt. NELC)
    
      ! We always try to handle nelementsTrial elements simultaneously.
      ! Of course, we will not handle more elements than the coarse
      ! mesh discretisation has.
      nelementsToDo = min(NELC-nelementsDone, nelementsTrial)
      
      ! Now comes the interesting part - we have to ensure that the DOF-mapping
      ! of the fine mesh discretisation fits into our DOF-array.
      ! If, for example, a coarse mesh quad was refined into more than 4 fine
      ! mesh quads, then it might happen that we cannot handle nelementsToDo
      ! coarse mesh elements at once, but we need to decrease nelementsToDo.
      
      NELF = 0
      do IELC = 1, nelementsToDo
      
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
        if((NELF+NELREF) .gt. nelementsTest) then
          nelementsToDo = IELC-1
          exit
        end if
        
        ! Copy the indices of the elements into the element list for the
        ! fine mesh discretisation
        do IELF = 1, NELREF
          p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IELIDX)+IELF-1)
        end do
        
        ! Add the number of refined elements to the counter
        NELF = NELF + NELREF
      
      end do
      
      ! If nelementsToDo is 0, then we have a serious problem...
      if (nelementsToDo .le. 0) then
        print *, "ERROR: mlop_create2LvlMatStruct9_uni"
        print *, "nelementsToDo = 0 !!!"
        call sys_halt()
      end if
      
      ! Call the DOF-mapping routine for the coarse and fine mesh
      call dof_locGlobMapping_mult(rdiscrCoarse, &
          p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
          p_IdofsTrial)
      call dof_locGlobMapping_mult(rdiscrFine, p_IelementRef(1:NELF), &
          p_IdofsTest)
      
      ! Reset the counter
      NELF = 0
      
      ! Loop through all the elements of the coarse mesh set
      do IELC = 1, nelementsToDo
      
        ! Get the index of the currently processed coarse mesh element
        IELIDX = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that are refined from the
        ! current coarse mesh element
        NELREF = p_IrefPatchIdx(IELIDX+1) - p_IrefPatchIdx(IELIDX)

        ! And loop through all elements of the current refinement patch
        do IELF = 1, NELREF
        
          ! For building the local matrices, we have first to
          ! loop through the test functions (the "O"'s), as these
          ! define the rows in the matrix.
          do IDOFE=1,indofTest

            ! The DOF IDOFE is now our "O".
            ! This global DOF gives us the row we have to build.
            IROW = p_IdofsTest(IDOFE,NELF+IELF)
            
            ! Now we loop through the other DOF's on the current element
            ! (the "X"'s).
            ! All these have common support with our current basis function
            ! and will therefore give an additive value to the global
            ! matrix.

            do JDOFE=1,indofTrial
              
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
              if (icurrentblock .ne. 1) then
                icurrentblock = 1
                istartidx = 0
                p_Icol => Rmemblock(1)%p_Icol
                p_Iindx => Rmemblock(1)%p_Iindx
              end if
              
              ! Is the list empty?

              if (p_Icol(IROW).eq.0) then

                ! Yes, row IROW is empty at the moment. Add the column as
                ! head of the list of row IROW.

                p_Icol(IROW) = JCOL
                
              else

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
                
                searchloop: do while ( (p_Icol(IPOS-istartIdx)) .ne. JCOL)
              
                  ! Did we reach the end of the list? Then we have to insert
                  ! a new element...
                  if (p_Iindx(IPOS-istartIdx) .eq. 0) then
                  
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
                    
                    if (NA .gt. iallocated) then
                     
                      ! Hmmm, we have to allocate more memory.
                      ! Do we have enough pointers left or do we have
                      ! to enlarge our list?
                      
                      if (iblocks .ge. size(p_Ihcol)) then 
                      
                        ! Not enough blocks, we have to reallocate the pointer lists!
                        allocate (p_IhTmp(iblocks+NmemBlkCount))
                        p_IhTmp(1:iblocks) = p_Ihcol(1:iblocks)
                        deallocate(p_Ihcol)
                        p_Ihcol => p_IhTmp

                        allocate (p_IhTmp(iblocks+NmemBlkCount))
                        p_IhTmp(1:iblocks) = p_Ihindx(1:iblocks)
                        deallocate(p_Ihindx)
                        p_Ihindx => p_IhTmp
                      
                        allocate (p_IsizeTmp(iblocks+NmemBlkCount))
                        p_IsizeTmp(1:iblocks) = p_Isize(1:iblocks)
                        deallocate(p_Isize)
                        p_Isize => p_IsizeTmp

                        allocate (RmemblockTmp(iblocks+NmemBlkCount))
                        RmemblockTmp(1:iblocks) = Rmemblock(1:iblocks)
                        deallocate(Rmemblock)
                        Rmemblock => RmemblockTmp
                        
                        ! Now we have enough blocks again.
                      end if

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

                      call storage_new1D ('mlop_create2LvlMatStruct9_uni', 'Ihicol', &
                                          p_Isize (iblocks), ST_INT, p_Ihcol(iblocks), &
                                          ST_NEWBLOCK_NOINIT)
                      call storage_getbase_int (p_Ihcol(iblocks),p_Icol)

                      ! The new index array must be filled with 0 - otherwise
                      ! the search routine below won't work!
                      call storage_new1D ('mlop_create2LvlMatStruct9_uni', 'p_Ihindx', &
                                          p_Isize (iblocks), ST_INT, p_Ihindx(iblocks), &
                                          ST_NEWBLOCK_ZERO)
                      call storage_getbase_int (p_Ihindx(iblocks),p_Iindx)
                      
                      Rmemblock(iblocks)%p_Icol => p_Icol
                      Rmemblock(iblocks)%p_Iindx => p_Iindx

                      iallocated = iallocated + p_Isize (iblocks)
                      
                    else
                      ! Be careful when leaving the current memory block
                      ! for insertion of an element at position NA!
                      !
                      ! If the new position is not in the current block,
                      ! it's in the last block... and so set the pointer
                      ! and indices appropriately!
                      
                      if ( NA .gt. (istartidx+p_Isize (icurrentblock))) then
                        istartidx = iallocated-p_Isize(iblocks)
                        icurrentblock = iblocks
                        p_Icol => Rmemblock(iblocks)%p_Icol
                        p_Iindx => Rmemblock(iblocks)%p_Iindx
                      end if
                      
                    end if
                  
                    ! Append JCOL to p_Icol
                    p_Icol(NA-istartIdx) = JCOL
                    
                    ! We have to make sure that p_Indx(NA)=0 to indicate the end of
                    ! the list. Ok, this is trivial because we allocated it with 
                    ! storage_new, configured to fill the memory with 0, so it is 0.
                    !
                    ! The searchloop ends here, continue with next JDOFE
                    
                    exit 
                  
                  else
                  
                    ! No, the list does not end here.
                    ! Take the next element in the list
                    IPOS = p_Iindx(IPOS-istartidx)
                    
                    ! Be careful when leaving the current memory block
                    do while ( IPOS .gt. (istartidx+p_Isize (icurrentblock)) )
                    
                      ! go to the next memory block and search there
                      istartidx = istartidx+p_Isize(icurrentblock)
                      icurrentblock = icurrentblock+1
                      p_Icol => Rmemblock(icurrentblock)%p_Icol
                      p_Iindx => Rmemblock(icurrentblock)%p_Iindx
                      
                    end do ! IPOS .GT. (istartidx+p_Isize (iblocks))
                  
                  end if ! p_Iindx(IPOS) = 0
                  
                end do searchloop
              
              end if ! p_Icol(IROW) = 0
                  
            end do ! JDOFE
          
          end do ! IDOFE
        
        end do ! IELF

        ! Add the number of refined elements to the counter
        NELF = NELF + NELREF
      
      end do ! IELC
      
      ! Add the number of processed coarse mesh elements
      nelementsDone = nelementsDone + nelementsToDo
    
    end do ! WHILE(nelementsDone .LT. NEL)
    
    ! Release the fine mesh element list
    deallocate(p_IelementRef)

    ! Clean up the DOF's arrays    
    deallocate(p_IdofsTest)
    deallocate(p_IdofsTrial)
    
    ! Ok, p_Icol is built. The hardest part is done!
    ! Now build KCOL by collecting the entries in the linear lists of 
    ! each row.
    !
    ! At first, as we now NA, we can allocate the real KCOL now!
    
    call storage_new1D ('mlop_create2LvlMatStruct9_uni', 'KCOL', &
                        NA, ST_INT, rmatrixScalar%h_KCOL, &
                        ST_NEWBLOCK_NOINIT)
    ! This must be a storage_getbase, no lsyssc_getbase, since this is the
    ! matrix construction routine!
    call storage_getbase_int (rmatrixScalar%h_Kcol,p_KCOL)
    
    ! Save NA in the matrix structure
    rmatrixScalar%NA = NA
    
    ! Set back NA to 0 at first.
    NA=0
        
    ! Loop through all of the NEQ linear lists:
    do IEQ=1,NEQ
    
      ! We are at the head of the list, now we have to walk
      ! through it to append the entries to KCOL.
      ! We always start in the first memory block.
      if (icurrentblock .ne. 1) then
        icurrentblock = 1
        istartidx = 0
        p_Icol => Rmemblock(1)%p_Icol
        p_Iindx => Rmemblock(1)%p_Iindx
      end if
      
      ! Add the head of the list to KCOL:
      NA=NA+1
      p_KCOL(NA)=p_Icol(IEQ)

      ! Set KLD appropriately:
      p_KLD(IEQ) = NA
      
      IPOS = IEQ
      
      do while (p_Iindx(IPOS-istartidx).ne.0)
      
        ! Get the position of the next entry in p_Icol:
        IPOS = p_Iindx(IPOS-istartidx)
        
        ! Be careful when leaving the current memory block
        do while ( IPOS .gt. (istartidx+p_Isize (icurrentblock)) )
        
          ! go to the next memory block and search there
          istartidx = istartidx+p_Isize(icurrentblock)
          icurrentblock = icurrentblock+1
          p_Icol => Rmemblock(icurrentblock)%p_Icol
          p_Iindx => Rmemblock(icurrentblock)%p_Iindx
          
        end do ! IPOS .GT. (istartidx+p_Isize (iblocks))
        
        ! Add the column number to the row in KCOL:
        NA=NA+1
        p_KCOL(NA)=p_Icol(IPOS-istartidx)
      
      end do ! KINDX(IPOS) <> 0

    end do ! IEQ
    
    ! Append the final entry to KLD:
    p_KLD(NEQ+1)=NA+1
    
    ! Sort entries on KCOL separately for each row.
    ! This is a small bubble-sort...
    !
    ! Loop through all rows:
    nmaxCol = 0

    do IEQ=1,NEQ

      ! Repeat until everything is sorted.
      
      BSORT=.false.
      do while (.not. BSORT)
      
        BSORT=.true.

        !  Loop through the line 

        do JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-2
        
          ! If the next element is larger...
        
          if (p_KCOL(JCOL) .gt. p_KCOL(JCOL+1)) then
          
            ! Change position of the current and next element
          
            IHELP=p_KCOL(JCOL)
            p_KCOL(JCOL)=p_KCOL(JCOL+1)
            p_KCOL(JCOL+1)=IHELP
            
            ! And repeat the sorting of that line
            
            BSORT=.false.
            
          end if
          
        end do ! JCOL
        
      end do ! (not BSORT)      

      ! Grab the largest column number. As the current line is sorted,
      ! we can find this using the end of the line.
      nmaxCol = max(nmaxCol,p_Kcol(p_Kld(IEQ+1)-1))

    end do ! IEQ
    
    ! HOORAY, THAT'S IT!
    ! Deallocate all temporary memory...
    
    do i=iblocks,1,-1
      call storage_free(p_Ihcol(i))
      call storage_free(p_Ihindx(i))
    end do
    
    deallocate(Rmemblock)
    deallocate(p_Isize)
    deallocate(p_Ihindx)
    deallocate(p_Ihcol)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine mlop_build2LvlMass9_uni (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar)
  
!<description>
  ! This routine calculates the entries of a 2-Level mass matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(IN) :: bclear
  
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>


  ! local variables
  integer :: i,k,JDFG, ICUBP, NELC,NELF
  integer(I32) :: IELC,IELF, IDXC, NELREF, IDOFE, JDOFE
  integer(PREC_DOFIDX) :: JCOL0,JCOL
  real(DP) :: OM, DB
  
  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! Cubature point coordinates on the reference element
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp,ncubpc
  
  ! Pointer to KLD, KCOL, DA
  integer(I32), dimension(:), pointer :: p_KLD, p_KCOL
  real(DP), dimension(:), pointer :: p_DA
  
  ! An allocateable array accepting the DOF's of a set of elements.
  integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsCoarse, IdofsFine
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: DbasCoarse,DbasFine
  
  ! Number of entries in the matrix - for quicker access
  integer(PREC_DOFIDX) :: NA
  integer(I32) :: NEQ
  
  ! Type of transformation from the reference to the real element 
  integer :: ctrafoCoarse, ctrafoFine
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevalTagCoarse, cevalTagFine
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofCoarse, indofFine
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  integer(I32), dimension(:), pointer :: p_IelementList, p_IelementRef
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer(PREC_DOFIDX), dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry

  ! An array that takes coordinates of the cubature formula on the reference element
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRefFine, p_DcubPtsRefCoarse
  
  ! Pointer to the jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relemDistCoarse, p_relemDistFine
  
  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  integer :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  integer :: nelementsCoarse, nelementsFine
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  integer(I32), dimension(:), pointer :: p_IrefPatchIdx, p_IrefPatch
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines as well as element evaluation routines.
  type(t_domainIntSubset) :: rintSubsetCoarse, rintSubsetFine

    ! We only need the function values as we want to assemble a mass matrix.
    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get information about the matrix:
    NA = rmatrixScalar%NA
    NEQ = rmatrixScalar%NEQ
    
    ! We need KCOL/KLD of our matrix
    if ((rmatrixScalar%h_KCOL .eq. ST_NOHANDLE) .or. &
        (rmatrixScalar%h_KLD .eq. ST_NOHANDLE)) then
      call output_line ('No discretisation structure! Cannot assemble matrix!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMass9_uni')
      call sys_halt()
    end if
    
    call lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    
    ! Check if the matrix entries exist. If not, allocate the matrix.
    if (rmatrixScalar%h_DA .eq. ST_NOHANDLE) then

      ! Clear the entries in the matrix - we need to start with zero
      ! when assembling a new matrix!
      call storage_new1D ('mlop_build2LvlMass9_uni', 'DA', &
                          NA, ST_DOUBLE, rmatrixScalar%h_DA, &
                          ST_NEWBLOCK_ZERO)
      call lsyssc_getbase_double (rmatrixScalar,p_DA)

    else
    
      call lsyssc_getbase_double (rmatrixScalar,p_DA)

      ! If desired, clear the matrix before assembling.
      if (bclear) then
        call lalg_clearVectorDble (p_DA)
      end if
      
    end if
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriaCoarse => rdiscretisationCoarse%p_rtriangulation
    p_rtriaFine => rdiscretisationFine%p_rtriangulation

    ! Get the refinement patch arrays from the fine triangulation
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Activate the current element distributions
    p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(1)
    p_relemDistFine => rdiscretisationFine%RelementDistr(1)
  
    ! Get the number of local DOF's for trial and test functions
    indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
    indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
      
    ! Calculate the number of coarse mesh elements we want to process
    ! in one run.
    nelementsCoarse = min(MLOP_NELEMSIM,p_relemDistCoarse%NEL) 
    
    ! Now calculate the number of fine mesh elements we want to process
    ! in one run. 
    select case(p_rtriaFine%ndim)
    case (1)
      nelementsFine = 2*nelementsCoarse
    case (2)
      nelementsFine = 4*nelementsCoarse
    case (3)
      nelementsFine = 8*nelementsCoarse
    end select
    
    ! Allocate an array saving a couple of DOF's for trial and test functions
    allocate(IdofsCoarse(indofCoarse,nelementsCoarse))
    allocate(IdofsFine(indofFine,nelementsFine))
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that the trial functions
    call storage_getbase_int (p_relemDistCoarse%h_IelementList, p_IelementList)
    
    ! And allocate the refinemed element list for the test functions
    allocate(p_IelementRef(nelementsFine))
    
    ! Get the number of coarse mesh elements there.
    NELC = p_relemDistCoarse%NEL
      
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
    ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
    
    ! Initialise the cubature formula, get cubature weights and point
    ! coordinates on the reference element of the fine mesh
    call cub_getCubPoints(p_relemDistFine%ccubTypeBilForm, ncubp, Dxi, Domega)
    
    ! Allocate some memory to hold the cubature points on the fine mesh
    allocate(p_DcubPtsRefFine(trafo_igetReferenceDimension(ctrafoFine),ncubp))

    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,ncubp
      do k=1,ubound(p_DcubPtsRefFine,1)
        p_DcubPtsRefFine(k,i) = Dxi(i,k)
      end do
    end do
    
    ! Now we need to transform the points from the fine mesh into the coarse mesh
    ! Please note that the following trick does only work for 2-level ordered
    ! meshes!
    select case(p_rtriaFine%ndim)
    case (1)
      ncubpc = 2*ncubp
      allocate(p_DcubPtsRefCoarse(trafo_igetReferenceDimension(ctrafoCoarse),ncubpc))
      call trafo_mapCubPtsRef2LvlEdge1D(ncubp,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
      
    case (2)
      ncubpc = 4*ncubp
      allocate(p_DcubPtsRefCoarse(trafo_igetReferenceDimension(ctrafoCoarse),ncubpc))
      call trafo_mapCubPtsRef2LvlQuad2D(ncubp,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

    case (3)
      ncubpc = 8*ncubp
      allocate(p_DcubPtsRefCoarse(trafo_igetReferenceDimension(ctrafoCoarse),ncubpc))
      call trafo_mapCubPtsRef2LvlHexa3D(ncubp,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
      
    end select
    
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    allocate(DbasFine(indofFine,&
             elem_getMaxDerivative(p_relemDistFine%celement),&
             ncubp,nelementsFine))
    allocate(DbasCoarse(indofCoarse,&
             elem_getMaxDerivative(p_relemDistCoarse%celement), &
             ncubpc,nelementsCoarse))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    allocate(Kentry(indofCoarse,indofFine,nelementsFine))
    allocate(Dentry(indofCoarse,indofFine,nelementsFine))
    
    ! Loop over the elements - blockwise.
    nelementsDone = 0
    do while(nelementsDone .lt. NELC)
    
      ! We always try to handle nelementsTrial elements simultaneously.
      ! Of course, we will not handle more elements than the coarse
      ! mesh discretisation has.
      nelementsToDo = min(NELC-nelementsDone, nelementsCoarse)
      
      ! Now comes the interesting part - we have to ensure that the DOF-mapping
      ! of the fine mesh discretisation fits into our DOF-array.
      ! If, for example, a coarse mesh quad was refined into more than 4 fine
      ! mesh quads, then it might happen that we cannot handle nelementsToDo
      ! coarse mesh elements at once, but we need to decrease nelementsToDo.
      NELF = 0
      do IELC = 1, nelementsToDo
      
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
        if((NELF+NELREF) .gt. nelementsFine) then
          nelementsToDo = IELC-1
          exit
        end if
        
        ! Copy the indices of the elements into the element list for the
        ! fine mesh discretisation
        do IELF = 1, NELREF
          p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IDXC)+IELF-1)
        end do
        
        ! Add the number of refined elements to the counter
        NELF = NELF + NELREF
      
      end do
      
      ! If nelementsToDo is 0, then we have a serious problem...
      if (nelementsToDo .le. 0) then
        print *, "ERROR: mlop_build2LvlMass9_uni"
        print *, "nelementsToDo = 0 !!!"
        call sys_halt()
      end if
      
      ! Call the DOF-mapping routine for the coarse and fine mesh
      call dof_locGlobMapping_mult(rdiscretisationCoarse, &
          p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
          IdofsCoarse)
      call dof_locGlobMapping_mult(rdiscretisationFine, p_IelementRef(1:NELF), &
          IdofsFine)
      
      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
      NELF = 0
      do IELC = 1, nelementsToDo
      
        ! Get the index of the currently processed coarse mesh element
        IDXC = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that are refined from the
        ! current coarse mesh element
        NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

        ! And loop through all elements of the current refinement patch
        do IELF = 1, NELREF
        
          ! Get the index of the currently processed fine mesh element
          !IDXF = p_IelementRef(NELF+IELF)
      
          ! For building the local matrices, we have first to
          ! loop through the test functions (the "O"'s), as these
          ! define the rows in the matrix.
          do IDOFE=1,indofFine
          
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
            do JDOFE=1,indofCoarse
              
              ! Get the global DOF of the "X" which interacts with 
              ! our "O".
              JDFG = IdofsCoarse(JDOFE,IELC)
              
              ! Starting in JCOL0 (which points to the beginning of
              ! the line initially), loop through the elements in
              ! the row to find the position of column IDFG.
              ! Jump out of the DO loop if we find the column.
              do JCOL=JCOL0,NA
                if (p_KCOL(JCOL) .eq. JDFG) exit
              end do

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
              
            end do ! IDOFE
            
          end do ! JDOFE
        
        end do ! IELF
        
        NELF = NELF + NELREF
        
      end do ! IELC
      
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
                      
      cevalTagCoarse = ior(cevalTagCoarse,EL_EVLTAG_REFPOINTS)
      cevalTagFine   = ior(cevalTagFine,EL_EVLTAG_REFPOINTS)

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (rintSubsetCoarse%revalElementSet,&
          cevalTagCoarse, p_rtriaCoarse, &
          p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
          ctrafoCoarse, p_DcubPtsRefCoarse(:,1:ncubpc))

      call elprep_prepareSetForEvaluation (rintSubsetFine%revalElementSet,&
          cevalTagFine, p_rtriaFine, p_IelementRef(1:NELF), &
          ctrafoFine, p_DcubPtsRefFine(:,1:ncubp))
      p_Ddetj => rintSubsetFine%revalElementSet%p_Ddetj
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relemDistCoarse%celement, &
          rintSubsetCoarse%revalElementSet, Bder, DbasCoarse)
      call elem_generic_sim2 (p_relemDistFine%celement, &
          rintSubsetFine%revalElementSet, Bder, DbasFine)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      
      ! Clear the local matrix
      Dentry = 0.0_DP

      ! Loop over the elements in the current set.
      NELF = 0
      do IELC=1,nelementsToDo
      
        ! Get the index of the currently processed coarse mesh element
        IDXC = p_IelementList(nelementsDone+IELC)
        
        ! Get the number of fine mesh elements that are refined from the
        ! current coarse mesh element
        NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

        ! And loop through all elements of the current refinement patch
        do IELF = 1, NELREF
          
          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,NELF+IELF))

            ! Now loop through all possible combinations of DOF's
            ! in the current cubature point. The outer loop
            ! loops through the "O"'s in the above picture,
            ! the test functions:
            do IDOFE=1,indofFine
            
              ! Get the value of the (test) basis function 
              ! phi_i (our "O") in the cubature point:
              DB = DbasFine(IDOFE,DER_FUNC,ICUBP,NELF+IELF)*OM
              
              ! Perform an inner loop through the other DOF's
              ! (the "X"). 
              do JDOFE=1,indofCoarse
              
                Dentry(JDOFE,IDOFE,NELF+IELF) = Dentry(JDOFE,IDOFE,NELF+IELF) + &
                         DB*DbasCoarse(JDOFE,DER_FUNC,&
                         ICUBP + (IELF-1)*ncubp,IELC)
              
              end do ! JDOFE
            
            end do ! IDOFE

          end do ! ICUBP 
          
        end do ! IELF
        
        NELF = NELF + NELREF

      end do ! IELC

      ! Incorporate the local matrix into the global one.
      ! Kentry gives the position of the additive contributions in Dentry.
      do IELF = 1, NELF
        do IDOFE=1,indofFine
          do JDOFE=1,indofCoarse
            p_DA(Kentry(JDOFE,IDOFE,IELF)) = &
              p_DA(Kentry(JDOFE,IDOFE,IELF)) + Dentry(JDOFE,IDOFE,IELF)
          end do
        end do
      end do

      ! Release the element sets here
      call elprep_releaseElementSet(rintSubsetFine%revalElementSet)
      call elprep_releaseElementSet(rintSubsetCoarse%revalElementSet)
      
      ! Increase the number of done elements
      nelementsDone = nelementsDone + nelementsToDo
      
    end do ! WHILE(nelementsDone .LE. NELC)
    !%OMP END DO
    
    ! Release memory
    deallocate(p_DcubPtsRefCoarse)
    deallocate(p_DcubPtsRefFine)
    deallocate(IdofsCoarse)
    deallocate(IdofsFine)
    deallocate(DbasCoarse)
    deallocate(DbasFine)
    deallocate(Kentry)
    deallocate(Dentry)
    deallocate(p_IelementRef)

    !%OMP END PARALLEL

  ! That's it
  
  end subroutine

end module
