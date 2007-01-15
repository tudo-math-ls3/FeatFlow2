!##############################################################################
!# ****************************************************************************
!# <name> bilinearformevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the discretisation of bilinear forms,
!# i.e. the creation of matrices and matrix structures. It contains the
!# following set of routines:
!#
!# 1.) bilf_createMatrixStructure
!#     -> Creates a 'scalar' matrix structure in a specified matrix
!#        format according to a discretisation.
!#
!# 2.) bilf_buildMatrixScalar
!#     -> Assembles the entries of a matrix, which structure was build
!#        with bilf_createMatrixStructure before.
!#
!# </purpose>
!##############################################################################

MODULE bilinearformevaluation

  USE fsystem
  USE linearsystemscalar
  USE spatialdiscretisation
  USE scalarpde
  USE derivatives
  USE cubature
  USE collection
  USE domainintegration
  
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

!<constantblock description="Method identifiers for construction of matrix structure.">

  ! Element-based matrix construction. This is the standard matrix construction method. 
  INTEGER, PARAMETER :: BILF_MATC_ELEMENTBASED = 0
  
  ! Edge-based matrix construction. The matrix stencil is extended in such a way,
  ! that the DOF's of one element may interact with the DOF's of all other elements
  ! that are adjacent via one of the edges.
  INTEGER, PARAMETER :: BILF_MATC_EDGEBASED    = 1

  ! Vertex-based matrix construction. The matrix stencil is extended in such a way,
  ! that the DOF's of one element may interact with the DOF's of all other elements
  ! that are adjacent via one of the corner vertices.
  INTEGER, PARAMETER :: BILF_MATC_VERTEXBASED  = 2

!</constantblock>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building matrices
  INTEGER, PARAMETER :: BILF_NELEMSIM   = 1000
  
!</constantblock>
!</constants>

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE bilf_createMatrixStructure (rdiscretisation,iformat,rmatrixScalar, &
                                         cconstrType,imemguess)
  
!<description>
  ! This routine allows to calculate the structure of a finite-element matrix
  ! on the heap. The size of the matrix is determined dynamically.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! Format of the matrix structure to be created. One of the LSYSSC_xxxx
  ! constants.
  INTEGER, INTENT(IN) :: iformat

  ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to specify
  ! the matrix construction method. If not specified,
  ! BILF_MATC_ELEMENTBASED is used.
  INTEGER, INTENT(IN), OPTIONAL      :: cconstrType

  ! OPTIONAL: An initial guess about how much memory the matrix needs. If set 
  ! to 0 or not given, an initial guess of 16*NEQ (but at least 10000 matrix 
  ! entries) is assumed.
  INTEGER(I32), INTENT(IN), OPTIONAL :: imemGuess
  
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation.
  ! Memory fo the structure is allocated dynamically on the heap.
  TYPE(t_matrixScalar), INTENT(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  INTEGER(I32) :: imem
  INTEGER :: ccType
  
  imem = 0
  IF (PRESENT(imemguess)) THEN
    imem = MAX(0,imemguess)
  END IF
  
  ccType = BILF_MATC_ELEMENTBASED
  IF (PRESENT(cconstrType)) ccType = cconstrType
  
  ! Do we have a not too complex triangulation? Would simplify a lot...
  IF ( (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) .OR. &
       (rdiscretisation%ccomplexity .EQ. SPDISC_CONFORMAL) ) THEN
  
    ! Which matrix structure do we have to create?
    SELECT CASE (iformat) 
    
    CASE (LSYSSC_MATRIX9)
    
      SELECT CASE (ccType)
      
      CASE (BILF_MATC_ELEMENTBASED)
        ! Call the creation routine for structure 9:
        CALL bilf_createMatStructure9_conf (rdiscretisation,rmatrixScalar,imem)
        
      CASE (BILF_MATC_EDGEBASED)
      
        IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
          CALL bilf_createMatStructure9eb_uni (rdiscretisation,rmatrixScalar,imem)
        ELSE
          PRINT *,'bilf_createMatrixStructure: Edge-based matrix constrution only for'//&
                  ' uniform discr., supported.'
          STOP
        END IF
        
      CASE DEFAULT
        PRINT *,'Invalid matrix construction method.'
        STOP
      END SELECT
      
    CASE (LSYSSC_MATRIX7)
    
      SELECT CASE (ccType)
      
      CASE (BILF_MATC_ELEMENTBASED)
      
        ! Call the creation routine for structure 9:
        CALL bilf_createMatStructure9_conf (rdiscretisation,rmatrixScalar,imem)

      CASE (BILF_MATC_EDGEBASED)
      
        IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
          CALL bilf_createMatStructure9eb_uni (rdiscretisation,rmatrixScalar,imem)
        ELSE
          PRINT *,'bilf_createMatrixStructure: Edge-based matrix constrution only for'//&
                  ' uniform discr., supported.'
          STOP
        END IF
        
      CASE DEFAULT
        PRINT *,'bilf_createMatrixStructure: Invalid matrix construction method.'
        STOP
      END SELECT
        
      ! Translate to matrix structure 7:
      CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
      
    CASE DEFAULT
      PRINT *,'bilf_createMatrixStructure: Not supported matrix structure!'
      STOP
    END SELECT
  
  ELSE
    PRINT *,'bilf_createMatrixStructure: General discretisation &
            & not implemented!'
    STOP
  END IF

  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE bilf_buildMatrixScalar (rform,bclear,rmatrixScalar,&
                                     fcoeff_buildMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix.
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! The matrix must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The bilinear form specifying the underlying PDE of the discretisation.
  TYPE(t_bilinearForm), INTENT(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  TYPE(t_collection), INTENT(IN), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  INCLUDE 'intf_coefficientMatrixSc.inc'
  OPTIONAL :: fcoeff_buildMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  TYPE(t_collection), POINTER :: p_rcollection
  
  ! Let p_rcollection point to rcollection - or NULL if it's not
  ! given.
  IF (PRESENT(rcollection)) THEN
    p_rcollection => rcollection
  ELSE
    p_rcollection => NULL()
  END IF

  ! The matrix must be unsorted, otherwise we can't set up the matrix.
  ! Note that we cannot switch off the sorting as easy as in the case
  ! of a vector, since there's a structure behind the matrix! So the caller
  ! has to make sure, the matrix is unsorted when this routine is called.
  IF (rmatrixScalar%isortStrategy .GT. 0) THEN
    PRINT *,'linf_buildMatrixScalar: Vector must be unsorted!'
    STOP
  END IF

  IF (.NOT. ASSOCIATED(rmatrixScalar%p_rspatialDiscretisation)) THEN
    PRINT *,'bilf_buildMatrixScalar: No discretisation associated!'
    STOP
  END IF

  ! Do we have a uniform triangulation? Would simplify a lot...
  SELECT CASE (rmatrixScalar%p_rspatialDiscretisation%ccomplexity)
  CASE (SPDISC_UNIFORM) 
    ! Uniform discretisation; only one type of elements, e.g. P1 or Q1
    SELECT CASE (rmatrixScalar%cdataType)
    CASE (ST_DOUBLE) 
      ! Which matrix structure do we have?
      SELECT CASE (rmatrixScalar%cmatrixFormat) 
      CASE (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,&  
                                         fcoeff_buildMatrixSc_sim,rcollection)
        !ELSE
        !  CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar)
        !END IF
      CASE (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9.
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,&  
                                       fcoeff_buildMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
                                       
      CASE DEFAULT
        PRINT *,'bilf_buildMatrix: Not supported matrix structure!'
        STOP
      END SELECT
    CASE DEFAULT
      PRINT *,'bilf_buildMatrix: Single precision matrices currently not supported!'
      STOP
    END SELECT
    
  CASE (SPDISC_CONFORMAL) 
    
    ! Conformal discretisation; may have mixed P1/Q1 elements e.g.
    SELECT CASE (rmatrixScalar%cdataType)
    CASE (ST_DOUBLE) 
      ! Which matrix structure do we have?
      SELECT CASE (rmatrixScalar%cmatrixFormat) 
      CASE (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,&  
                                         fcoeff_buildMatrixSc_sim,rcollection)
        !ELSE
        !  CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar)
        !END IF
        
      CASE (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,&  
                                       fcoeff_buildMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)

      CASE DEFAULT
        PRINT *,'bilf_buildMatrix: Not supported matrix structure!'
        STOP
      END SELECT
    CASE DEFAULT
      PRINT *,'bilf_buildMatrix: Single precision matrices currently not supported!'
      STOP
    END SELECT
  CASE DEFAULT
    PRINT *,'bilf_buildMatrix: General discretisation not implemented!'
    STOP
  END SELECT

  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE bilf_createMatStructure9_conf (rdiscretisation,rmatrixScalar,imemGuess)
  
!<description>
  ! This routine creates according to a given discretisation the matrix 
  ! structure of a structure-9 matrix. The discretisation is assumed to be
  ! conformal, i.e. the DOF's of different FE spaces in the trial space
  ! fit together. The function space for trial and test functions 
  ! may be different.
!</description>

!<input>
  
  ! The underlying discretisation structure which is to be used to
  ! create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
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
  INTEGER :: IDOFE, JDOFE, i, IHELP,NVE
  INTEGER(PREC_ELEMENTIDX) :: IEL, IELmax, IELset
  LOGICAL :: BSORT, bIdenticalTrialAndTest
  
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
  
  ! Number of currently active element distribution
  INTEGER :: icurrentElementDistr
  
  ! Currently active memory block
  INTEGER :: icurrentblock
  
  ! Size of memory blocks
  INTEGER(PREC_DOFIDX) :: imemblkSize
  
  ! Blocksize in terms of NEQ for guessing memory.
  ! The initial guess for memory is iblkSize*iblkSize*NEQ and every time
  ! memory is needed, another iblkSize*NEQ elements are added.
  INTEGER, PARAMETER :: iblkSize = 4
  
  ! Number of memory blocks to allocate
  INTEGER, PARAMETER :: NmemBlkCount = 5

  ! Pointer to KLD, KCOL, diagonal
  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL, p_Kdiagonal
  
  ! Size of memory currently allocated
  INTEGER(PREC_DOFIDX) :: iallocated
  
  ! An allocateable array accepting the DOF's of a set of elements.
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock

  ! The algorithm is: Test every DOF on one element against each other
  ! DOF on the same element and save the combination into a matrix
  ! in structure 9!
  !
  ! At first, initialise the structure-9 matrix:
  
  rmatrixScalar%p_rspatialDiscretisation => rdiscretisation
  rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9
  
  ! Get the #DOF's of the test space - as #DOF's of the test space is
  ! the number of equations in our matrix. The #DOF's in the trial space
  ! gives the number of columns of our matrix.
  rmatrixScalar%NCOLS         = dof_igetNDofGlob(rdiscretisation,.FALSE.)
  rmatrixScalar%NEQ           = dof_igetNDofGlob(rdiscretisation,.TRUE.)
  
  ! and get a pointer to the triangulation.
  p_rtriangulation => rdiscretisation%p_rtriangulation
  
  ! Get NEQ - we need it for guessing memory...
  NEQ = rmatrixScalar%NEQ
  
  IF (NEQ .EQ. 0) THEN
    PRINT *,'bilf_createMatrixStructure9_uni: Empty matrix!'
    STOP
  END IF
  
  ! Allocate KLD...
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'KLD', &
                      NEQ+1_I32, ST_INT, rmatrixScalar%h_KLD, &
                      ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rmatrixScalar%h_KLD,p_KLD)
  
  ! Allocate h_Kdiagonal
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'Kdiagonal', &
                      NEQ, ST_INT, rmatrixScalar%h_Kdiagonal, &
                      ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rmatrixScalar%h_Kdiagonal,p_Kdiagonal)
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocaing some arrays.
  nelementsPerBlock = MIN(BILF_NELEMSIM,p_rtriangulation%NEL)

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

  CALL storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                      p_Isize(1), ST_INT, p_Ihcol(1), ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (p_Ihcol(1),p_Icol)

  ! The new index array must be filled with 0 - otherwise
  ! the search routine below won't work!
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
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
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  
  DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
  
    ! Activate the current element distribution
    p_elementDistribution => rdiscretisation%RelementDistribution(icurrentElementDistr)

    ! Get the number of local DOF's for trial and test functions
    indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
    indofTest = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
      PRINT *,'bilf_createMatStructure9_conf: element spaces incompatible!'
      STOP
    END IF
    
    ! Allocate an array saving a couple of DOF's for trial and test functions
    !ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
    !ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))
    
    ! Test if trial/test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = &
      p_elementDistribution%itrialElement .EQ. p_elementDistribution%itestElement
      
    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    IF (bIdenticalTrialAndTest) THEN
      p_IdofsTrial => IdofsTest
    ELSE
      p_IdofsTrial => IdofsTrial
    END IF
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)
    

    ! Set the pointers/indices to the initial position. During the
    ! search for new DOF's, these might be changed if there's not enough
    ! memory in the first block.    
    icurrentblock = 1
    istartidx = 0
    p_Icol => Rmemblock(1)%p_Icol
    p_Iindx => Rmemblock(1)%p_Iindx
    
    ! Loop over the elements. 
    DO IELset = 1, p_rtriangulation%NEL, BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = MIN(p_rtriangulation%NEL,IELset-1+BILF_NELEMSIM)
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !        
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Call dof_locGlobMapping to get the global DOF's on our current
      ! element (the "X" and all "O"'s). 
      ! We don't need the local DOF's, so by setting IPAR=0,
      ! the call will only fill KDFG.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      ! Calculate the DOF's of the test functions:
      CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                  .TRUE.,IdofsTest)
                                   
      ! If the DOF's for the test functions are different, calculate them, too.
      IF (.NOT.bIdenticalTrialAndTest) THEN
        CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                    .FALSE.,IdofsTrial)
      END IF
      
      ! Loop through all the elements in the current set
      DO IEL=1,IELmax-IELset+1
        
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        DO IDOFE=1,indofTest

          ! The DOF IDOFE is now our "O".
          ! This global DOF gives us the row we have to build.
          
          IROW = IdofsTest(IDOFE,IEL)
          
          ! Now we loop through the other DOF's on the current element
          ! (the "X"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          DO JDOFE=1,indofTrial
            
            ! Get the global DOF - our "X". This gives the column number
            ! in the matrix where an entry occurs in row IROW (the line of 
            ! the current global DOF "O").
              
            JCOL = p_IdofsTrial(JDOFE,IEL)

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

                    CALL storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                                        p_Isize (iblocks), ST_INT, p_Ihcol(iblocks), &
                                        ST_NEWBLOCK_NOINIT)
                    CALL storage_getbase_int (p_Ihcol(iblocks),p_Icol)

                    ! The new index array must be filled with 0 - otherwise
                    ! the search routine below won't work!
                    CALL storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
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
      
      END DO ! IEL
    
    END DO ! IELset

    ! Clean up the DOF's arrays    
    !DEALLOCATE(IdofsTest)
    !DEALLOCATE(IdofsTrial)
    
  END DO ! icurrentElementDistr
  
  ! Ok, p_Icol is built. The hardest part is done!
  ! Now build KCOL by collecting the entries in the linear lists of 
  ! each row.
  !
  ! At first, as we now NA, we can allocate the real KCOL now!
  
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'KCOL', &
                      NA, ST_INT, rmatrixScalar%h_KCOL, &
                      ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rmatrixScalar%h_KCOL,p_KCOL)
  
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

    ! Grab the diagonal
    DO JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-1
      IF (p_KCOL(JCOL) .GE. IEQ) THEN
        p_Kdiagonal(IEQ) = JCOL
        EXIT
      END IF
    END DO   
    
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
  
  SUBROUTINE bilf_createMatStructure9eb_uni (rdiscretisation,rmatrixScalar,imemGuess)
  
!<description>
  ! This routine creates according to a given discretisation the matrix 
  ! structure of a structure-9 matrix. The discretisation is assumed to be
  ! uniform, i.e. there is only one combination of test- and trial functions
  ! allowed. The matrix is created by an edge-based approach, which increases
  ! the matrix stencil in such a way, that the DOF's of one element may
  ! interact with the DOF's of all elements that are adjacent via the
  ! edges.
!</description>

!<input>
  
  ! The underlying discretisation structure which is to be used to
  ! create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
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
  INTEGER :: IDOFE, JDOFE, i, IHELP,NVE, nelemBlockCount, IELidx
  INTEGER :: IELneighIdxJ
  INTEGER(PREC_ELEMENTIDX) :: IEL, IELmax, IELset
  LOGICAL :: BSORT, bIdenticalTrialAndTest
  
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
  
  ! Size of memory blocks
  INTEGER(PREC_DOFIDX) :: imemblkSize
  
  ! Blocksize in terms of NEQ for guessing memory.
  ! The initial guess for memory is iblkSize*iblkSize*NEQ and every time
  ! memory is needed, another iblkSize*NEQ elements are added.
  INTEGER, PARAMETER :: iblkSize = 4
  
  ! Number of memory blocks to allocate
  INTEGER, PARAMETER :: NmemBlkCount = 5

  ! Pointer to KLD, KCOL, diagonal
  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL, p_Kdiagonal
  
  ! Size of memory currently allocated
  INTEGER(PREC_DOFIDX) :: iallocated
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: IadjPtr, IadjElem
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock
  
  ! Adjacent elements
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_Kadj

  ! The algorithm is: Test every DOF on one element against each other
  ! DOF on the same element and save the combination into a matrix
  ! in structure 9!
  !
  ! At first, initialise the structure-9 matrix:
  
  rmatrixScalar%p_rspatialDiscretisation => rdiscretisation
  rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9
  
  ! Get the #DOF's of the test space - as #DOF's of the test space is
  ! the number of equations in our matrix. The #DOF's in the trial space
  ! gives the number of columns of our matrix.
  rmatrixScalar%NCOLS         = dof_igetNDofGlob(rdiscretisation,.FALSE.)
  rmatrixScalar%NEQ           = dof_igetNDofGlob(rdiscretisation,.TRUE.)
  
  ! and get a pointer to the triangulation.
  p_rtriangulation => rdiscretisation%p_rtriangulation
  
  CALL storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,p_Kadj)
  
  ! Get NEQ - we need it for guessing memory...
  NEQ = rmatrixScalar%NEQ
  
  IF (NEQ .EQ. 0) THEN
    PRINT *,'bilf_createMatrixStructure9_uni: Empty matrix!'
    STOP
  END IF
  
  ! Allocate KLD...
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'KLD', &
                      NEQ+1_I32, ST_INT, rmatrixScalar%h_KLD, &
                      ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rmatrixScalar%h_KLD,p_KLD)
  
  ! Allocate h_Kdiagonal
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'Kdiagonal', &
                      NEQ, ST_INT, rmatrixScalar%h_Kdiagonal, &
                      ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rmatrixScalar%h_Kdiagonal,p_Kdiagonal)
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocaing some arrays.
  nelementsPerBlock = MIN(BILF_NELEMSIM,p_rtriangulation%NEL)

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

  CALL storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                      p_Isize(1), ST_INT, p_Ihcol(1), ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (p_Ihcol(1),p_Icol)

  ! The new index array must be filled with 0 - otherwise
  ! the search routine below won't work!
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
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
  
  ! Activate the one and only element distribution
  p_elementDistribution => rdiscretisation%RelementDistribution(1)

  ! Get the number of local DOF's for trial and test functions
  indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
  indofTest = elem_igetNDofLoc(p_elementDistribution%itestElement)
  
  ! Get the number of corner vertices of the element
  NVE = elem_igetNVE(p_elementDistribution%itrialElement)
  IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
    PRINT *,'bilf_createMatStructure9_conf: element spaces incompatible!'
    STOP
  END IF
  
  ! Allocate the IadjCount array. This array counts for every element,
  ! how many elements are adjacent to that.
  ALLOCATE(IadjPtr(nelementsPerBlock+1))
  
  ! Allocate the IadjElem array. This collects for every element 
  ! the element number itself as well as all the adjacent elements. 
  ! (It's like a KLD-array...)
  ! IadjCount is a pointer into this array, so we can access directly 
  ! the numbers of the elements adjacent to one element.
  ! There is
  !   IadjElem(IadjPtr(IEL)) = IEL
  !   IadjElem(IadjPtr(IEL)+1..IadjPtr(IEL+1)-1) = adjacent elements
  ! As we know the maximum number of edges, we know the maximum size
  ! this array may need.
  ALLOCATE(IadjElem(nelementsPerBlock*(UBOUND(p_Kadj,1)+1)))

  ! Allocate an array saving a couple of DOF's for trial and test functions
  ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock*(UBOUND(p_Kadj,1)+1)))
  ALLOCATE(IdofsTest(indofTest,nelementsPerBlock*(UBOUND(p_Kadj,1)+1)))
  
  ! Test if trial/test functions are identical.
  ! We don't rely on bidenticalTrialAndTest purely, as this does not
  ! indicate whether there are identical trial and test functions
  ! in one block!
  bIdenticalTrialAndTest = &
    p_elementDistribution%itrialElement .EQ. p_elementDistribution%itestElement
    
  ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
  ! space IdofTest (if both spaces are identical). 
  ! We create a pointer for the trial space and not for the test space to
  ! prevent pointer-arithmetic in the innerst loop below!
  IF (bIdenticalTrialAndTest) THEN
    p_IdofsTrial => IdofsTest
  ELSE
    p_IdofsTrial => IdofsTrial
  END IF
  
  ! p_IelementList must point to our set of elements in the discretisation
  ! with that combination of trial/test functions
  CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                            p_IelementList)
  

  ! Set the pointers/indices to the initial position. During the
  ! search for new DOF's, these might be changed if there's not enough
  ! memory in the first block.    
  icurrentblock = 1
  istartidx = 0
  p_Icol => Rmemblock(1)%p_Icol
  p_Iindx => Rmemblock(1)%p_Iindx
  
  ! Loop over the elements. 
  DO IELset = 1, p_rtriangulation%NEL, BILF_NELEMSIM
  
    ! We always handle BILF_NELEMSIM elements simultaneously.
    ! How many elements have we actually here?
    ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
    ! elements simultaneously.
    
    IELmax = MIN(p_rtriangulation%NEL,IELset-1+BILF_NELEMSIM)
    
    ! --------------------- DOF SEARCH PHASE ------------------------
    
    ! In a first step, we search all the DOF's on one element
    ! and those on the elements adjacent via an edge to this.
    ! We want to get all the DOF numbers simultaneously!
    ! For this, we first set up an array that holds all the element
    ! numbers in our set as well as their neighbours.
    !
    ! Set up the IadjPtr array: Count 1 for the current element
    ! and count the number of adjacent elements to each element.
    ! Store the element number of one element to IadjElem and
    ! behind that the numbers of the adjacent elements.
    nelemBlockCount = 0
    DO IELidx=1,IELmax-IELset+1
      IEL = IELidx+IELset-1      ! actual element number
      
      nelemBlockCount = nelemBlockCount+1
      IadjPtr(IELidx) = nelemBlockCount
      IadjElem(nelemBlockCount) = IEL
      
      DO i=1,UBOUND(p_Kadj,1)
        IF (p_Kadj(i,IEL) .NE. 0) THEN
          nelemBlockCount = nelemBlockCount+1
          IadjElem(nelemBlockCount) = p_Kadj(i,IEL)
        END IF
      END DO
      
    END DO
    IadjPtr(IELmax-IELset+1+1) = nelemBlockCount+1
    
    ! nelemBlockCount is now the number of elements in the current
    ! block, consisting of the IELmax elements and their "edge-neighbours".
  
    ! The outstanding feature with finite elements is: A basis
    ! function for a DOF on one element has common support only
    ! with the DOF's on the same element! E.g. for Q1:
    !
    !        #. . .#. . .#. . .#
    !        .     .     .     .
    !        .  *  .  *  .  *  .
    !        #-----O-----O. . .#
    !        |     |     |     .
    !        |     | IEL |  *  .
    !        #-----X-----O. . .#
    !        |     |     |     .
    !        |     |     |  *  .
    !        #-----#-----#. . .#
    !        
    ! --> On element IEL, the basis function at "X" only interacts
    !     with the basis functions in "O". Elements in the 
    !     neighbourhood ("*") have no support, therefore we only have
    !     to collect all "O" DOF's.
    !
    ! Call dof_locGlobMapping to get the global DOF's on our current
    ! element (the "X" and all "O"'s). 
    ! We don't need the local DOF's, so by setting IPAR=0,
    ! the call will only fill KDFG.
    !
    ! More exactly, we call dof_locGlobMapping_mult to calculate all the
    ! global DOF's of our nelemBlockCount elements simultaneously.
    ! Calculate the DOF's of the test functions:
    CALL dof_locGlobMapping_mult(rdiscretisation, IadjElem(1:nelemBlockCount), &
                                .TRUE.,IdofsTest)
                                 
    ! If the DOF's for the test functions are different, calculate them, too.
    IF (.NOT.bIdenticalTrialAndTest) THEN
      CALL dof_locGlobMapping_mult(rdiscretisation, IadjElem(1:nelemBlockCount), &
                                  .FALSE.,IdofsTrial)
    END IF
  
    ! --------------------- DOF COMBINATION PHASE ------------------------
    
    ! Loop through all the elements in the current set
    DO IEL=1,IELmax-IELset+1
    
      ! For building the local matrices, we have first to
      ! loop through the test functions (the "O"'s), as these
      ! define the rows in the matrix.
      DO IDOFE=1,indofTest

        ! The DOF IDOFE is now our "O".
        ! This global DOF gives us the row we have to build.
        !
        ! The DOF's of element IEL start at position IadjPtr(IEL) in
        ! the IdofsTest array. 
        IROW = IdofsTest(IDOFE,IadjPtr(IEL)) 
        
        ! Loop through the "element-sets" in IadjElem, consisting
        ! of the element IEL itself as well as its neighbours - for the
        ! trial functions.
        DO IELneighIdxJ = IadjPtr(IEL),IadjPtr(IEL+1)-1

          ! Now we loop through the other DOF's on the current element
          ! (the "X"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          DO JDOFE=1,indofTrial
            
            ! Get the global DOF - our "X". This gives the column number
            ! in the matrix where an entry occurs in row IROW (the line of 
            ! the current global DOF "O").
              
            JCOL = p_IdofsTrial(JDOFE,IELneighIdxJ)

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
            ! p_Iindx(IPOS) is =0 if we reach the end of the linked list.
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

                    CALL storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                                        p_Isize (iblocks), ST_INT, p_Ihcol(iblocks), &
                                        ST_NEWBLOCK_NOINIT)
                    CALL storage_getbase_int (p_Ihcol(iblocks),p_Icol)

                    ! The new index array must be filled with 0 - otherwise
                    ! the search routine below won't work!
                    CALL storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
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
         
        END DO ! IELneighIdxJ
      
      END DO ! IDOFE
        
    END DO ! IEL
  
  END DO ! IELset

  ! Clean up the DOF's arrays    
  DEALLOCATE(IdofsTest)
  DEALLOCATE(IdofsTrial)

  ! --------------------- DOF COLLECTION PHASE ------------------------
    
  ! Ok, p_Icol is built. The hardest part is done!
  ! Now build KCOL by collecting the entries in the linear lists of 
  ! each row.
  !
  ! At first, as we now NA, we can allocate the real KCOL now!
  
  CALL storage_new1D ('bilf_createMatStructure9_conf', 'KCOL', &
                      NA, ST_INT, rmatrixScalar%h_KCOL, &
                      ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rmatrixScalar%h_KCOL,p_KCOL)
  
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

    ! Grab the diagonal
    DO JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-1
      IF (p_KCOL(JCOL) .GE. IEQ) THEN
        p_Kdiagonal(IEQ) = JCOL
        EXIT
      END IF
    END DO   
    
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
  
  DEALLOCATE(IadjElem)
  DEALLOCATE(IadjPtr)
  DEALLOCATE(Rmemblock)
  DEALLOCATE(p_Isize)
  DEALLOCATE(p_Ihindx)
  DEALLOCATE(p_Ihcol)
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE bilf_buildMatrix9d_conf (rdiscretisation,rform,bclear,rmatrixScalar,&
                                      fcoeff_buildMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix.
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance. The discretisation is assumed to be conformal, i.e. the DOF's
  ! of all finite elements must 'match'. Trial and test functions may be
  ! different.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! Double-precision version.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the matrix.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! The bilinear form specifying the underlying PDE of the discretisation.
  TYPE(t_bilinearForm), INTENT(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  TYPE(t_collection), INTENT(IN), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  INCLUDE 'intf_coefficientMatrixSc2.inc'
  OPTIONAL :: fcoeff_buildMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,i1,j,icurrentElementDistr,IDFG, ICUBP, IALBET, IA, IB, NVE
  LOGICAL :: bIdenticalTrialAndTest, bnonparTest, bnonparTrial
  INTEGER(I32) :: IEL, IELmax, IELset, IDOFE, JDOFE
  INTEGER(PREC_DOFIDX) :: JCOL0,JCOL
  REAL(DP) :: OM,AUX, DB
  
  ! Array to tell the element which derivatives to calculate
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrialTempl, BderTestTempl, BderTrial, BderTest
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp
  
  ! Pointer to KLD, KCOL, DA
  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest,DbasTrial
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasTrial
  
  ! Number of entries in the matrix - for quicker access
  INTEGER(PREC_DOFIDX) :: NA
  INTEGER(I32) :: NEQ
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dentry
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsRef

  ! An array receiving the coordinates of cubature points on
  ! the real element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsReal

  ! Pointer to the point coordinates to pass to the element function.
  ! Point either to DcubPtsRef or to DcubPtsReal, depending on whether
  ! the trial/test element is parametric or not.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTrial
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
  
  ! Array with coordinates of the corners that form the real element.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoords
  
  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Ddetj
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Djac
  
  ! Pointer to KVERT of the triangulation
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  
  ! Pointer to DCORVG of the triangulation
  REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
  
  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock
  
  ! Some variables to support nonconstant coefficients in the matrix.
  
  ! Pointer to the collection structure or to NULL()
  TYPE(t_collection), POINTER :: p_rcollection
  
  ! Pointer to the coefficients that are computed by the callback routine.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
  
  !REAL(DP), DIMENSION(11) :: DT
  
  !CHARACTER(LEN=20) :: CFILE
  
  IF (.NOT. ASSOCIATED(rmatrixScalar%p_rspatialDiscretisation)) THEN
    PRINT *,'bilf_buildMatrix9d_conf: No discretisation associated!'
    STOP
  END IF

  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDERxxxx
  ! according to these.

  !CALL ZTIME(DT(1))

  BderTrialTempl = .FALSE.
  BderTestTempl = .FALSE.
  
  ! Loop through the additive terms
  DO i=1,rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed! Build template's for BDER.
    ! We don't compute the actual BDER here, as there might be some special
    ! processing if trial/test functions are identical!
    !
    ! At first build the descriptors for the trial functions
    I1=rform%Idescriptors(1,I)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'bilf_buildMatrix9d_conf: Invalid descriptor'
      STOP
    ENDIF
    
    BderTrialTempl(I1)=.TRUE.

    ! Then those of the test functions
    I1=rform%Idescriptors(2,I)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'bilf_buildMatrix9d_conf: Invalid descriptor'
      STOP
    ENDIF
    
    BderTestTempl(I1)=.TRUE.
  END DO
  
  ! Get information about the matrix:
  NA = rmatrixScalar%NA
  NEQ = rmatrixScalar%NEQ
  
  ! We need KCOL/KLD of our matric
  CALL storage_getbase_int (rmatrixScalar%h_KCOL,p_KCOL)
  CALL storage_getbase_int (rmatrixScalar%h_KLD,p_KLD)
  
  ! Check if the matrix entries exist. If not, allocate the matrix.
  IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN

    ! Clear the entries in the matrix - we need to start with zero
    ! when assembling a new matrix!
    CALL storage_new1D ('bilf_buildMatrix9d_conf', 'DA', &
                        NA, ST_DOUBLE, rmatrixScalar%h_DA, &
                        ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double (rmatrixScalar%h_DA,p_DA)

  ELSE
  
    CALL storage_getbase_double (rmatrixScalar%h_DA,p_DA)

    ! If desired, clear the matrix before assembling.
    IF (bclear) THEN
      CALL lalg_clearVectorDble (p_DA)
    END IF
    
  END IF
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => rdiscretisation%p_rtriangulation
  
  ! Let p_rcollection point to rcollection - or NULL if it's not
  ! given.
  IF (PRESENT(rcollection)) THEN
    p_rcollection => rcollection
  ELSE
    p_rcollection => NULL()
  END IF

  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocaing some arrays.
  nelementsPerBlock = MIN(BILF_NELEMSIM,p_rtriangulation%NEL)
  
  ! Get a pointer to the KVERT and DCORVG array
  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
                             p_IverticesAtElement)
  CALL storage_getbase_double2D(p_rtriangulation%h_DcornerCoordinates, &
                             p_DcornerCoordinates)

  ! Allocate memory for corner coordinates
  ALLOCATE(DCoords(2,TRIA_MAXNVE2D,nelementsPerBlock))
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  !CALL ZTIME(DT(2))

  DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
  
    ! Activate the current element distribution
    p_elementDistribution => rdiscretisation%RelementDistribution(icurrentElementDistr)
  
    ! Get the number of local DOF's for trial and test functions
    indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
    indofTest = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
      PRINT *,'bilf_buildMatrix9d_conf: element spaces incompatible!'
      STOP
    END IF
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    CALL cub_getCubPoints(p_elementDistribution%ccubType, ncubp, Dxi, Domega)
    
    ! Allocate arrays accepting cubature point coordinates.
    ! It's at most as large as number of elements or length
    ! of the element set.
    ALLOCATE(DcubPtsRef(NDIM2D,ncubp,nelementsPerBlock))
    ALLOCATE(DcubPtsReal(NDIM2D,ncubp,nelementsPerBlock))
    
    ! Put the cubature point coordinates in the right format to the
    ! cubature-point array.
    ! Initialise all entries in DcubPtsRef with the same coordinates -
    ! as the cubature point coordinates are identical on all elements
    DO j=1,SIZE(DcubPtsRef,3)
      DO i=1,ncubp
        DcubPtsRef(1,i,j) = Dxi(i,1)
        DcubPtsRef(2,i,j) = Dxi(i,2)
      END DO
    END DO
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    ALLOCATE(Djac(NDIM2D*NDIM2D,ncubp,nelementsPerBlock))
    ALLOCATE(Ddetj(ncubp,nelementsPerBlock))
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    ALLOCATE(DbasTest(indofTest,elem_getMaxDerivative(p_elementDistribution%itestElement),&
             ncubp,nelementsPerBlock))
    ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(p_elementDistribution%itrialElement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))
    ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

    ! Check if one of the trial/test elements is nonparametric
    bnonparTrial = elem_isNonparametric(p_elementDistribution%itrialElement)
    bnonparTest  = elem_isNonparametric(p_elementDistribution%itestElement)
                    
    ! Let p_DcubPtsTrial / p_DcubPtsTest point either to DcubPtsReal or
    ! DcubPtsRef - depending on whether the space is parametric or not.
    IF (bnonparTrial) THEN
      p_DcubPtsTrial => DcubPtsReal
    ELSE
      p_DcubPtsTrial => DcubPtsRef
    END IF
    
    IF (bnonparTest) THEN
      p_DcubPtsTest => DcubPtsReal
    ELSE
      p_DcubPtsTest => DcubPtsRef
    END IF
    
    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    ALLOCATE(Kentry(indofTest,indofTrial,nelementsPerBlock))
    ALLOCATE(Dentry(indofTest,indofTrial))
    
    ! In case of nonconstant coefficients in that part of the matrix, we
    ! need an additional array to save all the coefficients:
    IF (.NOT. rform%BconstantCoeff(icurrentElementDistr)) THEN
      IF (rform%ballCoeffConstant) THEN
        PRINT *,'Error in bilf_buildMatrix9d_conf: Some oefficients are not constant &
                &although thy should be!'
        STOP
      END IF
      IF (.NOT. PRESENT(fcoeff_buildMatrixSc_sim)) THEN
        PRINT *,'Error in bilf_buildMatrix9d_conf: coefficient function not given!'
        STOP
      END IF
      ALLOCATE(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
    END IF
                    
    ! p_IdofsTest points either to the just created array or to the
    ! array with the DOF's of the trial functions - when trial and
    ! test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = p_elementDistribution%itrialElement .EQ. &
                             p_elementDistribution%itestElement

    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    IF (bIdenticalTrialAndTest) THEN
      p_IdofsTrial => IdofsTest
      p_DbasTrial  => DbasTest
      ! Build the actual combination of what the element should calculate.
      ! As we evaluate only once, what the element must calculate is an
      ! OR combination of the BDER from trial and test functions.
      BderTrial = BderTrialTempl .OR. BderTestTempl
      BderTest = BderTestTempl
    ELSE
      p_IdofsTrial => IdofsTrial
      p_DbasTrial  => DbasTrial
      
      ! Build the actual combination of what the element should calculate.
      ! Copy BDERxxxx to BDERxxxxAct
      BderTrial = BderTrialTempl
      BderTest = BderTestTempl
    END IF
    !CALL ZTIME(DT(3))
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)
                              
    ! Loop over the elements - blockwise.
    DO IELset = 1, p_rtriangulation%NEL, BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = MIN(p_rtriangulation%NEL,IELset-1+BILF_NELEMSIM)
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                  .TRUE.,IdofsTest)
                                   
      ! If the DOF's for the trial functions are different, calculate them, too.
      IF (.NOT.bIdenticalTrialAndTest) THEN
        CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                    .FALSE.,IdofsTrial)
      END IF
      !CALL ZTIME(DT(4))
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF's per element and
      ! indofTest test DOF's per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
      ! "active" (i.e. have common support) on our current element, each 
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTest,1..indofTrial) receives the position 
      !   in the global system matrix, where the corresponding value 
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately, 
      !  but we directly add them to the global matrix in this 
      !  approach.)
      !
      ! We build local matrices for all our elements 
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      DO IEL=1,IELmax-IELset+1

        ! Loop through the trial functions
        DO JDOFE=1,indofTrial
        
          ! Row JDOFE of the local matrix corresponds 
          ! to row=global DOF KDFG(JDOFE) in the global matrix.
          ! This is the "X" in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(p_IdofsTrial(JDOFE,IEL))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          DO IDOFE=1,indofTest
            
            ! Get the global DOF of the "O" which interacts with 
            ! our "X".
            
            IDFG=IdofsTest(IDOFE,IEL)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.
            
            DO JCOL=JCOL0,NA
              IF (p_KCOL(JCOL) .EQ. IDFG) EXIT
            END DO

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            
            Kentry(IDOFE,JDOFE,IEL)=JCOL
            
          END DO ! IDOFE
          
        END DO ! JDOFE
        
      END DO ! IEL
      !CALL ZTIME(DT(5))
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.
      !
      ! We have the coordinates of the cubature points saved in the
      ! coordinate array from above. Unfortunately for nonparametric
      ! elements, we need the real coordinate.
      ! Furthermore, we anyway need the coordinates of the element
      ! corners and the Jacobian determinants corresponding to
      ! all the points.
      !
      ! At first, get the coordinates of the corners of all the
      ! elements in the current set. 
      
!      DO IEL=1,IELmax-IELset+1
!        DCoords(:,:,IEL) = p_DcornerCoordinates(:, &
!                            p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!      END DO
      DO IEL=1,IELmax-IELset+1
        DO J = 1,NVE
          DO I = 1,NDIM2D
            DCoords(I,J,IEL) = p_DcornerCoordinates(I, &
                               p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
          END DO
        END DO
      END DO
      !CALL ZTIME(DT(6))
      
      ! Depending on the type of transformation, we must now choose
      ! the mapping between the reference and the real element.
      ! In case we use a nonparametric element or a nonconstant coefficient function,
      ! we need the coordinates of the points on the real element, too.
      IF (bnonparTrial .OR. bnonparTest .OR. (.NOT. rform%ballCoeffConstant)) THEN
      
        CALL trafo_calctrafo_sim (&
             rdiscretisation%RelementDistribution(icurrentElementDistr)%ctrafoType,&
             IELmax-IELset+1,ncubp,Dcoords,&
             DcubPtsRef,Djac(:,:,1:IELmax-IELset+1),Ddetj(:,1:IELmax-IELset+1),DcubPtsReal)
      
      ELSE
      
        CALL trafo_calctrafo_sim (p_elementDistribution%ctrafoType,&
             IELmax-IELset+1,ncubp,Dcoords,&
             DcubPtsRef,Djac(:,:,1:IELmax-IELset+1),Ddetj(:,1:IELmax-IELset+1))
             
      END IF
      
      !CALL ZTIME(DT(7))
      
      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
      IF (.NOT. rform%ballCoeffConstant) THEN
        CALL fcoeff_buildMatrixSc_sim (rdiscretisation,icurrentElementDistr, rform, &
                  IELset,IELmax-IELset+1,ncubp,p_IelementList(IELset:IELmax),Dcoords, &
                  DcubPtsRef,DcubPtsReal,p_IdofsTrial,IdofsTest,Djac,Ddetj,p_rcollection, &
                  Dcoefficients)
      END IF
      
      !CALL ZTIME(DT(8))                              
      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
      CALL elem_generic_sim (p_elementDistribution%itestElement, Dcoords, &
            Djac(:,:,1:IELmax-IELset+1), Ddetj(:,1:IELmax-IELset+1), &
            BderTest, DbasTest, ncubp, IELmax-IELset+1, p_DcubPtsTest)
            
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      IF (.NOT. bidenticalTrialAndTest) THEN
        CALL elem_generic_sim (p_elementDistribution%itrialElement, Dcoords, &
            Djac(:,:,1:IELmax-IELset+1), Ddetj(:,1:IELmax-IELset+1), &
            BderTrial, DbasTrial, ncubp, IELmax-IELset+1, p_DcubPtsTrial)
      END IF
      !CALL ZTIME(DT(9))
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! We have two different versions for the integration - one
      ! with constant coefficients and one with nonconstant coefficients.
      !
      ! Check the bilinear form which one to use:
      
      IF (rform%ballCoeffConstant) THEN
      
        ! Constant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable of the form.
        !
        ! Loop over the elements in the current set.

        DO IEL=1,IELmax-IELset+1
          
          ! Clear the local matrix
          Dentry = 0.0_DP
          
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.

            OM = Domega(ICUBP)*Ddetj(ICUBP,IEL)

            ! Loop over the additive factors in the bilinear form.
            DO IALBET = 1,rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be:
              !
              ! int_... ( phi_i )_IA  *  ( psi_j )_IB
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, 
              !      =2: 2nd derivative,...
              !    as defined in the module 'derivative'.
              
              IA = rform%Idescriptors(1,IALBET)
              IB = rform%Idescriptors(2,IALBET)
              
              ! Multiply OM with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              AUX = OM * rform%Dcoefficients(IALBET)
            
              ! Now loop through all possible combinations of DOF's
              ! in the current cubature point. The outer loop
              ! loops through the "X" in the above picture,
              ! the trial functions:

              DO JDOFE=1,indofTrial
              
                ! Get the value of the (trial) basis function 
                ! phi_i (our "X") in the cubature point:
                DB = p_DbasTrial(JDOFE,IA,ICUBP,IEL)
                
                !Perform an inner loop through the other DOF's
                ! (the "X"'s). 

                DO IDOFE=1,indofTest
                
                  ! Get the value of the basis function 
                  ! psi_j (our "O") in the cubature point. 
                  ! Them multiply:
                  !    DB * DBAS(..) * AUX
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up DB * DBAS(..) * AUX would give
                  ! the coefficient of the local matrix. We save this
                  ! contriobution in the local matrix.

                  !JCOLB = Kentry(IDOFE,JDOFE,IEL)
                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
                  Dentry(IDOFE,JDOFE) = Dentry(IDOFE,JDOFE)+DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
                
                END DO
              
              END DO ! JDOFE
              
            END DO ! IALBET

          END DO ! ICUBP 
          
          ! Incorporate the local matrices into the global one.
          ! Kentry gives the position of the additive contributions in Dentry.
          DO JDOFE=1,indofTrial
            DO IDOFE=1,indofTest
              p_DA(Kentry(IDOFE,JDOFE,IEL)) = p_DA(Kentry(IDOFE,JDOFE,IEL)) + Dentry(IDOFE,JDOFE)
            END DO
          END DO

        END DO ! IEL
        
      ELSE
      
        ! Nonconstant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable as computed above.
        !
        ! Loop over the elements in the current set.

        DO IEL=1,IELmax-IELset+1
          
          ! Clear the local matrix
          Dentry = 0.0_DP
          
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.

            OM = Domega(ICUBP)*Ddetj(ICUBP,IEL)

            ! Loop over the additive factors in the bilinear form.
            DO IALBET = 1,rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be:
              !
              ! int_... ( phi_i )_IA  *  ( psi_j )_IB
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, 
              !      =2: 2nd derivative,...
              !    as defined in the module 'derivative'.
              
              IA = rform%Idescriptors(1,IALBET)
              IB = rform%Idescriptors(2,IALBET)
              
              ! Multiply OM with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              ! Get the precalculated coefficient from the coefficient array.
              AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
            
              ! Now loop through all possible combinations of DOF's
              ! in the current cubature point. The outer loop
              ! loops through the "X" in the above picture,
              ! the trial functions:

              DO JDOFE=1,indofTrial
              
                ! Get the value of the (trial) basis function 
                ! phi_i (our "X") in the cubature point:
                DB = p_DbasTrial(JDOFE,IA,ICUBP,IEL)
                
                !Perform an inner loop through the other DOF's
                ! (the "O"'s). 

                DO IDOFE=1,indofTest
                
                  ! Get the value of the basis function 
                  ! psi_j (our "O") in the cubature point. This is
                  ! DBAS(KDFL(IDOFE),IB).
                  ! Them multiply:
                  !    DB * DBAS(..) * AUX
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up DB * DBAS(..) * AUX would give
                  ! the coefficient of the local matrix. We save this
                  ! contriobution in the local matrix of element IEL.

                  !JCOLB = Kentry(IDOFE,JDOFE,IEL)
                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
                  Dentry(IDOFE,JDOFE) = Dentry(IDOFE,JDOFE)+DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
                
                END DO
              
              END DO ! JDOFE
              
            END DO ! IALBET

          END DO ! ICUBP 
          
          ! Incorporate the local matrices into the global one.
          ! Kentry gives the position of the additive contributions in Dentry.
          DO JDOFE=1,indofTrial
            DO IDOFE=1,indofTest
              p_DA(Kentry(IDOFE,JDOFE,IEL)) = p_DA(Kentry(IDOFE,JDOFE,IEL)) + Dentry(IDOFE,JDOFE)
            END DO
          END DO

        END DO ! IEL

      END IF ! rform%ballCoeffConstant

      !CALL ZTIME(DT(10))
    END DO ! IELset
    
    IF (.NOT. rform%ballCoeffConstant) THEN
      DEALLOCATE(Dcoefficients)
    END IF
    DEALLOCATE(IdofsTest)
    DEALLOCATE(DbasTrial)
    DEALLOCATE(DbasTest)
    DEALLOCATE(Kentry)
    DEALLOCATE(Dentry)
    DEALLOCATE(Ddetj)
    DEALLOCATE(Djac)
    DEALLOCATE(DcubPtsReal)
    DEALLOCATE(DcubPtsRef)

  END DO ! icurrentElementDistr

  ! Clean up memory, finish

  DEALLOCATE(Dcoords)
  !CALL ZTIME(DT(11))
  
  !DO i=2,11
  !  PRINT *,'Time for assembly part ',i,': ',DT(i)-DT(i-1)
  !END DO
  
  !CFILE = 'MATRIX2.TXT'
  !CALL OWM17(p_DA,p_KCOL,p_KLD,&
  !           NEQ,NEQ,.TRUE.,0,'MAT1  ',CFILE,'(D20.10)')

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,&
                                       fcoeff_buildMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix.
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance. The discretisation is assumed to be conformal, i.e. the DOF's
  ! of all finite elements must 'match'. Trial and test functions may be
  ! different.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! Double-precision version.
!</description>

!<input>
  ! The bilinear form specifying the underlying PDE of the discretisation.
  TYPE(t_bilinearForm), INTENT(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  TYPE(t_collection), INTENT(IN), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  INCLUDE 'intf_coefficientMatrixSc.inc'
  OPTIONAL :: fcoeff_buildMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,i1,j,k,icurrentElementDistr,JDFG, ICUBP, IALBET, IA, IB
  LOGICAL :: bIdenticalTrialAndTest, bnonparTest, bnonparTrial
  INTEGER(I32) :: IEL, IELmax, IELset, IDOFE, JDOFE
  INTEGER(PREC_DOFIDX) :: JCOL0,JCOL
  REAL(DP) :: OM,AUX, DB
  
  ! Array to tell the element which derivatives to calculate
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrialTempl, BderTestTempl, BderTrial, BderTest
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp
  
  ! Pointer to KLD, KCOL, DA
  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest,DbasTrial
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasTrial
  
  ! Number of entries in the matrix - for quicker access
  INTEGER(PREC_DOFIDX) :: NA,NVE
  INTEGER(I32) :: NEQ
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dentry
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef

  ! An array receiving the coordinates of cubature points on
  ! the real element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal

  ! Pointer to the point coordinates to pass to the element function.
  ! Point either to p_DcubPtsRef or to p_DcubPtsReal, depending on whether
  ! the trial/test element is parametric or not.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTrial
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
  
  ! Array with coordinates of the corners that form the real element.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Dcoords
  
  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
  
  ! Pointer to KVERT of the triangulation
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  
  ! Pointer to DCORVG of the triangulation
  REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
  
  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock
  
  ! Some variables to support nonconstant coefficients in the matrix.
  
  ! Pointer to the collection structure or to NULL()
  TYPE(t_collection), POINTER :: p_rcollection
  
  ! Pointer to the coefficients that are computed by the callback routine.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  TYPE(t_domainIntSubset) :: rintSubset
  
  ! The discretisation - for easier access
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
  
  !REAL(DP), DIMENSION(11) :: DT
  
  !CHARACTER(LEN=20) :: CFILE
  
  IF (.NOT. ASSOCIATED(rmatrixScalar%p_rspatialDiscretisation)) THEN
    PRINT *,'bilf_buildMatrix9d_conf2: No discretisation associated!'
    STOP
  END IF

  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDERxxxx
  ! according to these.

  !CALL ZTIME(DT(1))

  BderTrialTempl = .FALSE.
  BderTestTempl = .FALSE.
  
  ! Loop through the additive terms
  DO i=1,rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed! Build template's for BDER.
    ! We don't compute the actual BDER here, as there might be some special
    ! processing if trial/test functions are identical!
    !
    ! At first build the descriptors for the trial functions
    I1=rform%Idescriptors(1,I)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'bilf_buildMatrix9d_conf: Invalid descriptor'
      STOP
    ENDIF
    
    BderTrialTempl(I1)=.TRUE.

    ! Then those of the test functions
    I1=rform%Idescriptors(2,I)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'bilf_buildMatrix9d_conf: Invalid descriptor'
      STOP
    ENDIF
    
    BderTestTempl(I1)=.TRUE.
  END DO
  
  ! Get information about the matrix:
  NA = rmatrixScalar%NA
  NEQ = rmatrixScalar%NEQ
  
  ! We need KCOL/KLD of our matric
  CALL storage_getbase_int (rmatrixScalar%h_KCOL,p_KCOL)
  CALL storage_getbase_int (rmatrixScalar%h_KLD,p_KLD)
  
  ! Check if the matrix entries exist. If not, allocate the matrix.
  IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN

    ! Clear the entries in the matrix - we need to start with zero
    ! when assembling a new matrix!
    CALL storage_new1D ('bilf_buildMatrix9d_conf', 'DA', &
                        NA, ST_DOUBLE, rmatrixScalar%h_DA, &
                        ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double (rmatrixScalar%h_DA,p_DA)

  ELSE
  
    CALL storage_getbase_double (rmatrixScalar%h_DA,p_DA)

    ! If desired, clear the matrix before assembling.
    IF (bclear) THEN
      CALL lalg_clearVectorDble (p_DA)
    END IF
    
  END IF
  
  ! Get the discretisation
  p_rdiscretisation => rmatrixScalar%p_rspatialDiscretisation
  
  IF (.NOT. ASSOCIATED(p_rdiscretisation)) THEN
    PRINT *,'bilf_buildMatrix9d_conf2 error: No discretisation attached to the matrix!'
    STOP
  END IF
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => p_rdiscretisation%p_rtriangulation
  
  ! Let p_rcollection point to rcollection - or NULL if it's not
  ! given.
  IF (PRESENT(rcollection)) THEN
    p_rcollection => rcollection
  ELSE
    p_rcollection => NULL()
  END IF

  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocaing some arrays.
  nelementsPerBlock = MIN(BILF_NELEMSIM,p_rtriangulation%NEL)
  
  ! Get a pointer to the KVERT and DCORVG array
  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
                             p_IverticesAtElement)
  CALL storage_getbase_double2D(p_rtriangulation%h_DcornerCoordinates, &
                             p_DcornerCoordinates)

  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  !CALL ZTIME(DT(2))

  DO icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
  
    ! Activate the current element distribution
    p_elementDistribution => p_rdiscretisation%RelementDistribution(icurrentElementDistr)
  
    ! Get the number of local DOF's for trial and test functions
    indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
    indofTest = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
      PRINT *,'bilf_buildMatrix9d_conf2: element spaces incompatible!'
      STOP
    END IF
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    CALL cub_getCubPoints(p_elementDistribution%ccubType, ncubp, Dxi, Domega)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    j = elem_igetCoordSystem(p_elementDistribution%itrialElement)
    
    ! Allocate memory and get local references to it.
    CALL domint_initIntegration (rintSubset,nelementsPerBlock,ncubp,j,NDIM2D)
    p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
    p_DcubPtsReal => rintSubset%p_DcubPtsReal
    p_Djac =>        rintSubset%p_Djac
    p_Ddetj =>       rintSubset%p_Ddetj
    p_Dcoords =>     rintSubset%p_DCoords
    
    ! Put the cubature point coordinates in the right format to the
    ! cubature-point array.
    ! Initialise all entries in p_DcubPtsRef with the same coordinates -
    ! as the cubature point coordinates are identical on all elements
    DO j=1,SIZE(p_DcubPtsRef,3)
      DO i=1,ncubp
        DO k=1,SIZE(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i,j) = Dxi(i,k)
        END DO
      END DO
    END DO
    
    ! Quickly check if one of the specified derivatives is out of the allowed range:
    DO IALBET = 1,rform%itermcount
      IA = rform%Idescriptors(1,IALBET)
      IB = rform%Idescriptors(2,IALBET)      
      IF ((IA.LT.0) .OR. &
          (IA .GT. elem_getMaxDerivative(p_elementDistribution%itrialElement))) THEN
        PRINT *,'bilf_buildMatrix9d_conf2: Specified trial-derivative',IA,&
                ' not available'
        STOP
      END IF

      IF ((IB.LT.0) .OR. &
          (IB .GT. elem_getMaxDerivative(p_elementDistribution%itestElement))) THEN
        PRINT *,'bilf_buildMatrix9d_conf2: Specified test-derivative',IA,&
                ' not available'
        STOP
      END IF
    END DO
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    ALLOCATE(DbasTest(indofTest,elem_getMaxDerivative(p_elementDistribution%itestElement),&
             ncubp,nelementsPerBlock))
    ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(p_elementDistribution%itrialElement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))
    ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

    ! Check if one of the trial/test elements is nonparametric
    bnonparTrial = elem_isNonparametric(p_elementDistribution%itrialElement)
    bnonparTest  = elem_isNonparametric(p_elementDistribution%itestElement)
                    
    ! Let p_DcubPtsTrial / p_DcubPtsTest point either to p_DcubPtsReal or
    ! p_DcubPtsRef - depending on whether the space is parametric or not.
    IF (bnonparTrial) THEN
      p_DcubPtsTrial => p_DcubPtsReal
    ELSE
      p_DcubPtsTrial => p_DcubPtsRef
    END IF
    
    IF (bnonparTest) THEN
      p_DcubPtsTest => p_DcubPtsReal
    ELSE
      p_DcubPtsTest => p_DcubPtsRef
    END IF
    
    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    ALLOCATE(Kentry(indofTrial,indofTest,nelementsPerBlock))
    ALLOCATE(Dentry(indofTrial,indofTest))
    
    ! In case of nonconstant coefficients in that part of the matrix, we
    ! need an additional array to save all the coefficients:
    IF (.NOT. rform%BconstantCoeff(icurrentElementDistr)) THEN
      IF (rform%ballCoeffConstant) THEN
        PRINT *,'Error in bilf_buildMatrix9d_conf: Some oefficients are not constant &
                &although thy should be!'
        STOP
      END IF
      IF (.NOT. PRESENT(fcoeff_buildMatrixSc_sim)) THEN
        PRINT *,'Error in bilf_buildMatrix9d_conf: coefficient function not given!'
        STOP
      END IF
      ALLOCATE(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
    END IF
                    
    ! p_IdofsTest points either to the just created array or to the
    ! array with the DOF's of the trial functions - when trial and
    ! test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = p_elementDistribution%itrialElement .EQ. &
                             p_elementDistribution%itestElement

    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    IF (bIdenticalTrialAndTest) THEN
      p_IdofsTrial => IdofsTest
      p_DbasTrial  => DbasTest
      ! Build the actual combination of what the element should calculate.
      ! As we evaluate only once, what the element must calculate is an
      ! OR combination of the BDER from trial and test functions.
      BderTrial = BderTrialTempl .OR. BderTestTempl
      BderTest = BderTrial
    ELSE
      p_IdofsTrial => IdofsTrial
      p_DbasTrial  => DbasTrial
      
      ! Build the actual combination of what the element should calculate.
      ! Copy BDERxxxx to BDERxxxxAct
      BderTrial = BderTrialTempl
      BderTest = BderTestTempl
    END IF
    !CALL ZTIME(DT(3))
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)
                              
    ! Loop over the elements - blockwise.
    DO IELset = 1, p_rtriangulation%NEL, BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = MIN(p_rtriangulation%NEL,IELset-1+BILF_NELEMSIM)
    
      ! --------------------- DOF SEARCH PHASE ------------------------
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  .TRUE.,IdofsTest)
                                   
      ! If the DOF's for the test functions are different, calculate them, too.
      IF (.NOT.bIdenticalTrialAndTest) THEN
        CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                    .FALSE.,IdofsTrial)
      END IF
      !CALL ZTIME(DT(4))
      
      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
      
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF's per element and
      ! indofTest test DOF's per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
      ! "active" (i.e. have common support) on our current element, each 
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position 
      !   in the global system matrix, where the corresponding value 
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately, 
      !  but we directly add them to the global matrix in this 
      !  approach.)
      !
      ! We build local matrices for all our elements 
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      DO IEL=1,IELmax-IELset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        DO IDOFE=1,indofTest
        
          ! Row IDOFE of the local matrix corresponds 
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(IdofsTest(IDOFE,IEL))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          DO JDOFE=1,indofTrial
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            JDFG=p_IdofsTrial(JDOFE,IEL)
            
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
            
            Kentry(JDOFE,IDOFE,IEL)=JCOL
            
          END DO ! IDOFE
          
        END DO ! JDOFE
        
      END DO ! IEL
      !CALL ZTIME(DT(5))
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.
      !
      ! We have the coordinates of the cubature points saved in the
      ! coordinate array from above. Unfortunately for nonparametric
      ! elements, we need the real coordinate.
      ! Furthermore, we anyway need the coordinates of the element
      ! corners and the Jacobian determinants corresponding to
      ! all the points.
      !
      ! At first, get the coordinates of the corners of all the
      ! elements in the current set. 
      
!      DO IEL=1,IELmax-IELset+1
!        p_Dcoords(:,:,IEL) = p_DcornerCoordinates(:, &
!                            p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!      END DO
      DO IEL=1,IELmax-IELset+1
        DO J = 1,NVE
          DO I = 1,NDIM2D
            p_Dcoords(I,J,IEL) = p_DcornerCoordinates(I, &
                               p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
          END DO
        END DO
      END DO
      !CALL ZTIME(DT(6))
      
      ! Depending on the type of transformation, we must now choose
      ! the mapping between the reference and the real element.
      ! In case we use a nonparametric element or a nonconstant coefficient function,
      ! we need the coordinates of the points on the real element, too.
      IF (bnonparTrial .OR. bnonparTest .OR. (.NOT. rform%ballCoeffConstant)) THEN
      
        CALL trafo_calctrafo_sim (&
             p_rdiscretisation%RelementDistribution(icurrentElementDistr)%ctrafoType,&
             IELmax-IELset+1,ncubp,p_Dcoords,&
             p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1),p_DcubPtsReal)
      
      ELSE
      
        CALL trafo_calctrafo_sim (p_elementDistribution%ctrafoType,&
             IELmax-IELset+1,ncubp,p_Dcoords,&
             p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1))
             
      END IF
      
      !CALL ZTIME(DT(7))
      
      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
      IF (.NOT. rform%ballCoeffConstant) THEN
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        CALL fcoeff_buildMatrixSc_sim (p_rdiscretisation,rform, &
                  IELmax-IELset+1_I32,ncubp,&
                  p_DcubPtsReal,p_IdofsTrial,IdofsTest,rintSubset,p_rcollection, &
                  Dcoefficients)
      END IF
      
      !CALL ZTIME(DT(8))                              
      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
      CALL elem_generic_sim (p_elementDistribution%itestElement, p_Dcoords, &
            p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
            BderTest, DbasTest, ncubp, IELmax-IELset+1, p_DcubPtsTest)
            
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      IF (.NOT. bidenticalTrialAndTest) THEN
        CALL elem_generic_sim (p_elementDistribution%itrialElement, p_Dcoords, &
            p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
            BderTrial, DbasTrial, ncubp, IELmax-IELset+1, p_DcubPtsTrial)
      END IF
      !CALL ZTIME(DT(9))
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! We have two different versions for the integration - one
      ! with constant coefficients and one with nonconstant coefficients.
      !
      ! Check the bilinear form which one to use:
      
      IF (rform%ballCoeffConstant) THEN
      
        ! Constant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable of the form.
        !
        ! Loop over the elements in the current set.

        DO IEL=1,IELmax-IELset+1
          
          ! Clear the local matrix
          Dentry = 0.0_DP
          
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.

            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

            ! Loop over the additive factors in the bilinear form.
            DO IALBET = 1,rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be added to the matrix entry:
              !
              ! a_ij  =  int_... ( psi_j )_IB  *  ( phi_i )_IA
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, ...
              !    as defined in the module 'derivative'.
              
              IA = rform%Idescriptors(1,IALBET)
              IB = rform%Idescriptors(2,IALBET)
              
              ! Multiply OM with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              AUX = OM * rform%Dcoefficients(IALBET)
            
              ! Now loop through all possible combinations of DOF's
              ! in the current cubature point. The outer loop
              ! loops through the "O"'s in the above picture,
              ! the test functions:

              DO IDOFE=1,indofTest
              
                ! Get the value of the (test) basis function 
                ! phi_i (our "O") in the cubature point:
                DB = DbasTest(IDOFE,IB,ICUBP,IEL)
                
                ! Perform an inner loop through the other DOF's
                ! (the "X"). 

                DO JDOFE=1,indofTrial
                
                  ! Get the value of the basis function 
                  ! psi_j (our "X") in the cubature point. 
                  ! Them multiply:
                  !    DB * DBAS(..) * AUX
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up DB * DBAS(..) * AUX would give
                  ! the coefficient of the local matrix. We save this
                  ! contribution in the local matrix.

                  !JCOLB = Kentry(JDOFE,IDOFE,IEL)
                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                  Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE) + &
                                        DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                
                END DO ! JDOFE
              
              END DO ! IDOFE
              
            END DO ! IALBET

          END DO ! ICUBP 
          
          ! Incorporate the local matrices into the global one.
          ! Kentry gives the position of the additive contributions in Dentry.
          DO IDOFE=1,indofTest
            DO JDOFE=1,indofTrial
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
                                              Dentry(JDOFE,IDOFE)
            END DO
          END DO

        END DO ! IEL
        
      ELSE
      
        ! Nonconstant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable as computed above.
        !
        ! Loop over the elements in the current set.

        DO IEL=1,IELmax-IELset+1
          
          ! Clear the local matrix
          Dentry = 0.0_DP
          
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.

            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

            ! Loop over the additive factors in the bilinear form.
            DO IALBET = 1,rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be added to the matrix entry:
              !
              ! a_ij  =  int_... ( psi_j )_IA  *  ( phi_i )_IB
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, ...
              !    as defined in the module 'derivative'.
              
              IA = rform%Idescriptors(1,IALBET)
              IB = rform%Idescriptors(2,IALBET)
              
              ! Multiply OM with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              ! Get the precalculated coefficient from the coefficient array.
              AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
            
              ! Now loop through all possible combinations of DOF's
              ! in the current cubature point. The outer loop
              ! loops through the "O" in the above picture,
              ! the test functions:

              DO IDOFE=1,indofTest
                
                ! Get the value of the (test) basis function 
                ! phi_i (our "O") in the cubature point:
                DB = DbasTest(IDOFE,IB,ICUBP,IEL)
                
                ! Perform an inner loop through the other DOF's
                ! (the "X"). 

                DO JDOFE=1,indofTrial
              
                  ! Get the value of the basis function 
                  ! psi_j (our "X") in the cubature point. 
                  ! Them multiply:
                  !    DB * DBAS(..) * AUX
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up DB * DBAS(..) * AUX would give
                  ! the coefficient of the local matrix. We save this
                  ! contribution in the local matrix of element IEL.

                  !JCOLB = Kentry(JDOFE,IDOFE,IEL)
                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                  Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                
                END DO
              
              END DO ! JDOFE
              
            END DO ! IALBET

          END DO ! ICUBP 
          
          ! Incorporate the local matrices into the global one.
          ! Kentry gives the position of the additive contributions in Dentry.
          DO IDOFE=1,indofTest
            DO JDOFE=1,indofTrial
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + Dentry(JDOFE,IDOFE)
            END DO
          END DO

        END DO ! IEL

      END IF ! rform%ballCoeffConstant

      !CALL ZTIME(DT(10))
    END DO ! IELset
    
    ! Release memory
    CALL domint_doneIntegration(rintSubset)

    IF (.NOT. rform%ballCoeffConstant) THEN
      DEALLOCATE(Dcoefficients)
    END IF
    DEALLOCATE(IdofsTrial)
    DEALLOCATE(IdofsTest)
    DEALLOCATE(DbasTrial)
    DEALLOCATE(DbasTest)
    DEALLOCATE(Kentry)
    DEALLOCATE(Dentry)

  END DO ! icurrentElementDistr

  ! Finish
  !CALL ZTIME(DT(11))
  
  !DO i=2,11
  !  PRINT *,'Time for assembly part ',i,': ',DT(i)-DT(i-1)
  !END DO
  
  !CFILE = 'MATRIX2.TXT'
  !CALL OWM17(p_DA,p_KCOL,p_KLD,&
  !           NEQ,NEQ,.TRUE.,0,'MAT1  ',CFILE,'(D20.10)')

  END SUBROUTINE
  
END MODULE
