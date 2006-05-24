!##############################################################################
!# ****************************************************************************
!# <name> linearsystemscalar </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic definitions and routines to maintain scalar
!# linear systems. A linear system is realised by matrices and vectors, where
!# the matrix does not decompose into smaller 'block-matrices'. Therefore,
!# we denote the matrices/vectors here as 'scalar', thus containing no 
!# ''blocks''.
!#
!# For storing matrices, we use the CSR format, realised as ''Format-7''
!# and ''Format-9'' matrices. Format-9 matrices describe the usual CSR-format
!# with the extension, that an index is stored for each line to the first
!# element on/above the diagonal. Format-7 is a variant of the CSR-format
!# used in old FEAT programs, where the diagonal element is stored first in
!# each line.
!#
!# The following routines can be found in this module:
!#
!# 1.) lsyssc_scalarProduct
!#     -> Calculate the scalar product of two vectors
!#
!# 2.) lsyssc_scalarMatVec
!#     -> Multiply a scalar matrix with a scalar vector
!#
!# 3.) lsyssc_releaseMatrix
!#     -> Release a scalar matrix from memory.
!#
!# 4.) lsyssc_releaseVector
!#     -> Release a scalar vector from memory.
!#
!# 5.) lsyssc_duplicateMatrix
!#     -> Create a duplicate of a given matrix or matrix-structure
!#
!# 6.) lsyssc_duplicateVector
!#     -> Create a duplicate of a given vector
!#
!# 7.) lsyssc_sortVector
!#     -> Resort the entries of a vector or unsort them
!#
!# 8.) lsyssc_sortMatrix
!#     -> Resort the entries of a matrix or unsort them
!# </purpose>
!##############################################################################

MODULE linearsystemscalar

  USE fsystem
  USE storage
  USE spatialdiscretisation
  USE discretebc
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Global format flags for matrices">

  ! Unidentified matrix format
  INTEGER, PARAMETER :: LSYSSC_MATRIXUNDEFINED = 0
  
  ! Identifier for matrix format 1 - full matrix
  INTEGER, PARAMETER :: LSYSSC_MATRIX1 = 1
  
  ! Identifier for matrix format 9 - CSR
  INTEGER, PARAMETER :: LSYSSC_MATRIX9 = 9

  ! Identifier for matrix format 7 - CSR with diagonal element in front
  INTEGER, PARAMETER :: LSYSSC_MATRIX7 = 7

  !</constantblock>

  !<constantblock description="Flags for the matrix specification bitfield">

  ! Standard matrix
  INTEGER, PARAMETER :: LSYSSC_MSPEC_STANDARD =        0
  
  ! Matrix structure is a copy of another matrix, shared via the same
  ! handles. 
  INTEGER, PARAMETER :: LSYSSC_MSPEC_STRUCTUREISCOPY = 2**0

  ! Matrix content is a copy of another matrix, shared via the same
  ! handles. 
  INTEGER, PARAMETER :: LSYSSC_MSPEC_CONTENTISCOPY   = 2**1
  
  ! Complete matrix is duplicate of another matrix and shares structure
  ! and entries via the same pointers
  INTEGER, PARAMETER :: LSYSSC_MSPEC_ISCOPY = LSYSSC_MSPEC_STRUCTUREISCOPY +&
                                              LSYSSC_MSPEC_CONTENTISCOPY

  ! Matrix is saved transposed
  INTEGER, PARAMETER :: LSYSSC_MSPEC_TRANSPOSED =      2**2

  ! Matrix not present in memory
  INTEGER, PARAMETER :: LSYSSC_MSPEC_NOTINMEMORY =     2**3
  
!</constantblock>

!<constantblock description="KIND values for matrix/vector data">
  
  ! kind value for indices in matrices
  INTEGER, PARAMETER :: PREC_MATIDX = I32

  ! kind value for indices in vectors
  INTEGER, PARAMETER :: PREC_VECIDX = I32

  ! kind value for precision that should be used in matrices
  INTEGER, PARAMETER :: PREC_MATRIX = DP

!</constantblock>

!<constantblock description="Constants for duplicating a matrix">
  
  ! Duplicate by ownership. What belongs to the matrix is duplicated,
  ! what belongs to another matrix is left as-is.
  INTEGER, PARAMETER :: LSYSSC_DUP_ASIS      = 0

  ! Duplicate the content, share the structure
  INTEGER, PARAMETER :: LSYSSC_DUP_CONTENT   = 2**0

  ! Duplicate the structure, share the content
  INTEGER, PARAMETER :: LSYSSC_DUP_STRUCTURE = 2**1

  ! Duplicate both, structure and content
  INTEGER, PARAMETER :: LSYSSC_DUP_ALL       = LSYSSC_DUP_STRUCTURE + &
                                               LSYSSC_DUP_CONTENT

!</constantblock>
  
!</constants>


!<types>
  
!<typeblock>
  
  ! A scalar vector that can be used by scalar linear algebra routines.
  
  TYPE t_vectorScalar
  
    ! Length of the vector; not necessarily = SIZE(p_Ddata), as this
    ! scalar vector might be a subvector of a larger vector allocated
    ! on the heap!
    INTEGER(PREC_VECIDX) :: NEQ       = 0
    
    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE.
    INTEGER         :: cdataType = ST_DOUBLE
    
    ! Handle identifying the vector entries.
    ! = ST_NOHANDLE if not allocated on the global heap. 
    INTEGER         :: h_Ddata   = ST_NOHANDLE
    
    ! Flag whether or not the vector is resorted.
    ! <=0: Vector is unsorted.
    !  >0: Vector is resorted according to a sorting strategy.
    ! The value identifies the sorting strategy used; this is
    ! usually one of the SSTRAT_xxxx constants from the module
    ! 'sortstrategy'.
    ! In most applications, the absolute value 
    !               |isortStrategy| = SSTRAT_xxxx
    ! indicates the sorting trategy to use, while the sign
    ! indicates whether the sorting strategy is active on the
    ! vector (+) or not (-).
    INTEGER         :: isortStrategy   = 0
    
    ! Handle to renumbering strategy for resorting the vector.
    ! The renumbering strategy is a vector
    !   array [1..2*NEQ] of integer
    ! The first NEQ entries (1..NEQ) represent the permutation how to
    ! sort an unsorted vector. The second NEQ entries (NEQ+1..2*NEQ)
    ! represent the inverse permutation.
    ! Whether or not the vector is actually sorted depends on the
    ! flag isortStrategy!
    INTEGER         :: h_IsortPermutation
    
    ! Start position of the vector data in the array identified by
    ! h_Ddata. Normally = 1. Can be set to > 1 if the vector is a subvector
    ! in a larger memory block allocated on the heap.
    INTEGER(PREC_VECIDX) :: iidxFirstEntry = 1
    
    ! Is set to true, if the handle h_Ddata belongs to another vector,
    ! i.e. when this vector shares data with another vector.
    ! This is usually used for block vectors containing a couple of
    ! scalar subvectors; in this case, the block vector is the actual
    ! 'owner' of the handle!
    LOGICAL              :: bisCopy        = .FALSE.
    
    ! A pointer to the spatial discretisation or NULL(), if the vector is
    ! just an array, not belonging to any discretisation.
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscretisation => NULL()
    
    ! A pointer to discretised boundary conditions for real boundary components.
    ! For every discrete BC that is valid for this scalar vector, there
    ! is an entry in the following list.
    TYPE(t_discreteBC), POINTER  :: p_rdiscreteBC => NULL()
    
    ! A pointer to discretised boundary conditions for fictitious boundary
    ! components
    TYPE(t_discreteBC), POINTER  :: p_rdiscreteBCfict => NULL()
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A scalar matrix that can be used by scalar linear algebra routines.
  
  TYPE t_matrixScalar
    
    ! Format-tag. Identifies the format of the matrix. 
    ! Can tage one of the LSYSSC_MATRIXx format flags.
    INTEGER :: imatrixFormat = LSYSSC_MATRIXUNDEFINED
    
    ! Matrix specification tag. This is a bitfield coming from an OR 
    ! combination of different LSYSSC_MSPEC_xxxx constants and specifies
    ! various details of the matrix. If it is =LSYSSC_MSPEC_STANDARD,
    ! the matrix is a usual matrix that needs no special handling.
    INTEGER(I32) :: imatrixSpec = LSYSSC_MSPEC_STANDARD
    
    ! Number of elements in the matrix
    INTEGER(PREC_MATIDX) :: NA  = 0
    
    ! Number of equations in the matrix
    INTEGER(PREC_VECIDX) :: NEQ = 0
    
    ! Multiplier for matrix entries. All entries in the matrix are
    ! scaled by this multiplier when doing Matrix-vector multiplication.
    REAL(DP)         :: dScaleFactor = 1.0_DP
    
    ! Flag whether or not the matrix is resorted.
    ! <=0: Matrix is unsorted.
    !  >0: Matrix is resorted according to a sorting strategy.
    ! The value identifies the sorting strategy used; this is
    ! usually one of the SSTRAT_xxxx constants from the module
    ! 'sortstrategy'.
    INTEGER         :: isortStrategy   = 0
    
    ! Handle to renumbering strategy for resorting the matrix.
    ! The renumbering strategy is a vector
    !   array [1..2*NEQ] of integer
    ! The first NEQ entries (1..NEQ) represent the permutation how to
    ! sort an unsorted matrix. The second NEQ entries (NEQ+1..2*NEQ)
    ! represent the inverse permutation.
    ! Whether or not the matrix is actually sorted depends on the
    ! flag isortStrategy!
    INTEGER         :: h_IsortPermutation
    
    
    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE.
    INTEGER         :: cdataType = ST_DOUBLE
    
    ! Format-7 and Format-9: Handle identifying the elements in the matrix
    !REAL(PREC_MATRIX), DIMENSION(:), POINTER       :: DA         => NULL()
    INTEGER                    :: h_DA = ST_NOHANDLE
    
    ! Format-7 and Format-9: Handle identifying the column structure
    !INTEGER(PREC_MATIDX), DIMENSION(:), POINTER    :: KCOL       => NULL()
    INTEGER                    :: h_Kcol = ST_NOHANDLE
    
    ! Format-7 and Format-9: Handle identifying the row structure
    !INTEGER, DIMENSION(:), POINTER            :: KLD        => NULL()
    INTEGER                    :: h_Kld = ST_NOHANDLE
    
    ! Format-9: Similar to row structure. Handle top array of length NEQ.
    ! For each row, pointer to the first element on the upper
    ! triangular part of the matrix. For matrices with elements on the
    ! diagonal, this points to the diagonal element in each row.
    ! For matrices with no elements in the upper right part of the
    ! matrix, this points to the first element on the next line.
    !INTEGER, DIMENSION(:), POINTER            :: Kdiagonal  => NULL()
    INTEGER                    :: h_Kdiagonal = ST_NOHANDLE
    
    ! A pointer to the spatial discretisation
    TYPE(t_spatialDiscretisation), POINTER     :: p_rspatialDiscretisation => NULL()
    
    ! A pointer to discretised boundary conditions for real boundary components.
    ! For every discrete BC that is valid for this scalar vector, there
    ! is an entry in the following list.
    TYPE(t_discreteBC), POINTER  :: p_rdiscreteBC     => NULL()
    
    ! A pointer to discretised boundary conditions for fictitious boundary
    ! components
    TYPE(t_discreteBC), POINTER  :: p_rdiscreteBCfict => NULL()
    
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

  !****************************************************************************
!<subroutine>
  
  REAL(DP) FUNCTION lsyssc_scalarProduct (rx, ry)
  
!<description>
  ! Calculates a scalar product of two vectors.
!</description>
  
!<input>
  ! First vector
  TYPE(t_vectorScalar), INTENT(IN)                  :: rx

  ! Second vector
  TYPE(t_vectorScalar), INTENT(IN)                  :: ry

!</input>

!<result>
  ! The scalar product <rx,ry> of the two vectors.
!</result>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata1dp
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata2dp
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata1dp
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata2dp
  REAL(DP) :: res
  INTEGER(PREC_VECIDX) ioffsetx,ioffsety,i
  
  ! Is there data at all?
  res = 0.0_DP
  
  IF ( (rx%NEQ .EQ. 0) .OR. (ry%NEQ .EQ. 0) .OR. (rx%NEQ .NE. rx%NEQ)) THEN
    PRINT *,'Error in lsyssc_scalarProduct: Vector dimensions wrong!'
    STOP
  END IF
  
  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsyssc_scalarProduct: Data types different!'
    STOP
  END IF
  
  ! Get the offset positions from the vector structures
  ioffsetx = rx%iidxFirstEntry
  ioffsety = ry%iidxFirstEntry
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    
    CALL storage_getbase_double (rx%h_Ddata,p_Ddata1dp)
    CALL storage_getbase_double (ry%h_Ddata,p_Ddata2dp)
    
    ! Change the pointer to point to the subvector
    p_Ddata1dp => p_Ddata1dp(ioffsetx:ioffsetx+rx%NEQ-1)
    p_Ddata2dp => p_Ddata2dp(ioffsety:ioffsety+ry%NEQ-1)
    
    ! Perform the scalar product
    res = p_Ddata1dp(1)*p_Ddata2dp(1)
    DO i=2,rx%NEQ
      res = res + p_Ddata1dp(i)*p_Ddata2dp(i)
    END DO
    
  CASE (ST_SINGLE)
    
    CALL storage_getbase_single (rx%h_Ddata,p_Fdata1dp)
    CALL storage_getbase_single (ry%h_Ddata,p_Fdata2dp)
    
    ! Change the pointer to point to the subvector
    p_Fdata1dp => p_Fdata1dp(ioffsetx:ioffsetx+rx%NEQ-1)
    p_Fdata2dp => p_Fdata2dp(ioffsety:ioffsety+ry%NEQ-1)

    ! Perform the scalar product
    res = p_Fdata1dp(1)*p_Fdata2dp(1)
    DO i=2,rx%NEQ
      res = res + p_Fdata1dp(i)*p_Fdata2dp(i)
    END DO
    
  CASE DEFAULT
    PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
    STOP
  END SELECT
  
  ! Return the scalar product, finish
  lsyssc_scalarProduct = res

  END FUNCTION

  !****************************************************************************

!<subroutine>
  
  SUBROUTINE lsyssc_scalarMatVec (rmatrix, rx, ry, cx, cy)
  
!<description>
  ! Performs a matrix vector multiplicationwith a given scalar matrix:
  !    $$ Dy   =   cx * rMatrix * rx   +   cy * ry $$
!</description>
  
!<input>
  
  ! Scalar matrix
  TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix

  ! Vector to multiply with the matrix.
  TYPE(t_vectorScalar), INTENT(IN)                  :: rx
  
  ! Multiplicative factor for rx
  REAL(DP), INTENT(IN)                              :: cx

  ! Multiplicative factor for ry
  REAL(DP), INTENT(IN)                              :: cy
  
!</input>

!<inputoutput>
  ! Additive vector. Receives the result of the matrix-vector multiplication
  TYPE(t_vectorScalar), INTENT(INOUT)               :: ry
!</inputoutput>

!</subroutine>
  
  ! rx and ry must have at least the same data type!
  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'MV with different data types for rx and ry not supported!'
    STOP
  END IF
  
  ! rx, ry and the matrix must have proper dimensions
  IF ((rx%NEQ .NE. ry%NEQ) .OR. (rx%NEQ .NE. rmatrix%NEQ)) THEN
    PRINT *,'Error in MV: Vectors and matrix have different size!'
    STOP
  END IF
  
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
    PRINT *,'Multiplication with transposed matrix not yet supported!'
    STOP
  END IF

  ! Select the right MV multiplication routine from the matrix format
  SELECT CASE (rmatrix%imatrixFormat)
  
  CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
  
    ! Take care of the precision of the entries
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      ! Format 7 and Format 9 multiplication
      SELECT CASE (rx%cdataType)
      
      CASE (ST_DOUBLE)
        ! double precision matrix, double precision vectors
        CALL lsyssc_LAX79doubledouble (rmatrix,rx,ry,cx,cy)
      
      CASE DEFAULT
        PRINT *,'Only double precision vectors supported for now in MV!'
        STOP
        
      END SELECT
      
    CASE DEFAULT
      PRINT *,'Only double precision matrices supported for now in MV!'
      STOP
    END SELECT
    
  CASE DEFAULT
    PRINT *,'Unknown matrix format in MV-multiplication!'
    STOP
  END SELECT
  
  CONTAINS
  
    ! Now the real MV multiplication routines follow.
    ! We create them in the scoping unit of the procedure to prevent
    ! direct calls from outside.
    
    !**************************************************************
    ! Format 7 and Format 9 multiplication
    ! double precision matrix,
    ! double precision vectors
    
    SUBROUTINE lsyssc_LAX79doubledouble (rmatrix,rx,ry,cx,cy)

    ! Save arguments as above - given as parameters as some compilers
    ! might have problems with scoping units...
    TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
    TYPE(t_vectorScalar), INTENT(IN)                  :: rx
    REAL(DP), INTENT(IN)                              :: cx
    REAL(DP), INTENT(IN)                              :: cy
    TYPE(t_vectorScalar), INTENT(INOUT)               :: ry

    REAL(DP), DIMENSION(:), POINTER :: p_DA, p_Dx, p_Dy
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_MATIDX) :: irow,icol
    INTEGER(PREC_VECIDX) :: ioffsetx,ioffsety
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ

      ! Get the matrix and the two vectors
      CALL storage_getbase_double (rmatrix%h_DA,p_DA)
      CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
      CALL storage_getbase_double (rx%h_Ddata,p_Dx)
      CALL storage_getbase_double (ry%h_Ddata,p_Dy)
      NEQ = rx%NEQ

      ! Get the offset positions from the vector structures
      ioffsetx = rx%iidxFirstEntry
      ioffsety = rx%iidxFirstEntry

      ! Change the pointer to point to the subvector
      p_Dx => p_Dx(ioffsetx:ioffsetx+NEQ-1)
      p_Dy => p_Dy(ioffsetx:ioffsety+NEQ-1)
      
      ! Perform the multiplication
      IF (cx .NE. 0.0_DP) THEN
      
        IF (cy .EQ. 0.0_DP) THEN
        
          ! cy = 0. We have simply to make matrix*vector without adding ry.
          ! Multiply the first entry in each line of the matrix with the
          ! corresponding entry in rx and add it to ry.
          ! Don't multiply with cy, this comes later.
          !
          ! What is this complicated IF-THEN structure for?
          ! Well, to prevent an initialisation of rx with zero in case cy=0!
          
          DO irow=1,NEQ
            icol = p_Kcol(p_Kld(irow))
            p_Dy(irow) = p_Dx(icol) * p_DA(p_Kld(irow))
          END DO
          
          ! Now we have an initial ry where we can do a usual MV
          ! with the rest of the matrix...
          
        ELSE 
        
          ! cy <> 0. We have to perform matrix*vector + vector.
          ! What we actually calculate here is:
          !    ry  =  cx * A * x  +  cy * y
          !        =  cx * ( A * x  +  cy/cx * y).
          !
          ! Scale down y:
        
          dtmp = cy/cx
          IF (dtmp .NE. 1.0_DP) THEN
            CALL lalg_vectorScaleDble(p_Dy,dtmp)
          END IF
          
          ! Multiply the first entry in each line of the matrix with the
          ! corresponding entry in rx and add it to the (scaled) ry.
          
          DO irow=1,NEQ
            ICOL = p_Kcol(p_Kld(irow))
            p_Dy(irow) = p_Dx(icol)*p_DA(p_Kld(irow)) + p_Dy(irow) 
          END DO
          
        ENDIF
        
        ! Multiply the rest of rx with the matrix and add it to ry:
        
        DO irow=1,NEQ
          DO icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            p_Dy(irow) = p_Dy(irow) + p_DA(icol)*p_Dx(p_Kcol(icol))
          END DO
        END DO
        
        ! Scale by cy, finish.
        
        IF (cx .NE. 1.0_DP) THEN
          CALL lalg_vectorScaleDble (p_Dy,cx)
        END IF
        
      ELSE 
        ! cx = 0. The formula is just a scaling of the vector ry!
        CALL lalg_vectorScaleDble(p_Dy,cy)
      ENDIF
      
    END SUBROUTINE
   
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_duplicateVector (roldVector,rnewVector)
  
!<description>
  ! Duplicates an existing vector. All structural data from roldVector
  ! is assigned to rnewVector. A new handle for vector data will be allocated
  ! in rnewVector and the data of the vector in roldVector is copied
  ! to the new data array.
!</description>
  
!<input>
  ! Vector to copy
  TYPE(t_vectorScalar), INTENT(IN)                :: roldVector
!</input>

!<output>
  ! The new vector which will be a copy of roldvector
  TYPE(t_vectorScalar), INTENT(OUT)               :: rnewVector
!</output>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Fsource,p_Fdest
  INTEGER(PREC_VECIDX) :: NEQ

  ! At first, copy all 'local' data.
  rnewVector = roldVector
  
  ! Then allocate a new array for the content in the same data
  ! format as the vector and copy the data.
  ! We can't use storage_copy here, as the content of roldVector
  ! maybe > NEQ!
  
  NEQ = roldVector%NEQ
  CALL storage_new ('lsyssc_duplicateVector','vec-copy',NEQ,&
                    roldVector%cdataType, rnewVector%h_Ddata, &
                    ST_NEWBLOCK_NOINIT)

  ! Take care of the index of the first entry when copying the data!

  SELECT CASE (roldVector%cdataType)
  CASE (ST_DOUBLE)
    CALL storage_getbase_double (roldVector%h_Ddata,p_Dsource)
    CALL storage_getbase_double (rnewVector%h_Ddata,p_Ddest)
    CALL lalg_vectorCopyDble (&
      p_Dsource(roldVector%iidxFirstEntry:roldVector%iidxFirstEntry+NEQ-1),p_Ddest)
    
  CASE (ST_SINGLE)
    CALL storage_getbase_single (roldVector%h_Ddata,p_Fsource)
    CALL storage_getbase_single (rnewVector%h_Ddata,p_Fdest)
    CALL lalg_vectorCopySngl (&
      p_Fsource(roldVector%iidxFirstEntry:roldVector%iidxFirstEntry+NEQ-1),p_Fdest)

  CASE DEFAULT
    PRINT *,'lsyssc_duplicateVector: Unsupported data type!'
    STOP
  END SELECT
   
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_duplicateMatrix (rsourceMatrix,rdestMatrix,idupFlag)
  
!<description>
  ! Duplicates an existing matrix. All structural data from roldVector
  ! is assigned to rnewMatrix. Whether the entries and/or structure
  ! is duplicated is decided on idupFlag.
!</description>
  
!<input>
  ! Source matrix.
  TYPE(t_matrixScalar), INTENT(IN)               :: rsourceMatrix
  
  ! Duplication flag. Decides on what is to copy and whether some information
  ! is to be shared with the old matrix by using the same handles.
  ! One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_ASIS      : Duplicate by ownership. What belongs to the 
  !   matrix is duplicated, what belongs to another matrix is left as-is.
  ! LSYSSC_DUP_CONTENT   : Duplicate the content, share the structure
  ! LSYSSC_DUP_ALL       : Duplicate both, structure and content
  INTEGER, INTENT(IN)                            :: idupflag
  
!</input>

!<output>
  ! Destination matrix.
  TYPE(t_matrixScalar), INTENT(OUT)               :: rdestMatrix
!</output>  

!</subroutine>

  ! local variables
  LOGICAL :: bdupContent, bdupStructure

  ! At first, copy all 'local' data.
  rdestMatrix = rsourceMatrix
  
  ! Now decide what to replace by a copy.
  SELECT CASE (idupflag)
  CASE (LSYSSC_DUP_ASIS)
    bdupContent = IAND(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .EQ. 0
    bdupStructure = IAND(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0
  CASE (LSYSSC_DUP_CONTENT)
    bdupContent = .FALSE.
    bdupStructure = .TRUE.
  CASE (LSYSSC_DUP_ALL)
    bdupContent = .TRUE.
    bdupStructure = .TRUE.
  CASE DEFAULT
    PRINT *,'lsyssc_duplicateMatrix: invalid idubflag'
    STOP
  END SELECT
  
  ! Depending on the matrix format choose how to duplicate
  SELECT CASE (rsourceMatrix%imatrixFormat)
  CASE (LSYSSC_MATRIX9)
    IF (bdupContent) THEN
      ! Put the destination handle to ST_NOHANDLE, so storage_copy
      ! will create new ones.
      rdestMatrix%h_DA = ST_NOHANDLE
      CALL storage_copy(rsourceMatrix%h_DA, rdestMatrix%h_DA)
    END IF
    IF (bdupStructure) THEN
      ! Put the destination handles to ST_NOHANDLE, so storage_copy
      ! will create new ones.
      rdestMatrix%h_Kcol = ST_NOHANDLE
      rdestMatrix%h_Kld = ST_NOHANDLE
      rdestMatrix%h_Kdiagonal = ST_NOHANDLE
      CALL storage_copy(rsourceMatrix%h_Kcol, rdestMatrix%h_Kcol)
      CALL storage_copy(rsourceMatrix%h_Kld, rdestMatrix%h_Kld)
      CALL storage_copy(rsourceMatrix%h_Kdiagonal, rdestMatrix%h_Kdiagonal)
    END IF

  CASE (LSYSSC_MATRIX7)
    IF (bdupContent) THEN
      ! Put the destination handle to ST_NOHANDLE, so storage_copy
      ! will create new ones.
      rdestMatrix%h_DA = ST_NOHANDLE
      CALL storage_copy(rsourceMatrix%h_DA, rdestMatrix%h_DA)
    END IF
    IF (bdupStructure) THEN
      ! Put the destination handles to ST_NOHANDLE, so storage_copy
      ! will create new ones.
      rdestMatrix%h_Kcol = ST_NOHANDLE
      rdestMatrix%h_Kld = ST_NOHANDLE
      CALL storage_copy(rsourceMatrix%h_Kcol, rdestMatrix%h_Kcol)
      CALL storage_copy(rsourceMatrix%h_Kld, rdestMatrix%h_Kld)
    END IF
    
  CASE DEFAULT
    PRINT *,'lsyssc_duplicateMatrix: Unsupported Matrix format'
    STOP
  END SELECT
  
  ! Set the duplication flag correctly
  SELECT CASE (idupflag)
  CASE (LSYSSC_DUP_ASIS)
    ! Nothing to do
  CASE (LSYSSC_DUP_CONTENT)
    rdestMatrix%imatrixSpec = IOR(IAND(rsourceMatrix%imatrixSpec,&
                                NOT(LSYSSC_MSPEC_ISCOPY)),&
                              LSYSSC_MSPEC_STRUCTUREISCOPY)
  CASE (LSYSSC_DUP_STRUCTURE)
    rdestMatrix%imatrixSpec = IOR(IAND(rsourceMatrix%imatrixSpec,&
                                NOT(LSYSSC_MSPEC_ISCOPY)),&
                              LSYSSC_MSPEC_CONTENTISCOPY)
  CASE (LSYSSC_DUP_ALL)
    rdestMatrix%imatrixSpec = IAND(rsourceMatrix%imatrixSpec,&
                                NOT(LSYSSC_MSPEC_ISCOPY))
  CASE DEFAULT
    PRINT *,'lsyssc_duplicateMatrix: invalid idubflag'
    STOP
  END SELECT
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_releaseVector (rvector)
  
!<description>
  ! Releases a vector from memory. The vector structure is cleaned up.
!</description>
  
!<inputoutput>
  
  ! Vector to release.
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rvector
  
!</inputoutput>

!</subroutine>

  ! Clean up the data structure.
  ! Don't release the vector data if the handle belongs to another
  ! vector!
  IF ((.NOT. rvector%bisCopy) .AND. (rvector%h_Ddata .NE. ST_NOHANDLE)) THEN
    CALL storage_free(rvector%h_Ddata)
  ELSE
    rvector%h_Ddata = ST_NOHANDLE
  END IF
  rvector%NEQ = 0
  rvector%cdataType = ST_NOHANDLE
  rvector%iidxFirstEntry = 1
  rvector%bisCopy = .FALSE.
  rvector%p_rspatialDiscretisation => NULL()
  rvector%p_rdiscreteBC => NULL()
  rvector%p_rdiscreteBCfict => NULL()
   
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_releaseMatrix (rmatrix)
  
!<description>
  ! Releases a matrix from memory. The structure is cleaned up.
  ! All data vectors belonging to the structure are released;
  ! data structures belonging to other matrices are not.
!</description>
  
!<inputoutput>
  
  ! Matrix to release.
  TYPE(t_matrixScalar), INTENT(INOUT)               :: rmatrix
  
!</inputoutput>

!</subroutine>

  ! Which handles do we have to release?
  !
  ! Release the matrix data if the handle is not a copy of another matrix
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .EQ. 0) THEN
    SELECT CASE (rmatrix%imatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
      ! Release matrix data, structure 9,7
      IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_free(rmatrix%h_DA)
      END IF
    END SELECT
  END IF
  
  ! Release the structure if it doesn't belong to another vector
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0) THEN
    ! Release matrix structure
    SELECT CASE (rmatrix%imatrixFormat)
    CASE (LSYSSC_MATRIX9)
      IF (rmatrix%h_Kcol .NE. ST_NOHANDLE)      CALL storage_free(rmatrix%h_Kcol)
      IF (rmatrix%h_Kld .NE. ST_NOHANDLE)       CALL storage_free(rmatrix%h_Kld)
      IF (rmatrix%h_Kdiagonal .NE. ST_NOHANDLE) CALL storage_free(rmatrix%h_Kdiagonal)
    CASE (LSYSSC_MATRIX7)
      IF (rmatrix%h_Kcol .NE. ST_NOHANDLE) CALL storage_free(rmatrix%h_Kcol)
      IF (rmatrix%h_Kld .NE. ST_NOHANDLE)  CALL storage_free(rmatrix%h_Kld)
    END SELECT
  END IF
  
  ! Clean up the rest
  rmatrix%h_DA        = ST_NOHANDLE
  rmatrix%h_Kcol      = ST_NOHANDLE
  rmatrix%h_Kld       = ST_NOHANDLE
  rmatrix%h_Kdiagonal = ST_NOHANDLE
  rmatrix%cdataType   = ST_DOUBLE
  rmatrix%NA  = 0
  rmatrix%NEQ = 0
  rmatrix%dScaleFactor = 1.0_DP
  rmatrix%p_rspatialDiscretisation => NULL()
  rmatrix%p_rdiscreteBC     => NULL()
  rmatrix%p_rdiscreteBCfict => NULL()

  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_sortVector (rvector,rtemp,isortStrategy,h_IsortPermutation)
  
!<description>
  ! Resorts the entries of the given vector rvector or unsorts it.
  ! rvector is the vector to be resorted, rtemp is a temporary vector and
  ! isortStrategy is a type flag identifying the sorting algorithm.
  !
  ! This routine can also be used to assign an unsorted vector a sorting 
  ! strategy without actually resorting it. For this purpose, isortStrategy
  ! should be '- SSTRAT_xxxx' and h_IsortPermutation a handle to a
  ! sorting strategy permutation.
!</description>
  
!<inputoutput>
  ! Vector to resort
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rvector

  ! A temporary vector of the same size and data type as rvector
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rtemp
!</inputoutput>

!<input>
  ! Identifier for the sorting strategy to apply to the vector.
  ! This is usually one of the SSTRAT_xxxx constants from the module
  ! 'sortstrategy', although it's actually used here as follows:
  ! <=0: Calculate the unsorted vector
  !  >0: Resort the vector according to a permutation;
  !      this is either the permutation specified in the vector
  !      or that one identified by h_IsortPermutation
 INTEGER, INTENT(IN)                                :: isortStrategy
  
  ! OPTIONAL: Handle to permutation to use for resorting the vector.
  ! 
  ! The array must be of the form
  !    array [1..2*NEQ] of integer
  ! with entries (1..NEQ)       = permutation
  ! and  entries (NEQ+1..2*NEQ) = inverse permutation.
  !
  ! If not specified, the associated permutation rvector%h_IsortPermutation
  ! is used.
  ! If specified and isortStrategy>0, the vector is resorted according to 
  ! the permutation h_IsortPermutation (probably unsorted before if necessary).
  !
  ! In any case, the associated permutation of the vector
  ! rvector%h_IsortPermutation is changed to h_IsortPermutation.
  ! Remark: The memory associated to the previous sorting strategy
  !  in this case is not released automatically; this has to be done by the
  !  application!
  INTEGER, INTENT(IN), OPTIONAL                     :: h_IsortPermutation
!</input> 

!</subroutine>

  ! local variables
  INTEGER :: h_Iperm
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata,p_Fdata2
  INTEGER(PREC_VECIDX) :: NEQ
  
    ! Desired sorting strategy and currently active sorting strategy identical?
    IF (.NOT. PRESENT(h_IsortPermutation)) THEN
    
      IF (isortStrategy .EQ. rvector%isortStrategy) RETURN
      IF ((isortStrategy .LE. 0) .AND. (rvector%isortStrategy .LE. 0)) RETURN
    
    ELSE
    
      IF ((isortStrategy .LE. 0) .AND. (rvector%isortStrategy .LE. 0)) THEN
        IF (h_IsortPermutation .NE. rvector%h_IsortPermutation) THEN
          ! Vector is unsorted and should stay unsorted, but
          ! permutation should change.
          rvector%h_IsortPermutation = h_IsortPermutation
        END IF
        RETURN
      END IF
      IF ((isortStrategy .GT. 0) .AND. (rvector%isortStrategy .GT. 0) .AND.&
          (h_IsortPermutation .EQ. rvector%h_IsortPermutation)) RETURN
    
    END IF
    
    NEQ = rvector%NEQ
    
    ! Get pointers to the vector data
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(rvector%h_Ddata,p_Ddata)
      CALL storage_getbase_double(rtemp%h_Ddata,p_Ddata2)
      ! Don't forget to take care of where the first entry is!
      p_Ddata => p_Ddata(rvector%iidxFirstEntry:rvector%iidxFirstEntry+NEQ-1)
      p_Ddata2 => p_Ddata2(rtemp%iidxFirstEntry:rtemp%iidxFirstEntry+NEQ-1)
    CASE (ST_SINGLE)   
      CALL storage_getbase_single(rvector%h_Ddata,p_Fdata)
      CALL storage_getbase_single(rtemp%h_Ddata,p_Fdata2)
      ! Don't forget to take care of where the first entry is!
      p_Fdata => p_Fdata(rvector%iidxFirstEntry:rvector%iidxFirstEntry+NEQ-1)
      p_Fdata2 => p_Fdata2(rtemp%iidxFirstEntry:rtemp%iidxFirstEntry+NEQ-1)
    CASE DEFAULT
      PRINT *,'lsyssc_sortVector: unsuppported data type'
      STOP
    END SELECT

    
    ! Sort the vector back?
    IF (isortStrategy .LE. 0) THEN
      ! Do it - with the associated permutation.
      CALL storage_getbase_int(rvector%h_IsortPermutation,p_Iperm)
      
      SELECT CASE (rvector%cdataType)
      CASE (ST_DOUBLE)
        ! Copy the entries to the temp vector
        CALL lalg_vectorCopyDble (p_Ddata,p_Ddata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortDble (p_Ddata2,p_Ddata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE (ST_SINGLE)   
        ! Copy the entries to the temp vector
        CALL lalg_vectorCopySngl (p_Fdata,p_Fdata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortSngl (p_Fdata2,p_Fdata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE DEFAULT
        PRINT *,'lsyssc_sortVector: unsuppported data type'
        STOP
        
      END SELECT
        
      ! Inform the vector about which sorting strategy we now use.
      rvector%isortStrategy = isortStrategy
      RETURN
    END IF
    
    ! Get the actual sorting strategy.
    h_Iperm = rvector%h_IsortPermutation
    IF (PRESENT(h_IsortPermutation)) h_Iperm = h_IsortPermutation 
    
    ! Do we have to sort back before resorting?
    IF ((h_Iperm .NE. rvector%h_IsortPermutation) .AND. &
        (rvector%h_IsortPermutation .NE. ST_NOHANDLE)) THEN
        
      ! Sort back at first - with the associated permutation
      CALL storage_getbase_int(rvector%h_IsortPermutation,p_Iperm)
      
      SELECT CASE (rvector%cdataType)
      CASE (ST_DOUBLE)
        ! Copy the entries to the temp vector
        CALL lalg_vectorCopyDble (p_Ddata,p_Ddata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortDble (p_Ddata2,p_Ddata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE (ST_SINGLE)   
        ! Copy the entries to the temp vector
        CALL lalg_vectorCopySngl (p_Fdata,p_Fdata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortSngl (p_Fdata2,p_Fdata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE DEFAULT
        PRINT *,'lsyssc_sortVector: unsuppported data type'
        STOP
        
      END SELECT
      
      ! Change the sorting strategy in the vector to the one
      ! we are now going to use. Throw away the old handle.
      rvector%h_IsortPermutation = h_Iperm
      
    END IF
    
    ! Now sort the vector according to h_Iperm
    CALL storage_getbase_int(h_Iperm,p_Iperm)
    
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      ! Copy the entries to the temp vector
      CALL lalg_vectorCopyDble (p_Ddata,p_Ddata2)
      ! Then do the sorting with the given permutation.
      CALL lalg_vectorSortDble (p_Ddata2,p_Ddata,p_Iperm(1:NEQ))
      
    CASE (ST_SINGLE)   
      ! Copy the entries to the temp vector
      CALL lalg_vectorCopySngl (p_Fdata,p_Fdata2)
      ! Then do the sorting with the given permutation.
      CALL lalg_vectorSortSngl (p_Fdata2,p_Fdata,p_Iperm(1:NEQ))
      
    CASE DEFAULT
      PRINT *,'lsyssc_sortVector: unsuppported data type'
      STOP
      
    END SELECT
    
    ! Inform the vector about which sorting strategy we now use.
    rvector%isortStrategy = isortStrategy
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_sortMatrix (rmatrix,bsortEntries,&
                                isortStrategy,h_IsortPermutation)
  
!<description>
  ! Matrix sorting. Sorts the entries and/or structure of a given
  ! matrix or unsorts them according to a permutation.
  !
  ! This routine can also be used to assign an unsorted vector a sorting 
  ! strategy without actually resorting it. For this purpose, isortStrategy
  ! should be '- SSTRAT_xxxx' and h_IsortPermutation a handle to a
  ! sorting strategy permutation.
!</description>
  
!<inputoutput>
  ! Vector to resort
  TYPE(t_matrixScalar), INTENT(INOUT)               :: rmatrix
!</inputoutput>

!<input>
  ! Sort the entries or only the structure of the matrix.
  ! = FALSE: Only sort the structure of the matrix,
  ! = TRUE : Sort both, entries and structure of the matrix.
  LOGICAL bsortEntries

  ! Identifier for the sorting strategy to apply to the matrix.
  ! This is usually one of the SSTRAT_xxxx constants from the module
  ! 'sortstrategy', although it's actually used here as follows:
  ! <=0: Calculate the unsorted matrix
  !  >0: Resort the vector according to a permutation;
  !      this is either the permutation specified in the vector
  !      or that one identified by h_IsortPermutation
 INTEGER, INTENT(IN)                                :: isortStrategy
  
  ! OPTIONAL: Handle to permutation to use for resorting the matrix.
  ! 
  ! The array must be of the form
  !    array [1..2*NEQ] of integer
  ! with entries (1..NEQ)       = permutation
  ! and  entries (NEQ+1..2*NEQ) = inverse permutation.
  !
  ! If not specified, the associated permutation rmatrix%h_IsortPermutation
  ! is used.
  ! If specified and isortStrategy>0, the matrix is resorted according to 
  ! the permutation h_IsortPermutation (probably unsorted before if necessary).
  !
  ! In any case, the associated permutation of the matrix
  ! rmatrix%h_IsortPermutation is changed to h_IsortPermutation.
  ! Remark: The memory associated to the previous sorting strategy
  !  in this case is not released automatically; this has to be done by the
  !  application!
  INTEGER, INTENT(IN), OPTIONAL                     :: h_IsortPermutation
!</input> 

!</subroutine>

  ! local variables
  INTEGER :: h_Iperm
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  INTEGER(PREC_VECIDX) :: NEQ
  
    ! Desired sorting strategy and currently active sorting strategy identical?
    IF (.NOT. PRESENT(h_IsortPermutation)) THEN
    
      IF (isortStrategy .EQ. rmatrix%isortStrategy) RETURN
      IF ((isortStrategy .LE. 0) .AND. (rmatrix%isortStrategy .LE. 0)) RETURN
    
    ELSE
    
      IF ((isortStrategy .LE. 0) .AND. (rmatrix%isortStrategy .LE. 0)) THEN
        IF (h_IsortPermutation .NE. rmatrix%h_IsortPermutation) THEN
          ! Vector is unsorted and should stay unsorted, but
          ! permutation should change.
          rmatrix%h_IsortPermutation = h_IsortPermutation
        END IF
        RETURN
      END IF
      IF ((isortStrategy .GT. 0) .AND. (rmatrix%isortStrategy .GT. 0) .AND.&
          (h_IsortPermutation .EQ. rmatrix%h_IsortPermutation)) RETURN
    
    END IF
    
    NEQ = rmatrix%NEQ
    
    ! Sort the matrix back?
    IF (isortStrategy .LE. 0) THEN
      ! Get the permutation that describes how to resort the matrix:
      CALL storage_getbase_int(rmatrix%h_IsortPermutation,p_Iperm)
      
      ! Exchange the roles of the first and second part of p_Iperm.
      ! This makes the permutation to the inverse permutation and vice versa.
      ! Call the resort-subroutine with that to do the actual sorting.
      CALL do_matsort (rmatrix,p_Iperm(NEQ+1:NEQ*2),p_Iperm(1:NEQ),bsortEntries)
      
      ! Inform the vector about which sorting strategy we now use.
      rmatrix%isortStrategy = isortStrategy
      RETURN
    END IF
    
    ! Get the actual sorting strategy.
    h_Iperm = rmatrix%h_IsortPermutation
    IF (PRESENT(h_IsortPermutation)) h_Iperm = h_IsortPermutation 
    
    ! Do we have to sort back before resorting?
    IF ((h_Iperm .NE. rmatrix%h_IsortPermutation) .AND. &
        (rmatrix%h_IsortPermutation .NE. ST_NOHANDLE)) THEN

      ! Sort back at first - with the associated permutation
      CALL storage_getbase_int(rmatrix%h_IsortPermutation,p_Iperm)
      
      ! Exchange the roles of the first and second part of p_Iperm.
      ! This makes the permutation to the inverse permutation and vice versa.
      ! Call the resort-subroutine with that to do the actual sorting.
      CALL do_matsort (rmatrix,p_Iperm(NEQ+1:NEQ*2),p_Iperm(1:NEQ),bsortEntries)
    END IF
    
    ! Now sort the vector according to h_Iperm
    CALL storage_getbase_int(h_Iperm,p_Iperm)

    ! This time, we don't exchange the roles of the permutation and
    ! its inverse :-)
    CALL do_matsort (rmatrix,p_Iperm(1:NEQ),p_Iperm(NEQ+1:NEQ*2),bsortEntries)
    
    ! Inform the vector about which sorting strategy we now use.
    rmatrix%isortStrategy = isortStrategy
  
  CONTAINS
    
    !----------------------------------------------------------------
    ! Sort matrix or matrix entries.
    ! This calls the actual resorting routine, depending
    ! on the information tags in the matrix.
    
    SUBROUTINE do_matsort (rmatrix,Itr1,Itr2,bsortEntries)
    
    ! The matrix to be resorted
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! The transformation to use
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Itr1
    
    ! The inverse transformation
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Itr2
    
    ! TRUE  = sort matrix structure + entries
    ! FALSE = sort only matrix structure
    LOGICAL, INTENT(IN) :: bsortEntries
    
    ! local variables
    
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_KldTmp,p_Kdiag
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol,p_KcolTmp
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_DdataTmp
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata,p_FdataTmp
    TYPE(t_matrixScalar) :: rtempMatrix
    INTEGER(PREC_VECIDX) :: NEQ
    
      NEQ = rmatrix%NEQ
    
      ! Which matrix configuration do we have?
      SELECT CASE (rmatrix%imatrixFormat)
      CASE (LSYSSC_MATRIX9)
        ! Duplicate the matrix before resorting.
        ! We need either a copy only of the structure or of the full matrix.
        IF (.NOT. bsortEntries) THEN
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,LSYSSC_DUP_STRUCTURE)
        ELSE
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,LSYSSC_DUP_ALL)
        END IF
        
        ! Get the structure of the original and the temporary matrix
        CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
        CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
        CALL storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiag)
        CALL storage_getbase_int (rtempMatrix%h_Kcol,p_KcolTmp)
        CALL storage_getbase_int (rtempMatrix%h_Kld,p_KldTmp)
        
        IF (.NOT. bsortEntries) THEN
        
          ! Sort only the structure of the matrix, keep the entries
          ! unchanged.
          CALL lsyssc_sortMat9Struc (p_Kcol, p_KcolTmp, p_KldTmp, p_KldTmp, &
                                     p_Kdiag, Itr1, Itr2, NEQ)        
        
        ELSE
        
          ! Sort structure + entries
          SELECT CASE (rmatrix%cdataType)
          CASE (ST_DOUBLE)
            ! Double precision version
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata)
            CALL storage_getbase_double (rtempMatrix%h_Da,p_DdataTmp)
            ! Switch the two permutations. This gives the 'back-sorting'
            CALL lsyssc_sortMat9_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, p_Kdiag, &
                                         Itr1, Itr2, NEQ)        
          CASE (ST_SINGLE)
            ! Single precision version
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata)
            CALL storage_getbase_double (rtempMatrix%h_Da,p_DdataTmp)
            ! Switch the two permutations. This gives the 'back-sorting'
            CALL lsyssc_sortMat9_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, p_Kdiag, &
                                         Itr1, Itr2, NEQ)        
          CASE DEFAULT
            PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
            STOP
          END SELECT
        END IF

        ! Remove temp matrix
        CALL lsyssc_releaseMatrix(rtempMatrix)
        
      CASE (LSYSSC_MATRIX7)
        ! Duplicate the matrix before resorting.
        ! We need either a copy only of the structure or of the full matrix.
        IF (.NOT. bsortEntries) THEN
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,LSYSSC_DUP_STRUCTURE)
        ELSE
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,LSYSSC_DUP_ALL)
        END IF
        
        ! Get the structure of the original and the temporary matrix
        CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
        CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
        CALL storage_getbase_int (rtempMatrix%h_Kcol,p_KcolTmp)
        CALL storage_getbase_int (rtempMatrix%h_Kld,p_KldTmp)
        
        IF (.NOT. bsortEntries) THEN
        
          ! Sort only the structure of the matrix, keep the entries
          ! unchanged.
          CALL lsyssc_sortMat7Struc (p_Kcol, p_KcolTmp, p_KldTmp, p_KldTmp, &
                                     Itr1, Itr2, NEQ)        
        
        ELSE
        
          ! Sort structure + entries
          SELECT CASE (rmatrix%cdataType)
          CASE (ST_DOUBLE)
            ! Double precision version
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata)
            CALL storage_getbase_double (rtempMatrix%h_Da,p_DdataTmp)
            ! Switch the two permutations. This gives the 'back-sorting'
            CALL lsyssc_sortMat7_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, &
                                         Itr1, Itr2, NEQ)        
          CASE (ST_SINGLE)
            ! Single precision version
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata)
            CALL storage_getbase_double (rtempMatrix%h_Da,p_DdataTmp)
            ! Switch the two permutations. This gives the 'back-sorting'
            CALL lsyssc_sortMat7_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, &
                                         Itr1, Itr2, NEQ)        
          CASE DEFAULT
            PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
            STOP
          END SELECT
        END IF

        ! Remove temp matrix
        CALL lsyssc_releaseMatrix(rtempMatrix)
    
      CASE DEFAULT
        PRINT *,'lsyssc_sortMatrix: Unsupported matrix format!'
        STOP
        
      END SELECT

    END SUBROUTINE
    
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7_double (Da, DaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 7 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(DP), DIMENSION(:), INTENT(IN) :: DaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Da

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq

      ! Which row should be moved to here?
      iidx = Itr1(i)

      ! Copy the diagonal element.
      ! The new row starts at position ildIdx. Save this to Ild.
      Da(ildIdx) = DaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1

      ! The row has length...
      isize = ih2Idx - ih1Idx + 1

      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j=1, isize
        ! Get the matrix entry of the current column and write it to DA
        Da(ildIdx) = DaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) !Achtung: Ist dies so richtig ? Ih1(j)->Ih2(j) ?

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7_single (Fa, FaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version.
    !
    ! Storage technique 7 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(SP), DIMENSION(:), INTENT(IN) :: FaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(SP), DIMENSION(:), INTENT(OUT) :: Fa

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    
    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
      ! Which row should be moved to here?
      iidx = Itr1(i)

      ! Copy the diagonal element.
      ! The new row starts at position ildIdx. Save this to Ild.
      Fa(ildIdx) = FaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1

      ! The row has length...
      isize = ih2Idx - ih1Idx + 1

      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j=1, isize
        ! Get the matrix entry of the current column and write it to DA
        Fa(ildIdx) = FaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)
    
    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7Struc (Icol, IcolH, Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the structure of the given matrix, corresponding to Itr1/Itr2.
    ! Only the structure of a given matrix is resorted, the routine 
    ! doesn't handle the entries.
    !
    ! Storage technique 7 version. 
  !</description>
    
  !<input>

    ! Number of equations
    INTEGER(PREC_VECIDX) , INTENT(IN) :: neq

    ! Permutation of 1..neq describing how to resort
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing how to sort
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<inputoutput>

    ! Column description of matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Icol

    ! Row description of matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(INOUT) :: Ild
    
    ! Column structure of source matrix -> resorted matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: IcolH
    
    ! Row positions of source matrix -> resorted matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(INOUT) :: IldH
    
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
    
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
      ! Which row should be moved to here?
      iidx = Itr1(i)

      ! Copy the diagonal element.
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx+1

      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1

      ! The row has length...
      isize = ih2Idx-ih1Idx+1

      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, Ih2Idx
        Ih1(j-ih1Idx+1)=Itr2(IcolH(j))
        Ih2(j-ih1Idx+1)=j
      END DO

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j=1, isize
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx+1
      END DO
    END DO
    Ild(neq+1) = IldH(neq+1) 
    
    DEALLOCATE(Ih1,Ih2)
  
  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9_double (Da, DaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 9 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(DP), DIMENSION(:), INTENT(IN) :: DaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH

    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the inverse permutation
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Da

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
    ! Positions of diagonal elements in the destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Idiag
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1
    
    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
    
      ! Which row should be moved to here?
      iidx = Itr1(i)
      
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)
      ih2Idx = IldH(iidx+1)-1
      
      ! The row has length...
      isize = ih2Idx - ih1Idx + 1
      
      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      
      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      DO idiagidx=1,isize
        IF (Ih1(idiagidx) .GE. i) EXIT
      END DO
      
      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        
        ! Get the matrix entry of the current column and write it to DA
        Da(ildIdx) = DaH(Ih2(j))
        
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9_dingle (Fa, FaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version.
    !
    ! Storage technique 9 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(SP), DIMENSION(:), INTENT(IN) :: FaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH

    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the inverse permutation
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(SP), DIMENSION(:), INTENT(OUT) :: Fa

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
    ! Positions of diagonal elements in the destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Idiag
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1
    
    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
    
      ! Which row should be moved to here?
      iidx = Itr1(i)
      
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)
      ih2Idx = IldH(iidx+1)-1
      
      ! The row has length...
      isize = ih2Idx - ih1Idx + 1
      
      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      
      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      DO idiagidx=1,isize
        IF (Ih1(idiagidx) .GE. i) EXIT
      END DO
      
      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        
        ! Get the matrix entry of the current column and write it to DA
        Fa(ildIdx) = FaH(Ih2(j))
        
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9Struc (Icol, IcolH, &
                                   Ild, IldH, Idiag, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the structure of the given matrix, corresponding to Itr1/Itr2.
    ! Only the structure of a given matrix is resorted, the routine 
    ! doesn't handle the entries.
    !
    ! Storage technique 9 version. 
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH

    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the inverse permutation
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
    ! Positions of diagonal elements in the destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Idiag
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1
    
    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
    
      ! Which row should be moved to here?
      iidx = Itr1(i)
      
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)
      ih2Idx = IldH(iidx+1)-1
      
      ! The row has length...
      isize = ih2Idx - ih1Idx + 1
      
      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      
      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      DO idiagidx=1,isize
        IF (Ih1(idiagidx) .GE. i) EXIT
      END DO
      
      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortCR (Ih1, Ih2)
  
  !<description>
    ! Performs bubble sort for the vector Ih1 and does the
    ! same swaps also on vector Ih2
  !</description>
    
  !<inputoutput>
      
    ! column vector
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ih1
      
    ! row vector; same size as Ih1
    INTEGER(PREC_VECIDX), DIMENSION(SIZE(Ih1)), INTENT(INOUT) :: Ih2
      
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(PREC_VECIDX) :: iaux, icomp
    LOGICAL :: bmore
    
    bmore = .TRUE.
    ! While there are swaps necessary...
    DO WHILE (bmore)
      bmore = .FALSE.
      ! Test for all elements minus the last one
      DO icomp=1, SIZE(Ih1)-1
        ! If the order is correct, next entry; othewise swap entries
        IF (Ih1(icomp) .GT. Ih1(icomp+1)) THEN
          iaux = Ih1(icomp)
          Ih1(icomp) = Ih1(icomp+1)
          Ih1(icomp+1) = iaux
          iaux = Ih2(icomp)
          Ih2(icomp) = Ih2(icomp+1)
          Ih2(icomp+1) = iaux
          bmore = .TRUE.
        ENDIF
      END DO  
    END DO
  
  END SUBROUTINE 

END MODULE
