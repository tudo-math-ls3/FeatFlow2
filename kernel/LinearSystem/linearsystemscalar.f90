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
  
  ! Complete matrix is duplicate of another matrix and shares structure
  ! and entries via the same pointers
  INTEGER, PARAMETER :: LSYSSC_MSPEC_ISCOPY =          2**0

  ! Matrix structure is a copy of another matrix, shared via the same
  ! pointers. Matrix entries are different
  INTEGER, PARAMETER :: LSYSSC_MSPEC_STRUCTUREISCOPY = 2**1
  
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
  REAL(SP), DIMENSION(:), POINTER :: p_Sdata1dp
  REAL(SP), DIMENSION(:), POINTER :: p_Sdata2dp
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
    
    CALL storage_getbase_single (rx%h_Ddata,p_Sdata1dp)
    CALL storage_getbase_single (ry%h_Ddata,p_Sdata2dp)
    
    ! Change the pointer to point to the subvector
    p_Sdata1dp => p_Sdata1dp(ioffsetx:ioffsetx+rx%NEQ-1)
    p_Sdata2dp => p_Sdata2dp(ioffsety:ioffsety+ry%NEQ-1)

    ! Perform the scalar product
    res = p_Sdata1dp(1)*p_Sdata2dp(1)
    DO i=2,rx%NEQ
      res = res + p_Sdata1dp(i)*p_Sdata2dp(i)
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
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .EQ. 0) THEN
    SELECT CASE (rmatrix%imatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
      ! Release matrix data, structure 9,7
      IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_free(rmatrix%h_DA)
      END IF
    END SELECT
  END IF
  
  ! Release the structure if it doesn't belong to another vector
  IF ((IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .EQ. 0) .AND. &
      (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0)) THEN
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

END MODULE
