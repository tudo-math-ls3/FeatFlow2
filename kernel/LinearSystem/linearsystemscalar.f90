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
!# 1.) lsyssc_scalarMatVec
!#     -> Multiply a scalar matrix with a scalar vector
!#
!# 2.) lsyssc_releaseScalarMatrix
!#     -> Release a scalar matrix from memory.
!#
!# 3.) lsyssc_releaseScalarVector
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
    
    ! A pointer to the spatial discretisation or NULL(), if the vector is
    ! just an array, not belonging to any discretisation.
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscretisation => NULL()
    
    ! A pointer to discretised boundary conditions for real boundary components.
    ! For every discrete BC that is valid for this scalar vector, there
    ! is an entry in the following list.
    TYPE(t_discreteBC), DIMENSION(:), POINTER  :: p_RdiscreteBC => NULL()
    
    ! A pointer to discretised boundary conditions for fictitious boundary
    ! components
    TYPE(t_discreteBC), DIMENSION(:), POINTER  :: p_RdiscreteBCfict => NULL()
    
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
    TYPE(t_discreteBC), DIMENSION(:), POINTER  :: p_RdiscreteBC     => NULL()
    
    ! A pointer to discretised boundary conditions for fictitious boundary
    ! components
    TYPE(t_discreteBC), DIMENSION(:), POINTER  :: p_RdiscreteBCfict => NULL()
    
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

!<subroutine>
  
  SUBROUTINE lsyssc_scalarMatVec (rMatrix, rx, ry, cx, cy)
  
!<description>
  ! Performs a matrix vector multiplicationwith a given scalar matrix:
  !    $$ Dy   =   cx * rMatrix * rx   +   cy * ry $$
!</description>
  
!<input>
  
  ! Scalar matrix
  TYPE(t_matrixScalar), INTENT(IN)                  :: rMatrix

  ! Vector to multiply with the matrix.
  TYPE(t_vectorScalar), INTENT(IN)                  :: rx
  
  ! Multiplicative factor for rx
  REAL(DP), INTENT(IN)                              :: cx

  ! Multiplicative factor for ry
  REAL(DP), INTENT(IN)                              :: cy
  
!</input>

!<inputoutput>
  ! Additive vector. Receives the result of the matrix-vector multiplication
  TYPE(t_vectorScalar), INTENT(OUT)                 :: ry
!</inputoutput>

  !...
  PRINT *,'MV not implemented'
   
  END SUBROUTINE
  
!</subroutine>
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_releaseScalarVector (rvector)
  
!<description>
  ! Releases a vector from memory. The vector structure is cleaned up.
!</description>
  
!<inputoutput>
  
  ! Vector to release.
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rvector
  
!</inputoutput>

!</subroutine>

  ! Clean up the data structure
  IF (rvector%h_Ddata .NE. ST_NOHANDLE) THEN
    CALL storage_free(rvector%h_Ddata)
  END IF
  rvector%NEQ = 0
  rvector%cdataType = ST_NOHANDLE
  rvector%iidxFirstEntry = 1
  rvector%p_rspatialDiscretisation => NULL()
  rvector%p_RdiscreteBC => NULL()
  rvector%p_RdiscreteBCfict => NULL()
   
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_releaseScalarMatrix (rmatrix)
  
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
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .EQ. 0) THEN
    SELECT CASE (rmatrix%imatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
      ! Release matrix data, structure 9,7
      IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_free(rmatrix%h_DA)
      END IF
    END SELECT
  END IF
  
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
  rmatrix%cdataType = ST_DOUBLE
  rmatrix%NA  = 0
  rmatrix%NEQ = 0
  rmatrix%dScaleFactor = 1.0_DP
  rmatrix%p_rspatialDiscretisation => NULL()
  rmatrix%p_RdiscreteBC     => NULL()
  rmatrix%p_RdiscreteBCfict => NULL()

  END SUBROUTINE

END MODULE
