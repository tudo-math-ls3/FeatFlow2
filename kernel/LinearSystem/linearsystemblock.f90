!##############################################################################
!# ****************************************************************************
!# <name> linearsystemblock </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic definitions and routines to maintain block
!# matrices and block vectors. A block matrix is realised as 2D-array of
!# scalar matrices, while a block vector is realised as 1D-array of
!# scalar vectors.
!#
!# The following routines can be found here:
!#
!#  1.) lsysbl_createVecBlockDirect
!#      -> Create a block vector by specifying the size of the subblocks
!# 
!#  2.) lsysbl_createVecBlockIndirect
!#      -> Create a block vector by copying the structure of another block
!#         vector
!#
!#  3.) lsysbl_createVecBlockIndMat
!#      -> Create a block vector according to a block matrix.
!#
!#  4.) lsysbl_createMatFromScalar
!#      -> Creates a 1x1 block matrix from a scalar matrix
!#
!#  5.) lsysbl_createVecFromScalar
!#      -> Creates a 1-block vector from a scalar vector
!#
!#  6.) lsysbl_assignDiscretIndirect
!#      -> Assign discretisation related information of one vector
!#         to another
!#
!#  7.) lsysbl_releaseVector
!#      -> Release a block vector from memory
!#
!#  8.) lsysbl_blockMatVec
!#      -> Multiply a block matrix with a block vector
!#
!#  9.) lsysbl_vectorCopy
!#       -> Copy a block vector over to another one
!#
!# 10.) lsysbl_vectorScale
!#      -> Scale a block vector by a constant
!#
!# 11.) lsysbl_vectorClear
!#      -> Clear a block vector
!#
!# 12.) lsysbl_vectorLinearComb
!#      -> Linear combination of two block vectors
!#
!# 13.) lsysbl_scalarProduct
!#      -> Calculate a scalar product of two vectors
!#
!# </purpose>
!##############################################################################

MODULE linearsystemblock

  USE fsystem
  USE storage
  USE linearsystemscalar
  USE linearalgebra
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Global constants for block matrices">

  ! Maximum number of blocks that is allowed in 'usual' block matrices,
  ! supported by the routines in this file.
  INTEGER, PARAMETER :: LSYSBL_MAXBLOCKS = 16
  
!</constantblock>

!<constantblock description="Flags for the matrix specification bitfield">

  ! Standard matrix
  INTEGER, PARAMETER :: LSYSBS_MSPEC_GENERAL           =        0
  
  ! Block matrix is a scalar (i.e. 1x1) matrix.
  INTEGER, PARAMETER :: LSYSBS_MSPEC_SCALAR            =     2**0

  ! Block matrix is of saddle-point type:
  !  (A  B1)
  !  (B2 0 )
  INTEGER, PARAMETER :: LSYSBS_MSPEC_SADDLEPOINT       =     2**1

  ! Block matrix is nearly of saddle-point type 
  !  (A  B1)
  !  (B2 C )
  ! with C~0 being a stabilisation matrix
  INTEGER, PARAMETER :: LSYSBS_MSPEC_NEARLYSADDLEPOINT =     2**2

!</constantblock>

!</constants>

!<types>
  
!<typeblock>
  
  ! A block vector that can be used by block linear algebra routines.
  
  TYPE t_vectorBlock
    
    ! Total number of equations in the vector
    INTEGER                    :: NEQ
    
    ! Handle identifying the vector entries or = ST_NOHANDLE if not
    ! allocated.
    INTEGER                    :: h_Ddata
    
    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE. The subvectors are all of the same type.
    INTEGER                    :: cdataType = ST_DOUBLE

    ! Is set to true, if the handle h_Ddata belongs to another vector,
    ! i.e. when this vector shares data with another vector.
    LOGICAL                    :: bisCopy   = .FALSE.

    ! Number of blocks allocated in RvectorBlock
    INTEGER                    :: nblocks = 0
    
    ! A 1D array with scalar vectors for all the blocks.
    ! The handle identifier inside of these blocks are set to h_Ddata.
    ! The iidxFirstEntry index pointer in each subblock is set according
    ! to the position inside of the large vector.
    TYPE(t_vectorScalar), DIMENSION(LSYSBL_MAXBLOCKS) :: RvectorBlock
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A block matrix that can be used by scalar linear algebra routines.
  
  TYPE t_matrixBlock
    
    ! Total number of equations in the vector
    INTEGER                    :: NEQ         = 0
    
    ! Number of diagonal blocks in the matrix
    INTEGER                    :: ndiagBlocks = 0

    ! Matrix specification tag. This is a bitfield coming from an OR 
    ! combination of different LSYSBS_MSPEC_xxxx constants and specifies
    ! various details of the matrix. If it is =LSYSBS_MSPEC_GENERAL,
    ! the matrix is a usual matrix that needs no special handling.
    INTEGER(I32) :: imatrixSpec = LSYSBS_MSPEC_GENERAL

    ! A 2D array with scalar matrices for all the blocks
    TYPE(t_matrixScalar), &
      DIMENSION(LSYSBL_MAXBLOCKS,LSYSBL_MAXBLOCKS) :: RmatrixBlock
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! This structure encapsules a pointer to a block matrix.
  
  TYPE t_matrixBlockPointer
    
    ! A pointer to a block matrix
    TYPE(t_matrixBlock), POINTER        :: ptr
    
  END TYPE
  
!</typeblock>

!</types>

  INTERFACE lsysbl_createVectorBlock
    MODULE PROCEDURE lsysbl_createVecBlockDirect 
    MODULE PROCEDURE lsysbl_createVecBlockIndirect 
  END INTERFACE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockDirect (rx, Isize, bclear)
  
!<description>
  ! Initialises the vector block structure rx. Isize is an array
  ! of integers containing the length of the individual blocks.
  ! Memory is allocated on the heap to hold vectors of the size
  ! according to Isize.
!</description>

!<input>
  
  ! An array with length-tags for the different blocks
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Isize
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  LOGICAL, INTENT(IN), OPTIONAL             :: bclear
  
!</input>

!<output>
  
  ! Destination structure. Memory is allocated for each of the blocks.
  TYPE(t_vectorBlock),INTENT(OUT) :: rx
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i,n
  
  ! Allocate one large vector holding all data.
  CALL storage_new1D ('lsysbl_createVecBlockDirect', 'Vector', SUM(Isize), ST_DOUBLE, &
                      rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  
  ! Initialise the sub-blocks. Save a pointer to the starting address of
  ! each sub-block.
  ! Denote in the subvector that the handle belongs to us - not to
  ! the subvector.
  
  n=1
  DO i = 1,SIZE(Isize)
    IF (Isize(i) .GT. 0) THEN
      rx%RvectorBlock(i)%NEQ = Isize(i)
      rx%RvectorBlock(i)%iidxFirstEntry = n
      rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
      rx%RvectorBlock(i)%cdataType = rx%cdataType
      rx%RvectorBlock(i)%bisCopy = .TRUE.
      n = n+Isize(i)
    ELSE
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
      rx%RvectorBlock(i)%h_Ddata = ST_NOHANDLE
    END IF
  END DO
  
  rx%NEQ = n
  rx%nblocks = SIZE(Isize)
  
  ! The data of the vector belongs to us (we created the handle), 
  ! not to somebody else.
  rx%bisCopy = .FALSE.
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(bclear)) THEN
    IF (bclear) THEN
      CALL lsysbl_vectorClear (rx)
    END IF
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockIndirect (rtemplate,rx,bclear)
  
!<description>
  ! Initialises the vector block structure rx. rtemplate is an
  ! existing block vector structure.
  ! Memory is allocated on the heap for rx, the different components
  ! of rx will have the same size as the corresponding ones in rtemplate.
!</description>
  
!<input>
  
  ! A template vector structure
  TYPE(t_vectorBlock),INTENT(IN) :: rtemplate
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  LOGICAL, INTENT(IN), OPTIONAL             :: bclear
  
!</input>

!<output>
  
  ! Destination structure. Memory is allocated for each of the blocks.
  TYPE(t_vectorBlock),INTENT(OUT) :: rx
  
!</output>
  
!</subroutine>

  ! Copy the vector structure with all pointers.
  ! The handle identifying the vector data on the heap will be overwritten later.
  rx = rtemplate
  
  ! Allocate one large new vector holding all data.
  CALL storage_new1D ('lsysbl_createVecBlockDirect', 'Vector', rtemplate%NEQ, ST_DOUBLE, &
                      rx%h_Ddata, ST_NEWBLOCK_NOINIT)
                      
  ! Put the new handle to all subvectors
  rx%RvectorBlock(:)%h_Ddata = rx%h_Ddata
  
  ! Note in the subvectors, that the handle belongs to somewhere else... to us!
  rx%RvectorBlock(:)%bisCopy = .TRUE.
  
  ! Our handle belongs to us, so it's not a copy of another vector.
  rx%bisCopy = .FALSE.

  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(bclear)) THEN
    IF (bclear) THEN
      CALL lsysbl_vectorClear (rx)
    END IF
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockIndMat (rtemplateMat,rx, bclear)
  
!<description>
  ! Initialises the vector block structure rx. rtemplateMat is an
  ! existing block matrix structure. The vector rx will be created
  ! according to the size of the blocks of the submatrices.
  ! Memory is allocated on the heap for rx, the different components
  ! of rx will have the same size as the corresponding diagonal
  ! blocks in rtemplateMat.
!</description>
  
!<input>
  ! A template vector structure
  TYPE(t_matrixBlock), INTENT(IN) :: rtemplateMat
  
  ! OPTIONAL: If set to TRUE, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  LOGICAL, INTENT(IN), OPTIONAL   :: bclear
!</input>

!<output>
  ! Destination structure. Memory is allocated for each of the blocks.
  TYPE(t_vectorBlock),INTENT(OUT) :: rx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i,n
  
  ! Allocate one large vector holding all data.
  CALL storage_new1D ('lsysbl_createVecBlockDirect', 'Vector', rtemplateMat%NEQ, &
                      ST_DOUBLE, rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  
  ! Initialise the sub-blocks. Save a pointer to the starting address of
  ! each sub-block.
  
  n=1
  DO i = 1,rtemplateMat%ndiagBlocks
    IF (rtemplateMat%RmatrixBlock(i,i)%NEQ .GT. 0) THEN
      rx%RvectorBlock(i)%NEQ = rtemplateMat%RmatrixBlock(i,i)%NEQ
      
      ! Take the handle of the complete-solution vector, but set the index of
      ! the first entry to a value >= 1 - so that it points to the first
      ! entry in the global solution vector!
      rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
      rx%RvectorBlock(i)%iidxFirstEntry = n
      
      ! Give the vector the dicretisation of the matrix
      rx%RvectorBlock(i)%p_rspatialDiscretisation => &
        rtemplateMat%RmatrixBlock(i,i)%p_rspatialDiscretisation
      
      ! Denote in the subvector that the handle belongs to us - not to
      ! the subvector.
      rx%RvectorBlock(i)%bisCopy = .TRUE.
      
      n = n+rtemplateMat%RmatrixBlock(i,i)%NEQ
    ELSE
      ! Let's hope this situation (an empty equation) never occurs - 
      ! might produce some errors elsewhere :)
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
    END IF
  END DO
  
  rx%NEQ = rtemplateMat%NEQ
  rx%nblocks = rtemplateMat%ndiagBlocks
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(bclear)) THEN
    IF (bclear) THEN
      CALL lsysbl_vectorClear (rx)
    END IF
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_assignDiscretIndirect (rtemplate,rx)
  
!<description>
  ! Assigns discretisation-related information (spatial discretisation, 
  ! boundary conditions,...) to rx. The vector rtemplate is used as a
  ! template, so at the end of the routine, rx and rtemplate share the
  ! same discretisation.
!</description>
  
!<input>
  ! A template vector structure
  TYPE(t_vectorBlock),INTENT(IN) :: rtemplate
!</input>

!<inputoutput>
  ! Destination structure. Discretisation-related information of rtemplate
  ! is assigned to rx.
  TYPE(t_vectorBlock),INTENT(INOUT) :: rx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i

  ! Simply modify all pointers of all subvectors, that's it.
  DO i=1,rx%nblocks
    ! Spatial discretisation
    rx%RvectorBlock(i)%p_rspatialDiscretisation => &
      rtemplate%RvectorBlock(i)%p_rspatialDiscretisation 

    ! Discrete boundary conditions
    rx%RvectorBlock(i)%p_rdiscreteBC => &
      rtemplate%RvectorBlock(i)%p_rdiscreteBC

    ! Discrete fictitious boundary conditions
    rx%RvectorBlock(i)%p_rdiscreteBCfict => &
      rtemplate%RvectorBlock(i)%p_rdiscreteBCfict
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_releaseVector (rx)
  
!<description>
  ! Releases the memory that is reserved by rx.
!</description>
  
!<inputoutput>
  
  ! Vector structure that is to be released.
  TYPE(t_vectorBlock),INTENT(INOUT) :: rx
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! Release the data - if the handle is not a copy of another vector!
  IF (.NOT. rx%bisCopy) THEN
    CALL storage_free (rx%h_Ddata)
  END IF
  
  ! Clean up the structure
  DO i = 1,SIZE(rx%RvectorBlock)
    rx%RvectorBlock%NEQ = 0
    rx%RvectorBlock%h_Ddata = ST_NOHANDLE
    rx%RvectorBlock%iidxFirstEntry = 0
    rx%bisCopy = .FALSE.
  END DO
  rx%NEQ = 0
  rx%nblocks = 0
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE lsysbl_releaseMatrix (rmatrix)
 
!<description>
  ! Releases the memory that is reserved by rmatrix.
!</description>
  
!<inputoutput>
  ! Block matrix to be released
  TYPE(t_matrixBlock), INTENT(INOUT)                :: rMatrix
!</inputoutput>

!</subroutine>
  
  ! local variables
  INTEGER :: x,y
  
  ! loop through all the submatrices and release them
  DO x=1,rmatrix%ndiagBlocks
    DO y=1,rmatrix%ndiagBlocks
      
      ! Only release the matrix if there is one.
      IF (rmatrix%RmatrixBlock(y,x)%NA .NE. 0) THEN
        CALL lsyssc_releaseMatrix (rmatrix%RmatrixBlock(y,x))
      END IF
      
    END DO
  END DO
  
  ! Clean up the other variables, finish
  rmatrix%NEQ = 0
  rmatrix%ndiagBlocks = 0
   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE lsysbl_blockMatVec (rMatrix, rx, ry, cx, cy)
 
!<description>
  ! Performs a matrix vector multiplicationwith a given scalar matrix:
  !    $$ Dy   =   cx * rMatrix * rx   +   cy * ry $$
!</description>
  
!<input>
  
  ! Block matrix
  TYPE(t_matrixBlock), INTENT(IN)                   :: rMatrix

  ! Block vector to multiply with the matrix.
  TYPE(t_vectorBlock), INTENT(IN)                   :: rx
  
  ! Multiplicative factor for rx
  REAL(DP), INTENT(IN)                              :: cx

  ! Multiplicative factor for ry
  REAL(DP), INTENT(IN)                              :: cy
  
!</input>

!<inputoutput>
  
  ! Additive vector. Receives the result of the matrix-vector multiplication
  TYPE(t_vectorBlock), INTENT(INOUT)                :: ry
  
!</inputoutput>

!</subroutine>
  
  ! local variables
  INTEGER :: x,y,mvok
  REAL(DP)     :: cyact
  
  ! loop through all the sub matrices and multiply.
  DO y=1,rMatrix%ndiagBlocks
  
    ! The original vector ry must only be respected once in the first
    ! matrix-vector multiplication. The additive contributions of
    ! the other blocks must simply be added to ry without ry being multiplied
    ! by cy again!
  
    cyact = cy
    MVOK = NO
    DO x=1,rMatrix%ndiagBlocks
      
      ! Only call the MV when there is a scalar matrix that we can use!
      IF (rMatrix%RmatrixBlock(y,x)%NA .NE. 0) THEN
        CALL lsyssc_scalarMatVec (rMatrix%RmatrixBlock(y,x), rx%RvectorBlock(x), &
                                  ry%RvectorBlock(x), cx, cyact)
        cyact = 1.0_DP
        mvok = YES
      END IF
      
    END DO
    
    ! If mvok=NO here, there was no matrix in the current line!
    ! A very special case, but nevertheless may happen. In this case,
    ! simply scale the vector ry by cyact!
    
    IF (mvok .EQ.NO) THEN
      CALL lsysbl_vectorScale (ry,cy)
    END IF
    
  END DO
   
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_vectorCopy (rx,ry)
  
!<description>
  ! Copies vector dx: Dy = Dx
!<input>
  
  ! Source vector
  TYPE(t_vectorBlock),INTENT(IN) :: rx
  
!</input>

!<inputoutput>
  
  ! Destination vector
  TYPE(t_vectorBlock),INTENT(INOUT) :: ry
  
  !</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource, p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource, p_Sdest
  
  ! Take care of the data type
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointers and copy the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
    CALL storage_getbase_double(ry%h_Ddata,p_Ddest)
    
    CALL lalg_vectorCopyDBle (p_Dsource,p_Ddest)
  
  CASE (ST_SINGLE)
    ! Get the pointers and copy the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Ssource)
    CALL storage_getbase_double(ry%h_Ddata,p_Sdest)
    
    CALL lalg_vectorCopySngl (p_Ssource,p_Sdest)

  CASE DEFAULT
    PRINT *,'lsysbl_vectorCopy: unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_vectorScale (rx,c)
  
!<description>
  ! Scales a vector vector rx: rx = c * rx
!</description>
  
!<inputoutput>
  
  ! Source and destination vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rx
  
!</inputoutput>

!<input>
  
  ! Multiplication factor
  REAL(DP), INTENT(IN) :: c
  
!</input>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Sdata
  
  ! Taje care of the data type!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Ddata)
    CALL lalg_vectorScaleDble (p_Ddata,c)  

  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_single(rx%h_Ddata,p_Sdata)
    CALL lalg_vectorScaleSngl (p_Sdata,REAL(c,SP))  

  CASE DEFAULT
    PRINT *,'lsysbl_vectorScale: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_vectorClear (rx)
  
!<description>
  ! Clears the block vector dx: Dx = 0
!<inputoutput>
  
  ! Destination vector to be cleared
  TYPE(t_vectorBlock), INTENT(INOUT) :: rx
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource
  
  ! Take care of the data type
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
    CALL lalg_vectorClearDble (p_Dsource)
  
  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_single(rx%h_Ddata,p_Ssource)
    CALL lalg_vectorClearSngl (p_Ssource)

  CASE DEFAULT
    PRINT *,'lsysbl_vectorClear: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_vectorLinearComb (rx,ry,cx,cy)
  
!<description>
  ! Performs a linear combination: ry = cx * rx  +  cy * ry
!</description>  
  
!<input>
  
  ! First source vector
  TYPE(t_vectorBlock), INTENT(IN)    :: rx
  
  ! Scaling factor for Dx
  REAL(DP), INTENT(IN)               :: cx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)               :: cy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  TYPE(t_vectorBlock), INTENT(INOUT) :: ry
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource, p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource, p_Sdest
  
  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsysbl_vectorLinearComb: different data types not supported!'
    STOP
  END IF
  
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointers and copy the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
    CALL storage_getbase_double(ry%h_Ddata,p_Ddest)
    
    CALL lalg_vectorLinearCombDble (p_Dsource,p_Ddest,cx,cy)

  CASE (ST_SINGLE)
    ! Get the pointers and copy the whole data array.
    CALL storage_getbase_single(rx%h_Ddata,p_Ssource)
    CALL storage_getbase_single(ry%h_Ddata,p_Sdest)
    
    CALL lalg_vectorLinearCombSngl (p_Ssource,p_Sdest,cx,cy)
  
  CASE DEFAULT
    PRINT *,'lsysbl_vectorLinearComb: Unsupported data type!'
    STOP
  END SELECT

  END SUBROUTINE
  
  !****************************************************************************
!<subroutine>
  
  REAL(DP) FUNCTION lsysbl_scalarProduct (rx, ry)
  
!<description>
  ! Calculates a scalar product of two block vectors.
!</description>
  
!<input>
  ! First vector
  TYPE(t_vectorBlock), INTENT(IN)                  :: rx

  ! Second vector
  TYPE(t_vectorBlock), INTENT(IN)                  :: ry

!</input>

!<result>
  ! The scalar product <rx,ry> of the two block vectors.
!</result>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: h_Ddata1dp
  REAL(DP), DIMENSION(:), POINTER :: h_Ddata2dp
  INTEGER(PREC_VECIDX) :: i
  REAL(DP) :: res
  
  ! Is there data at all?
  res = 0.0_DP
  
  IF ( (rx%NEQ .EQ. 0) .OR. (ry%NEQ .EQ. 0) .OR. (rx%NEQ .NE. rx%NEQ)) THEN
    PRINT *,'Error in lsyssc_scalarProduct: Vector dimensions wrong!'
    STOP
  END IF
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    
    CALL storage_getbase_double (rx%h_Ddata,h_Ddata1dp)
    
    SELECT CASE (ry%cdataType)
    CASE (ST_DOUBLE)
      ! Get the data arrays
      CALL storage_getbase_double (ry%h_Ddata,h_Ddata2dp)
      
      ! Perform the scalar product
      res = 0.0_DP
      DO i=1,rx%NEQ
        res = res + h_Ddata1dp(i)*h_Ddata2dp(i)
      END DO
      
    CASE DEFAULT
      PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
      STOP
    END SELECT
    
  CASE DEFAULT
    PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
    STOP
  END SELECT
  
  ! Return the scalar product, finish
  lsysbl_scalarProduct = res

  END FUNCTION

  !****************************************************************************
!<subroutine>
  
  REAL(DP) FUNCTION lsysbl_vectorNorm (rx,cnorm,iposMax)
  
!<description>
  ! Calculates the norm of a vector. cnorm specifies the type of norm to
  ! calculate.
!</description>
  
!<input>
  ! Vector to calculate the norm of.
  TYPE(t_vectorBlock), INTENT(IN)                  :: rx

  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  INTEGER, INTENT(IN) :: cnorm
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  INTEGER(I32), INTENT(OUT), OPTIONAL :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0, if an error occurred (unknown norm).
!</result>

!<result>
  ! The scalar product <rx,ry> of the two block vectors.
!</result>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Sdata

  ! Is there data at all?
  IF (rx%h_Ddata .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in lsysbl_vectorNorm: Vector empty!'
    STOP
  END IF
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the array and calculate the norm
    CALL storage_getbase_double (rx%h_Ddata,p_Ddata)
    lsysbl_vectorNorm = lalg_normDble (p_Ddata,cnorm,iposMax) 
    
  CASE (ST_SINGLE)
    ! Get the array and calculate the norm
    CALL storage_getbase_single (rx%h_Ddata,p_Sdata)
    lsysbl_vectorNorm = lalg_normSngl (p_Sdata,cnorm,iposMax) 
    
  CASE DEFAULT
    PRINT *,'lsysbl_vectorNorm: Unsupported data type!'
    STOP
  END SELECT
  
  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_createMatFromScalar (rscalarMat,rmatrix)
  
!<description>
  ! This routine creates a 1x1 block matrix rmatrix from a scalar matrix
  ! rscalarMat. Both, rscalarMat and rmatrix will share the same handles,
  ! so changing the content of rmatrix will change rscalarMat, too.
  ! Therefore, the imatrixSpec flag of the submatrix in the block matrix
  ! will be set to LSYSSC_MSPEC_ISCOPY to indicate that the matrix is a
  ! copy of another one.
!</description>
  
!<input>
  ! The scalar matrix which should provide the data
  TYPE(t_matrixScalar), INTENT(IN) :: rscalarMat
!</input>

!<output>
  ! The 1x1 block matrix, created from rscalarMat.
  TYPE(t_matrixBlock), INTENT(OUT) :: rmatrix
!</output>
  
!</subroutine>

    ! Fill the rmatrix structure with data.
    rmatrix%NEQ         = rscalarMat%NEQ
    rmatrix%ndiagBlocks = 1
    rmatrix%imatrixSpec = LSYSBS_MSPEC_SCALAR
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block matrix
    rmatrix%RmatrixBlock(1,1)             = rscalarMat

    ! The matrix is a copy of another one. Note this!
    rmatrix%RmatrixBlock(1,1)%imatrixSpec = LSYSSC_MSPEC_ISCOPY
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_createVecFromScalar (rscalarVec,rvector)
  
!<description>
  ! This routine creates a 1 block vector rvector from a scalar vector
  ! rscalarVec. Both, rscalarVec and rvector will share the same handles,
  ! so changing the content of rvector will change rscalarVec, too.
  ! Therefore, the bisCopy flag of the subvector in the block vector
  ! will be set to TRUE.
!</description>
  
!<input>
  ! The scalar vector which should provide the data
  TYPE(t_vectorScalar), INTENT(IN) :: rscalarVec
!</input>

!<output>
  ! The 1x1 block matrix, created from rscalarMat.
  TYPE(t_vectorBlock), INTENT(OUT) :: rvector
!</output>
  
!</subroutine>

    ! Fill the rvector structure with data.
    rvector%NEQ         = rscalarVec%NEQ
    rvector%cdataType   = rscalarVec%cdataType
    rvector%h_Ddata     = rscalarVec%h_Ddata
    rvector%nblocks     = 1
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block vector
    rvector%RvectorBlock(1)             = rscalarVec

    ! The handle belongs to another vector - note this in the structure.
    rvector%RvectorBlock(1)%bisCopy     = .TRUE.
    
    ! The data of the vector actually belongs to another one - to the
    ! scalar one. Note this in the structure.
    rvector%bisCopy = .TRUE.

  END SUBROUTINE

END MODULE
