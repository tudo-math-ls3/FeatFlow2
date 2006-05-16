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
!# 1.) lsysbl_createVecBlockDirect
!#     -> Create a block vector by specifying the size of the subblocks
!# 
!# 2.) lsysbl_createVecBlockInirect
!#     -> Create a block vector by copying the structure of another block
!#        vector
!#
!# 3.) lsysbl_releaseVectorBlock
!#     -> Release a block vector from memory
!#
!# 4.) lsysbl_blockMatVec
!#     -> Multiply a block matrix with a block vector
!#
!# 5.) lsysbl_blockCopy
!#     -> Copy a block vector over to another one
!#
!# 6.) lsysbl_blockScale
!#     -> Scale a block vector by a constant
!#
!# 7.) lsysbl_blockClear
!#     -> Clear a block vector
!#
!# 8.) lsysbl_blockLinearComb
!#     -> Linear combination of two block vectors
!# </purpose>
!##############################################################################

MODULE linearsystemblock

  USE fsystem
  USE storage
  USE linearsystemscalar
  USE linearalgebra

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

  SUBROUTINE lsysbl_createVecBlockDirect (rx, Isize, iclear)
  
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
  INTEGER, INTENT(IN), OPTIONAL             :: iclear
  
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
  
  n=1
  DO i = 1,SIZE(Isize)
    IF (Isize(i) .GT. 0) THEN
      rx%RvectorBlock(i)%NEQ = Isize(i)
      rx%RvectorBlock(i)%iidxFirstEntry = n
      rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
      rx%RvectorBlock(i)%cdataType = rx%cdataType
      n = n+Isize(i)
    ELSE
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
      rx%RvectorBlock(i)%h_Ddata = ST_NOHANDLE
    END IF
  END DO
  
  rx%NEQ = n
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(iclear)) THEN
    IF (iclear .EQ. YES) THEN
      CALL lsysbl_blockClear (rx)
    END IF
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockIndirect (rx, rtemplate,iclear)
  
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
  INTEGER, INTENT(IN), OPTIONAL             :: iclear
  
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

  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(iclear)) THEN
    IF (iclear .EQ. YES) THEN
      CALL lsysbl_blockClear (rx)
    END IF
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockIndMat (rx, rtemplateMat,bclear)
  
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
  DO i = 1,rtemplateMat%NEQ
    IF (rtemplateMat%RmatrixBlock(i,i)%NEQ .GT. 0) THEN
      rx%RvectorBlock(i)%NEQ = rtemplateMat%RmatrixBlock(i,i)%NEQ
      rx%RvectorBlock(i)%iidxFirstEntry = n
      n = n+rtemplateMat%RmatrixBlock(i,i)%NEQ
    ELSE
      ! Let's hope this situation (an empty equation) never occurs - 
      ! might produce some errors elsewhere :)
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
    END IF
  END DO
  
  rx%NEQ = rtemplateMat%NEQ
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(bclear)) THEN
    IF (bclear) THEN
      CALL lsysbl_blockClear (rx)
    END IF
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_releaseVectorBlock (rx)
  
!<description>
  ! Releases the memory that is reserved by rx.
!</description>
  
!<inputoutput>
  
  ! Vector structure that is to be released.
  TYPE(t_vectorBlock),INTENT(OUT) :: rx
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! Release the data
  CALL storage_free (rx%RvectorBlock(i)%h_Ddata)
  
  ! Clean up the structure
  DO i = 1,SIZE(rx%RvectorBlock)
    rx%RvectorBlock%NEQ = 0
    rx%RvectorBlock%h_Ddata = ST_NOHANDLE
    rx%RvectorBlock%iidxFirstEntry = 0
  END DO
  rx%NEQ = 0
  
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
  TYPE(t_vectorBlock), INTENT(OUT)                 :: ry
  
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
      CALL lsysbl_blockScale (ry,cy)
    END IF
    
  END DO
   
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_blockCopy (rx,ry)
  
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
    PRINT *,'lsysbl_blockCopy: unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_blockScale (rx,c)
  
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
    PRINT *,'lsysbl_blockScale: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_blockClear (rx)
  
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
    PRINT *,'lsysbl_blockClear: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_blockLinearComb (rx,ry,cx,cy)
  
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
    PRINT *,'lsysbl_blockLinearComb: different data types not supported!'
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
    PRINT *,'lsysbl_blockLinearComb: Unsupported data type!'
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

END MODULE
