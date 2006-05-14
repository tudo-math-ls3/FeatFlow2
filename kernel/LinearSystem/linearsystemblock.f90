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
      n = n+Isize(i)
    ELSE
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
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
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  
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
      IF (ry%RvectorBlock(x)%h_Ddata .NE. ST_NOHANDLE) THEN
        CALL storage_getbase_double(ry%RvectorBlock(x)%h_Ddata,p_Ddata)
        CALL lalg_vectorScale (p_Ddata,cy)
      END IF
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
  
  ! Get the pointers and copy the whole data array.
  CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
  CALL storage_getbase_double(ry%h_Ddata,p_Ddest)
  
  CALL lalg_vectorCopy (p_Dsource,p_Ddest)
  
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
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource
  
  ! Get the pointer and scale the whole data array.
  CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
  CALL lalg_vectorScale (p_DSource,c)  
  
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
  
  ! Get the pointer and scale the whole data array.
  CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
  CALL lalg_vectorClear (p_Dsource)
  
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
  
  ! Get the pointers and copy the whole data array.
  CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
  CALL storage_getbase_double(ry%h_Ddata,p_Ddest)
  
  CALL lalg_vectorLinearComb (p_Dsource,p_Ddest,cx,cy)

  END SUBROUTINE
  
END MODULE
