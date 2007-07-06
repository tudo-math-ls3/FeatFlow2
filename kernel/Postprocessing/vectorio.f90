!#########################################################################
!# ***********************************************************************
!# <name> vectorio </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines and constant necessary to output 
!# vectors/arrays to files or read them from files.
!#
!# The following routines can be found in this module:
!#
!# 1.) vecio_writeArray_Dble
!#     -> Writes an array into a (text or binary) file
!#
!# 2.) vecio_readArray_Dble
!#     -> Reads an array from a (text or binary) file
!#
!# 3.) vecio_writeBlockVectorHR
!#     -> Writes a block vector into a (text or binary) file
!#
!# 4.) vecio_writeVectorHR
!#     -> Writes a scalar vector into a (text or binary) file
!#
!# 5.) vecio_readBlockVectorHR
!#     -> Reads a block vector from a (text or binary) file
!#
!# 6.) vecio_readVectorHR
!#     -> Reads a scalar vector from a (text or binary) file
!#
!# </purpose>
!#########################################################################

MODULE vectorio

  USE fsystem
  USE storage
  USE io
  USE linearsystemscalar
  USE linearsystemblock
  
  IMPLICIT NONE

  CONTAINS

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE vecio_writeArray_Dble (Ddata, ifile, sfile, sformat, Ipermutation)
  
  !<description>
    ! Write double precision vector into a text file.
    ! The array is written out 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>
    
  !<input>
    ! vector: array [:] of double
    REAL(DP), DIMENSION(:), INTENT(IN) :: Ddata
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: sformat

    ! OPTIONAL: Permutation for unsorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    INTEGER(PREC_VECIDX), DIMENSION(:), OPTIONAL :: Ipermutation
  !</input>
    
!</subroutine>
    
    !local variables
    INTEGER :: i, cf
    REAL(DP) :: dval
    
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE, bformatted=PRESENT(sformat))
      IF (cf .EQ. -1) THEN
        PRINT *, 'vecio_writeArray_Dble: Could not open file '// &
                 trim(sfile)
        CALL sys_halt()
      END IF
    ELSE
      cf = ifile
    END IF
    
    IF (SIZE(Ddata) .LE. 0) RETURN

    ! Write the vector.
    ! Unsort the vector on the fly if necessary.
    IF (PRESENT(sformat)) THEN
      IF (PRESENT(Ipermutation)) THEN
        DO i=1, SIZE(Ddata)
          dval = Ddata(Ipermutation(i))
          WRITE (cf,sformat) dval
        END DO
      ELSE
        DO i=1, SIZE(Ddata)
          dval = Ddata(i)
          WRITE (cf,sformat) dval
        END DO
      END IF
    ELSE
      IF (PRESENT(Ipermutation)) THEN
        DO i=1, SIZE(Ddata)
          dval = Ddata(Ipermutation(i))
          WRITE (cf) dval
        END DO
      ELSE
        DO i=1, SIZE(Ddata)
          dval = Ddata(i)
          WRITE (cf) dval
        END DO
      END IF
    END IF
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE vecio_readArray_Dble (Ddata, ifile, sfile, sformat, Ipermutation)
  
  !<description>
    ! Reads a double precision vector from a file.
    ! The array is read in 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>
    
  !<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! OPTIONAL: Format string to use for the input; e.g. '(E20.10)'.
    ! If not specified, data is read from the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    ! When reading an array written out by vecio_writeArray_Dble,
    ! the format string shall match the setting of the 
    ! format string used there.
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: sformat
    
    ! OPTIONAL: Permutation for sorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    INTEGER(PREC_VECIDX), DIMENSION(:), OPTIONAL :: Ipermutation
  !</input>
  
  !<output>
    ! Array where to write the data to.
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddata
  !</output>
    
!</subroutine>
    
    ! local variables
    INTEGER :: i, cf
    REAL(DP) :: dval
    
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE, bformatted=PRESENT(sformat))
      IF (cf .EQ. -1) THEN
        PRINT *, 'vecio_readArray_Dble: Could not open file '// &
                 trim(sfile)
        CALL sys_halt()
      END IF
    ELSE
      cf = ifile
    END IF
    
    ! Read the array.
    IF (PRESENT(sformat)) THEN
      ! Unsort the vector on the fly if necessary.
      IF (PRESENT(Ipermutation)) THEN
        DO i=1, SIZE(Ddata)
          READ (cf,sformat) dval
          Ddata(Ipermutation(i)) = dval
        END DO
      ELSE
        DO i=1, SIZE(Ddata)
          READ (cf,sformat) dval
          Ddata(i) = dval
        END DO
      END IF
    ELSE
      ! Unsort the vector on the fly if necessary.
      IF (PRESENT(Ipermutation)) THEN
        DO i=1, SIZE(Ddata)
          Ddata(Ipermutation(i)) = dval
        END DO
      ELSE
        DO i=1, SIZE(Ddata)
          READ (cf) dval
          Ddata(i) = dval
        END DO
      END IF
    END IF
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE vecio_writeVectorHR (rvector, sarray, bunsort,&
                                  ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a (scalar) vector into a (text or binary) file.
  !</description>
    
  !<input>
    ! The vector to be written out
    TYPE(t_vectorScalar), INTENT(IN) :: rvector
    
    ! Name of the vector
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it's unsorted on the fly.
    ! =FALSE: Write vector as it is.
    LOGICAL, INTENT(IN) :: bunsort

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: sformat
  !</input>
    
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Ipermutation
    INTEGER :: cf,nchar

    CHARACTER(len=128) :: S
    CHARACTER(len=15) :: sarrayname
    CHARACTER(len=6) :: sformatChar

    IF (rvector%NEQ .EQ. 0) RETURN ! nothing to do

    ! Open the file
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE,bformatted=PRESENT(sformat))
      IF (cf .EQ. -1) THEN
        PRINT *, 'vecio_writeVectorHR: Could not open file '// &
                 trim(sfile)
        CALL sys_halt()
      END IF
    ELSE
      cf = ifile
    END IF

    IF (PRESENT(sformat)) THEN
      ! Get length of output strings
      S(:) = ' '
      WRITE (S,sformat) 0.0_DP
      nchar = LEN(trim(S))
      
      ! Build array format string
      sformatChar = '(A'//TRIM(sys_i3(nchar))//')'
      
      ! Write all format strings into the file
      WRITE (cf,'(A,3A15,2I15)') '# ',sarray, sformat, sformatChar, &
                                  nchar, rvector%NEQ
    ELSE
      sarrayname = sarray
      WRITE (cf) sarrayname, rvector%NEQ
    END IF
    
    ! Vector precision?
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      ! Permuted?
      NULLIFY(p_Ipermutation)
      IF (bunsort .AND. (lsyssc_isVectorSorted (rvector))) THEN
        CALL storage_getbase_int (rvector%h_IsortPermutation,p_Ipermutation)
        ! We must use the inverse permutation
        p_Ipermutation => p_Ipermutation(rvector%NEQ+1:)
      END IF

      CALL lsyssc_getbase_double (rvector,p_Ddata)
      
      IF (PRESENT(sformat)) THEN
        IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
          CALL vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat)
        ELSE 
          CALL vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
        END IF
      ELSE
        IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
          CALL vecio_writeArray_Dble (p_Ddata, cf, sfile)
        ELSE 
          CALL vecio_writeArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
        END IF
      END IF
    CASE DEFAULT
      PRINT *,'vecio_writeVectorHR: Unsupported vector precision.'
      CALL sys_halt()
    END SELECT
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
    
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE vecio_readVectorHR (rvector, sarray, bunsort,&
                                 ifile, sfile, bformatted)
  
  !<description>
    ! This routine reads a (scalar) vector from a (text or binary) file.
    !
    ! Note: The input data in the file may be written out with vecio_writeVectorHR
    ! or vecio_writeBlockVectorHR; in the latter case, a preciously written
    ! out block vector is read in as scalar vector.
  !</description>
    
  !<input>
    ! Input channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Read from channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! Name of the file where to read from. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! Read unsorted vector.
    ! =TRUE:  Read data and sort it according to the sorting strategy
    !         in the vector (if the vector is sorted)
    ! =FALSE: Read vector as it is.
    LOGICAL, INTENT(IN) :: bunsort

    ! Whether to read data formatted or unformatted.
    ! TRUE  = Treat data in input file formatted, i.e. in human readable form.
    ! FALSE = Data in the input file is unformatted, i.e. in processor
    !         dependent form.
    ! A vector written out by vecio_writeVectorHR with a format specifier
    ! sformat being specified shall be read with bformatted=TRUE.
    ! A vector written out by vecio_writeVectorHR without a format specifier
    ! sformat being specified shall be read with bformatted=FALSE.
    LOGICAL, INTENT(IN) :: bformatted
  !</input>

  !<inputoutput>
    ! The vector to be read in.
    ! If the vector is not initialised, a new vector is automatically created
    ! with the correct size.
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector
  !</inputoutput>

  !<output>    
    ! Name of the vector
    CHARACTER(len=*), INTENT(OUT) :: sarray
  !</output>
    
    
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Ipermutation
    INTEGER :: cf,nchar
    INTEGER(PREC_VECIDX) :: NEQ

    CHARACTER(len=15) :: sarrayname,sformat
    CHARACTER(len=6) :: sformatChar,S

    IF (rvector%NEQ .EQ. 0) RETURN ! nothing to do

    ! Open the file
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForReading(sfile, cf, bformatted=bformatted)
      IF (cf .EQ. -1) THEN
        PRINT *, 'vecio_writeVectorHR: Could not open file '// &
                 trim(sfile)
        CALL sys_halt()
      END IF
    ELSE
      cf = ifile
    END IF

    IF (bformatted) THEN
      ! Get the format specification from the file
      READ (cf,'(A2,3A15,2I15)') S,sarrayname, sformat, sformatChar, &
                                nchar, NEQ
    ELSE
      ! Get array information from the file
      READ (cf) sarrayname, NEQ
    END IF

    sarray = sarrayname
    
    ! Does the vector exist? If not, we create a new one.
    IF (rvector%NEQ .EQ. 0) THEN
      CALL lsyssc_createVector (rvector,NEQ,.FALSE.,ST_DOUBLE)
    END IF
    
    IF (rvector%NEQ .NE. NEQ) THEN
      PRINT *,'vecio_readVectorHR: Vector has wrong size!'
      CALL sys_halt()
    END IF
    
    ! Vector precision?
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      ! Permuted?
      NULLIFY(p_Ipermutation)
      IF (bunsort .AND. (lsyssc_isVectorSorted (rvector))) THEN
        CALL storage_getbase_int (rvector%h_IsortPermutation,p_Ipermutation)
        ! We must use the inverse permutation
        p_Ipermutation => p_Ipermutation(NEQ+1:)
      END IF

      CALL lsyssc_getbase_double (rvector,p_Ddata)
      
      IF (bformatted) THEN
        IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
          CALL vecio_readArray_Dble (p_Ddata, cf, sfile, sformat)
        ELSE 
          CALL vecio_readArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
        END IF
      ELSE
        IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
          CALL vecio_readArray_Dble (p_Ddata, cf, sfile)
        ELSE 
          CALL vecio_readArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
        END IF
      END IF
    CASE DEFAULT
      PRINT *,'vecio_readVectorHR: Unsupported vector precision.'
      CALL sys_halt()
    END SELECT
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
    
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE vecio_writeBlockVectorHR (rvector, sarray, bunsort,&
                                       ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a block vector into a (text or binary) file.
    ! The output file can be read in with vecio_readBlockVectorHR as
    ! block vector or with vecio_readVectorHR as scalar vector.
  !</description>
    
  !<input>
    ! The vector to be written out
    TYPE(t_vectorBlock), INTENT(IN) :: rvector
    
    ! Name of the vector
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it's unsorted on the fly.
    ! =FALSE: Write vector as it is.
    LOGICAL, INTENT(IN) :: bunsort

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: sformat
  !</input>
    
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Ipermutation
    INTEGER :: cf,nchar,I

    CHARACTER(len=128) :: S
    CHARACTER(len=15) :: sarrayname
    CHARACTER(len=6) :: sformatChar

    IF (rvector%NEQ .EQ. 0) RETURN ! nothing to do

    ! Open the file
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE,bformatted=PRESENT(sformat))
      IF (cf .EQ. -1) THEN
        PRINT *, 'vecio_writeBlockVectorHR: Could not open file '// &
                 trim(sfile)
        CALL sys_halt()
      END IF
    ELSE
      cf = ifile
    END IF

    ! Write all format strings into the file
    IF (PRESENT(sformat)) THEN
      ! Get length of output strings
      S(:) = ' '
      WRITE (S,sformat) 0.0_DP
      nchar = LEN(trim(S))
      
      ! Build array format string
      sformatChar = '(A'//TRIM(sys_i3(nchar))//')'
      
      WRITE (cf,'(A,3A15,3I15)',ADVANCE="NO") '# ',sarray, sformat, sformatChar, &
                                               nchar, rvector%NEQ, rvector%nblocks
      ! Write block structure
      DO i=1,rvector%nblocks
        WRITE (cf,'(I15)',ADVANCE="NO") rvector%RvectorBlock(i)%NEQ
      END DO
      
      ! New line
      WRITE (cf,*)
    ELSE
      sarrayname = sarray
      WRITE (cf) sarrayname, rvector%NEQ,rvector%nblocks
      ! Write block structure
      WRITE (cf) (rvector%RvectorBlock(i)%NEQ,i=1,rvector%nblocks)
    END IF
    
    ! Vector precision?
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      DO i=1,rvector%nblocks
        ! Permuted?
        NULLIFY(p_Ipermutation)
        IF (bunsort .AND. (lsyssc_isVectorSorted (rvector%RvectorBlock(i)))) THEN
          CALL storage_getbase_int (&
              rvector%RvectorBlock(i)%h_IsortPermutation,p_Ipermutation)
          ! We must use the inverse permutation
          p_Ipermutation => p_Ipermutation(rvector%RvectorBlock(i)%NEQ+1:)
        END IF

        CALL lsyssc_getbase_double (rvector%RvectorBlock(i),p_Ddata)
        
        IF (PRESENT(sformat)) THEN
          IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
            CALL vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat)
          ELSE 
            CALL vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
          END IF
        ELSE
          IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
            CALL vecio_writeArray_Dble (p_Ddata, cf, sfile)
          ELSE 
            CALL vecio_writeArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
          END IF
        END IF
      END DO
    CASE DEFAULT
      PRINT *,'vecio_writeBlockVectorHR: Unsupported vector precision.'
      CALL sys_halt()
    END SELECT
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
    
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE vecio_readBlockVectorHR (rvector, sarray, bunsorted,&
                                      ifile, sfile, bformatted)
  
  !<description>
    ! This routine reads a block vector from a (text or binary) file.
    !
    ! rvector may not be initialised when calling this routine; in this
    ! case, if the source file specifies a block vector, a new block
    ! vector is created in the same structure as written out
    ! by vecio_writeBlockVectorHR.
    !
    ! If rvector is initialised, the data is read from the file according
    ! to the structure of rvector.
    !
    ! If rvector specifies a structure and bformatted=TRUE, this routine
    ! is compatible to the FEAT 1.0 format of formatted vectors!
  !</description>
    
  !<input>
    ! Input channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Read from channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! Name of the file where to read from. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! Read unsorted vector.
    ! =TRUE:  Read data and sort it according to the sorting strategy
    !         in the vector (if the vector is sorted)
    ! =FALSE: Read vector as it is.
    LOGICAL, INTENT(IN) :: bunsorted

    ! Whether to read data formatted or unformatted.
    ! TRUE  = Treat data in input file formatted, i.e. in human readable form.
    ! FALSE = Data in the input file is unformatted, i.e. in processor
    !         dependent form.
    ! A vector written out by vecio_writeBlockVectorHR with a format specifier
    ! sformat being specified shall be read with bformatted=TRUE.
    ! A vector written out by vecio_writeBlockVectorHR without a format specifier
    ! sformat being specified shall be read with bformatted=FALSE.
    LOGICAL, INTENT(IN) :: bformatted
  !</input>

  !<inputoutput>
    ! The vector to be read in.
    ! If the vector is not initialised, a new vector is automatically created
    ! with the correct size.
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
  !</inputoutput>

  !<output>    
    ! Name of the vector
    CHARACTER(len=*), INTENT(OUT) :: sarray
  !</output>
    
    
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Ipermutation
    INTEGER :: cf,nchar,nblocks,i
    INTEGER(PREC_VECIDX) :: NEQ
    INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: IblockSize

    CHARACTER(len=15) :: sarrayname,sformat
    CHARACTER(len=6) :: sformatChar,S

    IF (rvector%NEQ .EQ. 0) RETURN ! nothing to do

    ! Open the file
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForReading(sfile, cf, bformatted=bformatted)
      IF (cf .EQ. -1) THEN
        PRINT *, 'vecio_writeVectorHR: Could not open file '// &
                 trim(sfile)
        CALL sys_halt()
      END IF
    ELSE
      cf = ifile
    END IF

    IF (bformatted) THEN
      ! Get the format specification from the file
      READ (cf,'(A2,3A15,3I15)',ADVANCE="NO") S,sarrayname, sformat, sformatChar, &
                                             nchar, NEQ
      ! If rvector does not specify the structure, read the vector structure
      ! from the file.
      IF (rvector%NEQ .EQ. 0) THEN
        READ (cf,'(I15)',ADVANCE="NO") nblocks
        
        ! Get the block structure
        ALLOCATE(IblockSize(nblocks))
        DO i=1,nblocks
          READ (cf,'(I15)',ADVANCE="NO") IblockSize(i)
        END DO
        
        ! New line
        READ (cf,*)
      ELSE
        ! Take the structure from rvector
        nblocks = rvector%nblocks
        ALLOCATE(IblockSize(nblocks))
        IblockSize(:) = rvector%RvectorBlock(1:nblocks)%NEQ

        ! New line
        READ (cf,*)
      END IF
    ELSE
      ! Get array information from the file
      READ (cf) sarrayname, NEQ, nblocks
      
      ! Get the block structure
      ALLOCATE(IblockSize(nblocks))
      READ (cf) (IblockSize(i),i=1,nblocks)
    END IF

    sarray = sarrayname
    
    ! Does the vector exist? If not, we create a new one.
    IF (rvector%NEQ .EQ. 0) THEN
      CALL lsysbl_createVecBlockDirect (rvector,IblockSize,.FALSE.,ST_DOUBLE)
    END IF
    
    ! Size of vector must match! Size of subvectors not -- it's not a bug,
    ! it's a feature ;-)
    IF (rvector%NEQ .NE. NEQ) THEN
      PRINT *,'vecio_readBlockVectorHR: Vector has wrong size!'
      CALL sys_halt()
    END IF

    ! We don't need the block size anymore.    
    DEALLOCATE (IblockSize)
    
    ! Vector precision?
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      DO i=1,nblocks
        ! Permuted?
        NULLIFY(p_Ipermutation)
        IF (bunsorted .AND. (lsyssc_isVectorSorted (rvector%RvectorBlock(i)))) THEN
          CALL storage_getbase_int (rvector%RvectorBlock(i)%h_IsortPermutation,&
              p_Ipermutation)
          ! We must use the inverse permutation
          p_Ipermutation => p_Ipermutation(rvector%RvectorBlock(i)%NEQ+1:)
        END IF

        CALL lsyssc_getbase_double (rvector%RvectorBlock(i),p_Ddata)
        
        IF (bformatted) THEN
          IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
            CALL vecio_readArray_Dble (p_Ddata, cf, sfile, sformat)
          ELSE 
            CALL vecio_readArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
          END IF
        ELSE
          IF (.NOT. ASSOCIATED(p_Ipermutation)) THEN
            CALL vecio_readArray_Dble (p_Ddata, cf, sfile)
          ELSE 
            CALL vecio_readArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
          END IF
        END IF
      END DO
    CASE DEFAULT
      PRINT *,'vecio_readBlockVectorHR: Unsupported vector precision.'
      CALL sys_halt()
    END SELECT
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
    
  END SUBROUTINE 

END MODULE
