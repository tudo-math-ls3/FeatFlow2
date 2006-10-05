!#########################################################################
!# ***********************************************************************
!# <name> vectorio </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines and constant necessary to output 
!# vectors to files.
!#
!# The following routines can be found in this module:
!#
!# 1.) vecio_writeBlockVectorHR
!#     -> Writes a vector in human readable form into a text file
!#
!# 2.) matio_writeArray_Dble
!#     -> Writes an array in human readable form into a text file
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
  SUBROUTINE vecio_writeBlockVectorHR (rvector, sarray, bunsort,&
                                       ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a block vector into a text file.
    ! The vector is written in human readable form.
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
    
    ! Format string to use for the output; e.g. '(E20.10)'
    CHARACTER(len=*), INTENT(IN) :: sformat
    
    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it's unsorted on the fly.
    ! =FALSE: Write vector as it is.
    LOGICAL, INTENT(IN) :: bunsort
  !</input>
    
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  !INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Ipermutation

  ! Vector precision?
  SELECT CASE (rvector%cdataType)
  CASE (ST_DOUBLE)
    ! Get the data arrays and write the matrix
    CALL lsysbl_getbase_double (rvector,p_Ddata)
    
    ! Vector sorted and we should unsort?
    IF (lsysbl_isVectorSorted (rvector) .AND. bunsort) THEN
      PRINT *,'Can''t write sorted block vectors at the moment. Sorry!'
      ! There must be a loop, every subvector must be unsorted and
      ! written out separately into the same file. Let's implement later...
      STOP
    ELSE
      CALL vecio_writeArray_Dble (p_Ddata, sarray, &
                                   ifile, sfile, sformat)
    END IF
  CASE DEFAULT
    PRINT *,'vecio_writeFullMatrix: Unsupported vector precision.'
    STOP
  END SELECT
    
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE vecio_writeArray_Dble (Ddata, sarray, &
                                    ifile, sfile, sformat, Ipermutation)
  
  !<description>
    ! Write double precision vector into a text file.
  !</description>
    
  !<input>
    ! vector: array [:] of double
    REAL(DP), DIMENSION(:), INTENT(IN) :: Ddata
    
    ! name of the vector
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    CHARACTER(len=*), INTENT(IN) :: sformat

    ! OPTIONAL: Permutation for unsorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    INTEGER(PREC_VECIDX), DIMENSION(:), OPTIONAL :: Ipermutation
  !</input>
    
!</subroutine>
    
    !local variables
    INTEGER :: i, cf, nchar
    REAL(DP) :: dval
    CHARACTER(len=128) :: S
    CHARACTER(len=6) :: sformatChar
    
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE)
      IF (cf .EQ. -1) THEN
        PRINT *, 'vecio_writeArray_Dble: Could not open file '// &
                 trim(sfile)
        STOP
      END IF
    ELSE
      cf = ifile
    END IF
    
    ! Get length of output strings
    S(:) = ' '
    WRITE (S,sformat) 0.0_DP
    nchar = LEN(trim(S))
    
    ! Build array format string
    sformatChar = '(A'//sys_i3(nchar)//')'
    
    ! Write all format strings into the file
    WRITE (cf,'(3A15,2I15)') sarray, sformat, sformatChar, &
                             nchar, SIZE(Ddata)
    
    IF (SIZE(Ddata) .LE. 0) RETURN

    ! Write the vector.
    ! Unsort the vector on the fly if necessary.
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
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE

END MODULE
