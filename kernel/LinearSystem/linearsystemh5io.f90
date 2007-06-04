!##############################################################################
!# ****************************************************************************
!# <name> linearsystemh5io </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines for handling the input/output
!# of matrices and vectors from/to the HDF subsystem.
!#
!# The following routines are available:
!#
!# 1.) lsysh5io_writeMatrix
!#     -> Writes a matrix into hdf5 file.
!#
!# 2.) lsysh5io_writeBlockMatrix
!#     -> Writes a block matrix into hdf5 file.
!#
!# 3.) lsysh5io_readMatrix
!#     -> Reads a matrix from hdf5 file.
!#
!# 4.) lsysh5io_readBlockMatrix
!#     -> Reads a block matrix from hdf5 file.
!#
!# 5.) lsysh5io_writeVector
!#     -> Writes a vector into hdf5 file.
!#
!# 6.) lsysh5io_writeBlockVector
!#     -> Write a block vector into hdf5 file
!#
!# 7.) lsysh5io_readVector
!#     -> Reads a vector from hdf5 file
!#
!# 8.) lsysh5io_readBlockVector
!#     -> Reads a block vector from hdf5 file.
!#
!# </purpose>
!#########################################################################

MODULE linearsystemh5io

  USE fsystem
  USE storage
  USE linearsystemscalar
  USE linearsystemblock
  USE h5lite
  USE hdf5

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: lsysh5io_writeMatrix
  PUBLIC :: lsysh5io_writeBlockMatrix
  PUBLIC :: lsysh5io_readMatrix
  PUBLIC :: lsysh5io_readBlockMatrix
  PUBLIC :: lsysh5io_writeVector
  PUBLIC :: lsysh5io_writeBlockVector
  PUBLIC :: lsysh5io_readVector
  PUBLIC :: lsysh5io_readBlockVector

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  INTERFACE lsysh5io_writeMatrix
    MODULE PROCEDURE lsysh5io_writeMatrix_Root
    MODULE PROCEDURE lsysh5io_writeMatrix_Location
  END INTERFACE

  INTERFACE lsysh5io_writeBlockMatrix
    MODULE PROCEDURE lsysh5io_writeBlockMatrix_Root
    MODULE PROCEDURE lsysh5io_writeBlockMatrix_Location
  END INTERFACE

  INTERFACE lsysh5io_readMatrix
    MODULE PROCEDURE lsysh5io_readMatrix_RootByID
    MODULE PROCEDURE lsysh5io_readMatrix_RootByName
    MODULE PROCEDURE lsysh5io_readMatrix_LocationByID
    MODULE PROCEDURE lsysh5io_readMatrix_LocationByName
  END INTERFACE

  INTERFACE lsysh5io_readBlockMatrix
    MODULE PROCEDURE lsysh5io_readBlockMatrix_RootByID
    MODULE PROCEDURE lsysh5io_readBlockMatrix_RootByName
    MODULE PROCEDURE lsysh5io_readBlockMatrix_LocationByID
    MODULE PROCEDURE lsysh5io_readBlockMatrix_LocationByName
  END INTERFACE

  INTERFACE lsysh5io_writeVector
    MODULE PROCEDURE lsysh5io_writeVector_Root
    MODULE PROCEDURE lsysh5io_writeVector_Location
  END INTERFACE

  INTERFACE lsysh5io_writeBlockVector
    MODULE PROCEDURE lsysh5io_writeBlockVector_Root
    MODULE PROCEDURE lsysh5io_writeBlockVector_Location
  END INTERFACE
  
  INTERFACE lsysh5io_readVector
    MODULE PROCEDURE lsysh5io_readVector_RootByID
    MODULE PROCEDURE lsysh5io_readVector_RootByName
    MODULE PROCEDURE lsysh5io_readVector_LocationByID
    MODULE PROCEDURE lsysh5io_readVector_LocationByName
  END INTERFACE

  INTERFACE lsysh5io_readBlockVector
    MODULE PROCEDURE lsysh5io_readBlockVector_RootByID
    MODULE PROCEDURE lsysh5io_readBlockVector_RootByName
    MODULE PROCEDURE lsysh5io_readBlockVector_LocationByID
    MODULE PROCEDURE lsysh5io_readBlockVector_LocationByName
  END INTERFACE

CONTAINS

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeMatrix_Root(rh5file,rmatrix,matrix_id,cmatrixName)

!<description>
    ! This subroutine writes the matrix to the root of the hdf5 file that
    ! is associated with the h5file descriptor. If the optional parameter
    ! cmatrixName is given, then this value is used to identify the
    ! matrix. Otherwise, t_matrixScalar.XXX is adopted, whereby XXX is the
    ! next free number that is not present in the file.
    ! The handle to the matrix object is returned as matrix_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)             :: rh5file

    ! The matrix to be written
    TYPE(t_matrixScalar), INTENT(IN)       :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_writeMatrix(rh5file%file_id,rmatrix,matrix_id,cmatrixName)
  END SUBROUTINE lsysh5io_writeMatrix_Root

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeMatrix_Location(loc_id,rmatrix,matrix_id,cmatrixName)

!<description>
    ! This subroutine writes the matrix to the location given by loc_id.
    ! If the optional parameter cmatrixName is given, then this value is 
    ! used to identify the matrix. Otherwise, t_matrixScalar.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the matrix object is returned as matrix_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)             :: loc_id

    ! The matrix to be written
    TYPE(t_matrixScalar), INTENT(IN)       :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    INTEGER,  DIMENSION(:), POINTER :: P_Idata

    ! Create new group for t_matrixScalar
    IF (PRESENT(cmatrixName)) THEN
      CALL h5lite_set_group(loc_id,cmatrixName,matrix_id)
    ELSE
      CALL h5lite_set_group(loc_id,"t_matrixScalar",matrix_id,.TRUE.)
    END IF

    ! Set attributes for scalar matrix
    CALL h5lite_set_attribute(matrix_id,"cmatrixFormat",          rmatrix%cmatrixFormat)
    CALL h5lite_set_attribute(matrix_id,"cinterleavematrixFormat",rmatrix%cinterleavematrixFormat)
    CALL h5lite_set_attribute(matrix_id,"imatrixSpec",            rmatrix%imatrixSpec)
    CALL h5lite_set_attribute(matrix_id,"NA",                     rmatrix%NA)
    CALL h5lite_set_attribute(matrix_id,"NEQ",                    rmatrix%NEQ)
    CALL h5lite_set_attribute(matrix_id,"NCOLS",                  rmatrix%NCOLS)
    CALL h5lite_set_attribute(matrix_id,"NVAR",                   rmatrix%NVAR)
    CALL h5lite_set_attribute(matrix_id,"dscaleFactor",           rmatrix%dscaleFactor)
    CALL h5lite_set_attribute(matrix_id,"isortStrategy",          rmatrix%isortStrategy)
    CALL h5lite_set_attribute(matrix_id,"cdataType",              rmatrix%cdataType)

    ! Set data for scalar matrix
    CALL h5lite_set_data(matrix_id,"ITags",rmatrix%ITags)
    CALL h5lite_set_data(matrix_id,"DTags",rmatrix%DTags)
    
    ! Write handle for sorting permutation
    IF (rmatrix%h_IsortPermutation /= ST_NOHANDLE) THEN
      CALL storage_getbase_int(rmatrix%h_IsortPermutation,p_Idata,2*rmatrix%NEQ)
      CALL h5lite_set_data(matrix_id,"h_IsortPermutation",p_Idata)
    END IF

    ! Write handle for matrix data
    IF (rmatrix%h_DA /= ST_NOHANDLE) THEN
      SELECT CASE(rmatrix%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double(rmatrix,p_Ddata)
        CALL h5lite_set_data(matrix_id,"h_DA",p_Ddata)
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single(rmatrix,p_Fdata)
        CALL h5lite_set_data(matrix_id,"h_DA",p_Fdata)
      CASE(ST_INT)
        CALL lsyssc_getbase_int(rmatrix,p_Idata)
        CALL h5lite_set_data(matrix_id,"h_DA",p_Idata)
      CASE DEFAULT
        PRINT *, "lsysh5io_writeMatrix_Location: Unsupported data format!"
        STOP
      END SELECT
    END IF

    ! Write handle for column structure
    IF (rmatrix%h_Kcol /= ST_NOHANDLE) THEN
      CALL lsyssc_getbase_Kcol(rmatrix,p_Idata)
      CALL h5lite_set_data(matrix_id,"h_Kcol",p_Idata)
    END IF

    ! Write handle for row structure
    IF (rmatrix%h_Kld /= ST_NOHANDLE) THEN
      CALL lsyssc_getbase_Kld(rmatrix,p_Idata)
      CALL h5lite_set_data(matrix_id,"h_Kld",p_Idata)
    END IF

    ! Write handle for diagonal structure
    IF (rmatrix%h_Kdiagonal /= ST_NOHANDLE) THEN
      CALL lsyssc_getbase_Kdiagonal(rmatrix,p_Idata)
      CALL h5lite_set_data(matrix_id,"h_Kdiagonal",p_Idata)
    END IF
  END SUBROUTINE lsysh5io_writeMatrix_Location

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeBlockMatrix_Root(rh5file,rmatrix,matrix_id,cmatrixName)

!<description>
    ! This subroutine writes the matrix to the root of the hdf5 file that
    ! is associated with the h5file descriptor. If the optional parameter
    ! cmatrixName is given, then this value is used to identify the
    ! matrix. Otherwise, t_matrixBlock.XXX is adopted, whereby XXX is the
    ! next free number that is not present in the file.
    ! The handle to the matrix object is returned as matrix_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)             :: rh5file

    ! The matrix to be written
    TYPE(t_matrixBlock), INTENT(IN)        :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_writeBlockMatrix(rh5file%file_id,rmatrix,matrix_id,cmatrixName)
  END SUBROUTINE lsysh5io_writeBlockMatrix_Root

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeBlockMatrix_Location(loc_id,rmatrix,matrix_id,cmatrixName)

!<description>
    ! This subroutine writes the matrix to the location given by loc_id.
    ! If the optional parameter cmatrixName is given, then this value is 
    ! used to identify the matrix. Otherwise, t_matrixBlock.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the matrix object is returned as matrix_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)             :: loc_id

    ! The matrix to be written
    TYPE(t_matrixBlock), INTENT(IN)        :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: matrixScalar_id
    INTEGER :: iblock,jblock

    ! Create new group for t_matrixScalar
    IF (PRESENT(cmatrixName)) THEN
      CALL h5lite_set_group(loc_id,cmatrixName,matrix_id)
    ELSE
      CALL h5lite_set_group(loc_id,"t_matrixBlock",matrix_id,.TRUE.)
    END IF

    ! Set attributes for block matrix
    CALL h5lite_set_attribute(matrix_id,"NEQ",        rmatrix%NEQ)
    CALL h5lite_set_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    CALL h5lite_set_attribute(matrix_id,"ndiagBlocks",rmatrix%ndiagBlocks)
    CALL h5lite_set_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    CALL h5lite_set_attribute(matrix_id,"imatrixSpec",rmatrix%imatrixSpec)
    
    ! Write each scalar submatrix
    IF (ASSOCIATED(rmatrix%RmatrixBlock)) THEN
      
      ! Loop over all blocks
      DO iblock=LBOUND(rmatrix%RmatrixBlock,1),UBOUND(rmatrix%RmatrixBlock,1)
        DO jblock=LBOUND(rmatrix%RmatrixBlock,2),UBOUND(rmatrix%RmatrixBlock,2)
          CALL lsysh5io_writeMatrix(matrix_id,rmatrix%RmatrixBlock(iblock,jblock),&
              matrixScalar_id,"RmatrixBlock("//TRIM(ADJUSTL(sys_si(iblock,6)))//","//&
              TRIM(ADJUSTL(sys_si(jblock,6)))//")")
        END DO
      END DO
      
    END IF
  END SUBROUTINE lsysh5io_writeBlockMatrix_Location

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysh5io_readMatrix_RootByID(rh5file,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)          :: rh5file

    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The matrix that should be read
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readMatrix(rh5file%file_id,matrix_id,rmatrix)
  END SUBROUTINE lsysh5io_readMatrix_RootByID

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readMatrix_RootByName(rh5file,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)          :: rh5file

    ! Name of the matrix
    CHARACTER(LEN=*), INTENT(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readMatrix(rh5file%file_id,cmatrixName,rmatrix)
  END SUBROUTINE lsysh5io_readMatrix_RootByName

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readMatrix_LocationByID(loc_id,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)          :: loc_id

    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    INTEGER,  DIMENSION(:), POINTER :: p_Idata
    INTEGER :: isize,sizeKld,sizeKdiagonal

    ! Get attributes
    CALL h5lite_get_attribute(matrix_id,"cmatrixFormat",          rmatrix%cmatrixFormat)
    CALL h5lite_get_attribute(matrix_id,"cinterleavematrixFormat",rmatrix%cinterleavematrixFormat)
    CALL h5lite_get_attribute(matrix_id,"imatrixSpec",            rmatrix%imatrixSpec)
    CALL h5lite_get_attribute(matrix_id,"NA",                     rmatrix%NA)
    CALL h5lite_get_attribute(matrix_id,"NEQ",                    rmatrix%NEQ)
    CALL h5lite_get_attribute(matrix_id,"NCOLS",                  rmatrix%NCOLS)
    CALL h5lite_get_attribute(matrix_id,"NVAR",                   rmatrix%NVAR)
    CALL h5lite_get_attribute(matrix_id,"dscaleFactor",           rmatrix%dscaleFactor)
    CALL h5lite_get_attribute(matrix_id,"isortStrategy",          rmatrix%isortStrategy)
    CALL h5lite_get_attribute(matrix_id,"cdataType",              rmatrix%cdataType)

    ! Get data
    CALL h5lite_get_data(matrix_id,"ITags",rmatrix%ITags)
    CALL h5lite_get_data(matrix_id,"DTags",rmatrix%DTags)

    ! Get handle for sorting permutation
    IF (h5lite_ispresent_data(matrix_id,"h_IsortPermutation")) THEN
      
      IF (rmatrix%h_IsortPermutation == ST_NOHANDLE) THEN
        ! If the permutation handle is not associated to some memory block, allocate new memory
        CALL storage_new1d ('lsysh5io_readMatrix_LocationByID','IsortPermutation',2*rmatrix%NEQ,&
            ST_INT,rmatrix%h_IsortPermutation,ST_NEWBLOCK_NOINIT)
      ELSE
        ! Otherwise, check that the provided memory is sufficient
        CALL storage_getsize(rmatrix%h_IsortPermutation,isize)
        IF (isize /= 2*rmatrix%NEQ) THEN
          PRINT *, "lsysh5io_readMatrix_LocationByID: Wrong size of permutation array!"
          STOP
        END IF
      END IF

      ! Get permutation handle
      CALL storage_getbase_int(rmatrix%h_IsortPermutation,p_Idata)
      CALL h5lite_get_data(matrix_id,"h_IsortPermutation",p_Idata)
    END IF

    ! Get handle for matrix data
    IF (h5lite_ispresent_data(matrix_id,"h_DA")) THEN

      IF (rmatrix%h_DA == ST_NOHANDLE) THEN
        ! If the data handle is not associated to some memory block, allocate new memory
        CALL storage_new1d ('lsysh5io_readMatrix_LocationByID','h_DA',rmatrix%NA,&
            rmatrix%cdataType,rmatrix%h_DA,ST_NEWBLOCK_NOINIT)
      ELSE
        ! Otherwise, check that the provided memory is sufficient
        CALL storage_getsize(rmatrix%h_DA,isize)
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIXUNDEFINED)
          IF (isize < rmatrix%NA) THEN
            PRINT *, "lsysh5io_readMatrix_LocationByID: Not enough memory!"
            STOP
          END IF

        CASE (LSYSSC_MATRIX1)
          IF (isize < rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR) THEN
            PRINT *, "lsysh5io_readMatrix_LocationByID: Not enough memory!"
            STOP
          END IF

        CASE (LSYSSC_MATRIXD)
          IF (isize < rmatrix%NA*rmatrix%NVAR) THEN
            PRINT *, "lsysh5io_readMatrix_LocationByID: Not enough memory!"
            STOP
          END IF

        CASE DEFAULT
          PRINT *, "lsysh5io_readMatrix_LocationByID: Unsupported interleave matrix format!"
          STOP
        END SELECT
      END IF

      ! Get data handle
      SELECT CASE(rmatrix%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double(rmatrix,p_Ddata)
        CALL h5lite_get_data(matrix_id,"h_DA",p_Ddata)
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single(rmatrix,p_Fdata)
        CALL h5lite_get_data(matrix_id,"h_DA",p_Fdata)
      CASE(ST_INT)
        CALL lsyssc_getbase_int(rmatrix,p_Idata)
        CALL h5lite_get_data(matrix_id,"h_DA",p_Idata)
      CASE DEFAULT
        PRINT *, "lsysh5io_readMatrix_Location: Unsupported data format!"
        STOP
      END SELECT
    END IF

    ! Get handle for column structure
    IF (h5lite_ispresent_data(matrix_id,"h_Kcol")) THEN

      IF (rmatrix%h_Kcol == ST_NOHANDLE) THEN
        ! If the column handle is not associated to some memory block, allocate new memory
        CALL storage_new1d ('lsysh5io_readMatrix_LocationByID','h_Kcol',rmatrix%NA,&
            ST_INT,rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
      ELSE
        ! Otherwise, check that the provided memory is sufficient
        CALL storage_getsize(rmatrix%h_Kcol,isize)
        IF (isize < rmatrix%NA) THEN
          PRINT *, "lsysh5io_readMatrix_Location: Not enough memory!"
          STOP
        END IF
      END IF

      ! Get column handle
      CALL lsyssc_getbase_Kcol(rmatrix,p_Idata)
      CALL h5lite_get_data(matrix_id,"h_Kcol",p_Idata)
    END IF

    ! Get handle for row structure
    IF (h5lite_ispresent_data(matrix_id,"h_Kld")) THEN
      
      ! Is the matrix transposed or not?
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1) THEN
        sizeKld=rmatrix%NCOLS+1
      ELSE
        sizeKld=rmatrix%NEQ+1
      END IF

      IF (rmatrix%h_Kld == ST_NOHANDLE) THEN
        ! If the column handle is not associated to some memory block, allocate new memory
        CALL storage_new1d ('lsysh5io_readMatrix_LocationByID','h_Kld',sizeKld,&
            ST_INT,rmatrix%h_Kld,ST_NEWBLOCK_NOINIT)
      ELSE
        ! Otherwise, check that the provided memory is sufficient
        CALL storage_getsize(rmatrix%h_Kld,isize)
        IF (isize < sizeKld) THEN
          PRINT *, "lsysh5io_readMatrix_Location: Not enough memory!"
          STOP
        END IF
      END IF

      ! Get row handle
      CALL lsyssc_getbase_Kld(rmatrix,p_Idata)
      CALL h5lite_get_data(matrix_id,"h_Kld",p_Idata)
    END IF

    ! Get handle for diagonal structure
    IF (h5lite_ispresent_data(matrix_id,"h_Kdiagonal")) THEN

      ! Is the matrix transposed or not?
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1) THEN
        sizeKdiagonal=rmatrix%NCOLS
      ELSE
        sizeKdiagonal=rmatrix%NEQ
      END IF

      IF (rmatrix%h_Kdiagonal == ST_NOHANDLE) THEN
        ! If the column handle is not associated to some memory block, allocate new memory
        CALL storage_new1d ('lsysh5io_readMatrix_LocationByID','h_Kdiagonal',sizeKdiagonal,&
            ST_INT,rmatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
      ELSE
        ! Otherwise, check that the provided memory is sufficient
        CALL storage_getsize(rmatrix%h_Kdiagonal,isize)
        IF (isize < sizeKdiagonal) THEN
          PRINT *, "lsysh5io_readMatrix_Location: Not enough memory!"
          STOP
        END IF
      END IF

      ! Get diagonal handle
      CALL lsyssc_getbase_Kdiagonal(rmatrix,p_Idata)
      CALL h5lite_get_data(matrix_id,"h_Kdiagonal",p_Idata)
    END IF
  END SUBROUTINE lsysh5io_readMatrix_LocationByID
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readMatrix_LocationByName(loc_id,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)          :: loc_id

    ! Name of the matrix
    CHARACTER(LEN=*), INTENT(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: matrix_id

    ! Determine the matrix_id and call the working routine
    CALL h5lite_get_group(loc_id,cmatrixName,matrix_id)
    IF (matrix_id < 0) THEN
      PRINT *, "lsysh5io_readMatrix_LocationByName: Unable to get matrix_id!"
      STOP
    END IF
    CALL lsysh5io_readMatrix(loc_id,matrix_id,rmatrix)
  END SUBROUTINE lsysh5io_readMatrix_LocationByName

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysh5io_readBlockMatrix_RootByID(rh5file,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)          :: rh5file

    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The matrix that should be read
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readBlockMatrix(rh5file%file_id,matrix_id,rmatrix)
  END SUBROUTINE lsysh5io_readBlockMatrix_RootByID

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readBlockMatrix_RootByName(rh5file,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)          :: rh5file

    ! Name of the matrix
    CHARACTER(LEN=*), INTENT(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readBlockMatrix(rh5file%file_id,cmatrixName,rmatrix)
  END SUBROUTINE lsysh5io_readBlockMatrix_RootByName

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readBlockMatrix_LocationByID(loc_id,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)          :: loc_id

    ! Handle to the matrix object
    INTEGER(HID_T), INTENT(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: iblock,jblock

    ! Get attributes for block matrix
    CALL h5lite_get_attribute(matrix_id,"NEQ",        rmatrix%NEQ)
    CALL h5lite_get_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    CALL h5lite_get_attribute(matrix_id,"ndiagBlocks",rmatrix%ndiagBlocks)
    CALL h5lite_get_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    CALL h5lite_get_attribute(matrix_id,"imatrixSpec",rmatrix%imatrixSpec)
    
    ! For the time being, we suppose that block matrices have quadratic shape,
    ! that is, ndiagBlocks is the maximum of column and row blocks
    IF (rmatrix%ndiagBlocks > 0) THEN
      ALLOCATE(rmatrix%RmatrixBlock(rmatrix%ndiagBlocks,rmatrix%ndiagBlocks))

      ! Loop over all blocks
      DO iblock=1,rmatrix%ndiagBlocks
        DO jblock=1,rmatrix%ndiagBlocks
          CALL lsysh5io_readMatrix(matrix_id,"RmatrixBlock("//&
              TRIM(ADJUSTL(sys_si(iblock,6)))//","//&
              TRIM(ADJUSTL(sys_si(jblock,6)))//")",rmatrix%RmatrixBlock(iblock,jblock))
        END DO
      END DO
    END IF
  END SUBROUTINE lsysh5io_readBlockMatrix_LocationByID

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readBlockMatrix_LocationByName(loc_id,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)          :: loc_id

    ! Name of the matrix
    CHARACTER(LEN=*), INTENT(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: matrix_id

    ! Determine the matrix_id and call the working routine
    CALL h5lite_get_group(loc_id,cmatrixName,matrix_id)
    IF (matrix_id < 0) THEN
      PRINT *, "lsysh5io_readBlockMatrix_LocationByName: Unable to get matrix_id!"
      STOP
    END IF
    CALL lsysh5io_readBlockMatrix(loc_id,matrix_id,rmatrix)
  END SUBROUTINE lsysh5io_readBlockMatrix_LocationByName

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeVector_Root(rh5file,rvector,vector_id,cvectorName)

!<description>
    ! This subroutine writes the vector to the root of the hdf5 file that
    ! is associated with the h5file descriptor. If the optional parameter
    ! cvectorName is given, then this value is used to identify the
    ! vector. Otherwise, t_vectorScalar.XXX is adopted, whereby XXX is the
    ! next free number that is not present in the file.
    ! The handle to the vector object is returned as vector_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)             :: rh5file

    ! The vector to be written
    TYPE(t_vectorScalar), INTENT(IN)       :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    INTEGER(HID_T), INTENT(OUT)            :: vector_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_writeVector(rh5file%file_id,rvector,vector_id,cvectorName)
  END SUBROUTINE lsysh5io_writeVector_Root

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeVector_Location(loc_id,rvector,vector_id,cvectorName)

!<description>
    ! This subroutine writes the vector to the location given by loc_id.
    ! If the optional parameter cvectorName is given, then this value is 
    ! used to identify the vector. Otherwise, t_vectorScalar.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the vector object is returned as vector_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)             :: loc_id

    ! The vector to be written
    TYPE(t_vectorScalar), INTENT(IN)       :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    INTEGER(HID_T), INTENT(OUT)            :: vector_id
!</output>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    INTEGER,  DIMENSION(:), POINTER :: P_Idata

    ! Create new group for t_vectorScalar
    IF (PRESENT(cvectorName)) THEN
      CALL h5lite_set_group(loc_id,cvectorName,vector_id)
    ELSE
      CALL h5lite_set_group(loc_id,"t_vectorScalar",vector_id,.TRUE.)
    END IF

    ! Set attributes for scalar vector
    CALL h5lite_set_attribute(vector_id,"NEQ",           rvector%NEQ)
    CALL h5lite_set_attribute(vector_id,"NVAR",          rvector%NVAR)
    CALL h5lite_set_attribute(vector_id,"cdataType",     rvector%cdataType)
    CALL h5lite_set_attribute(vector_id,"isortStrategy", rvector%isortStrategy)
    CALL h5lite_set_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    CALL h5lite_set_attribute(vector_id,"bisCopy",       rvector%bisCopy)

    ! Set data for scalar vector
    CALL h5lite_set_data(vector_id,"ITags",rvector%ITags)
    CALL h5lite_set_data(vector_id,"DTags",rvector%DTags)

    ! Write data of data handle
    IF (rvector%h_Ddata /= ST_NOHANDLE) THEN
      SELECT CASE(rvector%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double(rvector,p_Ddata)
        CALL h5lite_set_data(vector_id,"h_Ddata",p_Ddata)
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single(rvector,p_Fdata)
        CALL h5lite_set_data(vector_id,"h_Ddata",p_Fdata)
      CASE (ST_INT)
        CALL lsyssc_getbase_int(rvector,p_Idata)
        CALL h5lite_set_data(vector_id,"h_Ddata",p_Idata)
      CASE DEFAULT
        PRINT *, "lsysh5io_writeVector_Location: Unsupported data format!"
        STOP
      END SELECT
    END IF

    ! Write data of sorting array
    IF (rvector%h_IsortPermutation /= ST_NOHANDLE) THEN
      CALL storage_getbase_int(rvector%h_IsortPermutation,p_Idata,2*rvector%NEQ)
      CALL h5lite_set_data(vector_id,"h_IsortPermutation",p_Idata)
    END IF
  END SUBROUTINE lsysh5io_writeVector_Location

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeBlockVector_Root(rh5file,rvector,vector_id,cvectorName)

!<description>
    ! This subroutine writes the vector to the root of the hdf5 file that
    ! is associated with the h5file descriptor. If the optional parameter
    ! cvectorName is given, then this value is used to identify the
    ! vector. Otherwise, t_vectorBlock.XXX is adopted, whereby XXX is the
    ! next free number that is not present in the file.
    ! The handle to the vector object is returned as vector_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)             :: rh5file

    ! The vector to be written
    TYPE(t_vectorBlock), INTENT(IN)        :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    INTEGER(HID_T), INTENT(OUT)            :: vector_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_writeBlockVector(rh5file%file_id,rvector,vector_id,cvectorName)
  END SUBROUTINE lsysh5io_writeBlockVector_Root

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_writeBlockVector_Location(loc_id,rvector,vector_id,cvectorName)

!<description>
    ! This subroutine writes the vector to the location given by loc_id.
    ! If the optional parameter cvectorName is given, then this value is 
    ! used to identify the vector. Otherwise, t_vectorBlock.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the vector object is returned as vector_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)             :: loc_id

    ! The vector to be written
    TYPE(t_vectorBlock), INTENT(IN)        :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    INTEGER(HID_T), INTENT(OUT)            :: vector_id
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: vectorScalar_id
    INTEGER :: iblock

    ! Create new group for t_vectorScalar
    IF (PRESENT(cvectorName)) THEN
      CALL h5lite_set_group(loc_id,cvectorName,vector_id)
    ELSE
      CALL h5lite_set_group(loc_id,"t_vectorBlock",vector_id,.TRUE.)
    END IF
    
    ! Set attributes for block vector
    CALL h5lite_set_attribute(vector_id,"NEQ",           rvector%NEQ)
    CALL h5lite_set_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    CALL h5lite_set_attribute(vector_id,"cdataType",     rvector%cdataType)
    CALL h5lite_set_attribute(vector_id,"bisCopy",       rvector%bisCopy)
    CALL h5lite_set_attribute(vector_id,"nblocks",       rvector%nblocks)

    ! Write each scalar subvector
    IF (ASSOCIATED(rvector%RvectorBlock)) THEN
      
      ! Loop over all vector blocks
      DO iblock=1,rvector%nblocks
        CALL lsysh5io_writeVector(vector_id,rvector%RvectorBlock(iblock),&
            vectorScalar_id,"RvectorBlock("//TRIM(ADJUSTL(sys_si(iblock,6)))//")")
      END DO
      
    END IF
  END SUBROUTINE lsysh5io_writeBlockVector_Location

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readVector_RootByID(rh5file,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its unique vector_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)          :: rh5file

    ! Handle to the vector object
    INTEGER(HID_T), INTENT(IN)          :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readVector(rh5file%file_id,vector_id,rvector)
  END SUBROUTINE lsysh5io_readVector_RootByID

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readVector_RootByName(rh5file,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)          :: rh5file

    ! Name of the vector
    CHARACTER(LEN=*), INTENT(IN)        :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readVector(rh5file%file_id,cvectorName,rvector)
  END SUBROUTINE lsysh5io_readVector_RootByName
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readVector_LocationByID(loc_id,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its unique vector_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)          :: loc_id

    ! Handle to the vector object
    INTEGER(HID_T), INTENT(IN)          :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    INTEGER,  DIMENSION(:), POINTER :: p_Idata
    INTEGER :: isize

    ! Get attributes
    CALL h5lite_get_attribute(vector_id,"NEQ",           rvector%NEQ)
    CALL h5lite_get_attribute(vector_id,"NVAR",          rvector%NVAR)
    CALL h5lite_get_attribute(vector_id,"cdataType",     rvector%cdataType)
    CALL h5lite_get_attribute(vector_id,"isortStrategy", rvector%isortStrategy)
    CALL h5lite_get_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    CALL h5lite_get_attribute(vector_id,"bisCopy",       rvector%bisCopy)

    ! Get data
    CALL h5lite_get_data(vector_id,"ITags",rvector%ITags)
    CALL h5lite_get_data(vector_id,"DTags",rvector%DTags)

    ! Get data if present
    IF (h5lite_ispresent_data(vector_id,"h_Ddata")) THEN

      IF (rvector%h_Ddata == ST_NOHANDLE) THEN
        ! If the data handle is not associated to some memory block, allocate new memory
        CALL storage_new1d ('lsysh5io_readVector_LocationByID','Vector',rvector%NEQ,&
            rvector%cdataType,rvector%h_Ddata,ST_NEWBLOCK_NOINIT)
      ELSE
        ! Otherwise, check that the provided memory is sufficient
        CALL storage_getsize(rvector%h_Ddata,isize)
        IF (isize < rvector%iidxFirstEntry+rvector%NEQ-1) THEN
          PRINT *, "lsysh5io_readVector_LocationByID: Not enought memory!"
          STOP
        END IF
      END IF
      
      ! Get data handle
      SELECT CASE(rvector%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double(rvector,p_Ddata)
        CALL h5lite_get_data(vector_id,"h_Ddata",p_Ddata)
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single(rvector,p_Fdata)
        CALL h5lite_get_data(vector_id,"h_Ddata",p_Fdata)
        
      CASE (ST_INT)
        CALL lsyssc_getbase_int(rvector,p_Idata)
        CALL h5lite_get_data(vector_id,"h_Ddata",p_Idata)
        
      CASE DEFAULT
        PRINT *, "lsysh5io_readVector_LocationByID: Invalid data type!"
        STOP
      END SELECT
    END IF

    ! Get permutation array for sorting if present
    IF (h5lite_ispresent_data(vector_id,"h_IsortPermutation")) THEN

      IF (rvector%h_IsortPermutation == ST_NOHANDLE) THEN
        ! If the permutation handle is not associated to some memory block, allocate new memory
        CALL storage_new1d ('lsysh5io_readVector_LocationByID','IsortPermutation',2*rvector%NEQ,&
            ST_INT,rvector%h_IsortPermutation,ST_NEWBLOCK_NOINIT)
      ELSE
        ! Otherwise, check that the provided memory is sufficient
        CALL storage_getsize(rvector%h_IsortPermutation,isize)
        IF (isize /= 2*rvector%NEQ) THEN
          PRINT *, "lsysh5io_readVector_LocationByID: Wrong size of permutation array!"
          STOP
        END IF
      END IF

      ! Get permutation handle
      CALL storage_getbase_int(rvector%h_IsortPermutation,p_Idata)
      CALL h5lite_get_data(vector_id,"h_IsortPermutation",p_Idata)
    END IF
  END SUBROUTINE lsysh5io_readVector_LocationByID

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readVector_LocationByName(loc_id,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)          :: loc_id

    ! Name of the vector
    CHARACTER(LEN=*), INTENT(IN)        :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: vector_id

    ! Determine the vector_id and call the working routine
    CALL h5lite_get_group(loc_id,cvectorName,vector_id)
    IF (vector_id < 0) THEN
      PRINT *, "lsysh5io_readVector_LocationByName: Unable to get vector_id!"
      STOP
    END IF
    CALL lsysh5io_readVector(loc_id,vector_id,rvector)
  END SUBROUTINE lsysh5io_readVector_LocationByName

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readBlockVector_RootByID(rh5file,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its unique vector_id.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)         :: rh5file

    ! Handle to the vector object
    INTEGER(HID_T), INTENT(IN)         :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readBlockVector(rh5file%file_id,vector_id,rvector)
  END SUBROUTINE lsysh5io_readBlockVector_RootByID

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readBlockVector_RootByName(rh5file,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! The hdf file descriptor
    TYPE(t_h5file), INTENT(IN)         :: rh5file

    ! Name of the vector
    CHARACTER(LEN=*), INTENT(IN)       :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    CALL lsysh5io_readBlockVector(rh5file%file_id,cvectorName,rvector)
  END SUBROUTINE lsysh5io_readBlockVector_RootByName

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readBlockVector_LocationByID(loc_id,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its uniquad vector_id.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)         :: loc_id

    ! Handle to the vector object
    INTEGER(HID_T), INTENT(IN)         :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER :: iblock

    ! Get attributes
    CALL h5lite_get_attribute(vector_id,"NEQ",           rvector%NEQ)
    CALL h5lite_get_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    CALL h5lite_get_attribute(vector_id,"cdataType",     rvector%cdataType)
    CALL h5lite_set_attribute(vector_id,"bisCopy",       rvector%bisCopy)
    CALL h5lite_get_attribute(vector_id,"nblocks",       rvector%nblocks)
    
    ! If the data handle is not associated to some memory block, allocate new memory
    IF (rvector%h_Ddata == ST_NOHANDLE) &
        CALL storage_new1d ('lsysh5io_readBlockVector_LocationByID','Vector',rvector%NEQ,&
        rvector%cdataType,rvector%h_Ddata,ST_NEWBLOCK_NOINIT)

    ! Allocate 1D array  with scalar vectors
    ALLOCATE(rvector%RvectorBlock(rvector%nblocks))
    
    ! Read in scalar subvectors
    DO iblock=1,rvector%nblocks
      
      ! Set the handle for the scalar subvector
      rvector%RvectorBlock(iblock)%h_Ddata = rvector%h_Ddata
      
      ! Read in the scalar subvector
      CALL lsysh5io_readVector(vector_id,"RvectorBlock("//TRIM(ADJUSTL(sys_si(iblock,6)))//")",&
          rvector%RvectorBlock(iblock))
    END DO
  END SUBROUTINE lsysh5io_readBlockVector_LocationByID

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysh5io_readBlockVector_LocationByName(loc_id,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)         :: loc_id

    ! Name of the vector
    CHARACTER(LEN=*), INTENT(IN)       :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: vector_id

    ! Determine the vector_id and call the working routine
    CALL h5lite_get_group(loc_id,cvectorName,vector_id)
    IF (vector_id < 0) THEN
      PRINT *, "lsysh5io_readBlockVector_LocationByName: Unable to get vector_id!"
      STOP
    END IF
    CALL lsysh5io_readBlockVector(loc_id,vector_id,rvector)
  END SUBROUTINE lsysh5io_readBlockVector_LocationByName
END MODULE linearsystemh5io
