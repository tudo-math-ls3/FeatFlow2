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

module linearsystemh5io

  use fsystem
  use storage
  use linearsystemscalar
  use linearsystemblock
  use h5lite
  use hdf5

  implicit none
  
  private
  public :: lsysh5io_writeMatrix
  public :: lsysh5io_writeBlockMatrix
  public :: lsysh5io_readMatrix
  public :: lsysh5io_readBlockMatrix
  public :: lsysh5io_writeVector
  public :: lsysh5io_writeBlockVector
  public :: lsysh5io_readVector
  public :: lsysh5io_readBlockVector

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  interface lsysh5io_writeMatrix
    module procedure lsysh5io_writeMatrix_Root
    module procedure lsysh5io_writeMatrix_Location
  end interface

  interface lsysh5io_writeBlockMatrix
    module procedure lsysh5io_writeBlockMatrix_Root
    module procedure lsysh5io_writeBlockMatrix_Location
  end interface

  interface lsysh5io_readMatrix
    module procedure lsysh5io_readMatrix_RootByID
    module procedure lsysh5io_readMatrix_RootByName
    module procedure lsysh5io_readMatrix_LocationByID
    module procedure lsysh5io_readMatrix_LocationByName
  end interface

  interface lsysh5io_readBlockMatrix
    module procedure lsysh5io_readBlockMatrix_RootByID
    module procedure lsysh5io_readBlockMatrix_RootByName
    module procedure lsysh5io_readBlockMatrix_LocationByID
    module procedure lsysh5io_readBlockMatrix_LocationByName
  end interface

  interface lsysh5io_writeVector
    module procedure lsysh5io_writeVector_Root
    module procedure lsysh5io_writeVector_Location
  end interface

  interface lsysh5io_writeBlockVector
    module procedure lsysh5io_writeBlockVector_Root
    module procedure lsysh5io_writeBlockVector_Location
  end interface
  
  interface lsysh5io_readVector
    module procedure lsysh5io_readVector_RootByID
    module procedure lsysh5io_readVector_RootByName
    module procedure lsysh5io_readVector_LocationByID
    module procedure lsysh5io_readVector_LocationByName
  end interface

  interface lsysh5io_readBlockVector
    module procedure lsysh5io_readBlockVector_RootByID
    module procedure lsysh5io_readBlockVector_RootByName
    module procedure lsysh5io_readBlockVector_LocationByID
    module procedure lsysh5io_readBlockVector_LocationByName
  end interface

contains

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeMatrix_Root(rh5file,rmatrix,matrix_id,cmatrixName)

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
    type(t_h5file), intent(IN)             :: rh5file

    ! The matrix to be written
    type(t_matrixScalar), intent(IN)       :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    character(LEN=*), intent(IN), optional :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    integer(HID_T), intent(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_writeMatrix(rh5file%file_id,rmatrix,matrix_id,cmatrixName)
  end subroutine lsysh5io_writeMatrix_Root

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeMatrix_Location(loc_id,rmatrix,matrix_id,cmatrixName)

!<description>
    ! This subroutine writes the matrix to the location given by loc_id.
    ! If the optional parameter cmatrixName is given, then this value is 
    ! used to identify the matrix. Otherwise, t_matrixScalar.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the matrix object is returned as matrix_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)             :: loc_id

    ! The matrix to be written
    type(t_matrixScalar), intent(IN)       :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    character(LEN=*), intent(IN), optional :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    integer(HID_T), intent(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    integer,  dimension(:), pointer :: P_Idata

    ! Create new group for t_matrixScalar
    if (present(cmatrixName)) then
      call h5lite_set_group(loc_id,cmatrixName,matrix_id)
    else
      call h5lite_set_group(loc_id,"t_matrixScalar",matrix_id,.true.)
    end if

    ! Set attributes for scalar matrix
    call h5lite_set_attribute(matrix_id,"cmatrixFormat",          rmatrix%cmatrixFormat)
    call h5lite_set_attribute(matrix_id,"cinterleavematrixFormat",rmatrix%cinterleavematrixFormat)
    call h5lite_set_attribute(matrix_id,"imatrixSpec",            rmatrix%imatrixSpec)
    call h5lite_set_attribute(matrix_id,"NA",                     rmatrix%NA)
    call h5lite_set_attribute(matrix_id,"NEQ",                    rmatrix%NEQ)
    call h5lite_set_attribute(matrix_id,"NCOLS",                  rmatrix%NCOLS)
    call h5lite_set_attribute(matrix_id,"NVAR",                   rmatrix%NVAR)
    call h5lite_set_attribute(matrix_id,"dscaleFactor",           rmatrix%dscaleFactor)
    call h5lite_set_attribute(matrix_id,"isortStrategy",          rmatrix%isortStrategy)
    call h5lite_set_attribute(matrix_id,"cdataType",              rmatrix%cdataType)

    ! Set data for scalar matrix
    call h5lite_set_data(matrix_id,"ITags",rmatrix%ITags)
    call h5lite_set_data(matrix_id,"DTags",rmatrix%DTags)
    
    ! Write handle for sorting permutation
    if (rmatrix%h_IsortPermutation /= ST_NOHANDLE) then
      call storage_getbase_int(rmatrix%h_IsortPermutation,p_Idata,2*rmatrix%NEQ)
      call h5lite_set_data(matrix_id,"h_IsortPermutation",p_Idata)
    end if

    ! Write handle for matrix data
    if (rmatrix%h_DA /= ST_NOHANDLE) then
      select case(rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix,p_Ddata)
        call h5lite_set_data(matrix_id,"h_DA",p_Ddata)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix,p_Fdata)
        call h5lite_set_data(matrix_id,"h_DA",p_Fdata)
      case(ST_INT)
        call lsyssc_getbase_int(rmatrix,p_Idata)
        call h5lite_set_data(matrix_id,"h_DA",p_Idata)
      case DEFAULT
        print *, "lsysh5io_writeMatrix_Location: Unsupported data format!"
        stop
      end select
    end if

    ! Write handle for column structure
    if (rmatrix%h_Kcol /= ST_NOHANDLE) then
      call lsyssc_getbase_Kcol(rmatrix,p_Idata)
      call h5lite_set_data(matrix_id,"h_Kcol",p_Idata)
    end if

    ! Write handle for row structure
    if (rmatrix%h_Kld /= ST_NOHANDLE) then
      call lsyssc_getbase_Kld(rmatrix,p_Idata)
      call h5lite_set_data(matrix_id,"h_Kld",p_Idata)
    end if

    ! Write handle for diagonal structure
    if (rmatrix%h_Kdiagonal /= ST_NOHANDLE) then
      call lsyssc_getbase_Kdiagonal(rmatrix,p_Idata)
      call h5lite_set_data(matrix_id,"h_Kdiagonal",p_Idata)
    end if
  end subroutine lsysh5io_writeMatrix_Location

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeBlockMatrix_Root(rh5file,rmatrix,matrix_id,cmatrixName)

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
    type(t_h5file), intent(IN)             :: rh5file

    ! The matrix to be written
    type(t_matrixBlock), intent(IN)        :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    character(LEN=*), intent(IN), optional :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    integer(HID_T), intent(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_writeBlockMatrix(rh5file%file_id,rmatrix,matrix_id,cmatrixName)
  end subroutine lsysh5io_writeBlockMatrix_Root

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeBlockMatrix_Location(loc_id,rmatrix,matrix_id,cmatrixName)

!<description>
    ! This subroutine writes the matrix to the location given by loc_id.
    ! If the optional parameter cmatrixName is given, then this value is 
    ! used to identify the matrix. Otherwise, t_matrixBlock.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the matrix object is returned as matrix_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)             :: loc_id

    ! The matrix to be written
    type(t_matrixBlock), intent(IN)        :: rmatrix

    ! OPTIONAL: The name to be used to identify the matrix
    character(LEN=*), intent(IN), optional :: cmatrixName
!</input>

!<output>
    ! Handle to the matrix object
    integer(HID_T), intent(OUT)            :: matrix_id
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: matrixScalar_id
    integer :: iblock,jblock

    ! Create new group for t_matrixScalar
    if (present(cmatrixName)) then
      call h5lite_set_group(loc_id,cmatrixName,matrix_id)
    else
      call h5lite_set_group(loc_id,"t_matrixBlock",matrix_id,.true.)
    end if

    ! Set attributes for block matrix
    call h5lite_set_attribute(matrix_id,"NEQ",        rmatrix%NEQ)
    call h5lite_set_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    call h5lite_set_attribute(matrix_id,"ndiagBlocks",rmatrix%ndiagBlocks)
    call h5lite_set_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    call h5lite_set_attribute(matrix_id,"imatrixSpec",rmatrix%imatrixSpec)
    
    ! Write each scalar submatrix
    if (associated(rmatrix%RmatrixBlock)) then
      
      ! Loop over all blocks
      do iblock=lbound(rmatrix%RmatrixBlock,1),ubound(rmatrix%RmatrixBlock,1)
        do jblock=lbound(rmatrix%RmatrixBlock,2),ubound(rmatrix%RmatrixBlock,2)
          call lsysh5io_writeMatrix(matrix_id,rmatrix%RmatrixBlock(iblock,jblock),&
              matrixScalar_id,"RmatrixBlock("//trim(adjustl(sys_si(iblock,6)))//","//&
              trim(adjustl(sys_si(jblock,6)))//")")
        end do
      end do
      
    end if
  end subroutine lsysh5io_writeBlockMatrix_Location

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysh5io_readMatrix_RootByID(rh5file,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)          :: rh5file

    ! Handle to the matrix object
    integer(HID_T), intent(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The matrix that should be read
    type(t_matrixScalar), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readMatrix(rh5file%file_id,matrix_id,rmatrix)
  end subroutine lsysh5io_readMatrix_RootByID

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readMatrix_RootByName(rh5file,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)          :: rh5file

    ! Name of the matrix
    character(LEN=*), intent(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    type(t_matrixScalar), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readMatrix(rh5file%file_id,cmatrixName,rmatrix)
  end subroutine lsysh5io_readMatrix_RootByName

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readMatrix_LocationByID(loc_id,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)          :: loc_id

    ! Handle to the matrix object
    integer(HID_T), intent(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_matrixScalar), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    integer,  dimension(:), pointer :: p_Idata
    integer :: isize,sizeKld,sizeKdiagonal

    ! Get attributes
    call h5lite_get_attribute(matrix_id,"cmatrixFormat",          rmatrix%cmatrixFormat)
    call h5lite_get_attribute(matrix_id,"cinterleavematrixFormat",rmatrix%cinterleavematrixFormat)
    call h5lite_get_attribute(matrix_id,"imatrixSpec",            rmatrix%imatrixSpec)
    call h5lite_get_attribute(matrix_id,"NA",                     rmatrix%NA)
    call h5lite_get_attribute(matrix_id,"NEQ",                    rmatrix%NEQ)
    call h5lite_get_attribute(matrix_id,"NCOLS",                  rmatrix%NCOLS)
    call h5lite_get_attribute(matrix_id,"NVAR",                   rmatrix%NVAR)
    call h5lite_get_attribute(matrix_id,"dscaleFactor",           rmatrix%dscaleFactor)
    call h5lite_get_attribute(matrix_id,"isortStrategy",          rmatrix%isortStrategy)
    call h5lite_get_attribute(matrix_id,"cdataType",              rmatrix%cdataType)

    ! Get data
    call h5lite_get_data(matrix_id,"ITags",rmatrix%ITags)
    call h5lite_get_data(matrix_id,"DTags",rmatrix%DTags)

    ! Get handle for sorting permutation
    if (h5lite_ispresent_data(matrix_id,"h_IsortPermutation")) then
      
      if (rmatrix%h_IsortPermutation == ST_NOHANDLE) then
        ! If the permutation handle is not associated to some memory block, allocate new memory
        call storage_new1d ('lsysh5io_readMatrix_LocationByID','IsortPermutation',2*rmatrix%NEQ,&
            ST_INT,rmatrix%h_IsortPermutation,ST_NEWBLOCK_NOINIT)
      else
        ! Otherwise, check that the provided memory is sufficient
        call storage_getsize(rmatrix%h_IsortPermutation,isize)
        if (isize /= 2*rmatrix%NEQ) then
          print *, "lsysh5io_readMatrix_LocationByID: Wrong size of permutation array!"
          stop
        end if
      end if

      ! Get permutation handle
      call storage_getbase_int(rmatrix%h_IsortPermutation,p_Idata)
      call h5lite_get_data(matrix_id,"h_IsortPermutation",p_Idata)
    end if

    ! Get handle for matrix data
    if (h5lite_ispresent_data(matrix_id,"h_DA")) then

      if (rmatrix%h_DA == ST_NOHANDLE) then
        ! If the data handle is not associated to some memory block, allocate new memory
        call storage_new1d ('lsysh5io_readMatrix_LocationByID','h_DA',rmatrix%NA,&
            rmatrix%cdataType,rmatrix%h_DA,ST_NEWBLOCK_NOINIT)
      else
        ! Otherwise, check that the provided memory is sufficient
        call storage_getsize(rmatrix%h_DA,isize)
        select case(rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXUNDEFINED)
          if (isize < rmatrix%NA) then
            print *, "lsysh5io_readMatrix_LocationByID: Not enough memory!"
            stop
          end if

        case (LSYSSC_MATRIX1)
          if (isize < rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR) then
            print *, "lsysh5io_readMatrix_LocationByID: Not enough memory!"
            stop
          end if

        case (LSYSSC_MATRIXD)
          if (isize < rmatrix%NA*rmatrix%NVAR) then
            print *, "lsysh5io_readMatrix_LocationByID: Not enough memory!"
            stop
          end if

        case DEFAULT
          print *, "lsysh5io_readMatrix_LocationByID: Unsupported interleave matrix format!"
          stop
        end select
      end if

      ! Get data handle
      select case(rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix,p_Ddata)
        call h5lite_get_data(matrix_id,"h_DA",p_Ddata)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix,p_Fdata)
        call h5lite_get_data(matrix_id,"h_DA",p_Fdata)
      case(ST_INT)
        call lsyssc_getbase_int(rmatrix,p_Idata)
        call h5lite_get_data(matrix_id,"h_DA",p_Idata)
      case DEFAULT
        print *, "lsysh5io_readMatrix_Location: Unsupported data format!"
        stop
      end select
    end if

    ! Get handle for column structure
    if (h5lite_ispresent_data(matrix_id,"h_Kcol")) then

      if (rmatrix%h_Kcol == ST_NOHANDLE) then
        ! If the column handle is not associated to some memory block, allocate new memory
        call storage_new1d ('lsysh5io_readMatrix_LocationByID','h_Kcol',rmatrix%NA,&
            ST_INT,rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
      else
        ! Otherwise, check that the provided memory is sufficient
        call storage_getsize(rmatrix%h_Kcol,isize)
        if (isize < rmatrix%NA) then
          print *, "lsysh5io_readMatrix_Location: Not enough memory!"
          stop
        end if
      end if

      ! Get column handle
      call lsyssc_getbase_Kcol(rmatrix,p_Idata)
      call h5lite_get_data(matrix_id,"h_Kcol",p_Idata)
    end if

    ! Get handle for row structure
    if (h5lite_ispresent_data(matrix_id,"h_Kld")) then
      
      ! Is the matrix transposed or not?
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1) then
        sizeKld=rmatrix%NCOLS+1
      else
        sizeKld=rmatrix%NEQ+1
      end if

      if (rmatrix%h_Kld == ST_NOHANDLE) then
        ! If the column handle is not associated to some memory block, allocate new memory
        call storage_new1d ('lsysh5io_readMatrix_LocationByID','h_Kld',sizeKld,&
            ST_INT,rmatrix%h_Kld,ST_NEWBLOCK_NOINIT)
      else
        ! Otherwise, check that the provided memory is sufficient
        call storage_getsize(rmatrix%h_Kld,isize)
        if (isize < sizeKld) then
          print *, "lsysh5io_readMatrix_Location: Not enough memory!"
          stop
        end if
      end if

      ! Get row handle
      call lsyssc_getbase_Kld(rmatrix,p_Idata)
      call h5lite_get_data(matrix_id,"h_Kld",p_Idata)
    end if

    ! Get handle for diagonal structure
    if (h5lite_ispresent_data(matrix_id,"h_Kdiagonal")) then

      ! Is the matrix transposed or not?
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1) then
        sizeKdiagonal=rmatrix%NCOLS
      else
        sizeKdiagonal=rmatrix%NEQ
      end if

      if (rmatrix%h_Kdiagonal == ST_NOHANDLE) then
        ! If the column handle is not associated to some memory block, allocate new memory
        call storage_new1d ('lsysh5io_readMatrix_LocationByID','h_Kdiagonal',sizeKdiagonal,&
            ST_INT,rmatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
      else
        ! Otherwise, check that the provided memory is sufficient
        call storage_getsize(rmatrix%h_Kdiagonal,isize)
        if (isize < sizeKdiagonal) then
          print *, "lsysh5io_readMatrix_Location: Not enough memory!"
          stop
        end if
      end if

      ! Get diagonal handle
      call lsyssc_getbase_Kdiagonal(rmatrix,p_Idata)
      call h5lite_get_data(matrix_id,"h_Kdiagonal",p_Idata)
    end if
  end subroutine lsysh5io_readMatrix_LocationByID
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readMatrix_LocationByName(loc_id,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)          :: loc_id

    ! Name of the matrix
    character(LEN=*), intent(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    type(t_matrixScalar), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: matrix_id

    ! Determine the matrix_id and call the working routine
    call h5lite_get_group(loc_id,cmatrixName,matrix_id)
    if (matrix_id < 0) then
      print *, "lsysh5io_readMatrix_LocationByName: Unable to get matrix_id!"
      stop
    end if
    call lsysh5io_readMatrix(loc_id,matrix_id,rmatrix)
  end subroutine lsysh5io_readMatrix_LocationByName

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysh5io_readBlockMatrix_RootByID(rh5file,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)          :: rh5file

    ! Handle to the matrix object
    integer(HID_T), intent(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The matrix that should be read
    type(t_matrixBlock), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readBlockMatrix(rh5file%file_id,matrix_id,rmatrix)
  end subroutine lsysh5io_readBlockMatrix_RootByID

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readBlockMatrix_RootByName(rh5file,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)          :: rh5file

    ! Name of the matrix
    character(LEN=*), intent(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    type(t_matrixBlock), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readBlockMatrix(rh5file%file_id,cmatrixName,rmatrix)
  end subroutine lsysh5io_readBlockMatrix_RootByName

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readBlockMatrix_LocationByID(loc_id,matrix_id,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its unique matrix_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)          :: loc_id

    ! Handle to the matrix object
    integer(HID_T), intent(IN)          :: matrix_id
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_matrixBlock), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iblock,jblock

    ! Get attributes for block matrix
    call h5lite_get_attribute(matrix_id,"NEQ",        rmatrix%NEQ)
    call h5lite_get_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    call h5lite_get_attribute(matrix_id,"ndiagBlocks",rmatrix%ndiagBlocks)
    call h5lite_get_attribute(matrix_id,"NCOLS",      rmatrix%NCOLS)
    call h5lite_get_attribute(matrix_id,"imatrixSpec",rmatrix%imatrixSpec)
    
    ! For the time being, we suppose that block matrices have quadratic shape,
    ! that is, ndiagBlocks is the maximum of column and row blocks
    if (rmatrix%ndiagBlocks > 0) then
      allocate(rmatrix%RmatrixBlock(rmatrix%ndiagBlocks,rmatrix%ndiagBlocks))

      ! Loop over all blocks
      do iblock=1,rmatrix%ndiagBlocks
        do jblock=1,rmatrix%ndiagBlocks
          call lsysh5io_readMatrix(matrix_id,"RmatrixBlock("//&
              trim(adjustl(sys_si(iblock,6)))//","//&
              trim(adjustl(sys_si(jblock,6)))//")",rmatrix%RmatrixBlock(iblock,jblock))
        end do
      end do
    end if
  end subroutine lsysh5io_readBlockMatrix_LocationByID

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readBlockMatrix_LocationByName(loc_id,cmatrixName,rmatrix)

!<description>
    ! This subroutine reads the matrix from the location given by loc_id.
    ! The matrix is identified by its name given by cmatrixName.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)          :: loc_id

    ! Name of the matrix
    character(LEN=*), intent(IN)        :: cmatrixName
!</input>

!<inputoutput>
    ! The matrix that should be read
    type(t_matrixBlock), intent(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: matrix_id

    ! Determine the matrix_id and call the working routine
    call h5lite_get_group(loc_id,cmatrixName,matrix_id)
    if (matrix_id < 0) then
      print *, "lsysh5io_readBlockMatrix_LocationByName: Unable to get matrix_id!"
      stop
    end if
    call lsysh5io_readBlockMatrix(loc_id,matrix_id,rmatrix)
  end subroutine lsysh5io_readBlockMatrix_LocationByName

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeVector_Root(rh5file,rvector,vector_id,cvectorName)

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
    type(t_h5file), intent(IN)             :: rh5file

    ! The vector to be written
    type(t_vectorScalar), intent(IN)       :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    character(LEN=*), intent(IN), optional :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    integer(HID_T), intent(OUT)            :: vector_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_writeVector(rh5file%file_id,rvector,vector_id,cvectorName)
  end subroutine lsysh5io_writeVector_Root

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeVector_Location(loc_id,rvector,vector_id,cvectorName)

!<description>
    ! This subroutine writes the vector to the location given by loc_id.
    ! If the optional parameter cvectorName is given, then this value is 
    ! used to identify the vector. Otherwise, t_vectorScalar.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the vector object is returned as vector_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)             :: loc_id

    ! The vector to be written
    type(t_vectorScalar), intent(IN)       :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    character(LEN=*), intent(IN), optional :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    integer(HID_T), intent(OUT)            :: vector_id
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    integer,  dimension(:), pointer :: P_Idata

    ! Create new group for t_vectorScalar
    if (present(cvectorName)) then
      call h5lite_set_group(loc_id,cvectorName,vector_id)
    else
      call h5lite_set_group(loc_id,"t_vectorScalar",vector_id,.true.)
    end if

    ! Set attributes for scalar vector
    call h5lite_set_attribute(vector_id,"NEQ",           rvector%NEQ)
    call h5lite_set_attribute(vector_id,"NVAR",          rvector%NVAR)
    call h5lite_set_attribute(vector_id,"cdataType",     rvector%cdataType)
    call h5lite_set_attribute(vector_id,"isortStrategy", rvector%isortStrategy)
    call h5lite_set_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    call h5lite_set_attribute(vector_id,"bisCopy",       rvector%bisCopy)

    ! Set data for scalar vector
    call h5lite_set_data(vector_id,"ITags",rvector%ITags)
    call h5lite_set_data(vector_id,"DTags",rvector%DTags)

    ! Write data of data handle
    if (rvector%h_Ddata /= ST_NOHANDLE) then
      select case(rvector%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rvector,p_Ddata)
        call h5lite_set_data(vector_id,"h_Ddata",p_Ddata)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rvector,p_Fdata)
        call h5lite_set_data(vector_id,"h_Ddata",p_Fdata)
      case (ST_INT)
        call lsyssc_getbase_int(rvector,p_Idata)
        call h5lite_set_data(vector_id,"h_Ddata",p_Idata)
      case DEFAULT
        print *, "lsysh5io_writeVector_Location: Unsupported data format!"
        stop
      end select
    end if

    ! Write data of sorting array
    if (rvector%h_IsortPermutation /= ST_NOHANDLE) then
      call storage_getbase_int(rvector%h_IsortPermutation,p_Idata,2*rvector%NEQ)
      call h5lite_set_data(vector_id,"h_IsortPermutation",p_Idata)
    end if
  end subroutine lsysh5io_writeVector_Location

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeBlockVector_Root(rh5file,rvector,vector_id,cvectorName)

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
    type(t_h5file), intent(IN)             :: rh5file

    ! The vector to be written
    type(t_vectorBlock), intent(IN)        :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    character(LEN=*), intent(IN), optional :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    integer(HID_T), intent(OUT)            :: vector_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_writeBlockVector(rh5file%file_id,rvector,vector_id,cvectorName)
  end subroutine lsysh5io_writeBlockVector_Root

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_writeBlockVector_Location(loc_id,rvector,vector_id,cvectorName)

!<description>
    ! This subroutine writes the vector to the location given by loc_id.
    ! If the optional parameter cvectorName is given, then this value is 
    ! used to identify the vector. Otherwise, t_vectorBlock.XXX is adopted,
    ! whereby XXX is the next free number that is not present at the location
    ! loc_id. The handle to the vector object is returned as vector_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)             :: loc_id

    ! The vector to be written
    type(t_vectorBlock), intent(IN)        :: rvector

    ! OPTIONAL: The name to be used to identify the vector
    character(LEN=*), intent(IN), optional :: cvectorName
!</input>

!<output>
    ! Handle to the vector object
    integer(HID_T), intent(OUT)            :: vector_id
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: vectorScalar_id
    integer :: iblock

    ! Create new group for t_vectorScalar
    if (present(cvectorName)) then
      call h5lite_set_group(loc_id,cvectorName,vector_id)
    else
      call h5lite_set_group(loc_id,"t_vectorBlock",vector_id,.true.)
    end if
    
    ! Set attributes for block vector
    call h5lite_set_attribute(vector_id,"NEQ",           rvector%NEQ)
    call h5lite_set_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    call h5lite_set_attribute(vector_id,"cdataType",     rvector%cdataType)
    call h5lite_set_attribute(vector_id,"bisCopy",       rvector%bisCopy)
    call h5lite_set_attribute(vector_id,"nblocks",       rvector%nblocks)

    ! Write each scalar subvector
    if (associated(rvector%RvectorBlock)) then
      
      ! Loop over all vector blocks
      do iblock=1,rvector%nblocks
        call lsysh5io_writeVector(vector_id,rvector%RvectorBlock(iblock),&
            vectorScalar_id,"RvectorBlock("//trim(adjustl(sys_si(iblock,6)))//")")
      end do
      
    end if
  end subroutine lsysh5io_writeBlockVector_Location

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readVector_RootByID(rh5file,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its unique vector_id.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)          :: rh5file

    ! Handle to the vector object
    integer(HID_T), intent(IN)          :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorScalar), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readVector(rh5file%file_id,vector_id,rvector)
  end subroutine lsysh5io_readVector_RootByID

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readVector_RootByName(rh5file,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)          :: rh5file

    ! Name of the vector
    character(LEN=*), intent(IN)        :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorScalar), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readVector(rh5file%file_id,cvectorName,rvector)
  end subroutine lsysh5io_readVector_RootByName
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readVector_LocationByID(loc_id,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its unique vector_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)          :: loc_id

    ! Handle to the vector object
    integer(HID_T), intent(IN)          :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorScalar), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    integer,  dimension(:), pointer :: p_Idata
    integer :: isize

    ! Get attributes
    call h5lite_get_attribute(vector_id,"NEQ",           rvector%NEQ)
    call h5lite_get_attribute(vector_id,"NVAR",          rvector%NVAR)
    call h5lite_get_attribute(vector_id,"cdataType",     rvector%cdataType)
    call h5lite_get_attribute(vector_id,"isortStrategy", rvector%isortStrategy)
    call h5lite_get_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    call h5lite_get_attribute(vector_id,"bisCopy",       rvector%bisCopy)

    ! Get data
    call h5lite_get_data(vector_id,"ITags",rvector%ITags)
    call h5lite_get_data(vector_id,"DTags",rvector%DTags)

    ! Get data if present
    if (h5lite_ispresent_data(vector_id,"h_Ddata")) then

      if (rvector%h_Ddata == ST_NOHANDLE) then
        ! If the data handle is not associated to some memory block, allocate new memory
        call storage_new1d ('lsysh5io_readVector_LocationByID','Vector',rvector%NEQ,&
            rvector%cdataType,rvector%h_Ddata,ST_NEWBLOCK_NOINIT)
      else
        ! Otherwise, check that the provided memory is sufficient
        call storage_getsize(rvector%h_Ddata,isize)
        if (isize < rvector%iidxFirstEntry+rvector%NEQ-1) then
          print *, "lsysh5io_readVector_LocationByID: Not enought memory!"
          stop
        end if
      end if
      
      ! Get data handle
      select case(rvector%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rvector,p_Ddata)
        call h5lite_get_data(vector_id,"h_Ddata",p_Ddata)
        
      case (ST_SINGLE)
        call lsyssc_getbase_single(rvector,p_Fdata)
        call h5lite_get_data(vector_id,"h_Ddata",p_Fdata)
        
      case (ST_INT)
        call lsyssc_getbase_int(rvector,p_Idata)
        call h5lite_get_data(vector_id,"h_Ddata",p_Idata)
        
      case DEFAULT
        print *, "lsysh5io_readVector_LocationByID: Invalid data type!"
        stop
      end select
    end if

    ! Get permutation array for sorting if present
    if (h5lite_ispresent_data(vector_id,"h_IsortPermutation")) then

      if (rvector%h_IsortPermutation == ST_NOHANDLE) then
        ! If the permutation handle is not associated to some memory block, allocate new memory
        call storage_new1d ('lsysh5io_readVector_LocationByID','IsortPermutation',2*rvector%NEQ,&
            ST_INT,rvector%h_IsortPermutation,ST_NEWBLOCK_NOINIT)
      else
        ! Otherwise, check that the provided memory is sufficient
        call storage_getsize(rvector%h_IsortPermutation,isize)
        if (isize /= 2*rvector%NEQ) then
          print *, "lsysh5io_readVector_LocationByID: Wrong size of permutation array!"
          stop
        end if
      end if

      ! Get permutation handle
      call storage_getbase_int(rvector%h_IsortPermutation,p_Idata)
      call h5lite_get_data(vector_id,"h_IsortPermutation",p_Idata)
    end if
  end subroutine lsysh5io_readVector_LocationByID

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readVector_LocationByName(loc_id,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)          :: loc_id

    ! Name of the vector
    character(LEN=*), intent(IN)        :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorScalar), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: vector_id

    ! Determine the vector_id and call the working routine
    call h5lite_get_group(loc_id,cvectorName,vector_id)
    if (vector_id < 0) then
      print *, "lsysh5io_readVector_LocationByName: Unable to get vector_id!"
      stop
    end if
    call lsysh5io_readVector(loc_id,vector_id,rvector)
  end subroutine lsysh5io_readVector_LocationByName

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readBlockVector_RootByID(rh5file,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its unique vector_id.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)         :: rh5file

    ! Handle to the vector object
    integer(HID_T), intent(IN)         :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readBlockVector(rh5file%file_id,vector_id,rvector)
  end subroutine lsysh5io_readBlockVector_RootByID

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readBlockVector_RootByName(rh5file,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the root of the hdf5 file that
    ! is associated with the h5file descriptor. 
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! The hdf file descriptor
    type(t_h5file), intent(IN)         :: rh5file

    ! Name of the vector
    character(LEN=*), intent(IN)       :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>

    ! Call the working routine for the top-level
    call lsysh5io_readBlockVector(rh5file%file_id,cvectorName,rvector)
  end subroutine lsysh5io_readBlockVector_RootByName

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readBlockVector_LocationByID(loc_id,vector_id,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its uniquad vector_id.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)         :: loc_id

    ! Handle to the vector object
    integer(HID_T), intent(IN)         :: vector_id
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: iblock

    ! Get attributes
    call h5lite_get_attribute(vector_id,"NEQ",           rvector%NEQ)
    call h5lite_get_attribute(vector_id,"iidxFirstEntry",rvector%iidxFirstEntry)
    call h5lite_get_attribute(vector_id,"cdataType",     rvector%cdataType)
    call h5lite_set_attribute(vector_id,"bisCopy",       rvector%bisCopy)
    call h5lite_get_attribute(vector_id,"nblocks",       rvector%nblocks)
    
    ! If the data handle is not associated to some memory block, allocate new memory
    if (rvector%h_Ddata == ST_NOHANDLE) &
        call storage_new1d ('lsysh5io_readBlockVector_LocationByID','Vector',rvector%NEQ,&
        rvector%cdataType,rvector%h_Ddata,ST_NEWBLOCK_NOINIT)

    ! Allocate 1D array  with scalar vectors
    allocate(rvector%RvectorBlock(rvector%nblocks))
    
    ! Read in scalar subvectors
    do iblock=1,rvector%nblocks
      
      ! Set the handle for the scalar subvector
      rvector%RvectorBlock(iblock)%h_Ddata = rvector%h_Ddata
      
      ! Read in the scalar subvector
      call lsysh5io_readVector(vector_id,"RvectorBlock("//trim(adjustl(sys_si(iblock,6)))//")",&
          rvector%RvectorBlock(iblock))
    end do
  end subroutine lsysh5io_readBlockVector_LocationByID

  ! ***************************************************************************

!<subroutine>

  subroutine lsysh5io_readBlockVector_LocationByName(loc_id,cvectorName,rvector)

!<description>
    ! This subroutine reads the vector from the location given by loc_id.
    ! The vector is identified by its name given by cvectorName.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)         :: loc_id

    ! Name of the vector
    character(LEN=*), intent(IN)       :: cvectorName
!</input>

!<inputoutput>
    ! The vector that should be read
    type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: vector_id

    ! Determine the vector_id and call the working routine
    call h5lite_get_group(loc_id,cvectorName,vector_id)
    if (vector_id < 0) then
      print *, "lsysh5io_readBlockVector_LocationByName: Unable to get vector_id!"
      stop
    end if
    call lsysh5io_readBlockVector(loc_id,vector_id,rvector)
  end subroutine lsysh5io_readBlockVector_LocationByName
end module linearsystemh5io
