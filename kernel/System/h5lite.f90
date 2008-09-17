!##############################################################################
!# ****************************************************************************
!# <name> h5lite </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides basic routines for reading/writing simulation data in HDF5 format. 
!# This module makes use of the HDF5 core library which can be obtained from the HDF5 homepage
!# http://hdf.ncsa.uiuc.edu/HDF5/
!#
!# The following routines are available:
!#
!#
!# 1.) h5lite_init
!#     -> Initializes the HDF subsystem
!#
!# 2.) h5lite_done
!#     -> Finalizes the HDF subsystem
!#
!# 3.) h5lite_create
!#     -> Create a new h5file which is attached to a physical HDF5 file
!#
!# 4.) h5lite_release
!#     -> Release a h5file and close the attached physical HDF5 file.
!#
!#
!# 5.) h5lite_set_group = h5lite_set_group_Root /
!#                        h5lite_set_group_Location
!#     -> Set a new group with a given name
!#
!# 6.) h5lite_get_group = h5lite_get_group_Root /
!#                        h5lite_get_group_Location
!#     -> Returns the group_id for a given group name.
!#
!# 7.) h5lite_ispresent_attribute
!#     -> Returns true if the attribute is present in the file
!#
!# 8.) h5lite_set_attribute = h5lite_set_attribute_int    / 
!#                            h5lite_set_attribute_single /
!#                            h5lite_set_attribute_double / 
!#                            h5lite_set_attribute_char   /
!#                            h5lite_set_attribute_logical
!#    -> Sets an attribute for a given dataset
!#
!# 9.) h5lite_get_attribute = h5lite_get_attribute_integer / 
!#                            h5lite_get_attribute_single  /
!#                            h5lite_get_attribute_double  / 
!#                            h5lite_get_attribute_char    /
!#                            h5lite_get_attribute_logical
!#    -> Gets an attribute of a given dataset
!#
!# 10.) h5lite_ispresent_data
!#      -> Returns true if the data is present in the file
!#
!# 11.) h5lite_set_data = h5lite_set_data_handle /
!#                       h5lite_set_data_int    / h5lite_set_data_int2D    /
!#                       h5lite_set_data_single / h5lite_set_data_single2D /
!#                       h5lite_set_data_double / h5lite_set_data_double2D
!#
!# 12.) h5lite_get_data = h5lite_get_data_handle /
!#                       h5lite_get_data_int    / h5lite_get_data_int2D    / 
!#                       h5lite_get_data_single / h5lite_get_data_single2D /
!#                       h5lite_get_data_double / h5lite_get_data_double2D
!#
!# </purpose>
!##############################################################################
module h5lite
  use fsystem
  use storage
  use hdf5

  implicit none

  private
  public :: t_h5file
  public :: h5lite_init
  public :: h5lite_done
  public :: h5lite_create
  public :: h5lite_release
  public :: h5lite_set_group
  public :: h5lite_get_group
  public :: h5lite_ispresent_attribute
  public :: h5lite_set_attribute
  public :: h5lite_get_attribute
  public :: h5lite_ispresent_data
  public :: h5lite_set_data
  public :: h5lite_get_data
  
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>

!<constantblock description="Global file operations">

  ! Open an existing file readonly
  integer, parameter, public :: H5LITE_OPEN_READONLY    = 1

  ! Open an existing file 
  integer, parameter, public :: H5LITE_OPEN_READWRITE   = 2

  ! Create a new file which must not exist
  integer, parameter, public :: H5LITE_CREATE_NEW       = 3

  ! Create a new file and overwrite existing file
  integer, parameter, public :: H5LITE_CREATE_OVERWRITE = 4
!</constantblock>

!</constants>

  integer, parameter :: idcl=1

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface h5lite_set_group
    module procedure h5lite_set_group_Root
    module procedure h5lite_set_group_Location
  end interface

  interface h5lite_get_group
    module procedure h5lite_get_group_Root
    module procedure h5lite_get_group_Location
  end interface

  interface h5lite_set_attribute
    module procedure h5lite_set_attribute_int
    module procedure h5lite_set_attribute_single
    module procedure h5lite_set_attribute_double
    module procedure h5lite_set_attribute_char
    module procedure h5lite_set_attribute_logical
  end interface

  interface h5lite_get_attribute
    module procedure h5lite_get_attribute_int
    module procedure h5lite_get_attribute_single
    module procedure h5lite_get_attribute_double
    module procedure h5lite_get_attribute_char
    module procedure h5lite_get_attribute_logical
  end interface

  interface h5lite_set_data
    module procedure h5lite_set_data_handle
    module procedure h5lite_set_data_int
    module procedure h5lite_set_data_int2D
    module procedure h5lite_set_data_single
    module procedure h5lite_set_data_single2D
    module procedure h5lite_set_data_double
    module procedure h5lite_set_data_double2D
  end interface

  interface h5lite_get_data
    module procedure h5lite_get_data_handle
    module procedure h5lite_get_data_int
    module procedure h5lite_get_data_int2D
    module procedure h5lite_get_data_single
    module procedure h5lite_get_data_single2D
    module procedure h5lite_get_data_double
    module procedure h5lite_get_data_double2D
  end interface

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>

!<typeblock>

  ! This data structure contains settings/parameters for h5 file handling.
  
  type t_h5file

    ! Name of H5-file
    character(SYS_STRLEN) :: filename

    ! File identifier 
    integer(HID_T)        :: file_id

    ! File access mode
    ! Valid values are: 
    !   H5F_ACC_RDWR_F   - allow read/write access to file
    !   H5F_ACC_RDONLY_F - allow read-only access to file
    !   H5F_ACC_TRUNC_F  - overwrite existing file on creating 
    !   H5F_ACC_EXCL_F   - fail if file exists on creation
    integer               :: access_mode

    ! File access property list identifier
    integer(HID_T)        :: access_plist

    ! File creation property list identifier
    integer(HID_T)        :: creation_plist    
  end type t_h5file

contains

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_init()

!<description>
  ! This subroutine initializes the HDF subsystem
!</description>
!</subroutine>

    ! local variables
    integer :: hdferr

    ! Open HDF subsystem
    call h5open_f(hdferr)
    if (hdferr /= 0) then
      print *, "h5lite_init: Unable to open HDF subsystem!"
      stop
    end if
  end subroutine h5lite_init

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_done()

!<description>
  ! This subroutine finalizes the HDF subsystem
!</description>
!</subroutine>

    ! local variables
    integer :: hdferr

    ! Close HDF subsystem
    call h5close_f(hdferr)
    if (hdferr /= 0) then
      print *, "h5lite_done: Unable to close HDF subsystem!"
      stop
    end if
  end subroutine h5lite_done

  !*****************************************************************************

!<subroutine>
  
  subroutine h5lite_create(rh5file,filename,fileaccess)

!<description>
  ! This subroutine creates a new h5file object and attaches a physical
  ! HDF5 file to it. Depending on the optional parameter fileaccess, a
  ! new file is created or an existing one is opened.
!</description>

!<input>
    ! name of the physical HDF5 file without suffix '.h5'
    character(LEN=*), intent(IN) :: filename
  
    ! OPTIONAL: file access. If this parameter is not given, then the value
    ! H5LITE_CREATE_OVERWRITE is adopted.
    ! =H5LITE_OPEN_READONLY:    Try to open an existing file from disk and allow
    !                           only read-access to this file. If the file does
    !                           not exist, the routine stops with an error.
    ! =H5LITE_OPEN_READWRITE:   Try to open an existing file from disk and allow
    !                           read- and write-access to this file. If the file
    !                           does not exist, the routine stops with an error.
    ! =H5LITE_CREATE_NEW:       Create a new file. If another file with the same
    !                           filename already exists, then the routine stops
    !                           with an error.
    ! =H5LITE_CREATE_OVERWRITE: Create a new file and overwrite any existing
    !                           file with the same filename.
    integer, intent(IN), optional :: fileaccess
!</input>

!<output>
    ! h5file descriptor
    type(t_h5file), intent(OUT) :: rh5file
!</output>
!</subroutine>
  
    ! local variables
    integer :: fileacc,hdferr
    logical :: bisHdf5
    
    ! Set file access
    if (present(fileaccess)) then
      fileacc = fileaccess
    else
      fileacc = H5LITE_CREATE_OVERWRITE
    end if

    ! What kind of file access should be used
    select case(fileacc)
    case (H5LITE_OPEN_READONLY,H5LITE_OPEN_READWRITE)

      ! Check if file exists and is HDF5 file
      call h5fis_hdf5_f(trim(adjustl(filename))//'.h5',bisHdf5,hdferr)
      if (.not.bisHdf5) then
        print *, "h5lite_create: Unable to find HDF5 file with given filename!"
        stop
      end if

      ! Initialize file descriptor
      rh5file%filename = trim(adjustl(filename))//'.h5'
      if (fileacc == H5LITE_OPEN_READONLY) then
        rh5file%access_mode = H5F_ACC_RDONLY_F
      else
        rh5file%access_mode = H5F_ACC_RDWR_F
      end if
      
      ! Open the file
      call h5fopen_f(trim(adjustl(rh5file%filename)),rh5file%access_mode,rh5file%file_id,hdferr)
      if (hdferr /= 0) then
        print *, "h5lite_create: Unable to open HDF5 file!"
        stop
      end if

      ! Get file creation and access property lists
      call h5fget_access_plist_f(rh5file%file_id,rh5file%access_plist,hdferr)
      if (hdferr /= 0) then
        print *, "h5lite_create: An error occured while getting the accesss property list!"
        stop
      end if
      call h5fget_create_plist_f(rh5file%file_id,rh5file%creation_plist,hdferr)
      if (hdferr /= 0) then
        print *, "h5lite_create: An error occured while getting the creation property list!"
        stop
      end if

    case (H5LITE_CREATE_NEW,H5LITE_CREATE_OVERWRITE)

      ! Initialize file descriptor
      rh5file%filename       = trim(adjustl(filename))//'.h5'
      rh5file%creation_plist = H5P_DEFAULT_F
      rh5file%access_plist   = H5P_DEFAULT_F
      if (fileacc == H5LITE_CREATE_NEW) then
        rh5file%access_mode = H5F_ACC_EXCL_F
      else
        rh5file%access_mode = H5F_ACC_TRUNC_F
      end if
      
      ! Create file
      call h5fcreate_f(trim(adjustl(rh5file%filename)),rh5file%access_mode,rh5file%file_id,hdferr)
      if (hdferr /= 0) then
        print *, "h5lite_create: Unable to create new HDF5 file!"
        stop
      end if

    case DEFAULT
      print *, "h5lite_create: Unsupported file operation!"
      stop      
    end select
  end subroutine h5lite_create

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_release(rh5file)

!<description>
    ! This subroutine closes an existing h5file object
!</description>

!<inputoutput>
    ! h5file descriptor
    type(t_h5file), intent(INOUT) :: rh5file
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: hdferr
    
    ! Close file
    call h5fclose_f(rh5file%file_id,hdferr)
    if (hdferr /= 0) then
      print *, "h5lite_done: Unable to close HDF5 file!"
      stop
    end if
  end subroutine h5lite_release
  
  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_group_Root(rh5file,group_name,group_id,bappendNumber)

!<description>
    ! This subroutine creates a new group with a given name and returns
    ! the handle of the group. If the optional parameter bappendNumber 
    ! is present then the group_name can be changed to group_name.XXX,
    ! whereby XXX is the next free number in the h5file.
!</description>

!<input>
    ! h5file descriptor
    type(t_h5file), intent(IN) :: rh5file

    ! group name
    character(LEN=*), intent(IN) :: group_name

    ! OPTIONAL: append number to name if already present
    logical, intent(IN), optional :: bappendNumber
!</input>

!<output>
    ! group id
    integer(HID_T), intent(OUT) :: group_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    call h5lite_set_group(rh5file%file_id,group_name,group_id,bappendNumber)
  end subroutine h5lite_set_group_Root
  
  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_group_Location(loc_id,group_name,group_id,bappendNumber)

!<description>
    ! This subroutine creates a new group with a given name and returns
    ! the handle of the group. If the optional parameter bappendNumber 
    ! is present then the group_name can be changed to group_name.XXX, 
    ! whereby XXX is the next free number at the location loc_id. Otherwise,
    ! the routine stops with an error if an object with the same name already exists.
!</description>

!<input>
    ! location id
    integer(HID_T) :: loc_id

    ! group name
    character(LEN=*), intent(IN) :: group_name

    ! OPTIONAL: append number to name if already present
    logical, intent(IN), optional :: bappendNumber
!</input>

!<output>
    ! group id
    integer(HID_T), intent(OUT) :: group_id
!</output>
!</subroutine>

    ! local variable
    integer :: hdferr,iappend
    logical :: bappend

    ! Initialize appending of numbers
    bappend = .false.
    if (present(bappendNumber)) bappend=bappendNumber

    if (bappend) then

      ! Check if the name is already present
      call h5lite_get_group(loc_id,group_name,group_id)
      
      if (group_id < 0) then
        ! If it does not exist, then create it
        call h5gcreate_f(loc_id,group_name,group_id,hdferr)
        if (hdferr /= 0) then
          print *, "h5lite_set_group_Location: Unable to create group!"
          stop
        end if

      else
        
        ! If it does exist, then append '.XXXXXX'
        do iappend=1,999999

          ! Check if the new name is already present
          call h5lite_get_group(loc_id,group_name//'.'//sys_si0(iappend,6),group_id)
          
          ! ... and create it, if it does not exist
          if (group_id < 0) then
            call h5gcreate_f(loc_id,group_name//'.'//sys_si0(iappend,6),group_id,hdferr)
            if (hdferr /= 0) then
              print *, "h5lite_set_group_Location: Unable to create group!"
              stop
            end if
            exit
          end if
        end do
      end if

    else

      ! Create new group without checking if name exists
      call h5gcreate_f(loc_id,group_name,group_id,hdferr)
      if (hdferr /= 0) then
        print *, "h5lite_set_group_Location: Unable to create group!"
        stop
      end if
    end if
  end subroutine h5lite_set_group_Location

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_group_Root(rh5file,group_name,group_id)

!<description>
    ! This subroutine determines the group_id of an existing group at
    ! the root level of the h5file object. If a group with the given
    ! name does not exist, then a negative value is returned.
!</description>

!<input>
    ! h5file descriptor
    type(t_h5file), intent(IN) :: rh5file

    ! group name
    character(LEN=*), intent(IN) :: group_name
!</input>

!<output>
    ! group id
    integer(HID_T), intent(OUT) :: group_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    call h5lite_get_group(rh5file%file_id,group_name,group_id)
  end subroutine h5lite_get_group_Root

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_group_Location(loc_id,group_name,group_id)

!<description>
    ! This subroutine determines the group_id of an existing group at
    ! the location specified by loc_id. If a group with the given name
    ! does not exist, then a negative value is returned.
!</description>

!<input>
    ! location id
    integer(HID_T), intent(IN)   :: loc_id

    ! group name
    character(LEN=*), intent(IN) :: group_name
!</input>

!<output>
    ! group id
    integer(HID_T), intent(OUT)  :: group_id
!</output>
!</subroutine>

    ! local variables
    character(SYS_STRLEN) :: root_name,obj_name
    integer(SIZE_T) :: name_size
    integer :: obj_type,hdferr,imember,nmembers

    ! Initialize the group_id with negative value
    group_id = -1

    ! What type of object is loc_id?
    call h5iget_type_f(loc_id, obj_type, hdferr)
    if (hdferr /= 0) then
      print *, "h5lite_get_group_Location: Unable to determine type!"
      stop
    end if

    ! Determine name of location
    if (obj_type == H5I_FILE_F) then
      root_name="/"
    elseif (obj_type == H5I_GROUP_F) then
      call h5iget_name_f(loc_id,root_name,int(SYS_STRLEN,SIZE_T),name_size,hdferr)
      if (hdferr /= 0) then
        print *, "h5lite_get_group_Location: Unable to determine name!"
        stop
      end if
    else
      print *, "h5lite_get_group_Location: Invalid object ID passed!"
      stop
    end if
    
    ! Ok, we are now able to determine the number of objects stores in root
    call h5gn_members_f(loc_id,trim(adjustl(root_name)),nmembers,hdferr)
    if (hdferr /= 0) then
      print *, "h5lite_get_group_Location: Unable to determine number of objects!"
      stop
    end if

    ! Loop over all members (recall, HDF5 is written in C) 
    do imember=0,nmembers-1
      
      ! Retrieve information about current member
      call h5gget_obj_info_idx_f(loc_id,trim(adjustl(root_name)),imember,&
          obj_name,obj_type,hdferr)
      if (hdferr /= 0) then
        print *, "h5lite_get_group_Location: Unable to retrieve information!"
        stop
      end if

      ! Check if current object corresponds to a group object with correct name
      if (obj_type == H5G_GROUP_F .and. trim(adjustl(obj_name)) == group_name) then
        
        ! The desired group exists an can be opened
        call h5gopen_f(loc_id,group_name,group_id,hdferr)
        if (hdferr /= 0) then
          print *, "h5lite_get_group_Location: Unable to open group!"
          stop
        end if
      end if
    end do
  end subroutine h5lite_get_group_Location

  !*****************************************************************************

!<function>

  function h5lite_ispresent_attribute(set_id,attr_name) result(ispresent)

!<description>
    ! This function tests if a given attribute is present at the data set
!</description>


!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id
    
    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name
!</input>

!<result>

    ! presence of the attribute
    logical :: ispresent
!</result>
!</function>

    ! local variables
    integer(HID_T) :: attr_id
    integer :: iattr,nattrs,hdferr
    character(LEN=SYS_STRLEN) :: buf

    ! Initialization
    ispresent=.false.
    
    ! Get number of attributes
    call h5aget_num_attrs_f(set_id,nattrs,hdferr)
    if (hdferr /= 0) then
      print *, "h5lite_ispresent_attribute: Unable to determine number of arguments!"
      stop
    end if

    ! Iterate over attributes
    do iattr=0,nattrs-1
      call h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      call h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)

      ! Check if attribute is the desired one
      if (trim(adjustl(buf)) == trim(adjustl(attr_name))) then
        ispresent=.true.
        call h5aclose_f(attr_id,hdferr)
        return
      end if

      ! Close attribute
      call h5aclose_f(attr_id,hdferr)
    end do
  end function h5lite_ispresent_attribute

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_attribute_int(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets an integer attribute to a given data set
!</description>

!<input>
    
    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name

    ! value of attribut
    integer, intent(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id,space_id
    integer :: iattr,nattrs,hdferr
    character(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    call h5aget_num_attrs_f(set_id,nattrs,hdferr)
    
    ! Iterate over attributes
    do iattr=0,nattrs-1
      call h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      call h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)
      
      ! Check if attribute is the desired one
      if (trim(adjustl(buf))==trim(adjustl(attr_name))) then
        call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_val,(/1_HSIZE_T/),hdferr)
        call h5aclose_f(attr_id,hdferr)
        return
      end if
      
      ! Close attribute at index
      call h5aclose_f(attr_id,hdferr)
    end do
    
    ! Create new attribute
    call h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    call h5acreate_f(set_id,attr_name,H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_val,(/1_HSIZE_T/),hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)
  end subroutine h5lite_set_attribute_int

  !*****************************************************************************

!<subroutine>
  
  subroutine h5lite_set_attribute_single(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a single attribute to a given data set
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name

    ! value of attribute
    real(SP), intent(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id,space_id
    integer :: iattr,nattrs,hdferr
    character(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    call h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    do iattr=0,nattrs-1
      call h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      call h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)
      
      ! Check if attribute is the desired one
      if (trim(adjustl(buf))==trim(adjustl(attr_name))) then
        call h5awrite_f(attr_id,H5T_NATIVE_REAL,attr_val,(/1_HSIZE_T/),hdferr)
        call h5aclose_f(attr_id,hdferr)
        return
      end if

      ! Close attribute at index
      call h5aclose_f(attr_id,hdferr)
    end do

    ! Create new attribute
    call h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    call h5acreate_f(set_id,attr_name,H5T_NATIVE_REAL,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,H5T_NATIVE_REAL,attr_val,(/1_HSIZE_T/),hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)
  end subroutine h5lite_set_attribute_single
  
  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_attribute_double(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a double attribute to a given data set
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id
    
    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name

    ! value of attribute
    real(DP), intent(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id,space_id
    integer :: iattr,nattrs,hdferr
    character(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    call h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    do iattr=0,nattrs-1
      call h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      call h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)

      ! Check if attribute is the desired one
      if (trim(adjustl(buf))==trim(adjustl(attr_name))) then
        call h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,attr_val,(/1_HSIZE_T/),hdferr)
        call h5aclose_f(attr_id,hdferr)
        return
      end if

      ! Close attribute at index
      call h5aclose_f(attr_id,hdferr)
    end do
    
    ! Create new attribute
    call h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    call h5acreate_f(set_id,attr_name,H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,attr_val,(/1_HSIZE_T/),hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)
  end subroutine h5lite_set_attribute_double

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_attribute_char(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a character attribute to a given data set
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name

    ! value of attribute
    character(LEN=*), intent(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id,space_id,atype_id
    integer :: iattr,nattrs,hdferr
    character(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    call h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    do iattr=0,nattrs-1
      call h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      call h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)
      
      ! Check if attribute is the desired one
      if (trim(adjustl(buf))==trim(adjustl(attr_name))) then
        call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
        call h5tset_size_f(atype_id,len(attr_val),hdferr)
        call h5awrite_f(attr_id,atype_id,attr_val,(/1_HSIZE_T/),hdferr)
        call h5aclose_f(attr_id,hdferr)
        return
      end if

      ! Close attribute at index
      call h5aclose_f(attr_id,hdferr)
    end do
    
    ! Create new attribute
    call h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    call h5tset_size_f(atype_id,len(attr_val),hdferr)
    call h5acreate_f(set_id,attr_name,atype_id,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,atype_id,attr_val,(/1_HSIZE_T/),hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)
  end subroutine h5lite_set_attribute_char

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_attribute_logical(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a logical attribute to a given data set.
    ! Note that HDF5 does not support "logicals" by default. As a workaround,
    ! logicals are stored as a single character with values "T" or "F".
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name

    ! value of attribute
    logical, intent(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id,space_id,atype_id
    integer :: iattr,nattrs,hdferr
    character(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    call h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    do iattr=0,nattrs-1
      call h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      call h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)

      ! Check if attribute is the desired one
      if (trim(adjustl(buf))==trim(adjustl(attr_name))) then
        call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
        call h5tset_size_f(atype_id,1,hdferr)
        if (attr_val) then
          call h5awrite_f(attr_id,atype_id,"T",(/1_HSIZE_T/),hdferr)
        else
          call h5awrite_f(attr_id,atype_id,"F",(/1_HSIZE_T/),hdferr)
        end if
        call h5aclose_f(attr_id,hdferr)
        return
      end if

      ! Close attribute at index
      call h5aclose_f(attr_id,hdferr)
    end do
    
    ! Create new attribute
    call h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    call h5tset_size_f(atype_id,1,hdferr)
    call h5acreate_f(set_id,attr_name,atype_id,space_id,attr_id,hdferr)
    if (attr_val) then
      call h5awrite_f(attr_id,atype_id,"T",(/1_HSIZE_T/),hdferr)
    else
      call h5awrite_f(attr_id,atype_id,"F",(/1_HSIZE_T/),hdferr)
    end if
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)
  end subroutine h5lite_set_attribute_logical

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_attribute_int(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets an integer attribute from a given data set
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id
    
    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    integer, intent(OUT) :: attr_val
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id
    integer :: hdferr
    
    call h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    call h5aread_f(attr_id,H5T_NATIVE_INTEGER,attr_val,(/1_HSIZE_T/),hdferr) 
    call h5aclose_f(attr_id,hdferr)
  end subroutine h5lite_get_attribute_int

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_attribute_single(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a single attribute from a given data set
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    real(SP), intent(OUT) :: attr_val
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id
    integer :: hdferr
    
    call h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    call h5aread_f(attr_id,H5T_NATIVE_REAL,attr_val,(/1_HSIZE_T/),hdferr) 
    call h5aclose_f(attr_id,hdferr)
  end subroutine h5lite_get_attribute_single

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_attribute_double(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a double attribute from a given data set
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    real(DP), intent(OUT) :: attr_val
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: attr_id
    integer :: hdferr
    
    call h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    call h5aread_f(attr_id,H5T_NATIVE_DOUBLE,attr_val,(/1_HSIZE_T/),hdferr) 
    call h5aclose_f(attr_id,hdferr)
  end subroutine h5lite_get_attribute_double

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_attribute_char(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a character attribute from a given data set
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    character(LEN=*), intent(OUT) :: attr_val
!</output>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: attr_id,atype_id
    integer :: hdferr
    
    call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    call h5tset_size_f(atype_id,len(attr_val),hdferr)
    call h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    call h5aread_f(attr_id,atype_id,attr_val,(/1_HSIZE_T/),hdferr) 
    call h5aclose_f(attr_id,hdferr)
  end subroutine h5lite_get_attribute_char

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_attribute_logical(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a logical attribute from a given data set.
    ! Note that HDF5 does not support "logicals" by default. As a workaround,
    ! logicals are stored as a single character with values "T" or "F".
!</description>

!<input>

    ! data set
    integer(HID_T), intent(IN) :: set_id

    ! name of attribute
    character(LEN=*), intent(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    logical, intent(OUT) :: attr_val
!</output>
!</subroutine>
    
    ! local variables
    character(LEN=1) :: tmp_val
    integer(HID_T) :: attr_id,atype_id
    integer :: hdferr
    
    call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    call h5tset_size_f(atype_id,1,hdferr)
    call h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    call h5aread_f(attr_id,atype_id,tmp_val,(/1_HSIZE_T/),hdferr) 
    call h5aclose_f(attr_id,hdferr)
    
    select case(tmp_val)
    case("T","t")
      attr_val=.true.
    case("F","f")
      attr_val=.false.
    case DEFAULT
      print *, "h5lite_get_attribute_logical: Invalid attribute value!"
      stop
    end select
  end subroutine h5lite_get_attribute_logical

  !*****************************************************************************

!<function>

  function h5lite_ispresent_data(loc_id,name) result(ispresent)

!<description>
    ! This function tests if a given data set is present at location
!</description>


!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id
    
    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<result>

    ! presence of data
    logical :: ispresent
!</result>
!</function>

    ! local variables
    character(SYS_STRLEN) :: root_name,obj_name
    integer(HID_T) :: dataset_id,name_size
    integer :: obj_type,imember,nmembers,hdferr

    ! Initialization
    ispresent=.false.

    ! Determine name of location and number of members
    call h5iget_name_f(loc_id,root_name,int(SYS_STRLEN,SIZE_T),name_size,hdferr)
    call h5gn_members_f(loc_id,trim(adjustl(root_name)),nmembers,hdferr)

    ! Loop over members
    do imember=0,nmembers-1
      call h5gget_obj_info_idx_f(loc_id,trim(adjustl(root_name)),imember,&
          obj_name,obj_type,hdferr)
      
      if ((obj_type == H5G_DATASET_F .or. obj_type == H5G_LINK_F) .and.&
          trim(adjustl(obj_name)) == trim(adjustl(name))) then
        ispresent=.true.
        return
      end if
    end do
  end function h5lite_ispresent_data

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_data_handle(loc_id,name,ihandle)

!<description>
    ! This subroutine sets the data associated to a handle to a given group
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! handle to data
    integer, intent(IN) :: ihandle

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>
!</subroutine>

    ! local variables
    integer :: idatatype
    integer :: idimension
    integer, dimension(:,:), pointer  :: data_int2D
    integer, dimension(:), pointer    :: data_int
    real(DP), dimension(:,:), pointer :: data_double2D
    real(DP), dimension(:), pointer   :: data_double
    real(SP), dimension(:,:), pointer :: data_single2D
    real(SP), dimension(:), pointer   :: data_single
    
    ! Get dimension and type of data
    call storage_getdatatype(ihandle,idatatype)
    call storage_getdimension(ihandle,idimension)

    ! Which dimension are we?
    select case(idimension)
    case (1)
       
       ! Which datatype are we?
       select case (idatatype)
       case (ST_INT)
          call storage_getbase_int(ihandle,data_int)
          call h5lite_set_data_int(loc_id,name,data_int)

       case (ST_SINGLE)
          call storage_getbase_single(ihandle,data_single)
          call h5lite_set_data_single(loc_id,name,data_single)

       case (ST_DOUBLE)
          call storage_getbase_double(ihandle,data_double)
          call h5lite_set_data_double(loc_id,name,data_double)

       case DEFAULT
          print *, "(EE) h5lite_set_data_handle: invalid data type"
          stop
       end select

    case (2)

       ! Which datatype are we?
       select case (idatatype)
       case (ST_INT)
          call storage_getbase_int2D(ihandle,data_int2D)
          call h5lite_set_data_int2D(loc_id,name,data_int2D)

       case (ST_SINGLE)
          call storage_getbase_single2D(ihandle,data_single2D)
          call h5lite_set_data_single2D(loc_id,name,data_single2D)

       case (ST_DOUBLE)
          call storage_getbase_double2D(ihandle,data_double2D)
          call h5lite_set_data_double2D(loc_id,name,data_double2D)

       case DEFAULT
          print *, "(EE) h5lite_set_data_handle: invalid data type"
          stop
       end select

    case DEFAULT
       print *, "(EE) h5lite_set_data_handle: dimension exceeds 2"
       stop
    end select
  end subroutine h5lite_set_data_handle

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_data_int(loc_id,name,data)

!<description>
    ! This subroutine sets a 1D integer array to a given group
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id
    
    ! 1D integer array
    integer, dimension(:), intent(IN) :: data

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>
!</subroutine>

    ! local variables
    integer(HID_T) :: set_id,space_id,cpl
    integer(HSIZE_T) :: dim(1)
    integer :: hdferr
    
    dim=shape(data); if (sum(dim) == 0) return
    
    call h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    call h5sset_extent_simple_f(space_id,1,dim,dim,hdferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    call h5pset_chunk_f(cpl,1,dim,hdferr)
    call h5pset_deflate_f(cpl,idcl,hdferr)
    call h5dcreate_f(loc_id,name,H5T_NATIVE_INTEGER,space_id,set_id,hdferr,cpl)
    call h5dwrite_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
    call h5pclose_f(cpl,hdferr)
    call h5dclose_f(set_id,hdferr)
    call h5sclose_f(space_id,hdferr)   
  end subroutine h5lite_set_data_int

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_data_int2D(loc_id,name,data)

!<description>
    ! This subroutine sets a 2D integer array to a given group
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! 2D integer array
    integer, dimension(:,:), intent(IN) :: data

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: set_id,space_id,cpl
    integer(HSIZE_T) :: dim(2)
    integer :: hdferr

    dim=shape(data); if (sum(dim) == 0) return
    
    call h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    call h5sset_extent_simple_f(space_id,2,dim,dim,hdferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    call h5pset_chunk_f(cpl,2,dim,hdferr)
    call h5pset_deflate_f(cpl,idcl,hdferr)
    call h5dcreate_f(loc_id,name,H5T_NATIVE_INTEGER,space_id,set_id,hdferr,cpl)
    call h5dwrite_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
    call h5pclose_f(cpl,hdferr)
    call h5dclose_f(set_id,hdferr)
    call h5sclose_f(space_id,hdferr)   
  end subroutine h5lite_set_data_int2D

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_data_single(loc_id,name,data)

!<description>
    ! This subroutine sets a 1D single array to a given group
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! 1D single array
    real(SP), dimension(:), intent(IN) :: data

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: set_id,space_id,cpl
    integer(HSIZE_T) :: dim(1)
    integer :: hdferr
    
    dim=shape(data); if (sum(dim) == 0) return
    
    call h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    call h5sset_extent_simple_f(space_id,1,dim,dim,hdferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    call h5pset_chunk_f(cpl,1,dim,hdferr)
    call h5pset_deflate_f(cpl,idcl,hdferr)
    call h5dcreate_f(loc_id,name,H5T_NATIVE_REAL,space_id,set_id,hdferr,cpl)
    call h5dwrite_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
    call h5pclose_f(cpl,hdferr)
    call h5dclose_f(set_id,hdferr)
    call h5sclose_f(space_id,hdferr)   
  end subroutine h5lite_set_data_single

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_data_single2D(loc_id,name,data)

!<description>
    ! This subroutine sets a 2D single array to a given group
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! 2D single array
    real(SP), dimension(:,:), intent(IN) :: data

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: set_id,space_id,cpl
    integer(HSIZE_T) :: dim(2)
    integer :: hdferr

    dim=shape(data); if (sum(dim) == 0) return

    call h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    call h5sset_extent_simple_f(space_id,2,dim,dim,hdferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    call h5pset_chunk_f(cpl,2,dim,hdferr)
    call h5pset_deflate_f(cpl,idcl,hdferr)
    call h5dcreate_f(loc_id,name,H5T_NATIVE_REAL,space_id,set_id,hdferr,cpl)
    call h5dwrite_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
    call h5pclose_f(cpl,hdferr)
    call h5dclose_f(set_id,hdferr)
    call h5sclose_f(space_id,hdferr)   
  end subroutine h5lite_set_data_single2D

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_data_double(loc_id,name,data)

!<description>
    ! This subroutine sets a 1D double array to a given group
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! 1D double array
    real(DP), dimension(:), intent(IN) :: data

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>
!</subroutine>

    ! local variables
    integer(HID_T) :: set_id,space_id,cpl
    integer(HSIZE_T) :: dim(1)
    integer :: hdferr

    dim=shape(data); if (sum(dim) == 0) return

    call h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    call h5sset_extent_simple_f(space_id,1,dim,dim,hdferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    call h5pset_chunk_f(cpl,1,dim,hdferr)
    call h5pset_deflate_f(cpl,idcl,hdferr)
    call h5dcreate_f(loc_id,name,H5T_NATIVE_DOUBLE,space_id,set_id,hdferr,cpl)
    call h5dwrite_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
    call h5pclose_f(cpl,hdferr)
    call h5dclose_f(set_id,hdferr)
    call h5sclose_f(space_id,hdferr)   
  end subroutine h5lite_set_data_double

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_set_data_double2D(loc_id,name,data)

!<description>
    ! This subroutine sets a 2D double array to a given group
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! 2D double array
    real(DP), dimension(:,:), intent(IN) :: data

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: set_id,space_id,cpl
    integer(HSIZE_T) :: dim(2)
    integer :: hdferr
    
    dim=shape(data); if (sum(dim) == 0) return
    
    call h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    call h5sset_extent_simple_f(space_id,2,dim,dim,hdferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    call h5pset_chunk_f(cpl,2,dim,hdferr)
    call h5pset_deflate_f(cpl,idcl,hdferr)
    call h5dcreate_f(loc_id,name,H5T_NATIVE_DOUBLE,space_id,set_id,hdferr,cpl)
    call h5dwrite_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
    call h5pclose_f(cpl,hdferr)
    call h5dclose_f(set_id,hdferr)
    call h5sclose_f(space_id,hdferr)   
  end subroutine h5lite_set_data_double2D

  !*****************************************************************************

!<subroutine>
  
  subroutine h5lite_get_data_handle(loc_id,name,ihandle)

!<description>
    ! This subroutine gets a stored data and associates it to a handle
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id
    
    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<output>
    
    ! handle to data
    integer, intent(OUT) :: ihandle
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: idata
    real(SP), dimension(:), pointer :: sdata
    real(DP), dimension(:), pointer :: ddata
    integer, dimension(:,:), pointer :: idata2D
    real(SP), dimension(:,:), pointer :: sdata2D
    real(DP), dimension(:,:), pointer :: ddata2D
    
    integer(HID_T) :: dataset_id,datatype_id,dataspace_id
    integer(HSIZE_T), dimension(1) :: dim1D,maxdim1D
    integer(HSIZE_T), dimension(2) :: dim2D,maxdim2D
    integer :: ndims,hdferr
    logical :: istype
    
    ! Get the data set
    call h5dopen_f(loc_id,name,dataset_id,hdferr)
    
    ! Get the data type
    call h5dget_type_f(dataset_id,datatype_id,hdferr)
    
    ! Get the dataspace
    call h5dget_space_f(dataset_id,dataspace_id,hdferr)
    
    ! Get number of dimensions
    call h5sget_simple_extent_ndims_f(dataspace_id,ndims,hdferr)
    
    select case(ndims)
    case (1)
      ! Get dimensions
      call h5sget_simple_extent_dims_f(dataspace_id,dim1D,maxdim1D,hdferr)

      ! Are we integer?
      call h5tequal_f(datatype_id,H5T_NATIVE_INTEGER,istype,hdferr)
      if (istype) then
        call h5dclose_f(dataset_id,hdferr)
        call storage_new('h5lite_get_data_handle',trim(adjustl(name)), &
             int(dim1D(1)),ST_INT,ihandle,ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(ihandle,idata)
        call h5lite_get_data_int(loc_id,name,idata)
        return
      end if
      
      ! Are we single?
      call h5tequal_f(datatype_id,H5T_NATIVE_REAL,istype,hdferr)
      if (istype) then
        call h5dclose_f(dataset_id,hdferr)
        call storage_new('h5lite_get_data_handle',trim(adjustl(name)), &
             int(dim1D(1)),ST_SINGLE,ihandle,ST_NEWBLOCK_NOINIT)
        call storage_getbase_single(ihandle,sdata)
        call h5lite_get_data_single(loc_id,name,sdata)
        return
      end if
      
      ! Are we double?
      call h5tequal_f(datatype_id,H5T_NATIVE_DOUBLE,istype,hdferr)
      if (istype) then
        call h5dclose_f(dataset_id,hdferr)
        call storage_new('h5lite_get_data_handle',trim(adjustl(name)), &
             int(dim1D(1)),ST_DOUBLE,ihandle,ST_NEWBLOCK_NOINIT)
        call storage_getbase_double(ihandle,ddata)
        call h5lite_get_data_double(loc_id,name,ddata)
        return
      end if
      
      ! We are the wrong datatype
      print *, "(EE) h5lite_get_data_handle: unknown datatype"
      stop
      
    case (2)
      ! Get dimensions
      call h5sget_simple_extent_dims_f(dataspace_id,dim2D,maxdim2D,hdferr)

      ! Are we integer?
      call h5tequal_f(datatype_id,H5T_NATIVE_INTEGER,istype,hdferr)
      if (istype) then
        call h5dclose_f(dataset_id,hdferr)
        call storage_new('h5lite_get_data_handle',trim(adjustl(name)), &
             int(dim2D),ST_INT,ihandle,ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2D(ihandle,idata2D)
        call h5lite_get_data_int2D(loc_id,name,idata2D)
        return
      end if
      
      ! Are we single?
      call h5tequal_f(datatype_id,H5T_NATIVE_REAL,istype,hdferr)
      if (istype) then
        call h5dclose_f(dataset_id,hdferr)
        call storage_new('h5lite_get_data_handle',trim(adjustl(name)), &
             int(dim2D),ST_SINGLE,ihandle,ST_NEWBLOCK_NOINIT)
        call storage_getbase_single2D(ihandle,sdata2D)
        call h5lite_get_data_single2D(loc_id,name,sdata2D)
        return
      end if
      
      ! Are we double?
      call h5tequal_f(datatype_id,H5T_NATIVE_DOUBLE,istype,hdferr)
      if (istype) then
        call h5dclose_f(dataset_id,hdferr)
        call storage_new('h5lite_get_data_handle',trim(adjustl(name)), &
             int(dim2D),ST_DOUBLE,ihandle,ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2D(ihandle,ddata2D)
        call h5lite_get_data_double2D(loc_id,name,ddata2D)
        return
      end if
      
      ! We are the wrong datatype
      print *, "(EE) h5lite_get_data_handle: unknown datatype"
      stop
      
    case DEFAULT
      print *, "(EE) h5lite_get_data_handle: number of dimensions exceeds 2"
      stop
    end select
  end subroutine h5lite_get_data_handle

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_data_int(loc_id,name,data)

!<description>
    ! This subroutine gets a 1D integer array
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<output>
    ! 1D integer array
    integer, dimension(:), intent(OUT) :: data
!</output>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: set_id
    integer(HSIZE_T) :: dim(1)
    integer :: hdferr
    
    dim=shape(data)
    call h5dopen_f(loc_id,name,set_id,hdferr)
    if (hdferr == 0) then
      call h5dread_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
      call h5dclose_f(set_id,hdferr)
    end if
  end subroutine h5lite_get_data_int
  
  !*****************************************************************************

!<subroutine>

   subroutine h5lite_get_data_int2D(loc_id,name,data)

!<description>
    ! This subroutine gets a 2D integer array
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id
    
    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<output>
    ! 2D integer array
    integer, dimension(:,:), intent(OUT) :: data
!</output>
!</subroutine>
    
    ! local variables
    integer(HID_T) :: set_id
    integer(HSIZE_T) :: dim(2)
    integer :: hdferr
    
    dim=shape(data)
    call h5dopen_f(loc_id,name,set_id,hdferr)
    if (hdferr == 0) then
      call h5dread_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
      call h5dclose_f(set_id,hdferr)
    end if
  end subroutine h5lite_get_data_int2D
  
  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_data_single(loc_id,name,data)

!<description>
    ! This subroutine gets a 1D single array
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<output>
    ! 1D single array
    real(SP), dimension(:), intent(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: set_id
    integer(HSIZE_T) :: dim(1)
    integer :: hdferr
    
    dim=shape(data)
    call h5dopen_f(loc_id,name,set_id,hdferr)
    if (hdferr == 0) then
      call h5dread_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
      call h5dclose_f(set_id,hdferr)
    end if
  end subroutine h5lite_get_data_single

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_data_single2D(loc_id,name,data)

!<description>
    ! This subroutine gets a 2D single array
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<output>
    ! 2D single array
    real(SP), dimension(:,:), intent(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: set_id
    integer(HSIZE_T) :: dim(2)
    integer :: hdferr
    
    dim=shape(data)
    call h5dopen_f(loc_id,name,set_id,hdferr)
    if (hdferr == 0) then
      call h5dread_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
      call h5dclose_f(set_id,hdferr)
    end if
  end subroutine h5lite_get_data_single2D
  
  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_data_double(loc_id,name,data)

!<description>
    ! This subroutine gets a 1D double array
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<output>
    real(DP), dimension(:), intent(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: set_id
    integer(HSIZE_T) :: dim(1)
    integer :: hdferr
    
    dim=shape(data)
    call h5dopen_f(loc_id,name,set_id,hdferr)
    if (hdferr == 0) then
      call h5dread_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
      call h5dclose_f(set_id,hdferr)
    end if
  end subroutine h5lite_get_data_double

  !*****************************************************************************

!<subroutine>

  subroutine h5lite_get_data_double2D(loc_id,name,data)

!<description>
    ! This subroutine gets a 2D double array
!</description>

!<input>

    ! location identifier
    integer(HID_T), intent(IN) :: loc_id

    ! name of data
    character(LEN=*), intent(IN) :: name
!</input>

!<output>
    real(DP), dimension(:,:), intent(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    integer(HID_T) :: set_id
    integer(HSIZE_T) :: dim(2)
    integer :: hdferr
    
    dim=shape(data)
    call h5dopen_f(loc_id,name,set_id,hdferr)
    if (hdferr == 0) then
      call h5dread_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
      call h5dclose_f(set_id,hdferr)
    end if
  end subroutine h5lite_get_data_double2D
end module h5lite
