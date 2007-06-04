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
MODULE h5lite
  USE fsystem
  USE storage
  USE hdf5

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_h5file
  PUBLIC :: h5lite_init
  PUBLIC :: h5lite_done
  PUBLIC :: h5lite_create
  PUBLIC :: h5lite_release
  PUBLIC :: h5lite_set_group
  PUBLIC :: h5lite_get_group
  PUBLIC :: h5lite_ispresent_attribute
  PUBLIC :: h5lite_set_attribute
  PUBLIC :: h5lite_get_attribute
  PUBLIC :: h5lite_ispresent_data
  PUBLIC :: h5lite_set_data
  PUBLIC :: h5lite_get_data
  
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>

!<constantblock description="Global file operations">

  ! Open an existing file readonly
  INTEGER, PARAMETER, PUBLIC :: H5LITE_OPEN_READONLY    = 1

  ! Open an existing file 
  INTEGER, PARAMETER, PUBLIC :: H5LITE_OPEN_READWRITE   = 2

  ! Create a new file which must not exist
  INTEGER, PARAMETER, PUBLIC :: H5LITE_CREATE_NEW       = 3

  ! Create a new file and overwrite existing file
  INTEGER, PARAMETER, PUBLIC :: H5LITE_CREATE_OVERWRITE = 4
!</constantblock>

!</constants>

  INTEGER, PARAMETER :: idcl=1

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  INTERFACE h5lite_set_group
    MODULE PROCEDURE h5lite_set_group_Root
    MODULE PROCEDURE h5lite_set_group_Location
  END INTERFACE

  INTERFACE h5lite_get_group
    MODULE PROCEDURE h5lite_get_group_Root
    MODULE PROCEDURE h5lite_get_group_Location
  END INTERFACE

  INTERFACE h5lite_set_attribute
    MODULE PROCEDURE h5lite_set_attribute_int
    MODULE PROCEDURE h5lite_set_attribute_single
    MODULE PROCEDURE h5lite_set_attribute_double
    MODULE PROCEDURE h5lite_set_attribute_char
    MODULE PROCEDURE h5lite_set_attribute_logical
  END INTERFACE

  INTERFACE h5lite_get_attribute
    MODULE PROCEDURE h5lite_get_attribute_int
    MODULE PROCEDURE h5lite_get_attribute_single
    MODULE PROCEDURE h5lite_get_attribute_double
    MODULE PROCEDURE h5lite_get_attribute_char
    MODULE PROCEDURE h5lite_get_attribute_logical
  END INTERFACE

  INTERFACE h5lite_set_data
    MODULE PROCEDURE h5lite_set_data_handle
    MODULE PROCEDURE h5lite_set_data_int
    MODULE PROCEDURE h5lite_set_data_int2D
    MODULE PROCEDURE h5lite_set_data_single
    MODULE PROCEDURE h5lite_set_data_single2D
    MODULE PROCEDURE h5lite_set_data_double
    MODULE PROCEDURE h5lite_set_data_double2D
  END INTERFACE

  INTERFACE h5lite_get_data
    MODULE PROCEDURE h5lite_get_data_handle
    MODULE PROCEDURE h5lite_get_data_int
    MODULE PROCEDURE h5lite_get_data_int2D
    MODULE PROCEDURE h5lite_get_data_single
    MODULE PROCEDURE h5lite_get_data_single2D
    MODULE PROCEDURE h5lite_get_data_double
    MODULE PROCEDURE h5lite_get_data_double2D
  END INTERFACE

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>

!<typeblock>

  ! This data structure contains settings/parameters for h5 file handling.
  
  TYPE t_h5file

    ! Name of H5-file
    CHARACTER(SYS_STRLEN) :: filename

    ! File identifier 
    INTEGER(HID_T)        :: file_id

    ! File access mode
    ! Valid values are: 
    !   H5F_ACC_RDWR_F   - allow read/write access to file
    !   H5F_ACC_RDONLY_F - allow read-only access to file
    !   H5F_ACC_TRUNC_F  - overwrite existing file on creating 
    !   H5F_ACC_EXCL_F   - fail if file exists on creation
    INTEGER               :: access_mode

    ! File access property list identifier
    INTEGER(HID_T)        :: access_plist

    ! File creation property list identifier
    INTEGER(HID_T)        :: creation_plist    
  END TYPE t_h5file

CONTAINS

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_init()

!<description>
  ! This subroutine initializes the HDF subsystem
!</description>
!</subroutine>

    ! local variables
    INTEGER :: hdferr

    ! Open HDF subsystem
    CALL h5open_f(hdferr)
    IF (hdferr /= 0) THEN
      PRINT *, "h5lite_init: Unable to open HDF subsystem!"
      STOP
    END IF
  END SUBROUTINE h5lite_init

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_done()

!<description>
  ! This subroutine finalizes the HDF subsystem
!</description>
!</subroutine>

    ! local variables
    INTEGER :: hdferr

    ! Close HDF subsystem
    CALL h5close_f(hdferr)
    IF (hdferr /= 0) THEN
      PRINT *, "h5lite_done: Unable to close HDF subsystem!"
      STOP
    END IF
  END SUBROUTINE h5lite_done

  !*****************************************************************************

!<subroutine>
  
  SUBROUTINE h5lite_create(rh5file,filename,fileaccess)

!<description>
  ! This subroutine creates a new h5file object and attaches a physical
  ! HDF5 file to it. Depending on the optional parameter fileaccess, a
  ! new file is created or an existing one is opened.
!</description>

!<input>
    ! name of the physical HDF5 file without suffix '.h5'
    CHARACTER(LEN=*), INTENT(IN) :: filename
  
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
    INTEGER, INTENT(IN), OPTIONAL :: fileaccess
!</input>

!<output>
    ! h5file descriptor
    TYPE(t_h5file), INTENT(OUT) :: rh5file
!</output>
!</subroutine>
  
    ! local variables
    INTEGER :: fileacc,hdferr
    LOGICAL :: bisHdf5
    
    ! Set file access
    IF (PRESENT(fileaccess)) THEN
      fileacc = fileaccess
    ELSE
      fileacc = H5LITE_CREATE_OVERWRITE
    END IF

    ! What kind of file access should be used
    SELECT CASE(fileacc)
    CASE (H5LITE_OPEN_READONLY,H5LITE_OPEN_READWRITE)

      ! Check if file exists and is HDF5 file
      CALL h5fis_hdf5_f(TRIM(ADJUSTL(filename))//'.h5',bisHdf5,hdferr)
      IF (.NOT.bisHdf5) THEN
        PRINT *, "h5lite_create: Unable to find HDF5 file with given filename!"
        STOP
      END IF

      ! Initialize file descriptor
      rh5file%filename = TRIM(ADJUSTL(filename))//'.h5'
      IF (fileacc == H5LITE_OPEN_READONLY) THEN
        rh5file%access_mode = H5F_ACC_RDONLY_F
      ELSE
        rh5file%access_mode = H5F_ACC_RDWR_F
      END IF
      
      ! Open the file
      CALL h5fopen_f(TRIM(ADJUSTL(rh5file%filename)),rh5file%access_mode,rh5file%file_id,hdferr)
      IF (hdferr /= 0) THEN
        PRINT *, "h5lite_create: Unable to open HDF5 file!"
        STOP
      END IF

      ! Get file creation and access property lists
      CALL h5fget_access_plist_f(rh5file%file_id,rh5file%access_plist,hdferr)
      IF (hdferr /= 0) THEN
        PRINT *, "h5lite_create: An error occured while getting the accesss property list!"
        STOP
      END IF
      CALL h5fget_create_plist_f(rh5file%file_id,rh5file%creation_plist,hdferr)
      IF (hdferr /= 0) THEN
        PRINT *, "h5lite_create: An error occured while getting the creation property list!"
        STOP
      END IF

    CASE (H5LITE_CREATE_NEW,H5LITE_CREATE_OVERWRITE)

      ! Initialize file descriptor
      rh5file%filename       = TRIM(ADJUSTL(filename))//'.h5'
      rh5file%creation_plist = H5P_DEFAULT_F
      rh5file%access_plist   = H5P_DEFAULT_F
      IF (fileacc == H5LITE_CREATE_NEW) THEN
        rh5file%access_mode = H5F_ACC_EXCL_F
      ELSE
        rh5file%access_mode = H5F_ACC_TRUNC_F
      END IF
      
      ! Create file
      CALL h5fcreate_f(TRIM(ADJUSTL(rh5file%filename)),rh5file%access_mode,rh5file%file_id,hdferr)
      IF (hdferr /= 0) THEN
        PRINT *, "h5lite_create: Unable to create new HDF5 file!"
        STOP
      END IF

    CASE DEFAULT
      PRINT *, "h5lite_create: Unsupported file operation!"
      STOP      
    END SELECT
  END SUBROUTINE h5lite_create

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_release(rh5file)

!<description>
    ! This subroutine closes an existing h5file object
!</description>

!<inputoutput>
    ! h5file descriptor
    TYPE(t_h5file), INTENT(INOUT) :: rh5file
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: hdferr
    
    ! Close file
    CALL h5fclose_f(rh5file%file_id,hdferr)
    IF (hdferr /= 0) THEN
      PRINT *, "h5lite_done: Unable to close HDF5 file!"
      STOP
    END IF
  END SUBROUTINE h5lite_release
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_group_Root(rh5file,group_name,group_id,bappendNumber)

!<description>
    ! This subroutine creates a new group with a given name and returns
    ! the handle of the group. If the optional parameter bappendNumber 
    ! is present then the group_name can be changed to group_name.XXX,
    ! whereby XXX is the next free number in the h5file.
!</description>

!<input>
    ! h5file descriptor
    TYPE(t_h5file), INTENT(IN) :: rh5file

    ! group name
    CHARACTER(LEN=*), INTENT(IN) :: group_name

    ! OPTIONAL: append number to name if already present
    LOGICAL, INTENT(IN), OPTIONAL :: bappendNumber
!</input>

!<output>
    ! group id
    INTEGER(HID_T), INTENT(OUT) :: group_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    CALL h5lite_set_group(rh5file%file_id,group_name,group_id,bappendNumber)
  END SUBROUTINE h5lite_set_group_Root
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_group_Location(loc_id,group_name,group_id,bappendNumber)

!<description>
    ! This subroutine creates a new group with a given name and returns
    ! the handle of the group. If the optional parameter bappendNumber 
    ! is present then the group_name can be changed to group_name.XXX, 
    ! whereby XXX is the next free number at the location loc_id. Otherwise,
    ! the routine stops with an error if an object with the same name already exists.
!</description>

!<input>
    ! location id
    INTEGER(HID_T) :: loc_id

    ! group name
    CHARACTER(LEN=*), INTENT(IN) :: group_name

    ! OPTIONAL: append number to name if already present
    LOGICAL, INTENT(IN), OPTIONAL :: bappendNumber
!</input>

!<output>
    ! group id
    INTEGER(HID_T), INTENT(OUT) :: group_id
!</output>
!</subroutine>

    ! local variable
    INTEGER :: hdferr,iappend
    LOGICAL :: bappend

    ! Initialize appending of numbers
    bappend = .FALSE.
    IF (PRESENT(bappendNumber)) bappend=bappendNumber

    IF (bappend) THEN

      ! Check if the name is already present
      CALL h5lite_get_group(loc_id,group_name,group_id)
      
      IF (group_id < 0) THEN
        ! If it does not exist, then create it
        CALL h5gcreate_f(loc_id,group_name,group_id,hdferr)
        IF (hdferr /= 0) THEN
          PRINT *, "h5lite_set_group_Location: Unable to create group!"
          STOP
        END IF

      ELSE
        
        ! If it does exist, then append '.XXXXXX'
        DO iappend=1,999999

          ! Check if the new name is already present
          CALL h5lite_get_group(loc_id,group_name//'.'//sys_si0(iappend,6),group_id)
          
          ! ... and create it, if it does not exist
          IF (group_id < 0) THEN
            CALL h5gcreate_f(loc_id,group_name//'.'//sys_si0(iappend,6),group_id,hdferr)
            IF (hdferr /= 0) THEN
              PRINT *, "h5lite_set_group_Location: Unable to create group!"
              STOP
            END IF
            EXIT
          END IF
        END DO
      END IF

    ELSE

      ! Create new group without checking if name exists
      CALL h5gcreate_f(loc_id,group_name,group_id,hdferr)
      IF (hdferr /= 0) THEN
        PRINT *, "h5lite_set_group_Location: Unable to create group!"
        STOP
      END IF
    END IF
  END SUBROUTINE h5lite_set_group_Location

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_group_Root(rh5file,group_name,group_id)

!<description>
    ! This subroutine determines the group_id of an existing group at
    ! the root level of the h5file object. If a group with the given
    ! name does not exist, then a negative value is returned.
!</description>

!<input>
    ! h5file descriptor
    TYPE(t_h5file), INTENT(IN) :: rh5file

    ! group name
    CHARACTER(LEN=*), INTENT(IN) :: group_name
!</input>

!<output>
    ! group id
    INTEGER(HID_T), INTENT(OUT) :: group_id
!</output>
!</subroutine>

    ! Call the working routine for the top-level
    CALL h5lite_get_group(rh5file%file_id,group_name,group_id)
  END SUBROUTINE h5lite_get_group_Root

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_group_Location(loc_id,group_name,group_id)

!<description>
    ! This subroutine determines the group_id of an existing group at
    ! the location specified by loc_id. If a group with the given name
    ! does not exist, then a negative value is returned.
!</description>

!<input>
    ! location id
    INTEGER(HID_T), INTENT(IN)   :: loc_id

    ! group name
    CHARACTER(LEN=*), INTENT(IN) :: group_name
!</input>

!<output>
    ! group id
    INTEGER(HID_T), INTENT(OUT)  :: group_id
!</output>
!</subroutine>

    ! local variables
    CHARACTER(SYS_STRLEN) :: root_name,obj_name
    INTEGER(SIZE_T) :: name_size
    INTEGER :: obj_type,hdferr,imember,nmembers

    ! Initialize the group_id with negative value
    group_id = -1

    ! What type of object is loc_id?
    CALL h5iget_type_f(loc_id, obj_type, hdferr)
    IF (hdferr /= 0) THEN
      PRINT *, "h5lite_get_group_Location: Unable to determine type!"
      STOP
    END IF

    ! Determine name of location
    IF (obj_type == H5I_FILE_F) THEN
      root_name="/"
    ELSEIF (obj_type == H5I_GROUP_F) THEN
      CALL h5iget_name_f(loc_id,root_name,INT(SYS_STRLEN,SIZE_T),name_size,hdferr)
      IF (hdferr /= 0) THEN
        PRINT *, "h5lite_get_group_Location: Unable to determine name!"
        STOP
      END IF
    ELSE
      PRINT *, "h5lite_get_group_Location: Invalid object ID passed!"
      STOP
    END IF
    
    ! Ok, we are now able to determine the number of objects stores in root
    CALL h5gn_members_f(loc_id,TRIM(ADJUSTL(root_name)),nmembers,hdferr)
    IF (hdferr /= 0) THEN
      PRINT *, "h5lite_get_group_Location: Unable to determine number of objects!"
      STOP
    END IF

    ! Loop over all members (recall, HDF5 is written in C) 
    DO imember=0,nmembers-1
      
      ! Retrieve information about current member
      CALL h5gget_obj_info_idx_f(loc_id,TRIM(ADJUSTL(root_name)),imember,&
          obj_name,obj_type,hdferr)
      IF (hdferr /= 0) THEN
        PRINT *, "h5lite_get_group_Location: Unable to retrieve information!"
        STOP
      END IF

      ! Check if current object corresponds to a group object with correct name
      IF (obj_type == H5G_GROUP_F .AND. TRIM(ADJUSTL(obj_name)) == group_name) THEN
        
        ! The desired group exists an can be opened
        CALL h5gopen_f(loc_id,group_name,group_id,hdferr)
        IF (hdferr /= 0) THEN
          PRINT *, "h5lite_get_group_Location: Unable to open group!"
          STOP
        END IF
      END IF
    END DO
  END SUBROUTINE h5lite_get_group_Location

  !*****************************************************************************

!<function>

  FUNCTION h5lite_ispresent_attribute(set_id,attr_name) RESULT(ispresent)

!<description>
    ! This function tests if a given attribute is present at the data set
!</description>


!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id
    
    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name
!</input>

!<result>

    ! presence of the attribute
    LOGICAL :: ispresent
!</result>
!</function>

    ! local variables
    INTEGER(HID_T) :: attr_id
    INTEGER :: iattr,nattrs,hdferr
    CHARACTER(LEN=SYS_STRLEN) :: buf

    ! Initialization
    ispresent=.FALSE.
    
    ! Get number of attributes
    CALL h5aget_num_attrs_f(set_id,nattrs,hdferr)
    IF (hdferr /= 0) THEN
      PRINT *, "h5lite_ispresent_attribute: Unable to determine number of arguments!"
      STOP
    END IF

    ! Iterate over attributes
    DO iattr=0,nattrs-1
      CALL h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      CALL h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)

      ! Check if attribute is the desired one
      IF (TRIM(ADJUSTL(buf)) == TRIM(ADJUSTL(attr_name))) THEN
        ispresent=.TRUE.
        CALL h5aclose_f(attr_id,hdferr)
        RETURN
      END IF

      ! Close attribute
      CALL h5aclose_f(attr_id,hdferr)
    END DO
  END FUNCTION h5lite_ispresent_attribute

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_attribute_int(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets an integer attribute to a given data set
!</description>

!<input>
    
    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name

    ! value of attribut
    INTEGER, INTENT(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id,space_id
    INTEGER :: iattr,nattrs,hdferr
    CHARACTER(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    CALL h5aget_num_attrs_f(set_id,nattrs,hdferr)
    
    ! Iterate over attributes
    DO iattr=0,nattrs-1
      CALL h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      CALL h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)
      
      ! Check if attribute is the desired one
      IF (TRIM(ADJUSTL(buf))==TRIM(ADJUSTL(attr_name))) THEN
        CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_val,(/1_HSIZE_T/),hdferr)
        CALL h5aclose_f(attr_id,hdferr)
        RETURN
      END IF
      
      ! Close attribute at index
      CALL h5aclose_f(attr_id,hdferr)
    END DO
    
    ! Create new attribute
    CALL h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    CALL h5acreate_f(set_id,attr_name,H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
    CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr_val,(/1_HSIZE_T/),hdferr)
    CALL h5aclose_f(attr_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)
  END SUBROUTINE h5lite_set_attribute_int

  !*****************************************************************************

!<subroutine>
  
  SUBROUTINE h5lite_set_attribute_single(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a single attribute to a given data set
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name

    ! value of attribute
    REAL(SP), INTENT(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id,space_id
    INTEGER :: iattr,nattrs,hdferr
    CHARACTER(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    CALL h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    DO iattr=0,nattrs-1
      CALL h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      CALL h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)
      
      ! Check if attribute is the desired one
      IF (TRIM(ADJUSTL(buf))==TRIM(ADJUSTL(attr_name))) THEN
        CALL h5awrite_f(attr_id,H5T_NATIVE_REAL,attr_val,(/1_HSIZE_T/),hdferr)
        CALL h5aclose_f(attr_id,hdferr)
        RETURN
      END IF

      ! Close attribute at index
      CALL h5aclose_f(attr_id,hdferr)
    END DO

    ! Create new attribute
    CALL h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    CALL h5acreate_f(set_id,attr_name,H5T_NATIVE_REAL,space_id,attr_id,hdferr)
    CALL h5awrite_f(attr_id,H5T_NATIVE_REAL,attr_val,(/1_HSIZE_T/),hdferr)
    CALL h5aclose_f(attr_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)
  END SUBROUTINE h5lite_set_attribute_single
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_attribute_double(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a double attribute to a given data set
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id
    
    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name

    ! value of attribute
    REAL(DP), INTENT(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id,space_id
    INTEGER :: iattr,nattrs,hdferr
    CHARACTER(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    CALL h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    DO iattr=0,nattrs-1
      CALL h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      CALL h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)

      ! Check if attribute is the desired one
      IF (TRIM(ADJUSTL(buf))==TRIM(ADJUSTL(attr_name))) THEN
        CALL h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,attr_val,(/1_HSIZE_T/),hdferr)
        CALL h5aclose_f(attr_id,hdferr)
        RETURN
      END IF

      ! Close attribute at index
      CALL h5aclose_f(attr_id,hdferr)
    END DO
    
    ! Create new attribute
    CALL h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    CALL h5acreate_f(set_id,attr_name,H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
    CALL h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,attr_val,(/1_HSIZE_T/),hdferr)
    CALL h5aclose_f(attr_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)
  END SUBROUTINE h5lite_set_attribute_double

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_attribute_char(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a character attribute to a given data set
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name

    ! value of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id,space_id,atype_id
    INTEGER :: iattr,nattrs,hdferr
    CHARACTER(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    CALL h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    DO iattr=0,nattrs-1
      CALL h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      CALL h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)
      
      ! Check if attribute is the desired one
      IF (TRIM(ADJUSTL(buf))==TRIM(ADJUSTL(attr_name))) THEN
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
        CALL h5tset_size_f(atype_id,LEN(attr_val),hdferr)
        CALL h5awrite_f(attr_id,atype_id,attr_val,(/1_HSIZE_T/),hdferr)
        CALL h5aclose_f(attr_id,hdferr)
        RETURN
      END IF

      ! Close attribute at index
      CALL h5aclose_f(attr_id,hdferr)
    END DO
    
    ! Create new attribute
    CALL h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    CALL h5tset_size_f(atype_id,LEN(attr_val),hdferr)
    CALL h5acreate_f(set_id,attr_name,atype_id,space_id,attr_id,hdferr)
    CALL h5awrite_f(attr_id,atype_id,attr_val,(/1_HSIZE_T/),hdferr)
    CALL h5aclose_f(attr_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)
  END SUBROUTINE h5lite_set_attribute_char

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_attribute_logical(set_id,attr_name,attr_val)

!<description>
    ! This subroutine sets a logical attribute to a given data set.
    ! Note that HDF5 does not support "logicals" by default. As a workaround,
    ! logicals are stored as a single character with values "T" or "F".
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name

    ! value of attribute
    LOGICAL, INTENT(IN) :: attr_val
!</input>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id,space_id,atype_id
    INTEGER :: iattr,nattrs,hdferr
    CHARACTER(LEN=SYS_STRLEN) :: buf
    
    ! Check if attribute already exists
    CALL h5aget_num_attrs_f(set_id,nattrs,hdferr)

    ! Iterate over attributes
    DO iattr=0,nattrs-1
      CALL h5aopen_idx_f(set_id,iattr,attr_id,hdferr)
      CALL h5aget_name_f(attr_id,SYS_STRLEN,buf,hdferr)

      ! Check if attribute is the desired one
      IF (TRIM(ADJUSTL(buf))==TRIM(ADJUSTL(attr_name))) THEN
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
        CALL h5tset_size_f(atype_id,1,hdferr)
        IF (attr_val) THEN
          CALL h5awrite_f(attr_id,atype_id,"T",(/1_HSIZE_T/),hdferr)
        ELSE
          CALL h5awrite_f(attr_id,atype_id,"F",(/1_HSIZE_T/),hdferr)
        END IF
        CALL h5aclose_f(attr_id,hdferr)
        RETURN
      END IF

      ! Close attribute at index
      CALL h5aclose_f(attr_id,hdferr)
    END DO
    
    ! Create new attribute
    CALL h5screate_simple_f(1,(/1_HSIZE_T/),space_id,hdferr)
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    CALL h5tset_size_f(atype_id,1,hdferr)
    CALL h5acreate_f(set_id,attr_name,atype_id,space_id,attr_id,hdferr)
    IF (attr_val) THEN
      CALL h5awrite_f(attr_id,atype_id,"T",(/1_HSIZE_T/),hdferr)
    ELSE
      CALL h5awrite_f(attr_id,atype_id,"F",(/1_HSIZE_T/),hdferr)
    END IF
    CALL h5aclose_f(attr_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)
  END SUBROUTINE h5lite_set_attribute_logical

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_attribute_int(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets an integer attribute from a given data set
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id
    
    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    INTEGER, INTENT(OUT) :: attr_val
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id
    INTEGER :: hdferr
    
    CALL h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    CALL h5aread_f(attr_id,H5T_NATIVE_INTEGER,attr_val,(/1_HSIZE_T/),hdferr) 
    CALL h5aclose_f(attr_id,hdferr)
  END SUBROUTINE h5lite_get_attribute_int

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_attribute_single(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a single attribute from a given data set
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    REAL(SP), INTENT(OUT) :: attr_val
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id
    INTEGER :: hdferr
    
    CALL h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    CALL h5aread_f(attr_id,H5T_NATIVE_REAL,attr_val,(/1_HSIZE_T/),hdferr) 
    CALL h5aclose_f(attr_id,hdferr)
  END SUBROUTINE h5lite_get_attribute_single

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_attribute_double(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a double attribute from a given data set
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    REAL(DP), INTENT(OUT) :: attr_val
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: attr_id
    INTEGER :: hdferr
    
    CALL h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    CALL h5aread_f(attr_id,H5T_NATIVE_DOUBLE,attr_val,(/1_HSIZE_T/),hdferr) 
    CALL h5aclose_f(attr_id,hdferr)
  END SUBROUTINE h5lite_get_attribute_double

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_attribute_char(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a character attribute from a given data set
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    CHARACTER(LEN=*), INTENT(OUT) :: attr_val
!</output>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: attr_id,atype_id
    INTEGER :: hdferr
    
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    CALL h5tset_size_f(atype_id,LEN(attr_val),hdferr)
    CALL h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    CALL h5aread_f(attr_id,atype_id,attr_val,(/1_HSIZE_T/),hdferr) 
    CALL h5aclose_f(attr_id,hdferr)
  END SUBROUTINE h5lite_get_attribute_char

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_attribute_logical(set_id,attr_name,attr_val)

!<description>
    ! This subroutine gets a logical attribute from a given data set.
    ! Note that HDF5 does not support "logicals" by default. As a workaround,
    ! logicals are stored as a single character with values "T" or "F".
!</description>

!<input>

    ! data set
    INTEGER(HID_T), INTENT(IN) :: set_id

    ! name of attribute
    CHARACTER(LEN=*), INTENT(IN) :: attr_name
!</input>

!<output>

    ! value of attribute
    LOGICAL, INTENT(OUT) :: attr_val
!</output>
!</subroutine>
    
    ! local variables
    CHARACTER(LEN=1) :: tmp_val
    INTEGER(HID_T) :: attr_id,atype_id
    INTEGER :: hdferr
    
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdferr)
    CALL h5tset_size_f(atype_id,1,hdferr)
    CALL h5aopen_name_f(set_id,attr_name,attr_id,hdferr)
    CALL h5aread_f(attr_id,atype_id,tmp_val,(/1_HSIZE_T/),hdferr) 
    CALL h5aclose_f(attr_id,hdferr)
    
    SELECT CASE(tmp_val)
    CASE("T","t")
      attr_val=.TRUE.
    CASE("F","f")
      attr_val=.FALSE.
    CASE DEFAULT
      PRINT *, "h5lite_get_attribute_logical: Invalid attribute value!"
      STOP
    END SELECT
  END SUBROUTINE h5lite_get_attribute_logical

  !*****************************************************************************

!<function>

  FUNCTION h5lite_ispresent_data(loc_id,name) RESULT(ispresent)

!<description>
    ! This function tests if a given data set is present at location
!</description>


!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id
    
    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<result>

    ! presence of data
    LOGICAL :: ispresent
!</result>
!</function>

    ! local variables
    CHARACTER(SYS_STRLEN) :: root_name,obj_name
    INTEGER(HID_T) :: dataset_id,name_size
    INTEGER :: obj_type,imember,nmembers,hdferr

    ! Initialization
    ispresent=.FALSE.

    ! Determine name of location and number of members
    CALL h5iget_name_f(loc_id,root_name,INT(SYS_STRLEN,SIZE_T),name_size,hdferr)
    CALL h5gn_members_f(loc_id,TRIM(ADJUSTL(root_name)),nmembers,hdferr)

    ! Loop over members
    DO imember=0,nmembers-1
      CALL h5gget_obj_info_idx_f(loc_id,TRIM(ADJUSTL(root_name)),imember,&
          obj_name,obj_type,hdferr)
      
      IF ((obj_type == H5G_DATASET_F .OR. obj_type == H5G_LINK_F) .AND.&
          TRIM(ADJUSTL(obj_name)) == TRIM(ADJUSTL(name))) THEN
        ispresent=.TRUE.
        RETURN
      END IF
    END DO
  END FUNCTION h5lite_ispresent_data

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_data_handle(loc_id,name,ihandle)

!<description>
    ! This subroutine sets the data associated to a handle to a given group
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! handle to data
    INTEGER, INTENT(IN) :: ihandle

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>
!</subroutine>

    ! local variables
    INTEGER :: idatatype
    INTEGER :: idimension
    INTEGER, DIMENSION(:,:), POINTER  :: data_int2D
    INTEGER, DIMENSION(:), POINTER    :: data_int
    REAL(DP), DIMENSION(:,:), POINTER :: data_double2D
    REAL(DP), DIMENSION(:), POINTER   :: data_double
    REAL(SP), DIMENSION(:,:), POINTER :: data_single2D
    REAL(SP), DIMENSION(:), POINTER   :: data_single
    
    ! Get dimension and type of data
    CALL storage_getdatatype(ihandle,idatatype)
    CALL storage_getdimension(ihandle,idimension)

    ! Which dimension are we?
    SELECT CASE(idimension)
    CASE (1)
       
       ! Which datatype are we?
       SELECT CASE (idatatype)
       CASE (ST_INT)
          CALL storage_getbase_int(ihandle,data_int)
          CALL h5lite_set_data_int(loc_id,name,data_int)

       CASE (ST_SINGLE)
          CALL storage_getbase_single(ihandle,data_single)
          CALL h5lite_set_data_single(loc_id,name,data_single)

       CASE (ST_DOUBLE)
          CALL storage_getbase_double(ihandle,data_double)
          CALL h5lite_set_data_double(loc_id,name,data_double)

       CASE DEFAULT
          PRINT *, "(EE) h5lite_set_data_handle: invalid data type"
          STOP
       END SELECT

    CASE (2)

       ! Which datatype are we?
       SELECT CASE (idatatype)
       CASE (ST_INT)
          CALL storage_getbase_int2D(ihandle,data_int2D)
          CALL h5lite_set_data_int2D(loc_id,name,data_int2D)

       CASE (ST_SINGLE)
          CALL storage_getbase_single2D(ihandle,data_single2D)
          CALL h5lite_set_data_single2D(loc_id,name,data_single2D)

       CASE (ST_DOUBLE)
          CALL storage_getbase_double2D(ihandle,data_double2D)
          CALL h5lite_set_data_double2D(loc_id,name,data_double2D)

       CASE DEFAULT
          PRINT *, "(EE) h5lite_set_data_handle: invalid data type"
          STOP
       END SELECT

    CASE DEFAULT
       PRINT *, "(EE) h5lite_set_data_handle: dimension exceeds 2"
       STOP
    END SELECT
  END SUBROUTINE h5lite_set_data_handle

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_data_int(loc_id,name,data)

!<description>
    ! This subroutine sets a 1D integer array to a given group
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id
    
    ! 1D integer array
    INTEGER, DIMENSION(:), INTENT(IN) :: data

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: set_id,space_id,cpl
    INTEGER(HSIZE_T) :: dim(1)
    INTEGER :: hdferr
    
    dim=SHAPE(data); IF (SUM(dim) == 0) RETURN
    
    CALL h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    CALL h5sset_extent_simple_f(space_id,1,dim,dim,hdferr)
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    CALL h5pset_chunk_f(cpl,1,dim,hdferr)
    CALL h5pset_deflate_f(cpl,idcl,hdferr)
    CALL h5dcreate_f(loc_id,name,H5T_NATIVE_INTEGER,space_id,set_id,hdferr,cpl)
    CALL h5dwrite_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
    CALL h5pclose_f(cpl,hdferr)
    CALL h5dclose_f(set_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)   
  END SUBROUTINE h5lite_set_data_int

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_data_int2D(loc_id,name,data)

!<description>
    ! This subroutine sets a 2D integer array to a given group
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! 2D integer array
    INTEGER, DIMENSION(:,:), INTENT(IN) :: data

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: set_id,space_id,cpl
    INTEGER(HSIZE_T) :: dim(2)
    INTEGER :: hdferr

    dim=SHAPE(data); IF (SUM(dim) == 0) RETURN
    
    CALL h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    CALL h5sset_extent_simple_f(space_id,2,dim,dim,hdferr)
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    CALL h5pset_chunk_f(cpl,2,dim,hdferr)
    CALL h5pset_deflate_f(cpl,idcl,hdferr)
    CALL h5dcreate_f(loc_id,name,H5T_NATIVE_INTEGER,space_id,set_id,hdferr,cpl)
    CALL h5dwrite_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
    CALL h5pclose_f(cpl,hdferr)
    CALL h5dclose_f(set_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)   
  END SUBROUTINE h5lite_set_data_int2D

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_data_single(loc_id,name,data)

!<description>
    ! This subroutine sets a 1D single array to a given group
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! 1D single array
    REAL(SP), DIMENSION(:), INTENT(IN) :: data

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: set_id,space_id,cpl
    INTEGER(HSIZE_T) :: dim(1)
    INTEGER :: hdferr
    
    dim=SHAPE(data); IF (SUM(dim) == 0) RETURN
    
    CALL h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    CALL h5sset_extent_simple_f(space_id,1,dim,dim,hdferr)
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    CALL h5pset_chunk_f(cpl,1,dim,hdferr)
    CALL h5pset_deflate_f(cpl,idcl,hdferr)
    CALL h5dcreate_f(loc_id,name,H5T_NATIVE_REAL,space_id,set_id,hdferr,cpl)
    CALL h5dwrite_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
    CALL h5pclose_f(cpl,hdferr)
    CALL h5dclose_f(set_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)   
  END SUBROUTINE h5lite_set_data_single

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_data_single2D(loc_id,name,data)

!<description>
    ! This subroutine sets a 2D single array to a given group
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! 2D single array
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: data

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: set_id,space_id,cpl
    INTEGER(HSIZE_T) :: dim(2)
    INTEGER :: hdferr

    dim=SHAPE(data); IF (SUM(dim) == 0) RETURN

    CALL h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    CALL h5sset_extent_simple_f(space_id,2,dim,dim,hdferr)
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    CALL h5pset_chunk_f(cpl,2,dim,hdferr)
    CALL h5pset_deflate_f(cpl,idcl,hdferr)
    CALL h5dcreate_f(loc_id,name,H5T_NATIVE_REAL,space_id,set_id,hdferr,cpl)
    CALL h5dwrite_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
    CALL h5pclose_f(cpl,hdferr)
    CALL h5dclose_f(set_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)   
  END SUBROUTINE h5lite_set_data_single2D

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_data_double(loc_id,name,data)

!<description>
    ! This subroutine sets a 1D double array to a given group
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! 1D double array
    REAL(DP), DIMENSION(:), INTENT(IN) :: data

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: set_id,space_id,cpl
    INTEGER(HSIZE_T) :: dim(1)
    INTEGER :: hdferr

    dim=SHAPE(data); IF (SUM(dim) == 0) RETURN

    CALL h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    CALL h5sset_extent_simple_f(space_id,1,dim,dim,hdferr)
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    CALL h5pset_chunk_f(cpl,1,dim,hdferr)
    CALL h5pset_deflate_f(cpl,idcl,hdferr)
    CALL h5dcreate_f(loc_id,name,H5T_NATIVE_DOUBLE,space_id,set_id,hdferr,cpl)
    CALL h5dwrite_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
    CALL h5pclose_f(cpl,hdferr)
    CALL h5dclose_f(set_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)   
  END SUBROUTINE h5lite_set_data_double

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_set_data_double2D(loc_id,name,data)

!<description>
    ! This subroutine sets a 2D double array to a given group
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! 2D double array
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: data

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: set_id,space_id,cpl
    INTEGER(HSIZE_T) :: dim(2)
    INTEGER :: hdferr
    
    dim=SHAPE(data); IF (SUM(dim) == 0) RETURN
    
    CALL h5screate_f(H5S_SIMPLE_F,space_id,hdferr)
    CALL h5sset_extent_simple_f(space_id,2,dim,dim,hdferr)
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,cpl,hdferr)
    CALL h5pset_chunk_f(cpl,2,dim,hdferr)
    CALL h5pset_deflate_f(cpl,idcl,hdferr)
    CALL h5dcreate_f(loc_id,name,H5T_NATIVE_DOUBLE,space_id,set_id,hdferr,cpl)
    CALL h5dwrite_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
    CALL h5pclose_f(cpl,hdferr)
    CALL h5dclose_f(set_id,hdferr)
    CALL h5sclose_f(space_id,hdferr)   
  END SUBROUTINE h5lite_set_data_double2D

  !*****************************************************************************

!<subroutine>
  
  SUBROUTINE h5lite_get_data_handle(loc_id,name,ihandle)

!<description>
    ! This subroutine gets a stored data and associates it to a handle
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id
    
    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<output>
    
    ! handle to data
    INTEGER, INTENT(OUT) :: ihandle
!</output>
!</subroutine>

    ! local variables
    INTEGER, DIMENSION(:), POINTER :: idata
    REAL(SP), DIMENSION(:), POINTER :: sdata
    REAL(DP), DIMENSION(:), POINTER :: ddata
    INTEGER, DIMENSION(:,:), POINTER :: idata2D
    REAL(SP), DIMENSION(:,:), POINTER :: sdata2D
    REAL(DP), DIMENSION(:,:), POINTER :: ddata2D
    
    INTEGER(HID_T) :: dataset_id,datatype_id,dataspace_id
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1D,maxdim1D
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2D,maxdim2D
    INTEGER :: ndims,hdferr
    LOGICAL :: istype
    
    ! Get the data set
    CALL h5dopen_f(loc_id,name,dataset_id,hdferr)
    
    ! Get the data type
    CALL h5dget_type_f(dataset_id,datatype_id,hdferr)
    
    ! Get the dataspace
    CALL h5dget_space_f(dataset_id,dataspace_id,hdferr)
    
    ! Get number of dimensions
    CALL h5sget_simple_extent_ndims_f(dataspace_id,ndims,hdferr)
    
    SELECT CASE(ndims)
    CASE (1)
      ! Get dimensions
      CALL h5sget_simple_extent_dims_f(dataspace_id,dim1D,maxdim1D,hdferr)

      ! Are we integer?
      CALL h5tequal_f(datatype_id,H5T_NATIVE_INTEGER,istype,hdferr)
      IF (istype) THEN
        CALL h5dclose_f(dataset_id,hdferr)
        CALL storage_new('h5lite_get_data_handle',TRIM(ADJUSTL(name)), &
             INT(dim1D(1)),ST_INT,ihandle,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int(ihandle,idata)
        CALL h5lite_get_data_int(loc_id,name,idata)
        RETURN
      END IF
      
      ! Are we single?
      CALL h5tequal_f(datatype_id,H5T_NATIVE_REAL,istype,hdferr)
      IF (istype) THEN
        CALL h5dclose_f(dataset_id,hdferr)
        CALL storage_new('h5lite_get_data_handle',TRIM(ADJUSTL(name)), &
             INT(dim1D(1)),ST_SINGLE,ihandle,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single(ihandle,sdata)
        CALL h5lite_get_data_single(loc_id,name,sdata)
        RETURN
      END IF
      
      ! Are we double?
      CALL h5tequal_f(datatype_id,H5T_NATIVE_DOUBLE,istype,hdferr)
      IF (istype) THEN
        CALL h5dclose_f(dataset_id,hdferr)
        CALL storage_new('h5lite_get_data_handle',TRIM(ADJUSTL(name)), &
             INT(dim1D(1)),ST_DOUBLE,ihandle,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double(ihandle,ddata)
        CALL h5lite_get_data_double(loc_id,name,ddata)
        RETURN
      END IF
      
      ! We are the wrong datatype
      PRINT *, "(EE) h5lite_get_data_handle: unknown datatype"
      STOP
      
    CASE (2)
      ! Get dimensions
      CALL h5sget_simple_extent_dims_f(dataspace_id,dim2D,maxdim2D,hdferr)

      ! Are we integer?
      CALL h5tequal_f(datatype_id,H5T_NATIVE_INTEGER,istype,hdferr)
      IF (istype) THEN
        CALL h5dclose_f(dataset_id,hdferr)
        CALL storage_new('h5lite_get_data_handle',TRIM(ADJUSTL(name)), &
             INT(dim2D),ST_INT,ihandle,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int2D(ihandle,idata2D)
        CALL h5lite_get_data_int2D(loc_id,name,idata2D)
        RETURN
      END IF
      
      ! Are we single?
      CALL h5tequal_f(datatype_id,H5T_NATIVE_REAL,istype,hdferr)
      IF (istype) THEN
        CALL h5dclose_f(dataset_id,hdferr)
        CALL storage_new('h5lite_get_data_handle',TRIM(ADJUSTL(name)), &
             INT(dim2D),ST_SINGLE,ihandle,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single2D(ihandle,sdata2D)
        CALL h5lite_get_data_single2D(loc_id,name,sdata2D)
        RETURN
      END IF
      
      ! Are we double?
      CALL h5tequal_f(datatype_id,H5T_NATIVE_DOUBLE,istype,hdferr)
      IF (istype) THEN
        CALL h5dclose_f(dataset_id,hdferr)
        CALL storage_new('h5lite_get_data_handle',TRIM(ADJUSTL(name)), &
             INT(dim2D),ST_DOUBLE,ihandle,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2D(ihandle,ddata2D)
        CALL h5lite_get_data_double2D(loc_id,name,ddata2D)
        RETURN
      END IF
      
      ! We are the wrong datatype
      PRINT *, "(EE) h5lite_get_data_handle: unknown datatype"
      STOP
      
    CASE DEFAULT
      PRINT *, "(EE) h5lite_get_data_handle: number of dimensions exceeds 2"
      STOP
    END SELECT
  END SUBROUTINE h5lite_get_data_handle

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_data_int(loc_id,name,data)

!<description>
    ! This subroutine gets a 1D integer array
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<output>
    ! 1D integer array
    INTEGER, DIMENSION(:), INTENT(OUT) :: data
!</output>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: set_id
    INTEGER(HSIZE_T) :: dim(1)
    INTEGER :: hdferr
    
    dim=SHAPE(data)
    CALL h5dopen_f(loc_id,name,set_id,hdferr)
    IF (hdferr == 0) THEN
      CALL h5dread_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
      CALL h5dclose_f(set_id,hdferr)
    END IF
  END SUBROUTINE h5lite_get_data_int
  
  !*****************************************************************************

!<subroutine>

   SUBROUTINE h5lite_get_data_int2D(loc_id,name,data)

!<description>
    ! This subroutine gets a 2D integer array
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id
    
    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<output>
    ! 2D integer array
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: data
!</output>
!</subroutine>
    
    ! local variables
    INTEGER(HID_T) :: set_id
    INTEGER(HSIZE_T) :: dim(2)
    INTEGER :: hdferr
    
    dim=SHAPE(data)
    CALL h5dopen_f(loc_id,name,set_id,hdferr)
    IF (hdferr == 0) THEN
      CALL h5dread_f(set_id,H5T_NATIVE_INTEGER,data,dim,hdferr)
      CALL h5dclose_f(set_id,hdferr)
    END IF
  END SUBROUTINE h5lite_get_data_int2D
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_data_single(loc_id,name,data)

!<description>
    ! This subroutine gets a 1D single array
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<output>
    ! 1D single array
    REAL(SP), DIMENSION(:), INTENT(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: set_id
    INTEGER(HSIZE_T) :: dim(1)
    INTEGER :: hdferr
    
    dim=SHAPE(data)
    CALL h5dopen_f(loc_id,name,set_id,hdferr)
    IF (hdferr == 0) THEN
      CALL h5dread_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
      CALL h5dclose_f(set_id,hdferr)
    END IF
  END SUBROUTINE h5lite_get_data_single

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_data_single2D(loc_id,name,data)

!<description>
    ! This subroutine gets a 2D single array
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<output>
    ! 2D single array
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: set_id
    INTEGER(HSIZE_T) :: dim(2)
    INTEGER :: hdferr
    
    dim=SHAPE(data)
    CALL h5dopen_f(loc_id,name,set_id,hdferr)
    IF (hdferr == 0) THEN
      CALL h5dread_f(set_id,H5T_NATIVE_REAL,data,dim,hdferr)
      CALL h5dclose_f(set_id,hdferr)
    END IF
  END SUBROUTINE h5lite_get_data_single2D
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_data_double(loc_id,name,data)

!<description>
    ! This subroutine gets a 1D double array
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<output>
    REAL(DP), DIMENSION(:), INTENT(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: set_id
    INTEGER(HSIZE_T) :: dim(1)
    INTEGER :: hdferr
    
    dim=SHAPE(data)
    CALL h5dopen_f(loc_id,name,set_id,hdferr)
    IF (hdferr == 0) THEN
      CALL h5dread_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
      CALL h5dclose_f(set_id,hdferr)
    END IF
  END SUBROUTINE h5lite_get_data_double

  !*****************************************************************************

!<subroutine>

  SUBROUTINE h5lite_get_data_double2D(loc_id,name,data)

!<description>
    ! This subroutine gets a 2D double array
!</description>

!<input>

    ! location identifier
    INTEGER(HID_T), INTENT(IN) :: loc_id

    ! name of data
    CHARACTER(LEN=*), INTENT(IN) :: name
!</input>

!<output>
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: data
!</output>
!</subroutine>

    ! local variables
    INTEGER(HID_T) :: set_id
    INTEGER(HSIZE_T) :: dim(2)
    INTEGER :: hdferr
    
    dim=SHAPE(data)
    CALL h5dopen_f(loc_id,name,set_id,hdferr)
    IF (hdferr == 0) THEN
      CALL h5dread_f(set_id,H5T_NATIVE_DOUBLE,data,dim,hdferr)
      CALL h5dclose_f(set_id,hdferr)
    END IF
  END SUBROUTINE h5lite_get_data_double2D
END MODULE h5lite
