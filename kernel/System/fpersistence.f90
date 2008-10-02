!##############################################################################
!# ****************************************************************************
!# <name> fpersistence </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic data structures and routines to make all
!# data persistent, that is, store it to disk and retrieve it later.
!#
!# The basic concept of the persistence database is as follows:
!#
!#    t_fpdbObjectItem: Holds information about the abstract
!#                      object, e.g. a scalar matrix structure.
!#
!#    t_fpdbDataItem: Holds the data and its information, e.g.
!#                    the matrix entries and the sparsity pattern.
!#
!#
!# An abstract ObjectItem corresponds to a derived type that is
!# constructured in one of the kernel modules and/or the user's
!# application. Therefore, an ObjectItem can only provide very
!# rudimentary information about the derived type itself.
!# It is identified by the universally unique identifier (UUID)
!# which is unique for each ObjectItem. A derived type that should
!# be made persist has to implement the two following routines:
!#
!#    xxxx_createFpdbObject: Creates an ObjectItem from the given
!#                           derived type including its data.
!#
!#    xxx_restoreFpdbObject: Restores the derived type from the
!#                           given ObjectItem including its data.
!#
!#
!# In order to store the actual data, each ObjectItem provides an array
!# of DataItems which needs to be filled by the xxxx_createFpdbObject
!# object routine. In addition to the pure data, each DataItem provides
!# meta-information about the data such as size, type, dimension, etc.
!#
!# The derived type is restored from the ObjectItem and its associated
!# DataItems by the xxx_restoreFpdbObject routine which has to be
!# implemented in the particular kernel/application modules.
!#
!#
!# The ObjectItems are stored in the object table which is organised as
!# a binary tree using the UUIDs as unique keys. Typically, objects are
!# retrieved by their UUID so that an efficient search algorithm is 
!# essential. The complete object table is exported to disk by using
!# Fortran's direct-access facility to read and write files record-based.
!# A second file is created - the so-called data table. It stores all 
!# meta-information about the data such as data type, dimension and size.
!# In addition, it provides the UUID which is associated with the actual
!# data file (a standard binary file). For atomic data such as a single
!# integer value, the data is exported to the data table file directly.
!#
!# The following routines can be found here:
!#
!#  1.) fpdb_init
!#      -> Initialises the persistence database and imports
!#         all tables and data structures if required
!#
!#  2.) fpdb_done
!#      -> Finalises the persistence database
!#
!#  3.) fpdb_import
!#      -> Import the persistence database
!#
!#  4.) fpdb_export
!#      -> Export the persistence database
!#
!#  5.) fpdb_info
!#      -> Output information about persistence database
!#
!# The following auxiliary routines can be found here:
!#
!#  1.) fpdb_insertObject
!#      -> Inserts an item into the object table
!#
!#  2.) fpdb_removeObject
!#      -> Removes an item from the object table
!#
!#  3.) fpdb_retrieveObject = fpdb_retrieveObjectByUUID /
!#                            fpdb_retrieveObjectByName
!#      -> Retrieves an item from the object table
!#
!#  4.) fpdb_getObjectLength
!#      -> Returns the record length of the ObjectItem
!#
!#  5.) fpdb_getDataLength
!#      -> Returns the record length of the DataItem
!#
!#  6.) fpdb_getdata_single1d,
!#      fpdb_getdata_double1d,
!#      fpdb_getdata_int1d,
!#      fpdb_getdata_logical1d,
!#      fpdb_getdata_char1d
!#      -> Import the data item associated to a DataItem
!# </purpose>
!##############################################################################
module fpersistence

  use fsystem
  use io
  use uuid
  
  implicit none

!<constants>

!<constantblock description="Persistence database constants">

  ! Maximum number of dimensions supported by the persistence database
  integer, parameter :: FPDB_MAXDIM = 2
!</constantblock>


!<constantblock description="Flags for the persistence database specification bitfield">

  ! Database is readable
  integer, parameter :: FPDB_MSPEC_ISREADABLE  = 0

  ! Database is writeable
  integer, parameter :: FPDB_MSPEC_ISWRITABLE  = 2**1

  ! Database is imported from disk
  integer, parameter :: FPDB_MSPEC_IMPORTFIRST = 2**2

  ! Standard specifier for persistence database
  integer, parameter :: FPDB_MSPEC_STANDARD    = FPDB_MSPEC_ISREADABLE+&
                                                 FPDB_MSPEC_ISWRITABLE
!</constantblock>


!<constantblock description="Data item identifiers">

  ! defines an undefined data item
  ! This allows to create a nullified data item,
  ! e.g., to indicate that a pointer is nullified
  integer, parameter :: FPDB_NULL = 0
  
  ! defines an atomic ObjectItem as data item:
  ! This allows to create an ObjectItem for, say, a block matrix
  ! which has an array of scalar matrices as data. Note that the
  ! ObjectItem means that the corresponding item is physically
  ! created when it is imported from the persistence database.
  integer, parameter :: FPDB_OBJECT = 1
  
  ! defines an atomic Link as data item:
  ! This allows to create a link to an ObjectItem. As an example
  ! consider a scalar vectors which may feature a pointer to the
  ! spatial discretisation. Note that the spatial discretisation
  ! is not created when the vector is created but only a pointer
  ! to the spatial discretisation is associated.
  integer, parameter :: FPDB_LINK = 2

  ! defines an atomic single data item
  integer, parameter :: FPDB_SINGLE = 3

  ! defines an atomic double data item
  integer, parameter :: FPDB_DOUBLE = 4

  ! defines an atomic integer data item
  integer, parameter :: FPDB_INT = 5

  ! defines an atomic logical data item
  integer, parameter :: FPDB_LOGICAL = 6

  ! defines an atomic character data item
  integer, parameter :: FPDB_CHAR = 7

  ! defines a 1D single data item
  integer, parameter :: FPDB_SINGLE1D = 8

  ! defines a 1D double data item
  integer, parameter :: FPDB_DOUBLE1D = 9

  ! defines a 1D integer data item
  integer, parameter :: FPDB_INT1D = 10

  ! defines a 1D logical data item
  integer, parameter :: FPDB_LOGICAL1D = 11

  ! defines a 1D character data item
  integer, parameter :: FPDB_CHAR1D = 12

  ! defines a 2D single data item
  integer, parameter :: FPDB_SINGLE2D = 13

  ! defines a 2D double data item
  integer, parameter :: FPDB_DOUBLE2D = 14

  ! defines a 2D integer data item
  integer, parameter :: FPDB_INT2D = 15

  ! defines a 2D logical data item
  integer, parameter :: FPDB_LOGICAL2D = 16

  ! defines a 2D character data item
  integer, parameter :: FPDB_CHAR2D = 17
  
!</constantblock>

!</constants>

!************************************************************************

!<types>

!<typeblock>

  ! This type block specifies the persistence database.
  
  type t_fpdb
    
    private

    ! Name of the directory of the persistence database
    character(LEN=SYS_STRLEN) :: spath = './'

    ! Filename of the persistence database
    character(LEN=SYS_NAMELEN) :: sfilename = 'feat2pdb'
    
    ! Database specification flag. This is a bitfield coming from an OR 
    ! combination of different FPDB_MSPEC_xxxx constants and specifies
    ! various details of the database.
    integer :: idatabaseSpec = FPDB_MSPEC_STANDARD
    
    ! Unit number of object table
    integer :: iunitObjectTable = 0

    ! Unit number of data table
    integer :: iunitDataTable = 0
    
    ! Next free record number of object table
    integer :: irecnumObjectTable = 1

    ! Next free record number of data table
    integer :: irecnumDataTable = 1

    ! Root node of the object table of the persistence database
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectTable => null()
  end type t_fpdb

!</typeblock>

!<typeblock>

  ! This type block specifies a single item of the object table
  ! which is used in the  persistence database structure.
  
  type t_fpdbObjectItem
      
    ! Universally unique identifier of the item 
    type(t_uuid) :: ruuid

    ! Type descriptor of the item
    character(LEN=SYS_NAMELEN) :: stype = ''

    ! Full qualified name of the item
    character(LEN=SYS_NAMELEN) :: sname = ''

    ! First record number in data table which is associated with
    ! DataItems belonging to this ObjectItem
    integer :: ifirstDataRecnum = 0

    ! Number if records in data table which are associated with
    ! DataItems belonging to this ObjectItem
    integer :: ilengthDataRecnum = 0

    ! Array of items from the data table
    type(t_fpdbDataItem), dimension(:), pointer :: p_RfpdbDataItem => null()

    ! Father in the persistence database
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectFather => null()

    ! Left child in the persistence database
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectLeft => null()

    ! Right child in the persistence database
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectRight => null()
    
  end type t_fpdbObjectItem

!</typeblock>

!<typeblock>

  ! This type block specifies a single item of the data table
  ! which is used in the persistence database structure.
  ! Note that the multidimensional quantities are pointers which
  ! can be used either to 'point' to an existing quantity or
  ! they can be allocated to provide memory for the stored data.

  type t_fpdbDataItem
    
    ! Type descriptor of the data item
    integer :: ctype = FPDB_NULL
    
    ! Name of the data item
    character(LEN=SYS_NAMELEN) :: sname = ''

    ! Filename of the data item (if any)
    character(LEN=SYS_STRLEN) :: sfilename = './feat2pdb'
    
    ! Array of lower bounds of the associated data item
    integer, dimension(FPDB_MAXDIM) :: Ilbounds

    ! Array of upper bounds of the associated data item
    integer, dimension(FPDB_MAXDIM) :: Iubounds

    ! Pointer to an ObjectItem
    type(t_fpdbObjectItem), pointer :: p_fpdbObjectItem => null()

    ! Below, the implemented data types follow:

    ! -~-~-~ Auxiliary quantity for ObjectItems and links -~-~-~

    type(t_uuid) :: ruuid
    
    ! -~-~-~ Atomic quantities -~-~-~

    real(DP)  :: ddouble
    real(SP)  :: fsingle
    integer   :: iinteger
    logical   :: blogical
    character(len=SYS_NAMELEN) :: schar
    
    ! -~-~-~ 1D quantities -~-~-~
    
    real(DP), dimension(:), pointer  :: p_Ddouble1D
    real(SP), dimension(:), pointer  :: p_Fsingle1D
    integer, dimension(:), pointer   :: p_Iinteger1D
    logical, dimension(:), pointer   :: p_Blogical1D
    character, dimension(:), pointer :: p_Schar1D

    ! -~-~-~ 2D quantities -~-~-~
    
    real(DP), dimension(:,:), pointer  :: p_Ddouble2D
    real(SP), dimension(:,:), pointer  :: p_Fsingle2D
    integer, dimension(:,:), pointer   :: p_Iinteger2D
    logical, dimension(:,:), pointer   :: p_Blogical2D
    character, dimension(:,:), pointer :: p_Schar2D

  end type t_fpdbDataItem
  
!</typeblock>

!</types>

  interface fpdb_retrieveObject
    module procedure fpdb_retrieveObjectByUUID
    module procedure fpdb_retrieveObjectByName
  end interface

contains

!************************************************************************

!<subroutine>

  subroutine fpdb_init(rfpdb, spath, sfilename, idatabaseSpec)

!<description>
    ! This subroutine initialises the persistence database structure.
!</description>

!<input>
    ! Path of the persistence database
    character(LEN=*), intent(in) :: spath

    ! Filename of the persistence database
    character(LEN=*), intent(in) :: sfilename

    ! OPTIONAL: database specification bitfield
    integer, intent(in), optional :: idatabaseSpec
!</input>

!<output>
    ! Database structure to initialise.
    type(t_fpdb), intent(out) :: rfpdb    
!</output>
!</subroutine>
    
    ! local variables
    character(LEN=SYS_STRLEN) :: caction, cstatus
    integer :: irecordLength

    ! Initialise database structure
    if (sfilename .eq. '') then
      rfpdb%sfilename = 'feat2pdb'
    else
      rfpdb%sfilename = trim(sfilename)
    end if

    if (spath .eq. '') then
      rfpdb%spath = './'
    else
      rfpdb%spath = trim(spath)//'/'
    end if

    ! Set database specification (if required)
    if (present(idatabaseSpec)) then
      rfpdb%idatabaseSpec = idatabaseSpec
    else
      rfpdb%idatabaseSpec = FPDB_MSPEC_STANDARD
    end if

    ! Is database writeable?
    if (iand(rfpdb%idatabaseSpec, FPDB_MSPEC_ISWRITABLE) .ne. 0) then
      caction = 'readwrite'
    else
      caction = 'read'
    end if

    ! Is database present in disk?
    if (iand(rfpdb%idatabaseSpec, FPDB_MSPEC_IMPORTFIRST) .ne. 0) then
      cstatus = 'old'
    else
      cstatus = 'replace'
    end if

    ! Open tables for persistence database
    rfpdb%iunitObjectTable = sys_getFreeUnit()
    irecordLength          = fpdb_getObjectLength()
    open(unit=rfpdb%iunitObjectTable,&
         file=trim(rfpdb%spath)//trim(rfpdb%sfilename)//'_obj.tbl',&
         status=trim(cstatus), access='direct', recl=irecordLength,&
         form='unformatted', action=trim(caction), err=1)
    
    rfpdb%iunitDataTable = sys_getFreeUnit()
    irecordLength        = fpdb_getDataLength()
    open(unit=rfpdb%iunitDataTable,&
         file=trim(rfpdb%spath)//trim(rfpdb%sfilename)//'_data.tbl',&
         status=trim(cstatus), access='direct', recl=irecordLength,&
         form='unformatted', action=trim(caction), err=2)   
        
    ! That's it
    return


1   call output_line ('Unable to open the object table!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_init')
    call sys_halt()

2   call output_line ('Unable to open the data table!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_init')
    call sys_halt()

  end subroutine fpdb_init

!************************************************************************

!<subroutine>

  subroutine fpdb_done(rfpdb)

!<description>
    ! This subroutine finalises the persistence database structure.
!</description>

!<inputoutput>
    ! Database structure to finalise.
    type(t_fpdb), intent(inout) :: rfpdb
!</inputoutput>
!</subroutine>
    
    ! local variables
    logical :: bisOpened

    ! Close files for object, relation and storage tables
    inquire(rfpdb%iunitObjectTable, opened=bisOpened)
    if (bisOpened) close(rfpdb%iunitObjectTable)
    
    inquire(rfpdb%iunitDataTable, opened=bisOpened)
    if (bisOpened) close(rfpdb%iunitDataTable)

    ! Release the structure of ObjectItems
    call releaseObjectItem(rfpdb%p_rfpdbObjectTable)

  contains

    ! Here, the real cleaning routines follow.

    !**************************************************************
    ! Postorder traversal of the object table.

    recursive subroutine releaseObjectItem(p_rfpdbObjectItem)

      type(t_fpdbObjectItem), pointer :: p_rfpdbObjectItem

      integer :: i

      ! Proceed to left child
      if (associated(p_rfpdbObjectItem%p_rfpdbObjectLeft)) then
        call releaseObjectItem(p_rfpdbObjectItem%p_rfpdbObjectLeft)
      end if

      ! Proceed to right child
      if (associated(p_rfpdbObjectItem%p_rfpdbObjectRight)) then
        call releaseObjectItem(p_rfpdbObjectItem%p_rfpdbObjectRight)
      end if

      ! Remove associated data items
      if (associated(p_rfpdbObjectItem%p_RfpdbDataItem)) then
        do i = 1, size(p_rfpdbObjectItem%p_RfpdbDataItem)
          call releaseDataItem(p_rfpdbObjectItem%p_RfpdbDataItem(i))
        end do
        deallocate(p_rfpdbObjectItem%p_RfpdbDataItem)
      end if
      
      ! Release current ObjectItem
      deallocate(p_rfpdbObjectItem)

    end subroutine releaseObjectItem

    !**************************************************************
    ! Release the DataItem

    subroutine releaseDataItem(rfpdbDataItem)

      type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem

      rfpdbDataItem%ctype = FPDB_NULL
      rfpdbDataItem%sname = ''

      call uuid_createUUID(0, rfpdbDataItem%ruuid)

      nullify(rfpdbDataItem%p_fpdbObjectItem)

      nullify(rfpdbDataItem%p_Ddouble1D)
      nullify(rfpdbDataItem%p_Fsingle1D)
      nullify(rfpdbDataItem%p_Iinteger1D)
      nullify(rfpdbDataItem%p_Blogical1D)
      nullify(rfpdbDataItem%p_Schar1D)

      nullify(rfpdbDataItem%p_Ddouble2D)
      nullify(rfpdbDataItem%p_Fsingle2D)
      nullify(rfpdbDataItem%p_Iinteger2D)
      nullify(rfpdbDataItem%p_Blogical2D)
      nullify(rfpdbDataItem%p_Schar2D)

    end subroutine releaseDataItem

  end subroutine fpdb_done

!************************************************************************

!<subroutine>

  subroutine fpdb_import(rfpdb)

!<description>
    ! This subroutine imports the persistence database.
!</description>

!<inputoutput>
    ! The persistence database structure
    type(t_fpdb), intent(inout) :: rfpdb
!</inputoutput>

    ! Pointer to the ObjectItem
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectItem

    ! local variables
    integer :: irecordNumber
    
    ! Check that persistence database is completely empty
    if (associated(rfpdb%p_rfpdbObjectTable)) then
      call output_line ('Database must be completely empty!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'fpdb_import')
      call sys_halt()
    end if
    
    ! Phase1: Import persistence database starting at record number 1
    irecordNumber = 1
    import: do 

      ! Pseudo read to determine end-of-file
      read(rfpdb%iunitObjectTable, rec=irecordNumber, err=1)
      
      ! Allocate ObjectItem
      allocate(p_rfpdbObjectItem)

      ! Import ObjectItem from file
      call importObjectItem(p_rfpdbObjectItem, irecordNumber)
      
      ! Insert ObjectItem into object table
      call fpdb_insertObject(rfpdb, p_rfpdbObjectItem, .false.)

      ! Increase record number by one
      irecordNumber = irecordNumber+1
    end do import

    ! Phase2: Regenerate internal ObjectItem links
1   if (associated(rfpdb%p_rfpdbObjectTable)) then
      call relinkObjectItem(rfpdb%p_rfpdbObjectTable)
    end if


  contains
    
    ! Here, the real import routines follow.

    !**************************************************************
    ! Imports the nodes of the object database
    
    subroutine importObjectItem(rfpdbObjectItem, irecordNumber)

      ! The ObjectItem to be imported
      type(t_fpdbObjectItem), intent(out) :: rfpdbObjectItem

      ! Record number of the ObjectItem
      integer, intent(in) :: irecordNumber      
      
      ! local variables
      character(LEN=36) :: suuid
      integer :: i

      ! Import object from file using the first record number
      read(rfpdb%iunitObjectTable, rec=irecordNumber, err=1)&
            suuid,&
            rfpdbObjectItem%stype,&
            rfpdbObjectItem%sname,&
            rfpdbObjectItem%ifirstDataRecnum,&
            rfpdbObjectItem%ilengthDataRecnum

      ! Recreate UUID from string representation
      call uuid_createUUID(suuid, rfpdbObjectItem%ruuid)
    
      ! Allocate pointer for data
      allocate(rfpdbObjectItem%p_RfpdbDataItem(rfpdbObjectItem%ilengthdatarecnum))
      
      ! Import DataItems from file
      do i = 1, rfpdbObjectItem%ilengthDataRecnum
        call importDataItem(rfpdbObjectItem%p_RfpdbDataItem(i),&
                            rfpdbObjectItem%ifirstDataRecnum+i-1)
      end do

      ! That's it
      return
      

1     call output_line ('Unable to import object into table!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'importObjectItem')
      call sys_halt()

    end subroutine importObjectItem


    !**************************************************************
    ! Import of the data table
    ! Each DataItem is imported into the database.

    subroutine importDataItem(rfpdbDataItem, irecordNumber)

      integer, intent(in)               :: irecordNumber
      type(t_fpdbDataItem), intent(out) :: rfpdbDataItem

      ! local variables
      character(LEN=36) :: suuid
      integer :: iunit

      ! Read datatype from file
      read(rfpdb%iunitDataTable, rec=irecordNumber, err=1) rfpdbDataItem%ctype

      ! What data item are we?
      select case(rfpdbDataItem%ctype)

        ! ----------------------------------------------------------------------
        ! Below, we deal with data items which require special treatment
        ! ----------------------------------------------------------------------

      case (FPDB_NULL)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname
        
      case (FPDB_OBJECT,&
            FPDB_LINK)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname,&
             suuid
        call uuid_createUUID(suuid, rfpdbDataItem%ruuid)

        ! ----------------------------------------------------------------------
        ! Below, we deal with atomic data:
        !    The data item itself is read from the data table
        ! ----------------------------------------------------------------------

      case (FPDB_SINGLE)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname,&
             rfpdbDataItem%fsingle

      case (FPDB_DOUBLE)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname,&
             rfpdbDataItem%ddouble

      case (FPDB_INT)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname,&
             rfpdbDataItem%iinteger
      
      case (FPDB_LOGICAL)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname,&
             rfpdbDataItem%blogical

      case (FPDB_CHAR)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname,&
             rfpdbDataItem%schar
        
        ! ----------------------------------------------------------------------
        ! Below, we deal with 1D data:
        !    The UUID of the filename containing the data item is read from
        !    the data table. The actual data is read from a separate file.
        ! ----------------------------------------------------------------------

      case (FPDB_SINGLE1D,&
            FPDB_DOUBLE1D,&
            FPDB_INT1D,&
            FPDB_LOGICAL1D,&
            FPDB_CHAR1D)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname, suuid
        
        ! Set filename of data file
        rfpdbDataItem%sfilename = trim(rfpdb%spath)//suuid//'.fpd'

        ! Read content of UUID from data file
        call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                                   iunit, bformatted=.FALSE.)
        read(iunit, err=2) rfpdbDataItem%Ilbounds,&
                           rfpdbDataItem%Iubounds
        close(iunit)        

        ! ----------------------------------------------------------------------
        ! Below, we deal with 2D data:
        !    The UUID of the filename containing the data item is read from
        !    the data table. The actual data is read from a separate file.
        ! ----------------------------------------------------------------------

      case (FPDB_SINGLE2D,&
            FPDB_DOUBLE2D,&
            FPDB_INT2D,&
            FPDB_LOGICAL2D,&
            FPDB_CHAR2D)
        read(rfpdb%iunitDataTable, rec=irecordNumber, err=1)&
             rfpdbDataItem%ctype, rfpdbDataItem%sname, suuid

        ! Set filename of data file
        rfpdbDataItem%sfilename = trim(rfpdb%spath)//suuid//'.fpd'

        ! Read content of handle from data file
        call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                                   iunit, bformatted=.FALSE.)
        read(iunit, err=2) rfpdbDataItem%Ilbounds,&
                           rfpdbDataItem%Iubounds
        close(iunit)


      case DEFAULT
        call output_line ('Undefined data type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'importDataItem')
        call sys_halt()
      end select


      ! That's it
      return
      
      
1     call output_line ('Unable to import data into table!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'importDataItem')
      call sys_halt()

2     call output_line ('Unable to import data from file!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'importDataItem')
      call sys_halt()

    end subroutine importDataItem

    !**************************************************************
    ! Preorder traversal of the object table.
    ! Each DataItem which contains an ObjectItem is 
    ! relinked to the corresponding ObjectItem
    
    recursive subroutine relinkObjectItem(rfpdbObjectItem)
      
      type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem

      ! local variables
      type(t_fpdbDataItem), pointer :: p_rfpdbDataItem
      integer :: i

      ! Do we have associated DataItems?
      if (associated(rfpdbObjectItem%p_RfpdbDataItem)) then
        
        ! Loop over all DataItems
        do i = 1, size(rfpdbObjectItem%p_RfpdbDataItem)

          ! Set pointer
          p_rfpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(i)
          
          if ((p_rfpdbDataItem%ctype .eq. FPDB_OBJECT) .or.&
              (p_rfpdbDataItem%ctype .eq. FPDB_LINK)) then

            ! Retrieve corresponding ObjectItem by UUID
            call fpdb_retrieveObject(rfpdb, p_rfpdbDataItem%ruuid,&
                                     p_rfpdbDataItem%p_fpdbObjectItem)
            call uuid_createUUID(0, p_rfpdbDataItem%ruuid)
            if (.not.associated(p_rfpdbDataItem%p_fpdbObjectItem)) then
              call output_line ('Inconsistent database structure!', &
                                OU_CLASS_ERROR,OU_MODE_STD,'relinkObjectItem')
              call sys_halt()
            end if
          end if

        end do
      end if

      ! Proceed to left child
      if (associated(rfpdbObjectItem%p_rfpdbObjectLeft)) then
        call relinkObjectItem(rfpdbObjectItem%p_rfpdbObjectLeft)
      end if

      ! Proceed to right child
      if (associated(rfpdbObjectItem%p_rfpdbObjectRight)) then
        call relinkObjectItem(rfpdbObjectItem%p_rfpdbObjectRight)
      end if
      
    end subroutine relinkObjectItem

  end subroutine fpdb_import

!************************************************************************

!<subroutine>

  subroutine fpdb_export(rfpdb)

!<description>
    ! This subroutine exports the persistence database.
    ! A preorder traversal of the object database is performed and
    ! written to file. For
!</description>

!<inputoutput>
    ! The persistence database structure
    type(t_fpdb), intent(inout) :: rfpdb
!</inputoutput>

    ! Export object table to file
    if (associated(rfpdb%p_rfpdbObjectTable)) then
      call exportObjectItem(rfpdb%p_rfpdbObjectTable)
    end if
    
  contains

    ! Here, the real export routines follow.

    !**************************************************************
    ! Preorder traversal of the object table.
    ! Each ObjectItem is exported to the object database.
    
    recursive subroutine exportObjectItem(rfpdbObjectItem)
      
      type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem

      ! local variables
      integer :: i

      ! Get current record number of the data table
      rfpdbObjectItem%ifirstDataRecnum = rfpdb%irecnumDataTable

      ! Export associated data items to data table file
      if (associated(rfpdbObjectItem%p_RfpdbDataItem)) then
        do i = 1, size(rfpdbObjectItem%p_RfpdbDataItem)
          call exportDataItem(rfpdbObjectItem%p_RfpdbDataItem(i))
        end do
      end if

      ! Calculate the number of table entries for this ObjectItem
      rfpdbObjectItem%ilengthDataRecnum = rfpdb%irecnumDataTable- &
                                          rfpdbObjectItem%ifirstDataRecnum

      ! Export ObjectItem to file using the specified record number
      write(rfpdb%iunitObjectTable, rec=rfpdb%irecnumObjectTable, err=1)&
            uuid_conv2String(rfpdbObjectItem%ruuid),&
            rfpdbObjectItem%stype,&
            rfpdbObjectItem%sname,&
            rfpdbObjectItem%ifirstDataRecnum,&
            rfpdbObjectItem%ilengthDataRecnum

      ! Increase record number of object table
      rfpdb%irecnumObjectTable = rfpdb%irecnumObjectTable+1

      ! Proceed to left child
      if (associated(rfpdbObjectItem%p_rfpdbObjectLeft)) then
        call exportObjectItem(rfpdbObjectItem%p_rfpdbObjectLeft)
      end if

      ! Proceed to right child
      if (associated(rfpdbObjectItem%p_rfpdbObjectRight)) then
        call exportObjectItem(rfpdbObjectItem%p_rfpdbObjectRight)
      end if

      ! That's it
      return


1     call output_line ('Unable to export object to table!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exportObjectItem')
      call sys_halt()

    end subroutine exportObjectItem

    !**************************************************************
    ! Export the DataItem to the data table file.

    subroutine exportDataItem(rfpdbDataItem)

      type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem

      ! local variables
      integer      :: iunit

      ! What data item are we?
      select case(rfpdbDataItem%ctype)
        
        ! ----------------------------------------------------------------------
        ! Below, we deal with data items which require special treatment
        ! ----------------------------------------------------------------------

      case (FPDB_NULL)
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
        
      case (FPDB_OBJECT,&
            FPDB_LINK)
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%p_fpdbObjectItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
        
        ! ----------------------------------------------------------------------
        ! Below, we deal with atomic data:
        !    The data item itself is written to the data table
        ! ----------------------------------------------------------------------

      case (FPDB_SINGLE)
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              rfpdbDataItem%fsingle
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1

      case (FPDB_DOUBLE)
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              rfpdbDataItem%ddouble
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1

      case (FPDB_INT)
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              rfpdbDataItem%iinteger
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
      
      case (FPDB_LOGICAL)
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              rfpdbDataItem%blogical
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1

      case (FPDB_CHAR)
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              rfpdbDataItem%schar
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1

        ! ----------------------------------------------------------------------
        ! Below, we deal with 1D data:
        !    The UUID of the filename containing the data item is written to
        !    the data table. The actual data is written to a separate file.
        ! ----------------------------------------------------------------------

      case (FPDB_SINGLE1D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1) = lbound(rfpdbDataItem%p_Fsingle1D,1)
          rfpdbDataItem%Iubounds(1) = ubound(rfpdbDataItem%p_Fsingle1D,1)

          ! Write content of handle to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Fsingle1D
          close(iunit)
        end if

        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
        

      case (FPDB_DOUBLE1D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1) = lbound(rfpdbDataItem%p_Ddouble1D,1)
          rfpdbDataItem%Iubounds(1) = ubound(rfpdbDataItem%p_Ddouble1D,1)

          ! Write content of handle to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Ddouble1D
          close(iunit)
        end if

        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1



      case (FPDB_INT1D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1) = lbound(rfpdbDataItem%p_Iinteger1D,1)
          rfpdbDataItem%Iubounds(1) = ubound(rfpdbDataItem%p_Iinteger1D,1)

          ! Write content of handle to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Iinteger1D
          close(iunit)
        end if

        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
        

      case (FPDB_LOGICAL1D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1) = lbound(rfpdbDataItem%p_Blogical1D,1)
          rfpdbDataItem%Iubounds(1) = ubound(rfpdbDataItem%p_Blogical1D,1)

          ! Write content of handle to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Blogical1D
          close(iunit)
        end if

        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
       

      case (FPDB_CHAR1D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1) = lbound(rfpdbDataItem%p_Schar1D,1)
          rfpdbDataItem%Iubounds(1) = ubound(rfpdbDataItem%p_Schar1D,1)

          ! Write content of handle to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Schar1D
          close(iunit)
        end if

        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1


        ! ----------------------------------------------------------------------
        ! Below, we deal with 2D data:
        !    The UUID of the filename containing the data item is written to
        !    the data table. The actual data is written to a separate file.
        ! ----------------------------------------------------------------------
        
      case (FPDB_SINGLE2D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1:2) = lbound(rfpdbDataItem%p_Fsingle2D)
          rfpdbDataItem%Iubounds(1:2) = ubound(rfpdbDataItem%p_Fsingle2D)

          ! Write content of UUIDs to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Fsingle2D
          close(iunit)
        end if
          
        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1   
       

      case (FPDB_DOUBLE2D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1:2) = lbound(rfpdbDataItem%p_Ddouble2D)
          rfpdbDataItem%Iubounds(1:2) = ubound(rfpdbDataItem%p_Ddouble2D)

          ! Write content of UUIDs to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Ddouble2D
          close(iunit)
        end if
          
        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1


      case (FPDB_INT2D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1:2) = lbound(rfpdbDataItem%p_Iinteger2D)
          rfpdbDataItem%Iubounds(1:2) = ubound(rfpdbDataItem%p_Iinteger2D)

          ! Write content of UUIDs to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Iinteger2D
          close(iunit)
        end if
          
        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
       

      case (FPDB_LOGICAL2D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1:2) = lbound(rfpdbDataItem%p_Blogical2D)
          rfpdbDataItem%Iubounds(1:2) = ubound(rfpdbDataItem%p_Blogical2D)

          ! Write content of UUIDs to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Blogical2D
          close(iunit)
        end if
          
        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
        
        
      case (FPDB_CHAR2D)
        ! Check if this DataItem has been exported before
        if (uuid_isNil(rfpdbDataItem%ruuid)) then
          ! If not create a new UUID
          call uuid_createUUID(4, rfpdbDataItem%ruuid)

          ! Fill array of upper/lower bounds
          rfpdbDataItem%Ilbounds(1:2) = lbound(rfpdbDataItem%p_Schar2D)
          rfpdbDataItem%Iubounds(1:2) = ubound(rfpdbDataItem%p_Schar2D)

          ! Write content of UUIDs to data file
          call io_openFileForWriting(trim(rfpdb%spath)//&
                                     uuid_conv2String(rfpdbDataItem%ruuid)//'.fpd',&
                                     iunit, SYS_REPLACE, bformatted=.false.)
          write(iunit, err=2) rfpdbDataItem%Ilbounds,&
                              rfpdbDataItem%Iubounds,&
                              rfpdbDataItem%p_Schar2D
          close(iunit)
        end if
          
        ! Write the DataItem to data table file
        write(rfpdb%iunitDataTable, rec=rfpdb%irecnumDataTable, err=1)&
              rfpdbDataItem%ctype, rfpdbDataItem%sname,&
              uuid_conv2String(rfpdbDataItem%ruuid)
        rfpdb%irecnumDataTable = rfpdb%irecnumDataTable+1
        
        
      case DEFAULT
        call output_line ('Undefined data type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'exportDataItem')
        call sys_halt()
      end select

      ! That's it
      return


1     call output_line ('Unable to export data to table!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exportDataItem')
      call sys_halt()

2     call output_line ('Unable to export data to file!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'exportDataItem')
      call sys_halt()
        
    end subroutine exportDataItem
    
  end subroutine fpdb_export

!************************************************************************

!<subroutine>

  subroutine fpdb_info(rfpdb)

!<description>
    ! This subroutine prints statistics about the persistence database.
!</description>

!<input>
    ! Database structure
    type(t_fpdb), intent(in) :: rfpdb    
!</input>
!</subroutine>

    ! Print statistics
    if (associated(rfpdb%p_rfpdbObjectTable)) then
      call infoObjectItem(rfpdb%p_rfpdbObjectTable)
    end if

  contains
    
    ! Here, the statistics routines follow.
    
    !**************************************************************
    ! Inorder traversal of the object table.
    
    recursive subroutine infoObjectItem(rfpdbObjectItem)

      type(t_fpdbObjectItem), intent(in) :: rfpdbObjectItem

      ! local variables
      integer :: i

      ! Proceed to left child
      if (associated(rfpdbObjectItem%p_rfpdbObjectLeft)) then
        call infoObjectItem(rfpdbObjectItem%p_rfpdbObjectLeft)
      end if
  
      call output_line('ObjectItem: '//uuid_conv2String(rfpdbObjectItem%ruuid))
      call output_line(' stype: '//trim(rfpdbObjectItem%stype)//&
                       ' sname: '//trim(rfpdbObjectItem%sname))
      call output_line('------------------------------------------------')
      if (associated(rfpdbObjectItem%p_RfpdbDataItem)) then
        do i = 1, size(rfpdbObjectItem%p_RfpdbDataItem)
          call output_line('DataItem: '//trim(sys_siL(i,6))//&
                           ' ctype: '//trim(ctype2String(rfpdbObjectItem%p_RfpdbDataItem(i)%ctype))//&
                           ' sname: '//trim(rfpdbObjectItem%p_RfpdbDataItem(i)%sname))
        end do
        call output_line('------------------------------------------------')
      end if
      call output_lbrk()
    
      ! Proceed to right child
      if (associated(rfpdbObjectItem%p_rfpdbObjectRight)) then
        call infoObjectItem(rfpdbObjectItem%p_rfpdbObjectRight)
      end if
      
    end subroutine infoObjectItem

    function ctype2String(ctype) result(str)

      integer, intent(in) :: ctype
      character(len=SYS_STRLEN) :: str

      select case(ctype)
      case (FPDB_NULL)
        str = 'FPDB_NULL'
      case (FPDB_OBJECT)
        str = 'FPDB_OBJECT'
      case (FPDB_LINK)
        str = 'FPDB_LINK'
      case (FPDB_SINGLE)
        str = 'FPDB_SINGLE'
      case (FPDB_DOUBLE)
        str = 'FPDB_DOUBLE'
      case (FPDB_INT)
        str = 'FPDB_INT'
      case (FPDB_LOGICAL)
        str = 'FPDB_LOGICAL'
      case (FPDB_CHAR)
        str = 'FPDB_CHAR'
      case (FPDB_SINGLE1D)
        str = 'FPDB_SINGLE1D'
      case (FPDB_DOUBLE1D)
        str = 'FPDB_DOUBLE1D'
      case (FPDB_INT1D)
        str = 'FPDB_INT1D'
      case (FPDB_LOGICAL1D)
        str = 'FPDB_LOGICAL1D'
      case (FPDB_CHAR1D)
        str = 'FPDB_CHAR1D'
      case (FPDB_SINGLE2D)
        str = 'FPDB_SINGLE2D'
      case (FPDB_DOUBLE2D)
        str = 'FPDB_DOUBLE2D'
      case (FPDB_INT2D)
        str = 'FPDB_INT2D'
      case (FPDB_LOGICAL2D)
        str = 'FPDB_LOGICAL2D'
      case (FPDB_CHAR2D)
        str = 'FPDB_CHAR2D'
      case DEFAULT
        str = 'unknown'
      end select
    end function ctype2String

    end subroutine fpdb_info

!************************************************************************

!<subroutine>
  
  recursive subroutine fpdb_insertObject(rfpdb, rfpdbObjectItem, brecursive)

!<description>
    ! This subroutine inserts an item into the object table.
    ! The object table is stored as a binary tree, hence,
    ! standard search-insert algorithm is applied.
!</description>

!<input>
    ! Flag: if TRUE objects are inserted recursively
    logical, intent(in) :: brecursive
!</input>

!<inputoutput>
    ! The persistence database structure
    type(t_fpdb), intent(inout), target :: rfpdb

    ! The item that is inserted
    type(t_fpdbObjectItem), intent(inout), target :: rfpdbObjectItem
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectItem
    integer :: i
    
    ! Check if the tree is empty
    if (.not.associated(rfpdb%p_rfpdbObjectTable)) then

      ! Set ObjectItem as new root of the tree
      rfpdb%p_rfpdbObjectTable => rfpdbObjectItem
      
      ! Nullify all pointers
      nullify(rfpdbObjectItem%p_rfpdbObjectFather)
      nullify(rfpdbObjectItem%p_rfpdbObjectLeft)
      nullify(rfpdbObjectItem%p_rfpdbObjectRight)

      ! Do we have to insert objects recursively?
      if (.not.brecursive) return
      
      ! Check if DataItems contain ObjectItems
      if (associated(rfpdbObjectItem%p_RfpdbDataItem)) then
        do i = 1, size(rfpdbObjectItem%p_RfpdbDataItem)
          if ((rfpdbObjectItem%p_RfpdbDataItem(i)%ctype .eq. FPDB_OBJECT) .and.&
              (associated(rfpdbObjectItem%p_RfpdbDataItem(i)%p_fpdbObjectItem))) then
            call fpdb_insertObject(rfpdb,&
                 rfpdbObjectItem%p_RfpdbDataItem(i)%p_fpdbObjectItem,&
                 brecursive)
          end if
        end do
      end if

    else

      ! Set pointer to root of the tree
      p_rfpdbObjectItem => rfpdb%p_rfpdbObjectTable

      ! Walk down the tree until a leaf is found
      do while(associated(p_rfpdbObjectItem))
        
        if (p_rfpdbObjectItem%ruuid .eq. rfpdbObjectItem%ruuid) then
          call output_line ('Object is already present in database!', &
                            OU_CLASS_WARNING,OU_MODE_STD,'fpdb_insertObject')

          ! ObjectItem is already present in tree
          return

        elseif (p_rfpdbObjectItem%ruuid .lt. rfpdbObjectItem%ruuid) then
          
          if (associated(p_rfpdbObjectItem%p_rfpdbObjectRight)) then
            ! Proceed to right child
            p_rfpdbObjectItem => p_rfpdbObjectItem%p_rfpdbObjectRight
          else
            ! Attach object to right leaf and make item new leaf
            p_rfpdbObjectItem%p_rfpdbObjectRight => rfpdbObjectItem
            rfpdbObjectItem%p_rfpdbObjectFather => p_rfpdbObjectItem
            nullify(rfpdbObjectItem%p_rfpdbObjectLeft)
            nullify(rfpdbObjectItem%p_rfpdbObjectRight)

            ! Do we have to insert objects recursively?
            if (.not.brecursive) return
            
            ! Check if DataItems contain ObjectItems
            if (associated(rfpdbObjectItem%p_RfpdbDataItem)) then
              do i = 1, size(rfpdbObjectItem%p_RfpdbDataItem)
                if ((rfpdbObjectItem%p_RfpdbDataItem(i)%ctype .eq. FPDB_OBJECT) .and.&
                    (associated(rfpdbObjectItem%p_RfpdbDataItem(i)%p_fpdbObjectItem))) then
                  call fpdb_insertObject(rfpdb,&
                       rfpdbObjectItem%p_RfpdbDataItem(i)%p_fpdbObjectItem,&
                       brecursive)
                end if
              end do
            end if

            ! That's it
            return
          end if
          
        else
          
          if (associated(p_rfpdbObjectItem%p_rfpdbObjectLeft)) then
            ! Proceed to left child
            p_rfpdbObjectItem => p_rfpdbObjectItem%p_rfpdbObjectLeft
          else
            ! Attach object to left leaf and make item new leaf
            p_rfpdbObjectItem%p_rfpdbObjectLeft => rfpdbObjectItem
            rfpdbObjectItem%p_rfpdbObjectFather => p_rfpdbObjectItem
            nullify(rfpdbObjectItem%p_rfpdbObjectLeft)
            nullify(rfpdbObjectItem%p_rfpdbObjectRight)

            ! Do we have to insert objects recursively?
            if (.not.brecursive) return
            
            ! Check if DataItems contain ObjectItems
            if (associated(rfpdbObjectItem%p_RfpdbDataItem)) then
              do i = 1, size(rfpdbObjectItem%p_RfpdbDataItem)
                if ((rfpdbObjectItem%p_RfpdbDataItem(i)%ctype .eq. FPDB_OBJECT) .and.&
                    (associated(rfpdbObjectItem%p_RfpdbDataItem(i)%p_fpdbObjectItem))) then
                  call fpdb_insertObject(rfpdb,&
                       rfpdbObjectItem%p_RfpdbDataItem(i)%p_fpdbObjectItem,&
                       brecursive)
                end if
              end do
            end if
            
            ! That's it
            return
          end if
          
        end if
      end do
      
      ! If we end up here, some critical error occured
      call output_line ('Critical error occured!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_insertObject')
      call sys_halt()
    end if

  end subroutine fpdb_insertObject

!************************************************************************

!<subroutine>

  subroutine fpdb_removeObject(rfpdb, p_rfpdbObjectItem)

!<description>
    ! This subroutine removes an item from the object table.
    ! Note that this subroutine does not check if the given
    ! ObjectItem is actually stored in the tree. It only
    ! performs the standard remove algorithm, that is, the
    ! largest child (right-most leaf) in the left subtree is
    ! determined and used to replace the removed ObjectItem.    
!</description>

!<inputoutput>
    ! The persistence database structure
    type(t_fpdb), intent(inout) :: rfpdb

    ! The item that is removed
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectItem
!</inputoutput>
!</subroutine>

!!$    ! local variables
!!$    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectItemTmp
!!$
!!$    ! Check if the ObjectItem has a left subtree
!!$    if (associated(p_rfpdbObjectItem%p_rfpdbObjectLeft)) then
!!$
!!$      ! Set pointer to left subtree
!!$      p_rfpdbObjectItemTmp => p_rfpdbObjectItem%p_rfpdbObjectLeft
!!$
!!$      ! Walk down the tree until the right-most leaf is found
!!$      do while (associated(p_rfpdbObjectItemTmp%p_rfpdbObjectRight))
!!$        p_rfpdbObjectItemTmp => p_rfpdbObjectItemTmp%p_rfpdbObjectRight
!!$      end do
!!$
!!$      ! Update ObjectItem by the 
!!$      rfpdbObjectItem%ruuid             = p_rfpdbObjectItem%ruuid
!!$      rfpdbObjectItem%stype             = p_rfpdbObjectItem%stype
!!$      rfpdbObjectItem%sname             = p_rfpdbObjectItem%sname
!!$      rfpdbObjectItem%ifirstDataRecnum  = p_rfpdbObjectItem%ifirstDataRecnum
!!$      rfpdbObjectItem%ilengthDataRecnum = p_rfpdbObjectItem%ilengthDataRecnum
!!$
!!$    end if
!!$
!!$    ! Check if the ObjectItem has a father node
!!$    if (.not.associated(rfpdbObjectItem%p_rfpdbObjectFather)) then
!!$      
!!$      
!!$
!!$    else
!!$
!!$
!!$    end if
      
  end subroutine fpdb_removeObject

!************************************************************************

!<subroutine>

  subroutine fpdb_retrieveObjectByUUID(rfpdb, ruuid, p_rfpdbObjectItem)

!<description>
    ! This subroutine retrieves an item from the object table
    ! based on the Universally unique identifier ruuid.
    ! The object table is stored as a binary tree, hence,
    ! the standard search algorithm is applied.
!</description>

!<input>
    ! The persistence database structure
    type(t_fpdb), intent(in) :: rfpdb

    ! The universally unique identifier
    type(t_uuid), intent(in) :: ruuid
!</input>

!<output>
    ! The pointer associated to the object
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectItem
!</output>
!</subroutine>

    ! Set pointer to the root of the tree
    p_rfpdbObjectItem => rfpdb%p_rfpdbObjectTable
    
    ! Walk down the tree
    do while(associated(p_rfpdbObjectItem))
      
      if (p_rfpdbObjectItem%ruuid .eq. ruuid) then
        
        ! We found the ObjectItem
        return
        
      elseif (p_rfpdbObjectItem%ruuid .lt. ruuid) then
        
        ! Proceed to right child
        p_rfpdbObjectItem => p_rfpdbObjectItem%p_rfpdbObjectRight
        
      else
        
        ! Proceed to left child
        p_rfpdbObjectItem => p_rfpdbObjectItem%p_rfpdbObjectLeft
      end if
      
    end do

  end subroutine fpdb_retrieveObjectByUUID

!************************************************************************

!<subroutine>

  subroutine fpdb_retrieveObjectByName(rfpdb, sname, p_rfpdbObjectItem)

!<description>
    ! This subroutine retrieves an item from the object table
    ! based on the full qualified name given by parameter sname.
!</description>

!<input>
    ! The persistence database structure
    type(t_fpdb), intent(in) :: rfpdb

    ! The name of the ObjectItem to be retrieved
    character(len=*), intent(in) :: sname
!</input>

!<output>
    ! The pointer associated to the object
    type(t_fpdbObjectItem), pointer :: p_rfpdbObjectItem
!</output>
!</subroutine>

    ! local variables
    logical :: bfound

    ! Check if full qualified name is empty
    if (sname .eq. '') then
      call output_line ('Full qualified name must not be empty!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_retrieveObjectByName')
      call sys_halt()
    end if

    ! Initialisation
    nullify(p_rfpdbObjectItem)
    bfound = .false.

    ! Inorder traversal of the object table
    if (associated(rfpdb%p_rfpdbObjectTable)) then
      call searchObjectItem(rfpdb%p_rfpdbObjectTable)
    end if

  contains
    
    ! Here, the real searching routine follows

    !**************************************************************
    ! Preorder traversal of the object table.

    recursive subroutine searchObjectItem(rfpdbObjectItem)
      
      type(t_fpdbObjectItem), intent(in), target :: rfpdbObjectItem
      
      if (trim(rfpdbObjectItem%sname) .eq. sname) then
        p_rfpdbObjectItem => rfpdbObjectItem
        
        ! That's it
        bfound = .true.
        return
      end if

      ! Proceed to left child
      if (associated(rfpdbObjectItem%p_rfpdbObjectLeft)) then
        call searchObjectItem(rfpdbObjectItem%p_rfpdbObjectLeft)
        if (bfound) return
      end if

      ! Proceed to right child
      if (associated(rfpdbObjectItem%p_rfpdbObjectRight)) then
        call searchObjectItem(rfpdbObjectItem%p_rfpdbObjectRight)
        if (bfound) return
      end if
      
    end subroutine searchObjectItem

  end subroutine fpdb_retrieveObjectByName

!************************************************************************

!<function>
  
  function fpdb_getObjectLength() result(iolength)

!<description>
    ! This function computes the record length of the ObjectItem.
!</description>

!<result>
    integer :: iolength
!</result>
!</function>
    
    ! local variables
    type(t_fpdbObjectItem) :: rfpdbObjectItem

    inquire(iolength=iolength) uuid_conv2String(rfpdbObjectItem%ruuid),&
                               rfpdbObjectItem%stype,&
                               rfpdbObjectItem%sname,&
                               rfpdbObjectItem%ifirstDataRecnum,&
                               rfpdbObjectItem%ilengthDataRecnum

  end function fpdb_getObjectLength

!************************************************************************

!<function>
  
  function fpdb_getDataLength() result(iolength)

!<description>
    ! This function computes the record length of the DataItem.
!</description>

!<result>
    integer :: iolength
!</result>
!</function>
    
    ! local variables
    type(t_fpdbDataItem) :: rfpdbDataItem
    
    inquire(iolength=iolength) rfpdbDataItem%ctype,&
                               rfpdbDataItem%sname,&
                               uuid_conv2String(rfpdbDataItem%ruuid)
    
  end function fpdb_getDataLength

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_single1d(rfpdbDataItem, Farray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the single float array. If not given, 
    ! the pointer of the DataItem is used instead.
    real(SP), dimension(:), intent(out), target, optional :: Farray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    real(SP), dimension(:), pointer :: p_Fsingle
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_SINGLE1D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_single1d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Farray)) then
      p_Fsingle => Farray
    else
      if (.not.associated(rfpdbDataItem%p_Fsingle1D)) then
        allocate(rfpdbDataItem%p_Fsingle1D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1)))
      end if
      p_Fsingle => rfpdbDataItem%p_Fsingle1D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Fsingle,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Fsingle,1) .ne. rfpdbDataItem%Iubounds(1))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_single1d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Fsingle
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_single1d')
    call sys_halt()
  end subroutine fpdb_getdata_single1d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_double1d(rfpdbDataItem, Darray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the double float array. If not given, 
    ! the pointer of the DataItem is used instead.
    real(DP), dimension(:), intent(out), target, optional :: Darray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    real(DP), dimension(:), pointer :: p_Ddouble
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_DOUBLE1D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_double1d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Darray)) then
      p_Ddouble => Darray
    else
      if (.not.associated(rfpdbDataItem%p_Ddouble1D)) then
        allocate(rfpdbDataItem%p_Ddouble1D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1)))
      end if
      p_Ddouble => rfpdbDataItem%p_Ddouble1D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Ddouble,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Ddouble,1) .ne. rfpdbDataItem%Iubounds(1))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_double1d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Ddouble
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_double1d')
    call sys_halt()
  end subroutine fpdb_getdata_double1d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_int1d(rfpdbDataItem, Iarray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the integer array. If not given, 
    ! the pointer of the DataItem is used instead.
    integer, dimension(:), intent(out), target, optional :: Iarray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    integer, dimension(:), pointer :: p_Iinteger
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_INT1D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_int1d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Iarray)) then
      p_Iinteger => Iarray
    else
      if (.not.associated(rfpdbDataItem%p_Iinteger1D)) then
        allocate(rfpdbDataItem%p_Iinteger1D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1)))
      end if
      p_Iinteger => rfpdbDataItem%p_Iinteger1D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Iinteger,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Iinteger,1) .ne. rfpdbDataItem%Iubounds(1))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_int1d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Iinteger
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_int1d')
    call sys_halt()
  end subroutine fpdb_getdata_int1d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_logical1d(rfpdbDataItem, Larray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the logical array. If not given, 
    ! the pointer of the DataItem is used instead.
    logical, dimension(:), intent(out), target, optional :: Larray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    logical, dimension(:), pointer :: p_Blogical
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_LOGICAL1D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_logical1d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Larray)) then
      p_Blogical => Larray
    else
      if (.not.associated(rfpdbDataItem%p_Blogical1D)) then
        allocate(rfpdbDataItem%p_Blogical1D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1)))
      end if
      p_Blogical => rfpdbDataItem%p_Blogical1D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Blogical,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Blogical,1) .ne. rfpdbDataItem%Iubounds(1))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_logical1d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Blogical
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_blogical1d')
    call sys_halt()
  end subroutine fpdb_getdata_logical1d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_char1d(rfpdbDataItem, Sarray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the character array. If not given, 
    ! the pointer of the DataItem is used instead.
    character, dimension(:), intent(out), target, optional :: Sarray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    character, dimension(:), pointer :: p_Schar
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_CHAR1D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_char1d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Sarray)) then
      p_Schar => Sarray
    else
      if (.not.associated(rfpdbDataItem%p_Schar1D)) then
        allocate(rfpdbDataItem%p_Schar1D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1)))
      end if
      p_Schar => rfpdbDataItem%p_Schar1D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Schar,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Schar,1) .ne. rfpdbDataItem%Iubounds(1))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_char1d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Schar
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_char1d')
    call sys_halt()
  end subroutine fpdb_getdata_char1d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_single2d(rfpdbDataItem, Farray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the single float array. If not given, 
    ! the pointer of the DataItem is used instead.
    real(SP), dimension(:,:), intent(out), target, optional :: Farray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    real(SP), dimension(:,:), pointer :: p_Fsingle
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_SINGLE2D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_single2d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Farray)) then
      p_Fsingle => Farray
    else
      if (.not.associated(rfpdbDataItem%p_Fsingle2D)) then
        allocate(rfpdbDataItem%p_Fsingle2D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1),&
                 rfpdbDataItem%Ilbounds(2):rfpdbDataItem%Iubounds(2)))
      end if
      p_Fsingle => rfpdbDataItem%p_Fsingle2D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Fsingle,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Fsingle,1) .ne. rfpdbDataItem%Iubounds(1)) .or.&
        (lbound(p_Fsingle,2) .ne. rfpdbDataItem%Ilbounds(2)) .or.&
        (ubound(p_Fsingle,2) .ne. rfpdbDataItem%Iubounds(2))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_single2d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Fsingle
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_single2d')
    call sys_halt()
  end subroutine fpdb_getdata_single2d  

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_double2d(rfpdbDataItem, Darray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the double float array. If not given, 
    ! the pointer of the DataItem is used instead.
    real(DP), dimension(:,:), intent(out), target, optional :: Darray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    real(DP), dimension(:,:), pointer :: p_Ddouble
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_DOUBLE2D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_double2d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Darray)) then
      p_Ddouble => Darray
    else
      if (.not.associated(rfpdbDataItem%p_Ddouble2D)) then
        allocate(rfpdbDataItem%p_Ddouble2D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1),&
                 rfpdbDataItem%Ilbounds(2):rfpdbDataItem%Iubounds(2)))
      end if
      p_Ddouble => rfpdbDataItem%p_Ddouble2D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Ddouble,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Ddouble,1) .ne. rfpdbDataItem%Iubounds(1)) .or.&
        (lbound(p_Ddouble,2) .ne. rfpdbDataItem%Ilbounds(2)) .or.&
        (ubound(p_Ddouble,2) .ne. rfpdbDataItem%Iubounds(2))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_ddouble2d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Ddouble
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_ddouble2d')
    call sys_halt()
  end subroutine fpdb_getdata_double2d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_int2d(rfpdbDataItem, Iarray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the integer array. If not given, 
    ! the pointer of the DataItem is used instead.
    integer, dimension(:,:), intent(out), target, optional :: Iarray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    integer, dimension(:,:), pointer :: p_Iinteger
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_INT2D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_int2d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Iarray)) then
      p_Iinteger => Iarray
    else
      if (.not.associated(rfpdbDataItem%p_Iinteger2D)) then
        allocate(rfpdbDataItem%p_Iinteger2D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1),&
                 rfpdbDataItem%Ilbounds(2):rfpdbDataItem%Iubounds(2)))
      end if
      p_Iinteger => rfpdbDataItem%p_Iinteger2D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Iinteger,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Iinteger,1) .ne. rfpdbDataItem%Iubounds(1)) .or.&
        (lbound(p_Iinteger,2) .ne. rfpdbDataItem%Ilbounds(2)) .or.&
        (ubound(p_Iinteger,2) .ne. rfpdbDataItem%Iubounds(2))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_int2d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Iinteger
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_int2d')
    call sys_halt()
  end subroutine fpdb_getdata_int2d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_logical2d(rfpdbDataItem, Larray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the logical array. If not given, 
    ! the pointer of the DataItem is used instead.
    logical, dimension(:,:), intent(out), target, optional :: Larray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    logical, dimension(:,:), pointer :: p_Blogical
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_LOGICAL2D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_logical2d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Larray)) then
      p_Blogical => Larray
    else
      if (.not.associated(rfpdbDataItem%p_Blogical2D)) then
        allocate(rfpdbDataItem%p_Blogical2D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1),&
                 rfpdbDataItem%Ilbounds(2):rfpdbDataItem%Iubounds(2)))
      end if
      p_Blogical => rfpdbDataItem%p_Blogical2D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Blogical,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Blogical,1) .ne. rfpdbDataItem%Iubounds(1)) .or.&
        (lbound(p_Blogical,2) .ne. rfpdbDataItem%Ilbounds(2)) .or.&
        (ubound(p_Blogical,2) .ne. rfpdbDataItem%Iubounds(2))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_logical2d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Blogical
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_logical2d')
    call sys_halt()
  end subroutine fpdb_getdata_logical2d

!************************************************************************

!<subroutine>

  subroutine fpdb_getdata_char2d(rfpdbDataItem, Sarray)

!<description>
    ! This subroutine imports the data associated with the DataItem
!</description>

!<inputoutput>
    ! The DataItem to be imported
    type(t_fpdbDataItem), intent(inout) :: rfpdbDataItem
!</inputoutput>

!<output>
    ! OPTIONAL: the character array. If not given, 
    ! the pointer of the DataItem is used instead.
    character, dimension(:,:), intent(out), target, optional :: Sarray
!</output>
!</subroutine>

    ! local variable
    integer, dimension(FPDB_MAXDIM) :: Ilbounds,Iubounds
    character, dimension(:,:), pointer :: p_Schar
    integer :: iunit

    ! Check if data type is compatible
    if (rfpdbDataItem%ctype .ne. FPDB_CHAR2D) then
      call output_line ('Wrong data format!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_char2d')
      call sys_halt()
    end if
    
    ! Set the pointer to the data array
    if (present(Sarray)) then
      p_Schar => Sarray
    else
      if (.not.associated(rfpdbDataItem%p_Schar2D)) then
        allocate(rfpdbDataItem%p_Schar2D(&
                 rfpdbDataItem%Ilbounds(1):rfpdbDataItem%Iubounds(1),&
                 rfpdbDataItem%Ilbounds(2):rfpdbDataItem%Iubounds(2)))
      end if
      p_Schar => rfpdbDataItem%p_Schar2D
    end if
    
    ! Check if data array is compatible
    if ((lbound(p_Schar,1) .ne. rfpdbDataItem%Ilbounds(1)) .or.&
        (ubound(p_Schar,1) .ne. rfpdbDataItem%Iubounds(1)) .or.&
        (lbound(p_Schar,2) .ne. rfpdbDataItem%Ilbounds(2)) .or.&
        (ubound(p_Schar,2) .ne. rfpdbDataItem%Iubounds(2))) then
      call output_line ('Invalid structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_char2d')
      call sys_halt()
    end if

    ! Read content of UUID from data file
    call io_openFileForReading(trim(rfpdbDataItem%sfilename),&
                               iunit, bformatted=.FALSE.)
    read(iunit, err=1) Ilbounds, Iubounds, p_Schar
    close(iunit)

    ! That's it
    return

1   call output_line ('Unable to import data from file!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'fpdb_getdata_char2d')
    call sys_halt()
  end subroutine fpdb_getdata_char2d

end module fpersistence
