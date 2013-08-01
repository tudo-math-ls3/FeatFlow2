!##############################################################################
!# ****************************************************************************
!# <name> intgroupset </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a structure that encapsules groups of integers.
!# Such groups can be used, e.g., as column/row structure in CSR matrices
!# or to set up lists with triangle and quad cells in a mesh.
!#
!# The following subroutines can be found here:
!#
!#  1.) igs_initGroupSet
!#      -> Initialises an integer group set.
!#
!#  2.) igs_doneGroupSet
!#      -> Cleans up an integer group set
!#
!#  3.) igs_setItemCount
!#      -> Defines the number of items in each group
!#
!#  4.) igs_startItemCount
!#      -> Starts counting the number of items in each group
!#
!#  5.) igs_finishItemCount
!#      -> Stops/finishes counting the number of items in each group
!#
!#  6.) igs_allocItemGroupIDs / igs_deallocItemGroupIDs
!#      -> Adds/Removes information about the group of each item
!#
!#  7.) igs_calcItemGroupIDs
!#      -> Calculates information about the group of each item
!#
!#  8.) igs_allocItemTags / igs_deallocItemTags
!#      -> Adds/Removes a user defined tag to each item
!#
!#  9.) igs_allocGroupTags / igs_deallocGroupTags
!#      -> Adds/Removes a user defined tag to each group
!#
!#  10.) igs_sortItemsInGroups
!#      -> Sorts the items in all groups
!#
!#  11.) igs_getGroupTags
!#      -> Returns the list of user defined tags of each group
!#
!#  12.) igs_getGroup
!#      -> Retrieves a pointer to a group
!#
!#  13.) igs_getItemGroupIDs
!#      -> Returns a pointer to a list with group IDs for each item
!#
!#  14.) igs_initGroupAppend
!#      -> Initialises appending groups
!#
!#  15.) igs_appendGroup
!#      -> Appends a group to the group set
!#
!#  16.) igs_doneGroupAppend
!#      -> Finishes appending groups
!#
!#  17.) igs_copyGroupSet
!#      -> Copies a group to another one. Probably shares data.
!# </purpose>
!##############################################################################

module intgroupset

  use fsystem
  use genoutput
  use storage
  use sort
  use linearalgebra

  implicit none
  private

!<types>

!<typeblock>

  ! Defines a set of groups of integers.
  type t_intGroupSet
  
    ! Number of groups currently in this set
    integer :: ngroups = 0
    
    ! Maximum allowed number of groups supported by this structure.
    integer :: nmaxGroups = 0
  
    ! Total number of items in the set
    integer :: nitems = 0
    
    ! Total maximum number of items supported by this structure.
    integer :: nmaxItems = 0
    
    ! Maximum number of items per group
    integer :: nmaxItemsPerGroup = 0

    ! Handle to integer array p_Iitems with all items in this set. The items are
    ! sorted for the group and indexed by p_IitemIndex.
    integer :: h_Iitems = ST_NOHANDLE
    
    ! Index array that defines the starting (and possibly stopping) indices of the
    ! different groups in p_Iitems.
    integer :: h_IitemIndex = ST_NOHANDLE
    
    ! Optional: Array that defines for every item the group it is associated to.
    ! Not available if the handle is ST_NOHANDLE.
    integer :: h_IitemGroupIDs = ST_NOHANDLE
    
    ! Optional: Array that defines for every group a corresponding user defined tag.
    integer :: h_IgroupTags = ST_NOHANDLE
    
    ! Optional: Array that defines for every item a corresponding user defined tag.
    integer :: h_IitemTags = ST_NOHANDLE
    
    ! Defines if the group data is shared with another group.
    ! If set to TRUE, the group and item data is shared with another group.
    logical :: bsharedGroupData = .false.
    
    ! Defines if the group tag is shared with another group.
    ! If set to TRUE, the group tags data is shared with another group tag.
    logical :: bsharedGroupTags = .false.
    
    ! Defines if the group tag is shared with another group.
    ! If set to TRUE, the item tags are shared with another item tag array.
    logical :: bsharedItemTags = .false.

    ! Defines if h_IitemIndex defines a simple or double indexed array.
    ! =FALSE: p_IitemIndex has "ngroups+1" entries.
    !         "p_Iitems( p_IitemIndex(i) : p_IitemIndex(i+1)-1 )" contains the items from group i.
    !         There is p_IitemIndex(ngroups+1) = nitems+1.
    ! =TRUE:  p_IitemIndex has "2*ngroups" entries.
    !         "p_Iitems( p_IitemIndex(2*i-1) : p_IitemIndex(2*i) )" contains the items from group i.
    !         There may be "holes" in p_Iitems with reserved space for additional items
    !         added later.
    logical :: bdoubleIndex = .false.
    
  end type


!</typeblock>

  public :: t_intGroupSet

!</types>

  interface igs_initGroupSet
    module procedure igs_initGroupSetSimple
  end interface
  
  interface igs_getGroup
    module procedure igs_getGroupByPointer
    module procedure igs_getGroupByIndex
  end interface

  public :: igs_initGroupSet
  public :: igs_doneGroupSet
  public :: igs_allocItemGroupIDs
  public :: igs_deallocItemGroupIDs
  public :: igs_calcItemGroupIDs
  public :: igs_allocItemTags
  public :: igs_deallocItemTags
  public :: igs_allocGroupTags
  public :: igs_deallocGroupTags
  public :: igs_setItemCount
  public :: igs_startItemCount
  public :: igs_finishItemCount
  public :: igs_getGroup
  public :: igs_getGroupTags
  public :: igs_sortItemsInGroups
  public :: igs_printGroupSet
  public :: igs_copyGroupSet
  public :: igs_appendGroup
  public :: igs_initGroupAppend
  public :: igs_doneGroupAppend

contains

!****************************************************************************

!<subroutine>

  subroutine igs_initGroupSetSimple (rintGroupSet, ngroups, bhasGroupTag)

!<description>
  ! Initialises a new group set with ngroup groups and nitems items.
!</description>

!<input>
  ! Number of groups initially present in the group set.
  integer, intent(in) :: ngroups
  
  ! OPTIONAL: If set to TRUE, a user defined group tag will be associated
  ! to each group. Default is FALSE.
  logical, intent(in), optional :: bhasGroupTag
!</input>

!<output>
  ! Integer group set to be created.
  type(t_intGroupSet), intent(out) :: rintGroupSet
!</output>

!</subroutine>

    ! Define the content of the structure
    rintGroupSet%ngroups = ngroups
    rintGroupSet%nmaxGroups = ngroups
    rintGroupSet%nitems = 0
    rintGroupSet%nmaxItems = 0
    rintGroupSet%nmaxItemsPerGroup = 0

    ! Simply indexed group pointers
    rintGroupSet%bdoubleIndex = .false.

    ! Allocate main memory
    call storage_new ("igs_initGroupSetSimple", &
        "h_IitemIndex", ngroups+1, ST_INT, rintGroupSet%h_IitemIndex, ST_NEWBLOCK_NOINIT)
    
    if (present(bhasGroupTag)) then
      if (bhasGroupTag) then
        call igs_allocGroupTags (rintGroupSet)
      end if
    end if

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_doneGroupSet (rintGroupSet)

!<description>
  ! Cleans up an int group set, releases memory.
!</description>

!<inputoutput>
  ! Integer group set to be cleaned up.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! Define the content of the structure
    rintGroupSet%ngroups = 0
    rintGroupSet%nmaxGroups = 0
    rintGroupSet%nitems = 0
    rintGroupSet%nmaxItems = 0
    rintGroupSet%nmaxItemsPerGroup = 0

    ! Deallocate memory
    if ((.not. rintGroupSet%bsharedGroupData) .and. &
        (rintGroupSet%h_Iitems .ne. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_Iitems)
    end if
    
    if ((.not. rintGroupSet%bsharedGroupData) .and. &
        (rintGroupSet%h_IitemIndex .ne. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_IitemIndex)
    end if

    if ((.not. rintGroupSet%bsharedGroupData) .and. &
        (rintGroupSet%h_IitemGroupIDs .ne. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_IitemGroupIDs)
    end if

    if ((.not. rintGroupSet%bsharedGroupTags) .and. &
        (rintGroupSet%h_IitemTags .ne. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_IitemTags)
    end if

    if ((.not. rintGroupSet%bsharedItemTags) .and. &
        (rintGroupSet%h_IgroupTags .ne. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_IgroupTags)
    end if

    rintGroupSet%bsharedGroupData = .false.
    rintGroupSet%bsharedGroupTags = .false.
    rintGroupSet%bsharedItemTags = .false.
    rintGroupSet%bdoubleIndex = .false.
    
    rintGroupSet%h_Iitems = ST_NOHANDLE
    rintGroupSet%h_IitemIndex = ST_NOHANDLE
    rintGroupSet%h_IitemGroupIDs = ST_NOHANDLE
    rintGroupSet%h_IitemTags = ST_NOHANDLE
    rintGroupSet%h_IgroupTags = ST_NOHANDLE

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_allocItemGroupIDs (rintGroupSet)

!<description>
  ! If not present, information for determining the group of an
  ! item is allocated.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    if (rintGroupSet%h_IitemGroupIDs .eq. ST_NOHANDLE) then
      call storage_new ("igs_initGroupSetSimple", &
          "h_IitemIndex", rintGroupSet%nmaxItems, ST_INT, &
          rintGroupSet%h_IitemGroupIDs, ST_NEWBLOCK_ZERO)
    end if
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_deallocItemGroupIDs (rintGroupSet)

!<description>
  ! Removes information about the group of each item.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    if ((.not. rintGroupSet%bsharedGroupData) .and. &
        (rintGroupSet%h_IitemGroupIDs .ne. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_IitemGroupIDs)
    end if
    
    rintGroupSet%h_IitemGroupIDs = ST_NOHANDLE
    rintGroupSet%bsharedGroupData = .false.
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_calcItemGroupIDs (rintGroupSet)

!<description>
  ! Calculates for every element the group it is associated to.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: igroup,iitem
    integer, dimension(:), pointer :: p_Iitems, p_IitemGroupIDs

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_calcItemGroupIDs")
      call sys_halt()
    end if
    
    if (rintGroupSet%h_IitemGroupIDs .eq. ST_NOHANDLE) then
      ! Auto-allocate
      call igs_allocItemGroupIDs (rintGroupSet)
    end if
    
    ! Get a pointer to the item groups
    call storage_getbase_int (rintGroupSet%h_IitemGroupIDs,p_IitemGroupIDs)
    
    ! Loop through all groups
    do igroup = 1,rintGroupSet%ngroups
    
      ! Get the items
      call igs_getGroup (rintGroupSet,igroup,p_Iitems)
      
      ! Set the group id for each item
      do iitem = 1,ubound (p_Iitems,1)
        p_IitemGroupIDs(p_Iitems(iitem)) = igroup
      end do
    
    end do

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_allocItemTags (rintGroupSet)

!<description>
  ! If not present, information for user defined item tags is allocated.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    if (rintGroupSet%h_IitemGroupIDs .eq. ST_NOHANDLE) then
      call storage_new ("igs_allocItemTags", &
          "h_IitemTags", rintGroupSet%nmaxItems, ST_INT, &
          rintGroupSet%h_IitemTags, ST_NEWBLOCK_ZERO)
    end if
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_deallocItemTags (rintGroupSet)

!<description>
  ! Removes user defined item tags.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    if ((.not. rintGroupSet%bsharedItemTags) .and. &
        (rintGroupSet%h_IitemGroupIDs .ne. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_IitemTags)
    end if
    
    rintGroupSet%h_IitemGroupIDs = ST_NOHANDLE
    rintGroupSet%bsharedItemTags = .false.
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_allocGroupTags (rintGroupSet)

!<description>
  ! If not present, information for user defined item tags is allocated.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    if (rintGroupSet%h_IgroupTags .eq. ST_NOHANDLE) then
      call storage_new ("igs_allocGroupTags", &
          "h_IgroupTags", rintGroupSet%nmaxGroups, ST_INT, rintGroupSet%h_IgroupTags, ST_NEWBLOCK_ZERO)
    end if
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_deallocGroupTags (rintGroupSet)

!<description>
  ! Removes user defined group tags.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    if ((.not. rintGroupSet%bsharedGroupTags) .and. &
        (rintGroupSet%h_IgroupTags .eq. ST_NOHANDLE)) then
      call storage_free (rintGroupSet%h_IgroupTags)
    end if
    
    rintGroupSet%h_IgroupTags = ST_NOHANDLE
    rintGroupSet%bsharedGroupTags = .false.
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_setItemCount (rintGroupSet,IitemsInGroup)

!<description>
  ! Defines for each group in rintGroupSet the number of items in the group.
  ! Allocates memory for the items.
!</description>

!<input>
  ! Integer list that contains for each group the number of items in the group.
  integer, dimension(:), intent(in) :: IitemsInGroup
!</input>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,nmaxItemsPerGroup
    integer, dimension(:), pointer :: p_IitemIndex

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_setItemCount")
      call sys_halt()
    end if

    if (rintGroupSet%h_Iitems .ne. ST_NOHANDLE) then
      call output_line ("Items in rintGroupSet already initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_setItemCount")
      call sys_halt()
    end if

    if (.not. rintGroupSet%bdoubleIndex) then
    
      ! ---------------------
      ! Simple index. 
      ! ---------------------
      
      ! Sum up the number of items to get the item pointer.
      call storage_getbase_int (rintGroupSet%h_IitemIndex,p_IitemIndex)
      
      p_IitemIndex(1) = 1
      nmaxItemsPerGroup = IitemsInGroup(1)
      do i=1,rintGroupSet%ngroups
        p_IitemIndex(i+1) = p_IitemIndex(i) + IitemsInGroup(i)
        nmaxItemsPerGroup = max(IitemsInGroup(i),IitemsInGroup(i))
      end do
    
      ! Save maximum number of items per group
      rintGroupSet%nmaxItemsPerGroup = nmaxItemsPerGroup

      ! Total number of items
      rintGroupSet%nitems = p_IitemIndex(rintGroupSet%ngroups+1)-1
      rintGroupSet%nmaxItems = rintGroupSet%nitems

      ! Allocate memory for the items
      call storage_new ("igs_initGroupSetSimple", &
          "h_Iitems", rintGroupSet%nitems, ST_INT, rintGroupSet%h_Iitems, ST_NEWBLOCK_NOINIT)
      
    else
    
      ! ---------------------
      ! Double index. Not implemented.
      ! ---------------------
    
      ! ...
    end if
    
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_startItemCount (rintGroupSet,p_IitemsInGroup)

!<description>
  ! Initialises an item counter. Used for the initialisation of the elements
  ! in the group. Used as an alternative to igs_setItemCount.
  !
  ! Returns and integer array p_IitemsInGroup initialised
  ! with zero. The calling subroutine has to write into p_IitemsInGroup(i)
  ! the number of items in group i. As soon as p_IitemsInGroup is
  ! done writing to, the calling subroutine has to call
  ! igs_finishItemCount which sets up the rintGroupSet accordingly.
!</description>

!<input>
  ! Integer list. The caller has to write the number of elements
  ! in each group to this array.
  integer, dimension(:), pointer :: p_IitemsInGroup
!</input>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IitemIndex

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_startItemCount")
      call sys_halt()
    end if

    if (rintGroupSet%h_Iitems .ne. ST_NOHANDLE) then
      call output_line ("Items in rintGroupSet already initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_startItemCount")
      call sys_halt()
    end if

    ! Get the basic item-index list
    call storage_clear(rintGroupSet%h_IitemIndex)
    call storage_getbase_int (rintGroupSet%h_IitemIndex,p_IitemIndex)
    
    ! Only the first ngroups items.
    p_IitemsinGroup => p_IitemIndex(1:rintGroupSet%ngroups)

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_finishItemCount (rintGroupSet,p_IitemsInGroup)

!<description>
  ! Finishes counting the items in each group. Has to be called after
  ! igs_startItemCount and after teh caller has counted the number of
  ! items in each group.
  ! Initialises rintGroupSet according to the defined number of items in
  ! each group. The pointer p_IitemsInGroup gets invalid.
!</description>

!<inputoutput>
  ! Integer list with number of entries in the group
  integer, dimension(:), pointer :: p_IitemsInGroup

  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IitemIndex
    integer :: i,ncount,nitems,nmaxItemsPerGroup

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_finishItemCount")
      call sys_halt()
    end if

    if (rintGroupSet%h_Iitems .ne. ST_NOHANDLE) then
      call output_line ("Items in rintGroupSet already initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_finishItemCount")
      call sys_halt()
    end if

    ! Pointer gets invalid.
    nullify(p_IitemsInGroup)

    ! Get the basic item-index list
    call storage_getbase_int (rintGroupSet%h_IitemIndex,p_IitemIndex)

    if (.not. rintGroupSet%bdoubleIndex) then
    
      ! ---------------------
      ! Simple index. 
      ! ---------------------
      
      ! Sum up the number of items to get the item pointer.
      ncount = 1
      nitems = p_IitemIndex(1)
      nmaxItemsPerGroup = nitems
      p_IitemIndex(1) = 1
      
      do i=2,rintGroupSet%ngroups
        ncount = ncount + nitems
        nitems = p_IitemIndex(i)
        nmaxItemsPerGroup = max(nmaxItemsPerGroup,nitems)
        p_IitemIndex(i) = ncount
      end do

      ! Final element in the index array. Should be the number of items + 1.
      ncount = ncount + nitems
      p_IitemIndex(rintGroupSet%ngroups+1) = ncount
      
      ! Save maximum number of items per group
      rintGroupSet%nmaxItemsPerGroup = nmaxItemsPerGroup
    
      ! Total number of items
      rintGroupSet%nitems = ncount-1
      rintGroupSet%nmaxItems = rintGroupSet%nitems

      ! Allocate memory for the items
      call storage_new ("igs_initGroupSetSimple", &
          "h_Iitems", rintGroupSet%nitems, ST_INT, rintGroupSet%h_Iitems, ST_NEWBLOCK_NOINIT)

    else
    
      ! ---------------------
      ! Double index. Not implemented.
      ! ---------------------
    
      ! ...
    end if

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_initGroupAppend (rintGroupSet,nitems)

!<description>
  ! Starts appending groups. Initialises rintGroupSet for nitems items.
  ! Afterwards, all groups can be appended to rintGroupSet by 
  ! subsequent calls to igs_appendGroup.
!</description>

!<input>
  ! Number of items in the group set.
  integer, intent(in) :: nitems
!</input>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_startItemCount")
      call sys_halt()
    end if

    if (rintGroupSet%h_Iitems .ne. ST_NOHANDLE) then
      call output_line ("Items in rintGroupSet already initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_startItemCount")
      call sys_halt()
    end if

    ! Clear the item index counter
    call storage_clear(rintGroupSet%h_IitemIndex)

    ! Initialise index arrays
    rintGroupSet%nitems = 0
    rintGroupSet%ngroups = 0
    rintGroupSet%nmaxItems = nitems

    ! Allocate memory for the items
    call storage_new ("igs_initGroupSetSimple", &
        "h_Iitems", nitems, ST_INT, rintGroupSet%h_Iitems, ST_NEWBLOCK_NOINIT)

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_appendGroup (rintGroupSet,Iitems)

!<description>
  ! Appends the items Iitems to rintGroupSet. igs_startGroupAppend must have
  ! been called in advance.
!</description>

!<input>
  ! List of items to be appended
  integer, dimension(:), intent(in) :: Iitems
!</input>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IitemIndex
    integer, dimension(:), pointer :: p_Iitems
    integer :: i,istartindex

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_appendGroup")
      call sys_halt()
    end if

    if (rintGroupSet%h_Iitems .eq. ST_NOHANDLE) then
      call output_line ("Items in rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_appendGroup")
      call sys_halt()
    end if

    if (rintGroupSet%nitems + ubound(Iitems,1) .gt. rintGroupSet%nmaxItems) then
      call output_line ("Too many elements in Iitems.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_appendGroup")
      call sys_halt()
    end if

    call storage_getbase_int (rintGroupSet%h_IitemIndex,p_IitemIndex)
    call storage_getbase_int (rintGroupSet%h_Iitems,p_Iitems)

    if (.not. rintGroupSet%bdoubleIndex) then

      ! ---------------------
      ! Simple index. 
      ! ---------------------

      ! New group
      rintGroupSet%ngroups = rintGroupSet%ngroups + 1
      istartindex = rintGroupSet%nitems + 1
      p_IitemIndex (rintGroupSet%ngroups) = istartindex 
      
      ! Copy the items
      if (ubound(Iitems,1) .le. 1000) then
        do i=1,ubound(Iitems,1)
          p_Iitems(i-1+istartindex) = Iitems(i)
        end do
      else
        call lalg_copyVectorInt (Iitems(:),p_Iitems(istartindex:),ubound(Iitems,1))
      end if
      
      ! Increase number of items
      rintGroupSet%nitems = rintGroupSet%nitems + ubound(Iitems,1)

    else
    
      ! ---------------------
      ! Double index. Not implemented.
      ! ---------------------

      ! ...
    end if

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_doneGroupAppend (rintGroupSet)

!<description>
  ! Finishes appending groups. Must be called after all froups
  ! are appended to rintGroupSet.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IitemIndex

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_doneGroupAppend")
      call sys_halt()
    end if

    if (rintGroupSet%h_Iitems .eq. ST_NOHANDLE) then
      call output_line ("Items in rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_doneGroupAppend")
      call sys_halt()
    end if

    call storage_getbase_int (rintGroupSet%h_IitemIndex,p_IitemIndex)

    if (.not. rintGroupSet%bdoubleIndex) then

      ! ---------------------
      ! Simple index. 
      ! ---------------------

      ! End pointer
      p_IitemIndex (rintGroupSet%ngroups+1) = rintGroupSet%nitems + 1
      
    else
    
      ! ---------------------
      ! Double index. Not implemented.
      ! ---------------------

      ! ...

    end if

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_getGroupByPointer (rintGroupSet,igroup,p_Iitems,p_IitemTags)

!<description>
  ! Returns a pointer to the items in a group.
!</description>

!<input>
  ! Integer group set to be used.
  type(t_intGroupSet), intent(inout) :: rintGroupSet

  ! Number of the group to retrieve information for.
  integer, intent(in) :: igroup
!</input>

!<output>
  ! Pointer to the items of the group
  integer, dimension(:), pointer, optional :: p_Iitems
  
  ! OPTIONAL: Pointer to the associated item tags.
  integer, dimension(:), pointer, optional :: p_IitemTags
!</output>

!</subroutine>

    ! local variables
    integer :: istart,iend
    integer, dimension(:), pointer :: p_IitemIndex
    
    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_getGroupByPointer")
      call sys_halt()
    end if

    if ((igroup .lt. 1) .or. (igroup .gt. rintGroupSet%ngroups)) then
      call output_line ("igroup out of bounds.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_getGroupByPointer")
      call sys_halt()
    end if

    ! Get the index.
    call storage_getbase_int (rintGroupSet%h_IitemIndex,p_IitemIndex)
    
    if (.not. rintGroupSet%bdoubleIndex) then
    
      ! ---------------------
      ! Simple index
      ! ---------------------
      
      istart = p_IitemIndex(igroup)
      iend = p_IitemIndex(igroup+1)-1
    
    else
    
      ! ---------------------
      ! Double index.
      ! ---------------------
    
      istart = p_IitemIndex(2*igroup-1)
      iend = p_IitemIndex(2*igroup)

    end if
    
    ! Get the items
    if (present(p_Iitems)) then
      nullify(p_Iitems)
      
      if (rintGroupSet%h_Iitems .ne. ST_NOHANDLE) then
        call storage_getbase_int (rintGroupSet%h_Iitems,p_Iitems)
        
        ! Get a sub-pointer to the actual items.
        p_Iitems => p_Iitems(istart:iend)
      end if
    end if

    if (present(p_IitemTags)) then
      nullify(p_IitemTags) 
      
      if (rintGroupSet%h_IitemTags .ne. ST_NOHANDLE) then
        call storage_getbase_int (rintGroupSet%h_IitemTags,p_IitemTags)
        
        ! Get a sub-pointer to the actual items.
        p_IitemTags => p_IitemTags(istart:iend)
      end if
    end if
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_getGroupByIndex (rintGroupSet,igroup,istart,iend,p_Iitems,p_IitemTags)

!<description>
  ! Returns a pointer to all items plus the start/end indices of group
  ! igroup in the item list.
!</description>

!<input>
  ! Integer group set to be used.
  type(t_intGroupSet), intent(inout) :: rintGroupSet

  ! Number of the group to retrieve information for.
  integer, intent(in) :: igroup
!</input>

!<output>
  ! Starting index of group igroup
  integer, intent(out) :: istart
  
  ! Ending index of group igroup
  integer, intent(out) :: iend

  ! OPTIONAL: Pointer to all items.
  integer, dimension(:), pointer, optional :: p_Iitems
  
  ! OPTIONAL: Pointer to all item tags.
  integer, dimension(:), pointer, optional :: p_IitemTags
!</utput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IitemIndex
    
    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_getGroupByIndex")
      call sys_halt()
    end if

    if ((igroup .lt. 1) .or. (igroup .gt. rintGroupSet%ngroups)) then
      call output_line ("igroup out of bounds.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_getGroupByIndex")
      call sys_halt()
    end if

    ! Get the index.
    call storage_getbase_int (rintGroupSet%h_IitemIndex,p_IitemIndex)
    
    if (.not. rintGroupSet%bdoubleIndex) then
    
      ! ---------------------
      ! Simple index
      ! ---------------------
      
      istart = p_IitemIndex(igroup)
      iend = p_IitemIndex(igroup+1)-1
    
    else
    
      ! ---------------------
      ! Double index.
      ! ---------------------
    
      istart = p_IitemIndex(2*igroup-1)
      iend = p_IitemIndex(2*igroup)

    end if
    
    ! Get the items
    if (present(p_Iitems)) then
      nullify(p_Iitems)
      if (rintGroupSet%h_Iitems .ne. ST_NOHANDLE) then
        call storage_getbase_int (rintGroupSet%h_Iitems,p_Iitems)
      end if
    end if

    if (present(p_IitemTags)) then
      nullify(p_IitemTags) 
      if (rintGroupSet%h_IitemTags .ne. ST_NOHANDLE) then
        call storage_getbase_int (rintGroupSet%h_IitemTags,p_IitemTags)
      end if
    end if

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_getGroupTags (rintGroupSet,p_IgroupTags)

!<description>
  ! Returns a pointer to the list of group tags.
!</description>

!<input>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</input>

!<inputoutput>
  ! Pointer to a list of group tags
  integer, dimension(:), pointer :: p_IgroupTags
!</inputoutput>

!</subroutine>

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_getGroupTags")
      call sys_halt()
    end if

    ! Get the list
    nullify(p_IgroupTags)

    if (rintGroupSet%h_IgroupTags .ne. ST_NOHANDLE) then
      call storage_getbase_int (rintGroupSet%h_IgroupTags,p_IgroupTags)
    end if
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_getItemGroupIDs (rintGroupSet,p_IitemGroupIDs)

!<description>
  ! Returns a pointer to a list with group IDs for each item
!</description>

!<input>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</input>

!<inputoutput>
  ! Pointer to a list of group IDs for each item
  integer, dimension(:), pointer :: p_IitemGroupIDs
!</inputoutput>

!</subroutine>

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_getItemGroupIDs")
      call sys_halt()
    end if

    ! Get the list
    nullify(p_IitemGroupIDs)

    if (rintGroupSet%h_IitemGroupIDs .ne. ST_NOHANDLE) then
      call storage_getbase_int (rintGroupSet%h_IitemGroupIDs,p_IitemGroupIDs)
    end if
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_sortItemsInGroups (rintGroupSet)

!<description>
  ! Sorts the items in the groups for increasing order.
!</description>

!<inputoutput>
  ! Integer group set to be modified.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Iitems,p_IitemTags
    integer :: igroup

    ! Some basic checks
    if (rintGroupSet%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_sortItemsInGroups")
      call sys_halt()
    end if

    if (rintGroupSet%h_IitemTags .eq. ST_NOHANDLE) then
    
      ! Loop through all groups
      do igroup = 1,rintGroupSet%ngroups
        ! Get the items
        call igs_getGroup (rintGroupSet,igroup,p_Iitems)
      
        ! Sort the items
        call sort_int (p_Iitems,SORT_QUICK)
      end do
    
    else

      ! Loop through all groups
      do igroup = 1,rintGroupSet%ngroups
        ! Get the items
        call igs_getGroup (rintGroupSet,igroup,p_Iitems,p_IitemTags)
      
        ! Sort the items.
        ! Simultaneously, provide the item tags as "mapping".
        ! This is sorted in the same order as the items.
        call sort_int (p_Iitems,SORT_QUICK,p_IitemTags)
      end do
      
    end if
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_printGroupSet (rintGroupSet)

!<description>
  ! Prints the groups to the terminal
!</description>

!<inputoutput>
  ! Integer group set to be printed.
  type(t_intGroupSet), intent(inout) :: rintGroupSet
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: igroup, iitem
    integer, dimension(:), pointer :: p_Iitems, p_IitemGroupIDs, p_IitemTags

    ! Loop through all groups
    do igroup = 1,rintGroupSet%ngroups
      ! Get the items
      call igs_getGroup (rintGroupSet,igroup,p_Iitems,p_IitemTags)
    
      ! Print the items
      call output_line ("Group "//trim(sys_siL(igroup,10))//":",bnolinebreak=.true.)
      
      do iitem = 1,ubound(p_Iitems,1)
        call output_line (" "//trim(sys_siL(p_Iitems(iitem),10)),bnolinebreak=.true.)
      end do
      
      call output_lbrk()
      
      if (associated(p_IitemTags)) then
        ! Print the item tags
        call output_line ("Tags :",bnolinebreak=.true.)
        
        do iitem = 1,ubound(p_Iitems,1)
          call output_line (" "//trim(sys_siL(p_IitemTags(iitem),10)),bnolinebreak=.true.)
        end do
        
        call output_lbrk()
      end if
       
    end do
    
    call igs_getItemGroupIDs (rintGroupSet,p_IitemGroupIDs)
    if (associated(p_IitemGroupIDs)) then
      ! Print the items
      call output_line ("Group of each item:",bnolinebreak=.true.)
      
      do iitem = 1,ubound(p_IitemGroupIDs,1)
        call output_line (" "//trim(sys_siL(p_IitemGroupIDs(iitem),10)),bnolinebreak=.true.)
      end do
      
      call output_lbrk()
    end if

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine igs_copyGroupSet (rgroupSource,rgroupDest,&
      bshareGroupData,bshareGroupTags,bshareItemTags)

!<description>
  ! Copies the integer group set to another structure
!</description>

!<input>
  ! Source structure
  type(t_intGroupSet), intent(in) :: rgroupSource
  
  ! OPTIONAL: If set to TRUE, the group data is shared, no memory is allocated.
  ! Default is FALSE.
  logical, intent(in), optional :: bshareGroupData

  ! OPTIONAL: If set to TRUE, the group tags are shared, no memory is allocated.
  ! Default is FALSE.
  logical, intent(in), optional :: bshareGroupTags

  ! OPTIONAL: If set to TRUE, the item tags are shared, no memory is allocated.
  ! Default is FALSE.
  logical, intent(in), optional :: bshareItemTags
!</input>

!<output>
  ! Destination structure
  type(t_intGroupSet), intent(out) :: rgroupDest
!</output>

!</subroutine>

    ! Some basic checks
    if (rgroupSource%h_IitemIndex .eq. ST_NOHANDLE) then
      call output_line ("rintGroupSet not initialised.",&
          OU_CLASS_ERROR,OU_MODE_STD,"igs_sortItemsInGroups")
      call sys_halt()
    end if

    ! First, just copy the structure
    rgroupDest = rgroupSource
    
    ! What about sharing the group data?
    if (present(bshareGroupData)) then
      rgroupDest%bsharedGroupData = bshareGroupData
    else
      rgroupDest%bsharedGroupData = .false.
    end if
    
    ! Not shared? Copy the data
    if (.not. rgroupDest%bsharedGroupData) then
      rgroupDest%h_IitemIndex = ST_NOHANDLE
      call storage_copy (rgroupSource%h_IitemIndex,rgroupDest%h_IitemIndex)
      
      if (rgroupSource%h_Iitems .ne. ST_NOHANDLE) then
        rgroupDest%h_Iitems = ST_NOHANDLE
        call storage_copy (rgroupSource%h_Iitems,rgroupDest%h_Iitems)
      end if
    end if

    ! What about sharing the group tags
    if (present(bshareGroupTags)) then
      rgroupDest%bsharedGroupTags = bshareGroupTags
    else
      rgroupDest%bsharedGroupTags = .false.
    end if

    ! Not shared? Copy the data
    if (.not. rgroupDest%bsharedGroupTags) then
      if (rgroupSource%h_IgroupTags .ne. ST_NOHANDLE) then
        rgroupDest%h_IgroupTags = ST_NOHANDLE
        call storage_copy (rgroupSource%h_IgroupTags,rgroupDest%h_IgroupTags)
      end if
    end if

    ! What about sharing the item tags
    if (present(bshareItemTags)) then
      rgroupDest%bsharedItemTags = bshareItemTags
    else
      rgroupDest%bsharedItemTags = .false.
    end if

    ! Not shared? Copy the data
    if (.not. rgroupDest%bsharedItemTags) then
      if (rgroupSource%h_IitemTags .ne. ST_NOHANDLE) then
        rgroupDest%h_IitemTags = ST_NOHANDLE
        call storage_copy (rgroupSource%h_IitemTags,rgroupDest%h_IitemTags)
      end if
    end if

  end subroutine

end module
