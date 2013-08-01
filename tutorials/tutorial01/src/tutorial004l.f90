!##############################################################################
!# Tutorial 004l: Datastructures - integer groups
!##############################################################################

module tutorial004l

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use intgroupset

  implicit none
  private
  
  public :: start_tutorial004l

contains

  ! ***************************************************************************

  subroutine start_tutorial004l

    ! Declare some variables
    type(t_intGroupSet) :: rintGroupSet, rintGroupSetCopy
    integer, dimension(:), pointer :: p_Iitems,p_IitemCounter,p_IitemTags

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004l")
    call output_separator (OU_SEP_MINUS)
    
    ! ==================================================================
    ! Create an integer group set
    ! by directly adding the groups
    ! ==================================================================
    
    ! ---------------------------------
    ! Create a group set with 4 groups.
    ! ---------------------------------
    
    call igs_initGroupSet (rintGroupSet, 4)
    
    ! ---------------------------------
    ! Initialise the groups and the entries
    ! ---------------------------------

    ! Start appending groups. We have 10 itemes.
    call igs_initGroupAppend (rintGroupSet,10)
    
    ! Define the four groups
    call igs_appendGroup  (rintGroupSet,(/ 1, 5 /))
    call igs_appendGroup  (rintGroupSet,(/ 6, 2 /))
    call igs_appendGroup  (rintGroupSet,(/ 7, 3, 4 /))
    call igs_appendGroup  (rintGroupSet,(/ 10, 9, 8 /))
    
    ! Finish
    call igs_doneGroupAppend (rintGroupSet)

    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------

    call output_lbrk()
    call output_line ("Integer group, unsorted.")
    call igs_printGroupSet (rintGroupSet)
    
    ! ---------------------------------
    ! Calculate the group ID of each item
    ! ---------------------------------

    call igs_calcItemGroupIDs (rintGroupSet)

    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------

    call output_lbrk ()
    call output_line ("Integer group, direct creation, unsorted, with group ID.")
    call igs_printGroupSet (rintGroupSet)

    ! ---------------------------------
    ! Sort the groups
    ! ---------------------------------
    
    call igs_sortItemsInGroups (rintGroupSet)
    
    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------
    
    call output_lbrk ()
    call output_line ("Integer group, direct creation, sorted, with group ID.")
    call igs_printGroupSet (rintGroupSet)
    
    call output_lbrk ()
    call output_line ("Integer group, direct creation, sorted, with group ID.")
    call igs_printGroupSet (rintGroupSet)

    ! ---------------------------------
    ! Clean up.
    ! ---------------------------------
    
    call igs_doneGroupSet (rintGroupSet)

    ! ==================================================================
    ! Create an integer group set
    ! by directly specifying the number
    ! of elements in each group.
    ! ==================================================================
    
    ! ---------------------------------
    ! Create a group set with 4 groups.
    ! ---------------------------------
    
    call igs_initGroupSet (rintGroupSet, 4)
    
    ! ---------------------------------
    ! Initialise the groups and the entries
    ! ---------------------------------
    
    ! Two groups with 2 entries, two groups with 3 entries
    call igs_setItemCount (rintGroupSet,(/ 2, 2, 3, 3 /))
    
    ! Define the four groups
    call igs_getGroup (rintGroupSet,1,p_Iitems)
    p_Iitems = (/ 1, 5 /)

    call igs_getGroup (rintGroupSet,2,p_Iitems)
    p_Iitems = (/ 6, 2 /)
    
    call igs_getGroup (rintGroupSet,3,p_Iitems)
    p_Iitems = (/ 7, 3, 4 /)

    call igs_getGroup (rintGroupSet,4,p_Iitems)
    p_Iitems = (/ 10, 9, 8 /)
    
    ! ---------------------------------
    ! Print the content
    ! ---------------------------------
    
    call output_lbrk()
    call output_line ("Integer group, manual creation, unsorted.")
    call igs_printGroupSet (rintGroupSet)
    
    ! ---------------------------------
    ! Calculate the group ID of each item
    ! ---------------------------------
    
    call igs_calcItemGroupIDs (rintGroupSet)

    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------

    call output_lbrk ()
    call output_line ("Integer group, manual creation, unsorted, with group ID.")
    call igs_printGroupSet (rintGroupSet)
    
    ! ---------------------------------
    ! Sort the groups
    ! ---------------------------------

    call igs_sortItemsInGroups (rintGroupSet)
    
    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------

    call output_lbrk ()
    call output_line ("Integer group, manual creation, sorted, with group ID.")
    call igs_printGroupSet (rintGroupSet)
    
    ! ---------------------------------
    ! Clean up.
    ! ---------------------------------

    call igs_doneGroupSet (rintGroupSet)
    
    ! ==================================================================
    ! Create an integer group set
    ! by defining a counter for the
    ! number of items in each group
    ! ==================================================================
    
    ! ---------------------------------
    ! Create a group set with 4 groups.
    ! ---------------------------------
    
    call igs_initGroupSet (rintGroupSet, 4)
    
    ! ---------------------------------
    ! Initialise the groups and the entries
    ! ---------------------------------

    ! Start the item counter
    call igs_startItemCount (rintGroupSet,p_IitemCounter)
    
    ! Define the number of items / count the items
    p_IitemCounter(1) = p_IitemCounter(1) + 2
    p_IitemCounter(2) = p_IitemCounter(2) + 2
    p_IitemCounter(3) = p_IitemCounter(3) + 3
    p_IitemCounter(4) = p_IitemCounter(4) + 3
    
    ! Done. Initialise the structure appropriately
    call igs_finishItemCount (rintGroupSet,p_IitemCounter)
    
    ! Define the four groups
    call igs_getGroup (rintGroupSet,1,p_Iitems)
    p_Iitems = (/ 1, 5 /)

    call igs_getGroup (rintGroupSet,2,p_Iitems)
    p_Iitems = (/ 6, 2 /)
    
    call igs_getGroup (rintGroupSet,3,p_Iitems)
    p_Iitems = (/ 7, 3, 4 /)

    call igs_getGroup (rintGroupSet,4,p_Iitems)
    p_Iitems = (/ 10, 9, 8 /)
    
    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------

    call output_lbrk()
    call output_line ("Integer group, counter creation, unsorted.")
    call igs_printGroupSet (rintGroupSet)
    
    ! ---------------------------------
    ! Calculate the group ID of each item
    ! ---------------------------------

    call igs_calcItemGroupIDs (rintGroupSet)

    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------

    call output_lbrk ()
    call output_line ("Integer group, counter creation, unsorted, with group ID.")
    call igs_printGroupSet (rintGroupSet)

    ! ---------------------------------
    ! Sort the groups
    ! ---------------------------------
    
    call igs_sortItemsInGroups (rintGroupSet)
    
    ! ---------------------------------
    ! Print the groups
    ! ---------------------------------
    
    call output_lbrk ()
    call output_line ("Integer group, counter creation, sorted, with group ID.")
    call igs_printGroupSet (rintGroupSet)
    
    ! ---------------------------------
    ! Clean up.
    ! ---------------------------------
    
    call igs_doneGroupSet (rintGroupSet)

    ! ==================================================================
    ! Copy test for groups
    ! ==================================================================
    
    ! ---------------------------------
    ! Create a group set with 4 groups.
    ! ---------------------------------
    
    call igs_initGroupSet (rintGroupSet, 4)
    
    ! ---------------------------------
    ! Initialise the groups and the entries
    ! ---------------------------------

    ! Start appending groups. We have 10 itemes.
    call igs_initGroupAppend (rintGroupSet,10)
    
    ! Define the four groups
    call igs_appendGroup  (rintGroupSet,(/ 1, 5 /))
    call igs_appendGroup  (rintGroupSet,(/ 6, 2 /))
    call igs_appendGroup  (rintGroupSet,(/ 7, 3, 4 /))
    call igs_appendGroup  (rintGroupSet,(/ 10, 9, 8 /))
    
    ! Finish
    call igs_doneGroupAppend (rintGroupSet)

    ! ---------------------------------
    ! Copy the group, share data, print it
    ! ---------------------------------
    
    call igs_copyGroupSet (rintGroupSet,rintGroupSetCopy,.true.)

    call output_lbrk ()
    call output_line ("Integer group, copy, shared.")
    call igs_printGroupSet (rintGroupSetCopy)

    ! ---------------------------------
    ! Clean up copy.
    ! ---------------------------------
    
    call igs_doneGroupSet (rintGroupSetCopy)

    ! ---------------------------------
    ! Copy the group, do not share data, print it
    ! ---------------------------------
    
    call igs_copyGroupSet (rintGroupSet,rintGroupSetCopy,.false.)

    call output_lbrk ()
    call output_line ("Integer group, copy, not shared.")
    call igs_printGroupSet (rintGroupSetCopy)

    ! ---------------------------------
    ! Clean up copy.
    ! ---------------------------------
    
    call igs_doneGroupSet (rintGroupSetCopy)

    ! ---------------------------------
    ! Final cleanup
    ! ---------------------------------
    
    call igs_doneGroupSet (rintGroupSet)

    ! ==================================================================
    ! Associating tags
    ! ==================================================================
    
    ! ---------------------------------
    ! Create a group set with 4 groups.
    ! ---------------------------------
    
    call igs_initGroupSet (rintGroupSet, 4)
    
    ! ---------------------------------
    ! Initialise the groups and the entries
    ! ---------------------------------

    ! Start appending groups. We have 10 itemes.
    call igs_initGroupAppend (rintGroupSet,10)
    
    ! Define the four groups
    call igs_appendGroup  (rintGroupSet,(/ 1, 5 /))
    call igs_appendGroup  (rintGroupSet,(/ 6, 2 /))
    call igs_appendGroup  (rintGroupSet,(/ 7, 3, 4 /))
    call igs_appendGroup  (rintGroupSet,(/ 10, 9, 8 /))
    
    ! Finish
    call igs_doneGroupAppend (rintGroupSet)

    ! ---------------------------------
    ! Append tags
    ! ---------------------------------

    ! Append a user defined tag array
    call igs_allocItemTags (rintGroupSet)
    
    ! Define some tags -- let us say, the item number in the group
    call igs_getGroup (rintGroupSet,1,p_IitemTags=p_IitemTags)
    p_IitemTags = (/ 1, 2 /)

    call igs_getGroup (rintGroupSet,2,p_IitemTags=p_IitemTags)
    p_IitemTags = (/ 1, 2 /)
    
    call igs_getGroup (rintGroupSet,3,p_IitemTags=p_IitemTags)
    p_IitemTags = (/ 1, 2, 3 /)

    call igs_getGroup (rintGroupSet,4,p_IitemTags=p_IitemTags)
    p_IitemTags = (/ 1, 2, 3 /)
    
    ! ---------------------------------
    ! Print, sort, print again
    ! ---------------------------------
    call output_lbrk ()
    call output_line ("Integer group, with tags, unsorted.")
    call igs_printGroupSet (rintGroupSet)
    
    ! Sort
    call igs_sortItemsInGroups (rintGroupSet)

    ! Print again
    call output_lbrk ()
    call output_line ("Integer group, with tags, sorted.")
    call igs_printGroupSet (rintGroupSet)

    ! ---------------------------------
    ! Cleanup
    ! ---------------------------------
    
    call igs_doneGroupSet (rintGroupSet)
  end subroutine
  
end module
