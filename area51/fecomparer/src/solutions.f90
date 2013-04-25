!####################################################################
!# Comparison of FE functions.
!####################################################################

module solutions

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  
  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  
  use linearsystemscalar
  use linearsystemblock
  
  use paramlist
  use vectorio

  implicit none
  
  ! Encapsules a solution
  type t_solution
      
    ! Whether or not the solution is nonstationary.
    ! =0: Stationary solution

    integer :: cnonstationary = 0

    ! Underlying FE space and solution type

    integer :: cfespace = 0

    ! Filename of the file containing the solution

    character(LEN=SYS_STRLEN) :: sfilename = ""

    ! Refinement level w.r.t. the global mesh

    integer :: ireflevel = 0

    ! Component to analyse.
    ! =0: Analyse all components.

    integer :: icomponent = 0
    
    ! Underlying FE space hierarchy
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    
    ! Array of vectors for all levels
    type(t_vectorBlock), dimension(:), pointer :: p_Rvectors
    
  end type
  
  public :: t_solution
  
  public :: sol_readparams
  public :: sol_initSolution
  public :: sol_doneSolution
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine sol_readparams (rparlist,iid,rsolution)

!<description>
  ! Reads parameters for a solution.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! ID of the solution
  integer :: iid
!</input>

!<output>
  ! Solution structure
  type(t_solution), intent(out) :: rsolution
!</output>

!</subroutine>

    ! local variables
    character(LEN=SYS_STRLEN) :: ssection,sparent
    
    ! Get the section
    ssection = "SOLUTION"//trim(sys_siL(iid,10))
    
    do
      call parlst_getvalue_string (rparlist, ssection, &
          "sparent", sparent, "", bdequote=.true.)
      if (sparent .eq. "") exit
      
      ! Go to the parent.
      ssection = sparent
    end do   

    ! Get parameters
    call parlst_getvalue_int (rparlist, ssection, &
        "cnonstationary", rsolution%cnonstationary,rsolution%cnonstationary)

    call parlst_getvalue_int (rparlist, ssection, &
        "cfespace", rsolution%cfespace,rsolution%cfespace)
        
    call parlst_getvalue_int (rparlist, ssection, &
        "ireflevel", rsolution%ireflevel,rsolution%ireflevel)
        
    call parlst_getvalue_string (rparlist, ssection, &
        "sfilename", rsolution%sfilename,rsolution%sfilename, bdequote=.true.)
    
    call parlst_getvalue_int (rparlist, ssection, &
        "icomponent", rsolution%icomponent,rsolution%icomponent)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sol_initSolution (rsolution,rfeHierarchy)

!<description>
  ! Initialises memory.
!</description>

!<input>
  ! Underlying FE space hierarchy.
  type(t_feHierarchy), target :: rfeHierarchy
!</input>

!<inputoutput>
  ! Solution structure
  type(t_solution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    integer :: i
    
    ! Save the hierarchy.
    rsolution%p_rfeHierarchy => rfeHierarchy
    
    ! Create solution vectors.
    allocate (rsolution%p_Rvectors (rfeHierarchy%nlevels))
    
    ! Create a vector on each level
    do i=1,rfeHierarchy%nlevels
      call lsysbl_createVector (rfeHierarchy%p_RfeSpaces(i)%p_rdiscretisation,&
          rsolution%p_Rvectors(i),.true.)
    end do    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sol_doneSolution (rsolution)

!<description>
  ! Initialises memory.
!</description>

!<inputoutput>
  ! Solution structure
  type(t_solution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    integer :: i
    
    ! Release everything.
    do i=1,size(rsolution%p_Rvectors)
      call lsysbl_releaseVector (rsolution%p_Rvectors(i))
    end do
    
    deallocate (rsolution%p_Rvectors)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sol_readSolutionData (rsolution)

!<description>
  ! Reads the data of a solution from a file.
!</description>

!<inputoutput>
  ! Solution structure
  type(t_solution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    character(LEN=SYS_STRLEN) :: sarray

    ! If a filename is given, read the file.
    if (rsolution%sfilename .ne. "") then
      call output_line ("Reading solution.")
      call vecio_readBlockVectorHR (rsolution%p_Rvectors(rsolution%ireflevel), &
          sarray, .true.,0, rsolution%sfilename, .true.)
    else
      call lsysbl_clearVector (rsolution%p_Rvectors(rsolution%ireflevel))
    end if
    
  end subroutine

end module
