!##############################################################################
!# ****************************************************************************
!# <name> spacetimehierarchy </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module maintains realises a hierarchy of space-time meshes.
!#
!# A space-time level consists of a combination of a time-scale with a
!# spatial FEM space. The time-scale defines the discretisation in time,
!# the FEM-space thge discretisation in space.
!#
!# The following routines can be found here:
!#
!# 1.) sth_initHierarchy
!#     -> Initialises a space-time hierarchy
!#
!# 2.) sth_doneHierarchy
!#     -> Releases a space-time hierarchy
!#
!# 3.) sth_defineHierarchyByCoarsening
!#     -> Defines the refinement of a space-time hierarchy.
!#
!# 4.) sth_getLevel
!#     -> Returns pointers to the space and time level of a space-time level.
!# </purpose>
!##############################################################################

module spacetimehierarchy

  use fsystem
  use genoutput
  use boundary
  use basicgeometry
  use triangulation
  use collection
  use dofmapping

  use spatialdiscretisation
  use fespacehierarchy
  
  use timediscretisation
  use timescalehierarchy
  
  !use spacetimediscretisation
  
  implicit none
  
  private
  
  public :: t_spacetimeHierarchy
  public :: sth_initHierarchy
  public :: sth_doneHierarchy
  public :: sth_defineHierarchyByCoarsening
  public :: sth_getLevel
  public :: sth_printHierStatistics
  
!<types>

!<typeblock>

  ! A FE hierarchy that describes a hierarchy of FE spaces.
  type t_spacetimeHierarchy
  
    ! Number of levels available in this structure.
    integer :: nlevels = 0
    
    ! Pointer to a FE space hierarchy that defines space levels.
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    
    ! Pointer to a time scale hierarchy that defines time levels.
    type(t_timescaleHierarchy), pointer :: p_rtimeHierarchy
    
    ! Level definition of the space-time mesh.
    ! For each space-time level i, IlevelDef defines the corresponding
    ! space and time mesh:
    ! IlevelDef(1,i) = Space-level in rfeHierarchy.
    ! IlevelDef(2,i) = Time-level in rtimeHierarchy.
    integer, dimension(:,:), pointer :: p_IlevelDef
    
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sth_initHierarchy (rhierarchy,rfeHierarchy,rtimeHierarchy)

!<description>
  ! Basic initialisation of a space-time hierarchy.
  ! References to rfeHierarchy and rtimeHierarchy are saved to rhierarchy.
!</description>
 
!<input>
  ! Underlying space hierarchy.
  type(t_feHierarchy), intent(in), target :: rfeHierarchy

  ! Underlying time scale hierarchy.
  type(t_timescaleHierarchy), intent(in), target :: rtimeHierarchy
!</input>

!<output>
  ! A space-time hierarchy structure.
  type(t_spacetimeHierarchy), intent(out) :: rhierarchy
!</output>
  
!</subroutine>

    ! Remember the underlying hierarchies.
    rhierarchy%p_rfeHierarchy => rfeHierarchy
    rhierarchy%p_rtimeHierarchy => rtimeHierarchy
    
    ! Initialisation of the other variables.
    nullify(rhierarchy%p_IlevelDef)
    rhierarchy%nlevels = 0
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sth_doneHierarchy (rhierarchy)

!<description>
  ! Cleans up a space-time hierarchy.
!</description>
 
!<inputoutput>
  ! A space-time hierarchy structure to be cleaned up.
  type(t_spacetimeHierarchy), intent(inout) :: rhierarchy
!</inputoutput>
  
!</subroutine>

    ! Remove the pointers.
    nullify(rhierarchy%p_rfeHierarchy)
    nullify(rhierarchy%p_rtimeHierarchy)
    
    ! Release memory if necessary
    if (associated(rhierarchy%p_IlevelDef)) then
      deallocate(rhierarchy%p_IlevelDef)
    end if
    rhierarchy%nlevels = 0
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sth_defineHierarchyByCoarsening (rhierarchy,&
      nminSpace,nmaxSpace,nminTime,nmaxTime,dspacetimeRefFactor)

!<description>
  ! Defines the refinement of a space-time hierarchy.
  ! nminSpace/nmaxSpace define the lowest/highest space level to usem
  ! nmintime/nmaxTime the lowest/highest time-level.
  ! nmaxSpace/nmaxTime defines the finest space-time mesh.
  ! Lower levels are created by coarsening.
!</description>
 
!<input>
  ! Minimum allowed space level
  integer, intent(in) :: nminSpace
  
  ! Maximum allowed space level
  integer, intent(in) :: nmaxSpace
  
  ! Minimum allowed time level
  integer, intent(in) :: nminTime
  
  ! Maximum allowed time level
  integer, intent(in) :: nmaxTime

  ! OPTIONAL: Factor that specifies the number of refinement levels in 
  ! time before refining in space.
  ! Values e.g.:
  !   SYS_INFINITY = refinement only in time, always use maximum space level.
  !   2.0 = two refinements in time per refinement in space
  !   1.0 = each refinement in time results in a refinement in space. (Default)
  !   0.5 = two refinements in space per refinement in time
  !   0.0 = refinement only in space, always use maximum time level.
  real(dp), intent(in), optional :: dspacetimeRefFactor
!</input>

!<inputoutput>
  ! A space-time hierarchy structure.
  type(t_spacetimeHierarchy), intent(inout) :: rhierarchy
!</inputoutput>
  
!</subroutine>

    real(DP) :: drefFactor, dcurrentfactor
    integer :: ntime1,ntime2,nspace1,nspace2
    integer :: icurrent, i, j, nmax
    
    ! Cut off the levels if necessary.
    ntime1 = max(nminTime,1)
    nspace1 = max(nminSpace,1)
    
    ntime2 = min(nmaxTime,rhierarchy%p_rtimeHierarchy%nlevels)
    nspace2 = min(nmaxSpace,rhierarchy%p_rfeHierarchy%nlevels)
    
    ! Determine the refinement factor
    drefFactor = 1.0_DP
    if (present(dspacetimeRefFactor)) drefFactor = dspacetimeRefFactor

    ! Which type of refinement do we have?
    if (drefFactor .eq. SYS_INFINITY) then
      ! Refinement only in time.
      rhierarchy%nlevels = ntime2-ntime1+1
      allocate(rhierarchy%p_IlevelDef(2,rhierarchy%nlevels))
      rhierarchy%p_IlevelDef(1,:) = nspace2
      rhierarchy%p_IlevelDef(2,:) = (/(j,j=1,rhierarchy%nlevels)/)
      
    else if (drefFactor .eq. 0.0_DP) then
      ! Refinement only in space
      rhierarchy%nlevels = nspace2-nspace1+1
      allocate(rhierarchy%p_IlevelDef(2,rhierarchy%nlevels))
      rhierarchy%p_IlevelDef(1,:) = (/(j,j=1,rhierarchy%nlevels)/)
      rhierarchy%p_IlevelDef(2,:) = ntime2

    else if (drefFactor .ge. 1.0_DP) then
      ! Simultaneous refinement in space and time.
      ! Time refinement is the reference.
      !
      ! At first, determine how many space coarsenings we have left.
      nmax = int(real(nspace2-nspace1+1,dp) * drefFactor + 0.5_DP)

      ! Initialise the level counter depending on how much levels we have
      !rhierarchy%nlevels = min(ntime2-ntime1+1,nmax)
      rhierarchy%nlevels = ntime2-ntime1+1
      allocate(rhierarchy%p_IlevelDef(2,rhierarchy%nlevels))
      rhierarchy%p_IlevelDef(2,:) = (/(j,j=1,rhierarchy%nlevels)/)

      ! Initialise the space levels
      icurrent = nspace2
      dcurrentfactor = 0.0_DP
      do i=rhierarchy%nlevels, 1,-1
        rhierarchy%p_IlevelDef(1,i) = icurrent
        rhierarchy%p_IlevelDef(2,i) = ntime2-(rhierarchy%nlevels-i)

        ! If the current factor overflows, coarsen in space.   
        if (present(dspacetimeRefFactor)) then
          dcurrentfactor = dcurrentfactor + 1.0_DP/dspacetimeRefFactor
        else
          dcurrentfactor = dcurrentfactor + 1.0_DP
        end if
        if (dcurrentfactor .ge. 1.0_DP) then
          ! Coarsening in space if possible
          if (icurrent .gt. nminSpace) &
            icurrent = icurrent - 1
          dcurrentfactor = mod(dcurrentfactor,1.0_DP)
        end if
      end do
      
    else if (drefFactor .gt. 0.0_DP) then

      ! Simultaneous refinement in space and time.
      ! Space refinement is the reference.
      !
      ! At first, determine how many time coarsenings we have left.
      nmax = int(real(ntime2-ntime1+1,dp) / drefFactor + 0.5_DP)

      ! Initialise the level counter depending on how much levels we have
      rhierarchy%nlevels = min(nspace2-nspace1+1,nmax)
      allocate(rhierarchy%p_IlevelDef(2,rhierarchy%nlevels))

      ! Initialise the time levels
      icurrent = nspace2
      dcurrentfactor = 0.0_DP
      do i=rhierarchy%nlevels, 1,-1
        rhierarchy%p_IlevelDef(2,i) = nspace2-(rhierarchy%nlevels-i)
        rhierarchy%p_IlevelDef(2,i) = icurrent

        ! If the current factor overflows, coarsen in time.
        if (present(dspacetimeRefFactor)) then
          dcurrentfactor = dcurrentfactor + dspacetimeRefFactor
        else
          dcurrentfactor = dcurrentfactor + 1.0_DP
        end if
        if (dcurrentfactor .ge. 1.0_DP) then
          ! Coarsening in space
          icurrent = icurrent - 1
          dcurrentfactor = mod(dcurrentfactor,1.0_DP)
        end if
      end do
    
    else
      call output_line ('Invalid refinement!.', &
          OU_CLASS_ERROR,OU_MODE_STD,'sth_defineHierarchyByCoarsening')
      call sys_halt()
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sth_getLevel (rhierarchy,ilevel,p_rfeSpaceLevel,p_rtimeDiscr,&
      ispaceLevel,itimeLevel)

!<description>
  ! Returns pointer to the space and time discretisation on a specific
  ! space-time level.
!</description>
 
!<input>
  ! A space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in), target :: rhierarchy
  
  ! Space-time level.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! OPTIONAL: Returns a pointer to the corresponding space-level
  type(t_feSpaceLevel), pointer, optional :: p_rfeSpaceLevel
  
  ! OPTIONAL: Returns a pointer to the corresponding time level
  type(t_timeDiscretisation), pointer, optional :: p_rtimeDiscr
  
  ! Optional: Returns the space-level. Level 1 identifies the coarse mesh.
  integer, intent(out), optional :: ispaceLevel

  ! Optional: Returns the time-level. Level 1 identifies the coarse mesh.
  integer, intent(out), optional :: itimeLevel
  
!</inputoutout>
  
!</subroutine>

    if ((ilevel .lt. 1) .or. (ilevel .gt. rhierarchy%nlevels)) then
      call output_line ('Invalid level!.', &
          OU_CLASS_ERROR,OU_MODE_STD,'sth_getLevel')
      call sys_halt()
    end if
    
    ! Get the pointers.
    if (present(p_rfeSpaceLevel)) then
      p_rfeSpaceLevel => rhierarchy%p_rfeHierarchy%p_rfeSpaces( &
          rhierarchy%p_IlevelDef(1,ilevel) )
    end if

    if (present(p_rtimeDiscr)) then
      p_rtimeDiscr => rhierarchy%p_rtimeHierarchy%p_rtimeLevels( &
          rhierarchy%p_IlevelDef(2,ilevel) )
    end if
    
    if (present(ispaceLevel)) then
      ispaceLevel = rhierarchy%p_IlevelDef(1,ilevel)
    end if

    if (present(itimeLevel)) then
      itimeLevel = rhierarchy%p_IlevelDef(2,ilevel)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sth_printHierStatistics (rhierarchy)

!<description>
  ! Writes statistics about the mesh hierarchy to the terminal.
!</description>

!<inputoutput>
  ! A space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rhierarchy
!</inputoutput>
  
!</subroutine>
  
    integer :: i,ispacelevel,itimelevel
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
  
    do i=1,rhierarchy%nlevels
      if (i .eq. 1) then
        ! Print a headline
        call output_line("Lv. S-Lv. T-Lv. #dof(space)  #dof(time)")
        call output_line("---------------------------------------")
      end if
      
      ! Get information about the level
      call sth_getLevel (rhierarchy,i,p_rfeSpaceLevel,p_rtimeDiscr,&
          ispaceLevel,itimeLevel)
      
      ! Print statistics about that level
      call output_line ( &
          trim(sys_si(i,3))&
        //trim(sys_si(ispaceLevel,6))&
        //trim(sys_si(itimeLevel,6))&
        //trim(sys_si(dof_igetNDofGlobBlock(p_rfeSpaceLevel%p_rdiscretisation),12))&
        //trim(sys_si(tdiscr_igetNDofGlob(p_rtimeDiscr),12)) )
    end do
    
  end subroutine

end module
