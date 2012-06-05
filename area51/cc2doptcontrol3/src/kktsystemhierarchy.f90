!##############################################################################
!# ****************************************************************************
!# <name> kktsystemhierarchy </name>
!# ****************************************************************************
!#
!# <purpose>
!# Realises a hierarchy of solutions of a KKT system.
!# </purpose>
!##############################################################################

module kktsystemhierarchy

  use fsystem
  use genoutput
  use spatialdiscretisation
  use timediscretisation

  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral

  use fespacehierarchybase  
  use spacetimehierarchy

  use kktsystemspaces
  use kktsystem
  
  use structuresoperatorasm
  use spacematvecassembly

  implicit none
  
  private

!<types>

!<typeblock>
  
  ! A hierarchy of solutions of a KKT system
  type t_kktsystemHierarchy
  
    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierPrimal => null()
    
    ! A space-time hierarchy based on the dual space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierDual => null()

    ! A space-time hierarchy based on the control space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierControl => null()

    ! Number of levels in the hierarchy
    integer :: nlevels = 0

    ! Array of KKT system (solutions) for all the levels
    type(t_kktsystem), dimension(:), pointer :: p_RkktSystems => null()
  
  end type

!</typeblock>

  public :: t_kktsystemHierarchy

!<typeblock>
  
  ! A hierarchy of directional derivatives of a KKT system
  type t_kktsystemDirDerivHierarchy
  
    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierPrimal => null()
    
    ! A space-time hierarchy based on the dual space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierDual => null()

    ! A space-time hierarchy based on the control space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierControl => null()

    ! Number of levels in the hierarchy
    integer :: nlevels = 0

    ! Array of KKT system (solutions) for all the levels
    type(t_kktsystemDirDeriv), dimension(:), pointer :: p_RkktSysDirDeriv => null()
  
  end type

!</typeblock>

  public :: t_kktsystemDirDerivHierarchy

!</types>

  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
  public :: kkth_initPrimalSol
  
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
  public :: kkth_initDualSol

  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
  public :: kkth_initControl
  
  ! Initialises a KKT system structure for level ilevel of the hierarchy
  public :: kkth_initKKTSystem
  
  ! Initialises a hierarchy of KKT systems.
  public :: kkth_initHierarchy

  ! Cleans up a hierarchy of KKT systems.
  public :: kkth_doneHierarchy
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initPrimalSol (rvector,rkktsystemHierarchy,ilevel)
  
!<description>
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy
  
  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Vector to be created
  type(t_primalSpace), intent(out) :: rvector
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev
    type(t_feSpaceLevel), pointer :: p_rspaceDiscr
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierPrimal,ilev,&
        p_rspaceDiscr,p_rtimeDiscr)

    ! Create the vector
    call kktsp_initPrimalVector (rvector,&
        p_rspaceDiscr%p_rdiscretisation,p_rtimeDiscr)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initDualSol (rvector,rkktsystemHierarchy,ilevel)
  
!<description>
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy
  
  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Vector to be created
  type(t_dualSpace), intent(out) :: rvector
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev
    type(t_feSpaceLevel), pointer :: p_rspaceDiscr
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierDual,ilev,&
        p_rspaceDiscr,p_rtimeDiscr)

    ! Create the vector
    call kktsp_initDualVector (rvector,&
        p_rspaceDiscr%p_rdiscretisation,p_rtimeDiscr)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initControl (rvector,rkktsystemHierarchy,ilevel)
  
!<description>
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy
  
  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Vector to be created
  type(t_controlSpace), intent(out) :: rvector
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev
    type(t_feSpaceLevel), pointer :: p_rspaceDiscr
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierControl,ilev,&
        p_rspaceDiscr,p_rtimeDiscr)

    ! Create the vector
    call kktsp_initControlVector (rvector,&
        p_rspaceDiscr%p_rdiscretisation,p_rtimeDiscr)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initKKTSystem (rkktsystem,rkktsystemHierarchy,ilevel,roperatorAsmHier)
  
!<description>
  ! Initialises a KKT system structure for level ilevel of the hierarchy
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy

  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
  
  ! Hierarchy of all possible space-time operators
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
!</input>

!<output>
  ! KKT sytem to be created
  type(t_kktsystem), intent(out) :: rkktsystem
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev, ispacelevel, itimelevel
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierPrimal,ilev,&
        ispacelevel=ispacelevel, itimelevel=itimelevel)
    
    ! Initialise the KKT system on level i
    call kkt_initKKTsystem (rkktsystem,roperatorAsmHier,&
        ispacelevel, itimelevel)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initHierarchy (rkktsystemHierarchy,roperatorAsmHier,&
      rspaceTimeHierPrimal,rspaceTimeHierDual,rspaceTimeHierControl,ballocate)
  
!<description>
  ! Initialises a hierarchy of KKT systems.
!</description>

!<input>
  ! Hierarchy of space-time operator assembly structures.
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
  
  ! Space-time hierarchy in the primal space
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierPrimal

  ! Space-time hierarchy in the dual space
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierDual

  ! Space-time hierarchy in the control space
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierControl
  
  ! Whether or not to allocate memory for the solutions
  logical, intent(in) :: ballocate
!</input>

!<output>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(out) :: rkktsystemHierarchy
!</output>

!</subroutine>
    
    ! local variables
    integer :: i

    ! Remember the structures
    rkktsystemHierarchy%p_rspaceTimeHierPrimal => rspaceTimeHierPrimal
    rkktsystemHierarchy%p_rspaceTimeHierDual => rspaceTimeHierDual
    rkktsystemHierarchy%p_rspaceTimeHierControl => rspaceTimeHierControl
    
    ! Numer of levels in the hierarchy?
    rkktsystemHierarchy%nlevels = rspaceTimeHierPrimal%nlevels
    
    nullify(rkktsystemHierarchy%p_RkktSystems)
    
    if (ballocate) then
    
      ! Create the KKT system solutions
      allocate(rkktsystemHierarchy%p_RkktSystems(rkktsystemHierarchy%nlevels))
      do i=1,rkktsystemHierarchy%nlevels
        ! Initialise the KKT system on level i
        call kkth_initKKTSystem (&
            rkktsystemHierarchy%p_RkktSystems(i),rkktsystemHierarchy,i,&
            roperatorAsmHier)
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_doneHierarchy (rkktsystemHierarchy)
  
!<description>
  ! Cleans up a hierarchy of KKT systems.
!</description>

!<inputoutput>
  ! Hierarchy of KKT system solutions to be cleaned up
  type(t_kktsystemHierarchy), intent(inout) :: rkktsystemHierarchy
!</inputoutput>

!</subroutine>
    
    ! local variables
    integer :: i

    ! Remember the structures
    nullify(rkktsystemHierarchy%p_rspaceTimeHierPrimal)
    nullify(rkktsystemHierarchy%p_rspaceTimeHierDual)
    nullify(rkktsystemHierarchy%p_rspaceTimeHierControl)
    
    if (associated(rkktsystemHierarchy%p_RkktSystems)) then
      ! Clean up the solutions on every level
      do i=rkktsystemHierarchy%nlevels,1,-1
        call kkt_doneKKTsystem (rkktsystemHierarchy%p_RkktSystems(i))
      end do

      deallocate(rkktsystemHierarchy%p_RkktSystems)
    end if
    
    rkktsystemHierarchy%nlevels = 0

  end subroutine

end module
