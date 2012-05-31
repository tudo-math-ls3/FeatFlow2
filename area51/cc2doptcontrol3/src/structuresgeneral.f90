!##############################################################################
!# ****************************************************************************
!# <name> structuresgeneral </name>
!# ****************************************************************************
!#
!# <purpose>
!# General basic structures.
!#
!# Routines in this module:
!#
!# 1.) struc_getDebugFlags
!#     -> Reads debug flags
!#
!# 2.) struc_getRefinementParams
!#     -> Reads parameters about the space-time refinement.
!# </purpose>
!##############################################################################

module structuresgeneral

  use fsystem
  use storage
  use paramlist
  
  use analyticsolution
  
  use timediscretisation
  use structuresdiscretisation
  use structuresoptcontrol
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! Type block that specifies settings about the refinement of a mesh.
  type t_settings_refinement
  
    ! Type of refinement.
    ! =0: 2-level refinement
    integer :: crefType = 0

    ! Number of pre-refinements with 2-level refinement to calculate the
    ! coarse mesh in rmeshHierarchy from rtriaCoarse.
    integer :: npreref = 0
    
    ! Total number of levels.
    integer :: nlevels = 0
    
    ! Output level that defines which information to print
    ! during the refinement.
    ! =-1: no information.
    ! = 0: only errors/warnings.
    ! = 1: basic information.
    ! = 2: standard information.
    integer :: coutputlevel = 2
    
  end type

!</typeblock>

  public :: t_settings_refinement

!<typeblock>

  ! This type encapsules all debug flags that may be of use somewhere
  ! in the program.
  type t_optcDebugFlags
  
    ! If the following constant is set from 1.0 to 0.0, the primal system is
    ! decoupled from the dual system!
    real(dp) :: dprimalDualCoupling = 1.0_DP

    ! If the following constant is set from 1.0 to 0.0, the dual system is
    ! decoupled from the primal system!
    real(dp) :: ddualPrimalCoupling = 1.0_DP

    ! If the following parameter is set from 1.0 to 0.0, the terminal
    ! condition between the primal and dual equation is decoupled, i.e.
    ! the dual equation gets independent from the primal one.
    real(dp) :: dterminalCondDecoupled = 1.0_DP

    ! If the following parameter is set from 1.0 to 0.0, the time coupling
    ! is disabled, resulting in a stationary simulation in every timestep.
    real(dp) :: dtimeCoupling = 1.0_DP
    
    ! Modification to the discrete RHS.
    ! =0: No modification (standard).
    ! =1: Disturb the primal velocity RHS in all DOF's with a random value.
    ! =2: Disturb the primal velocity RHS in all DOF's except for the boundary DOF's
    !     with a random value.
    ! =3: Disturb the dual velocity RHS in all DOF's with a random value.
    ! =4: Disturb the dual velocity RHS in all DOF's except for the boundary DOF's
    !     with a random value.
    ! =5: Disturb the velocity RHS in all DOF's with a random value.
    ! =6: Disturb the velocity RHS in all DOF's except for the boundary DOF's
    !     with a random value.
    integer :: crhsmodification = 0

    ! Maximum error to be introduced to the RHS if crhsmodification=1/2.
    real(DP) :: drhsrandomMax = 1E-13

  end type

!</typeblock>

  public :: t_optcDebugFlags

!<typeblock>

  ! A collection of global data which is passen to the callback routines
  ! in user_callback.f90
  type t_globalData
    ! An application specific parameter list.
    type(t_parlist), pointer :: p_rparlist => null()
    
    ! Pointer to the coarse time discretisation.
    type(t_timeDiscretisation), pointer :: p_rtimeCoarse => null()
    
    ! Reference to the physics parameters
    type(t_settings_physics), pointer :: p_rphysics => null()

    ! Reference to the optimal control parameters.
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl => null()
    
    ! Pointer to the primal right hand side.
    type(t_anSolution), pointer :: p_rrhsPrimal => null()

    ! Pointer to the dual right hand side.
    type(t_anSolution), pointer :: p_rrhsDual => null()

    ! Pointer to the target function.
    type(t_anSolution), pointer :: p_rtargetFunction => null()
  end type

!</typeblock>

  public :: t_globalData

!</types>

  ! Reads all debug flags from the parameter list.
  public :: struc_getDebugFlags
  
  ! Reads in parameters that define the refinement.
  public :: struc_getRefinementParams

contains

  ! ***************************************************************************

!<subroutine>

  subroutine struc_getDebugFlags (rparlist,ssection,rdebugFlags)
  
!<description>
  ! Reads all debug flags from the parameter list.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the parameters.
  character(len=*), intent(in) :: ssection
!</input>

!<output>
  ! Structure with the parameters of the nonlinear soace-time solver.
  type(t_optcDebugFlags), intent(out) :: rdebugFlags
!</output>

!</subroutine>

    call parlst_getvalue_double (rparlist,ssection,&
        'dprimalDualCoupling',rdebugFlags%dprimalDualCoupling,1.0_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        'ddualPrimalCoupling',rdebugFlags%ddualPrimalCoupling,1.0_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        'dterminalCondDecoupled',rdebugFlags%dterminalCondDecoupled,1.0_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        'dtimeCoupling',rdebugFlags%dtimeCoupling,1.0_DP)

    call parlst_getvalue_int (rparlist,ssection,&
        'crhsmodification',rdebugFlags%crhsmodification,0)

    call parlst_getvalue_double (rparlist,ssection,&
        'drhsrandomMax',rdebugFlags%drhsrandomMax,0.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine struc_getRefinementParams (rparlist,rrefinementSpace,rrefinementTime,&
      ssectionSpace,ssectionTime)
  
!<description>
  ! Reads in parameters that define the refinement.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the spatial mesh/domain can be found.
  character(len=*), intent(in) :: ssectionSpace

  ! Section where the parameters of the time mesh can be found.
  character(len=*), intent(in) :: ssectionTime
!</input>

!<output>
  ! Description of the refinement in space
  type(t_settings_refinement), intent(out) :: rrefinementSpace

  ! Description of the refinement in time
  type(t_settings_refinement), intent(out) :: rrefinementTime
!</output>

!</subroutine>

    ! local variables
    integer :: nlmin,nlmax

    ! Get the level numbers.
    call parlst_getvalue_int (rparlist,ssectionSpace,'NLMIN',nlmin,1)
    call parlst_getvalue_int (rparlist,ssectionSpace,'NLMAX',nlmax,2)

    ! Do some correction to the level numbers
    if (nlmin .le. 0) nlmin = nlmax+nlmin
    nlmin = min(nlmin,nlmax)
    nlmax = max(nlmin,nlmax)
    
    ! Calculate the total number of levels
    rrefinementSpace%crefType = 0
    rrefinementSpace%nlevels = nlmax-nlmin+1
    rrefinementSpace%npreref = nlmin-1

    ! The same for the time mesh.
    call parlst_getvalue_int (rparlist,ssectionTime,'TIMENLMIN',nlmin,1)
    call parlst_getvalue_int (rparlist,ssectionTime,'TIMENLMAX',nlmax,2)

    ! Do some correction to the level numbers
    if (nlmin .le. 0) nlmin = nlmax+nlmin
    nlmin = min(nlmin,nlmax)
    nlmax = max(nlmin,nlmax)
    
    ! Calculate the total number of levels
    rrefinementTime%crefType = 0
    rrefinementTime%nlevels = nlmax-nlmin+1
    rrefinementTime%npreref = nlmin-1

  end subroutine

end module
