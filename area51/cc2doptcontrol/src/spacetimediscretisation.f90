!##############################################################################
!# ****************************************************************************
!# <name> spacetimediscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the space time discretisation.
!# The structure t_ccoptSpaceTimeDiscretisation defines the general
!# discretisation of the underlying problem in space and time and can be
!# initialised and cleaned up with the routines in this module. The following
!# routines can be found here:
!#
!# 1.) sptidis_initDiscretisation
!#     -> Initialise a space-time discretisation structure.
!#
!# 2.) sptidis_doneDiscretisation
!#     -> Clean up a space-time discretisation structure.
!##############################################################################

module spacetimediscretisation

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use paramlist
  use timestepping
  
  use collection
  use convection
    
  use basicstructures

  implicit none

!<types>

!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  type t_ccoptSpaceTimeDiscretisation
  
    ! Spatial refinement level of this matrix.
    integer :: ilevel = 0
  
    ! Number of time steps
    !INTEGER :: niterations         = 0

    ! Absolute start time of the simulation
    !REAL(DP) :: dtimeInit          = 0.0_DP
    
    ! Maximum time of the simulation
    !REAL(DP) :: dtimeMax           = 1.0_DP
    
    ! Time step length of the time discretisation
    !REAL(DP) :: dtstep

    ! Defines the discretisation in time
    type(t_timeDiscretisation) :: rtimeDiscr
    
    ! Number of time-DOF's. For Theta-schemes, this coincides with
    ! the number of intervals, but for more complicated schemes,
    ! this may differ.
    integer :: NEQtime = 0

    ! Time-stepping scheme;
    ! 0=one step FD scheme (Euler, CN)
    ! 2=dG(0)
    !INTEGER :: ctimeStepScheme     = 0
    
    ! Parameter for one step scheme (THETA) if itimeStepScheme=0;
    ! =0:Forward Euler(instable), =1: Backward Euler, =0.5: Crank-Nicolson
    !REAL(DP) :: dtimeStepTheta     = 1.0_DP
    
    ! Regularisation parameter for the control $\alpha$. Must be <> 0.
    ! A value of 0.0 disables the terminal condition.
    real(DP) :: dalphaC = 1.0_DP
    
    ! Regularisation parameter for the terminal condition
    ! $\gamma/2*||y(T)-z(T)||$.
    ! A value of 0.0 disables the terminal condition.
    real(DP) :: dgammaC = 0.0_DP
    
    ! Problem-related structure that provides the templates for
    ! matrices/vectors on the spatial level of the matrix.
    type(t_problem_lvl), pointer :: p_rlevelInfo

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine sptidis_initDiscretisation (rproblem,ilevelTime,ilevelSpace,&
      rspaceTimeDiscr)
  
!<description>
  ! Initialises a space time discretisation structure according to the
  ! parameters in the DAT file and the problem structure.
!</description>

!<input>
  ! The problem structure describing the whole problem.
  type(t_problem), intent(IN), target :: rproblem
  
  ! 'Refinement level in time'. =1: Initialise as described in
  ! rproblem. >1: Refine (irefineInTime-1) times regularly in time.
  ! (-> #timesteps * 2**(irefineInTime-1) )
  integer, intent(IN) :: ilevelTime

  ! 'Refinement level in space'. =1: Calculate on the coarse mesh
  ! >0: calculate on level ilevelSpace. Must be <= rproblem%NLMAX,
  ! >= rproblem%NLMIN.
  integer, intent(IN) :: ilevelSpace
!</input>

!<inputoutput>
  ! Supersystem-structure to be initialised.
  type(t_ccoptSpaceTimeDiscretisation), intent(OUT) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: niterations

    ! Copy most relevant data from the problem structure.
    rspaceTimeDiscr%ilevel          = ilevelSpace
    rspaceTimeDiscr%p_rlevelInfo    => rproblem%RlevelInfo(ilevelSpace)
    niterations = rproblem%rtimedependence%niterations * 2**max(0,ilevelTime-1)
    
    ! Initialise the time discretisation
    select case (rproblem%rtimedependence%ctimeStepScheme)
    case (0)
      call tdiscr_initOneStepTheta(rspaceTimeDiscr%rtimeDiscr,&
          rproblem%rtimedependence%dtimeInit,&
          rproblem%rtimedependence%dtimeMax,&
          niterations,rproblem%rtimedependence%dtimeStepTheta)
    case (2)
      call tdiscr_initdG0(rspaceTimeDiscr%rtimeDiscr,&
          rproblem%rtimedependence%dtimeInit,&
          rproblem%rtimedependence%dtimeMax,niterations)
    case DEFAULT
      print *,'cc_initParamsSupersystem: Unsupported time discretisation.'
      call sys_halt()
    end select
    
    ! Save the number of time-DOF's
    rspaceTimeDiscr%NEQtime = tdiscr_igetNDofGlob(rspaceTimeDiscr%rtimeDiscr)
    
    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',rspaceTimeDiscr%dalphaC,1.0_DP)
    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dgammaC',rspaceTimeDiscr%dgammaC,0.0_DP)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine sptidis_doneDiscretisation (rspaceTimeDiscr)
  
!<description>
  ! Cleans up a given space time discrtisation structure.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  type(t_ccoptSpaceTimeDiscretisation), intent(OUT) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>

    ! Currently there is nothing to do!

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine sptidis_infoDiscretisation (rspaceTimeDiscr)
  
!<description>
  ! Prints statistical information about a time discretisation.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>
    
    call output_line ('NEQ in time   : '//sys_siL(rspaceTimeDiscr%NEQTime,10))
    call output_line ('NEQ in space  : '//&
      sys_siL(dof_igetNDofGlobBlock(rspaceTimeDiscr%p_rlevelInfo%rdiscretisation),10))

  end subroutine

end module
