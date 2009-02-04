!##############################################################################
!# ****************************************************************************
!# <name> euler_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global variables
!# which are required to solve the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!#  1.) euler_getNVAR
!#      Returns the number of variables depending on the spatial dimension
!#
!# </purpose>
!##############################################################################

module euler_basic

  use afcstabilisation
  use fsystem
  use linearsystemscalar
  use paramlist
  use statistics
  use thermodynamics
  use triangulation
  use problem


  implicit none

  private
  public :: t_euler
  public :: euler_getNVAR

  interface euler_getNVAR
    module procedure euler_getNvarProblemLevel
    module procedure euler_getNvarDescriptor
  end interface

!<constants>

!<constantblock description="Global type of mass matrix">

  ! no mass matrix, i.e., no time derivative
  integer, parameter, public :: MASS_ZERO       = 0
  
  ! consistent mass matrix
  integer, parameter, public :: MASS_CONSISTENT = 1

  ! lumped mass matrix
  integer, parameter, public :: MASS_LUMPED     = 2

!</constantblock>


!<constantblock description="Global type of coupling approach">

  ! time-dependent flow
  integer, parameter, public :: FLOW_SEGREGATED = 0

  ! steady-state flow
  integer, parameter, public :: FLOW_ALLCOUPLED = 1

!</constantblock>


!<constantblock description="Global types of perturbation parameters">

  ! Perturbation parameter is chosen as in the NITSOL package
  integer, parameter, public :: PERTURB_NITSOL  = 1

  ! Perturbation parameter is chosen as SQRT(machine precision)
  integer, parameter, public :: PERTURB_SQRTEPS = 2

!</constantblock>


!<constantblock description="Global constants from Gas Dynamic">

  ! ratio of specific heats
  real(DP), parameter, public :: GAMMA = GAMMA_AIR
  real(DP), parameter, public :: RGAS  = R_AIR

  ! auxiliary parameters related to gamma
  real(DP), parameter, public :: G1  = GAMMA-1.0
  real(DP), parameter, public :: G2  = (GAMMA-1.0)/2.0
  real(DP), parameter, public :: G3  = 1.0/(GAMMA-1.0)
  real(DP), parameter, public :: G4  = 1.0/GAMMA
  real(DP), parameter, public :: G5  = GAMMA/(GAMMA-1.0)
  real(DP), parameter, public :: G6  = GAMMA-2.0
  real(DP), parameter, public :: G7  = (GAMMA-1.0)/(2.0*GAMMA)
  real(DP), parameter, public :: G8  = (GAMMA+1.0)/(2.0*GAMMA)
  real(DP), parameter, public :: G9  = 2.0/(GAMMA-1.0)
  real(DP), parameter, public :: G10 = 2.0/(GAMMA+1.0)
  real(DP), parameter, public :: G11 = 2.0*GAMMA/(GAMMA-1.0)
  real(DP), parameter, public :: G12 = (GAMMA-1.0)/(GAMMA+1.0)

!</constantblock>


!<constantblock description="Global constants for number of variables">

  ! number of solution components in 1D
  integer, parameter, public :: NVAR1D = 3

  ! number of solution components in 2D
  integer, parameter, public :: NVAR2D = 4

  ! number of solution components in 3D
  integer, parameter, public :: NVAR3D = 5

!</constantblock>

!</constants>

  !*****************************************************************************

!<types>
  
!<typeblock>

  ! This structure contains all required data to describe an instance
  ! of the compressible Euler flow benchmark application
  type t_euler

    ! Dimension of the problem
    integer :: ndimension

    ! Type of the mass matrix
    integer :: imasstype

    ! Type of mass antidiffusion
    integer :: imassantidiffusion

    ! Type of coupling
    integer :: icoupled

    ! Type of preconditioner
    integer :: iprecond

    ! Type of the element(s)
    integer :: ieltype

    ! Matrix format
    integer :: imatrixFormat

    ! System format
    integer :: isystemFormat
    
    ! Timer for the solution process
    type(t_timer) :: rtimerSolution

    ! Timer for the adaptation process
    type(t_timer) :: rtimerAdaptation

    ! Timer for the error estimation process
    type(t_timer) :: rtimerErrorEstimation

    ! Timer for the triangulation process
    type(t_timer) :: rtimerTriangulation

    ! Timer for the assembly of constant coefficient matrices
    type(t_timer) :: rtimerAssemblyCoeff
    
    ! Timer for the assembly of system matrices
    type(t_timer) :: rtimerAssemblyMatrix

    ! Timer for the assembly of residual/right-hand side vectors
    type(t_timer) :: rtimerAssemblyVector

    ! Timer for pre- and post-processing
    type(t_timer) :: rtimerPrePostprocess
    
  end type t_euler

!</typeblock>

contains

  !*****************************************************************************

!<function>

  pure function euler_getNvarProblemLevel(rproblemLevel) result(NVAR)

!<description>
    ! This function returns the number of flow variables
    ! using the information from the problem level structure
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel
!</input>

!<result>
    ! number of variables
    integer :: NVAR
!</result>
!</function>

    NVAR = rproblemLevel%rtriangulation%ndim + 2
    
  end function euler_getNvarProblemLevel

  !*****************************************************************************

!<function>

  pure function euler_getNvarDescriptor(rappDescriptor) result(NVAR)

!<description>
    ! This function returns the number of flow variables
    ! using the information from the application descriptor
!</description>

!<input>
    ! application descriptor
    type(t_euler), intent(IN) :: rappDescriptor
!</input>

!<result>
    ! number of variables
    integer :: NVAR
!</result>
!</function>

    NVAR = rappDescriptor%ndimension + 2
    
  end function euler_getNvarDescriptor

!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$  
!!$  subroutine euler_getVariableNodewise(NEQ, NVAR, Dsolution, ivariable, Dvalue)
!!$
!!$!<description>
!!$    ! This subroutine extracts an individual variable from the
!!$    ! global solution vector store in interleave format.
!!$!</description>
!!$
!!$!<input>
!!$    ! Number of equations
!!$    integer, intent(IN) :: NEQ
!!$
!!$    ! Number of variables
!!$    integer, intent(IN) :: NVAR
!!$    
!!$    ! Solution vector
!!$    real(DP), dimension(NVAR,NEQ), intent(IN) :: Dsolution
!!$
!!$    ! Variable to extract
!!$    integer, intent(IN) :: ivariable
!!$!</input>
!!$
!!$!<output>
!!$    ! Value of the extracted variable
!!$    real(DP), dimension(:), intent(OUT) :: Dvalue
!!$!</output>
!!$!</subroutine>
!!$
!!$    ! local variables
!!$    real(DP) :: p
!!$    integer :: ieq
!!$    
!!$    
!!$    ! Which variable should be extracted?
!!$    if (ivariable .eq. VAR_DENSITY) then
!!$      ! density
!!$      do ieq = 1, NEQ
!!$        Dvalue(ieq) = Dsolution(VAR_DENSITY, ieq)
!!$      end do
!!$      
!!$    elseif ((ivariable .ge. VAR_MOMENTUM_X) .and.&
!!$            (ivariable .le. NVAR)) then
!!$      ! velocities and energy (located at NVAR)
!!$      do ieq = 1, NEQ
!!$        Dvalue(ieq) = Dsolution(ivariable, ieq)/&
!!$                      Dsolution(VAR_DENSITY, ieq)
!!$      end do
!!$      
!!$    elseif (ivariable .eq. NVAR+1) then
!!$      ! pressure ...
!!$      select case(NVAR)
!!$      case (NVAR1D)
!!$        ! ... in 1D
!!$        do ieq = 1, NEQ
!!$          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(NVAR1D, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$        end do
!!$
!!$      case (NVAR2D)
!!$        ! ... in 2D
!!$        do ieq = 1, NEQ
!!$          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(NVAR2D, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$        end do
!!$      
!!$      case (NVAR3D)
!!$        ! ... in 3D
!!$        do ieq = 1, NEQ
!!$          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(NVAR2D, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Z, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$        end do
!!$      end select
!!$
!!$    elseif(ivariable .eq. NVAR+2) then
!!$      ! Mach number ...
!!$      select case(NVAR)
!!$      case (NVAR1D)
!!$        ! ... in 1D
!!$        do ieq = 1, NEQ
!!$          p = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(NVAR1D, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$        end do
!!$
!!$      case (NVAR2D)
!!$        ! ... in 2D
!!$        do ieq = 1, NEQ
!!$          p = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(NVAR2D, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$        end do
!!$        
!!$      case (NVAR3D)
!!$        ! ... in 3D
!!$        do ieq = 1, NEQ
!!$          p = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(NVAR3D, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Z, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
!!$              Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq),&
!!$              Dsolution(VAR_MOMENTUM_Z, ieq)/Dsolution(VAR_DENSITY, ieq))
!!$        end do
!!$      end select
!!$      
!!$    else
!!$      call output_line('Unsupported variable!',&
!!$                        OU_CLASS_ERROR,OU_MODE_STD,'euler_getVariableNodewise')
!!$      call sys_halt()
!!$    end if
!!$  end subroutine euler_getVariableNodewise
!!$  
!!$  !*****************************************************************************
!!$  
!!$!<subroutine>
!!$  
!!$  subroutine euler_getVariableBlockwise(NEQ, NVAR, Dsolution, ivariable, Dvalue)
!!$
!!$!<description>
!!$    ! This subroutine extracts an individual variable from the
!!$    ! global solution vector store in block format.
!!$!</description>
!!$
!!$!<input>
!!$    ! Number of equations
!!$    integer, intent(IN) :: NEQ
!!$
!!$    ! Number of variables
!!$    integer, intent(IN) :: NVAR
!!$    
!!$    ! Solution vector
!!$    real(DP), dimension(NEQ,NVAR), intent(IN) :: Dsolution
!!$
!!$    ! Variable to extract
!!$    integer, intent(IN) :: ivariable
!!$!</input>
!!$
!!$!<output>
!!$    ! Value of the extracted variable
!!$    real(DP), dimension(NEQ), intent(OUT) :: Dvalue
!!$!</output>
!!$!</subroutine>     
!!$    
!!$    ! local variables
!!$    real(DP) :: p
!!$    integer :: ieq
!!$    
!!$    ! Which variable should be extracted?
!!$    if (ivariable .eq. VAR_DENSITY) then
!!$      ! density
!!$      do ieq = 1, NEQ
!!$        Dvalue(ieq) = Dsolution(ieq, VAR_DENSITY)
!!$      end do
!!$      
!!$    elseif ((ivariable .ge. VAR_MOMENTUM_X) .and.&
!!$            (ivariable .le. NVAR)) then
!!$      ! velocities and energy (located at NVAR)
!!$      do ieq = 1, NEQ
!!$        Dvalue(ieq) = Dsolution(ieq,ivariable)/&
!!$                      Dsolution(ieq, VAR_DENSITY)
!!$      end do
!!$      
!!$    elseif(ivariable .eq. NVAR+1) then
!!$      ! pressure ...
!!$      select case(NVAR)
!!$      case (NVAR1D)
!!$        ! ... in 1D
!!$        do ieq = 1, NEQ
!!$          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(ieq, NVAR1D)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY))
!!$        end do
!!$
!!$      case (NVAR2D)
!!$        ! ... in 2D
!!$        do ieq = 1, NEQ
!!$          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(ieq, NVAR2D)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY))
!!$        end do
!!$        
!!$      case (NVAR3D)
!!$        ! ... in 3D
!!$        do ieq = 1, NEQ
!!$          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(ieq, NVAR3D)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Z)/Dsolution(ieq, VAR_DENSITY))
!!$        end do
!!$      end select
!!$      
!!$    elseif (ivariable .eq. NVAR+2) then
!!$      ! Mach number ...
!!$      select case(NVAR)
!!$      case (NVAR1D)
!!$        ! ... in 1D
!!$        do ieq = 1, NEQ
!!$          p = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(ieq, NVAR1D)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY))
!!$          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY))
!!$        end do
!!$
!!$      case (NVAR2D)
!!$        ! ... in 2D
!!$        do ieq = 1, NEQ
!!$          p = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(ieq, NVAR2D)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY))
!!$          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY))
!!$        end do
!!$
!!$      case (NVAR3D)
!!$        ! ... in 3D
!!$        do ieq = 1, NEQ
!!$          p = thdyn_pressure(GAMMA_AIR,&
!!$              Dsolution(ieq, NVAR3D)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Z)/Dsolution(ieq, VAR_DENSITY))
!!$          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
!!$              Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY),&
!!$              Dsolution(ieq, VAR_MOMENTUM_Z)/Dsolution(ieq, VAR_DENSITY))
!!$        end do
!!$      end select
!!$      
!!$    else
!!$      call output_line('Unsupported variable!',&
!!$                        OU_CLASS_ERROR,OU_MODE_STD,'euler_getVariableBlockwise')
!!$      call sys_halt()
!!$    end if
!!$  end subroutine euler_getVariableBlockwise

end module euler_basic
