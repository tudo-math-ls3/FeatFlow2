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
!#  1.) euler_getNVAR = euler_getNVARProblemLevel /
!#                      euler_getNVARDescriptor
!#      Returns the number of variables depending on the spatial dimension
!#
!# 2.) euler_getVariable
!#     -> Extracts a single variable from the vector of conservative
!#        variables stored in interleave or block format
!#
!# 3.) euler_getVarInterleaveFormat
!#     -> Extracts a single variable from the scalar vector of 
!#        conservative variables stored in interleave format
!#
!# 4.) euler_getVarBlockFormat
!#     -> Extracts a single variable from the vector of conservative
!#        variables stored in block format
!#
!# </purpose>
!##############################################################################

module euler_basic

  use fparser
  use fsystem
  use linearsystemblock
  use linearsystemscalar
  use problem
  use statistics
  use thermodynamics
  use triangulation

  implicit none

  private
  public :: t_euler
  public :: euler_getNVAR
  public :: euler_getVariable
  public :: euler_getVarInterleaveFormat
  public :: euler_getVarBlockFormat

  interface euler_getNVAR
    module procedure euler_getNVARProblemLevel
    module procedure euler_getNVARDescriptor
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


!<constantblock description="Global type of system coupling approach">

  ! time-dependent flow
  integer, parameter, public :: SYSTEM_SEGREGATED = 0

  ! steady-state flow
  integer, parameter, public :: SYSTEM_ALLCOUPLED = 1

!</constantblock>


!<constantblock description="Global type of dissipation">

  ! Scalar dissipation
  integer, parameter, public :: DISSIPATION_SCALAR = 1

  ! Tensorial dissipation
  integer, parameter, public :: DISSIPATION_TENSOR = 2

  ! Rusanov flux
  integer, parameter, public :: DISSIPATION_RUSANOV = 3

!</constantblock>


!<constantblock description="Global type of recovery-based error estimation">

  ! L2-projection
  integer, parameter, public :: ERREST_L2PROJECTION  = 1

  ! Superconvergent patch recovery (vertex-based)
  integer, parameter, public :: ERREST_SPR_VERTEX    = 2

  ! Superconvergent patch recovery (element-based)
  integer, parameter, public :: ERREST_SPR_ELEMENT   = 3

  ! Superconvergent patch recovery (face-based)
  integer, parameter, public :: ERREST_SPR_FACE      = 4
 
  ! Limited averaging gradient recovery
  integer, parameter, public :: ERREST_LIMAVR        = 5

  ! Second-difference indicator (by Loehner)
  integer, parameter, public :: ERREST_SECONDDIFF    = 6

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

  ! universal gas constant
  real(DP), parameter, public :: RGAS  = RGAS_AIR

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
  real(DP), parameter, public :: G13 = 3.0-GAMMA
  real(DP), parameter, public :: G14 = (GAMMA-3.0)/2.0

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
    integer :: imassantidiffusiontype

    ! Type of dissipation
    integer :: idissipationtype

    ! Type of coupling
    integer :: isystemCoupling

    ! Type of preconditioner
    integer :: isystemPrecond

    ! Type of the element(s)
    integer :: celement

    ! Matrix format
    integer :: imatrixFormat

    ! System format
    integer :: isystemFormat
    
    ! Timer for the solution process
    type(t_timer) :: rtimerSolution

    ! Timer for the nonlinear solver
    type(t_timer) :: rtimerNonlinearSolver

    ! Timer for the linear solver
    type(t_timer) :: rtimerLinearSolver

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

!</types>

contains

  !*****************************************************************************

!<function>

  pure function euler_getNVARProblemLevel(rproblemLevel) result(NVAR)

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
    
  end function euler_getNVARProblemLevel

  !*****************************************************************************

!<function>

  pure function euler_getNVARDescriptor(rappDescriptor) result(NVAR)

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
    
  end function euler_getNVARDescriptor

  !*****************************************************************************

!<subroutine>

  subroutine euler_getVariable(rvectorBlock, cvariable, rvectorScalar)

!<description>
    ! This subroutine extracts a single variable from the vector of
    ! conservative variables which is stored in interleave of block format
!</description>

!<input>
    ! Vector of conservative variables
    type(t_vectorBlock), intent(IN) :: rvectorBlock

    ! Identifier for the variable
    character(LEN=*), intent(IN) :: cvariable
!</input>

!<inputoutput>
    ! Extracted single variable
    type(t_vectorScalar), intent(INOUT) :: rvectorScalar
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata, p_Dvalue
    integer :: neq, nvar

    
    ! Check if we are in interleave of block format
    if (rvectorBlock%nblocks .eq. 1) then

      ! Set dimensions
      neq  = rvectorBlock%RvectorBlock(1)%NEQ
      nvar = rvectorBlock%RvectorBlock(1)%NVAR
      
      ! Check if vector is compatible
      if (rvectorScalar%NEQ .eq. 0) then
        call lsyssc_createVector(rvectorScalar, neq, .false.)
      elseif (rvectorScalar%NEQ .lt. neq) then
        call lsyssc_resizeVector(rvectorScalar, neq, .false., .false.)
      end if

      ! Set pointers
      call lsysbl_getbase_double(rvectorBlock, p_Ddata)
      call lsyssc_getbase_double(rvectorScalar, p_Dvalue)

      ! Fill the scalar vector with data from the variable
      call euler_getVarInterleaveFormat(neq, nvar, cvariable, p_Ddata, p_Dvalue)

    else

      ! Set dimensions
      neq  = int(rvectorBlock%NEQ/rvectorBlock%nblocks)
      nvar = rvectorBlock%nblocks
      
      ! Check if vector is compatible
      if (rvectorScalar%NEQ .eq. 0) then
        call lsyssc_createVector(rvectorScalar, neq, .false.)
      elseif (rvectorScalar%NEQ .lt. neq) then
        call lsyssc_resizeVector(rvectorScalar, neq, .false., .false.)
      end if

      ! Set pointers
      call lsysbl_getbase_double(rvectorBlock, p_Ddata)
      call lsyssc_getbase_double(rvectorScalar, p_Dvalue)

      ! Fill the scalar vector with data from the variable
      call euler_getVarBlockFormat(neq, nvar, cvariable, p_Ddata, p_Dvalue)

    end if
    
    ! Attach spatial discretization from first subvector
    rvectorScalar%p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
    
  end subroutine euler_getVariable

  !*****************************************************************************
  
!<subroutine>

  pure subroutine euler_getVarInterleaveFormat(neq, nvar, cvariable, Ddata, Dvalue)

!<description>
    ! This subroutine extracs a single variable from the vector of
    ! conservative variabels which is stored in interleave format
!</description>

!<input>
    ! Number of equations
    integer, intent(IN) :: neq

    ! Number of variables
    integer, intent(IN) :: nvar
    
    ! Identifier for the variable
    character(LEN=*), intent(IN) :: cvariable
    
    ! Vector of conservative variables
    real(DP), dimension(nvar,neq), intent(IN) :: Ddata
!</input>

!<output>
    ! Extracted single variable
    real(DP), dimension(:), intent(OUT) :: Dvalue
!</output>    
!</subroutine>
    
    ! local variables
    real(DP) :: p
    integer :: ieq,ivar
    
    
    select case (trim(adjustl(sys_upcase(cvariable))))
    case ('DENSITY')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(1, ieq)
      end do
      
    case ('VELOCITY_X')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(2, ieq)/Ddata(1, ieq)
      end do
      
    case ('VELOCITY_Y')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(3, ieq)/Ddata(1, ieq)
      end do
      
    case ('VELOCITY_Z')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(4, ieq)/Ddata(1, ieq)
      end do
      
    case ('ENERGY')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(nvar, ieq)/Ddata(1, ieq)
      end do

    case ('EFFECTIVE_ENERGY')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(nvar, ieq)
      end do
      
    case ('PRESSURE')
      select case (nvar)
      case (NVAR1D)
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq))
        end do
        
      case (NVAR2D)
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq))
        end do
        
      case (NVAR3D)
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq),&
              Ddata(4, ieq)/Ddata(1, ieq))
        end do
      end select
      
    case ('MACHNUMBER')
      select case (nvar)
      case (NVAR1D)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq))
        end do
        
      case (NVAR2D)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq))
        end do
        
      case (NVAR3D)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq),&
              Ddata(4, ieq)/Ddata(1, ieq))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq),&
              Ddata(4, ieq)/Ddata(1, ieq))
        end do
      end select

    case DEFAULT
      do ieq = 1, neq
        Dvalue(ieq) = 0.0_DP
      end do
      
    end select
    
  end subroutine euler_getVarInterleaveFormat
  
  !*****************************************************************************

!<subroutine>
  
   pure subroutine euler_getVarBlockformat(neq, nvar, cvariable, Ddata, Dvalue)

!<description>
    ! This subroutine extracs a single variable from the vector of
    ! conservative variabels which is stored in block format
!</description>

!<input>
    ! Number of equations
    integer, intent(IN) :: neq

    ! Number of variables
    integer, intent(IN) :: nvar
    
    ! Identifier for the variable
    character(LEN=*), intent(IN) :: cvariable
    
    ! Vector of conservative variables
    real(DP), dimension(neq,nvar), intent(IN) :: Ddata
!</input>

!<output>
    ! Extracted single variable
    real(DP), dimension(:), intent(OUT) :: Dvalue
!</output>    
!</subroutine>
      
    ! local variables
    real(DP) :: p
    integer :: ieq,ivar
    
    
    select case (trim(adjustl(sys_upcase(cvariable))))
    case ('DENSITY')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 1)
      end do
      
    case ('VELOCITY_X')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 2)/Ddata(ieq, 1)
      end do
      
    case ('VELOCITY_Y')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 3)/Ddata(ieq, 1)
      end do
      
    case ('VELOCITY_Z')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 4)/Ddata(ieq, 1)
      end do
      
    case ('ENERGY')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, nvar)/Ddata(ieq, 1)
      end do

    case ('EFFECTIVE_ENERGY')
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, nvar)
      end do
      
    case ('PRESSURE')
      select case (nvar)
      case (NVAR1D)
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1))
        end do
        
      case (NVAR2D)
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1))
        end do
        
      case (NVAR3D)
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1),&
              Ddata(ieq, 4)/Ddata(ieq, 1))
        end do
      end select
      
    case ('MACHNUMBER')
      select case (nvar)
      case (NVAR1D)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1))
        end do
        
      case (NVAR2D)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1))
        end do
        
      case (NVAR3D)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1),&
              Ddata(ieq, 4)/Ddata(ieq, 1))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1),&
              Ddata(ieq, 4)/Ddata(ieq, 1))
        end do
      end select

    case DEFAULT
      do ieq = 1, neq
        Dvalue(ieq) = 0.0_DP
      end do

    end select
    
  end subroutine euler_getVarBlockformat

end module euler_basic
