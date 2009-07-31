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
!# 5.) euler_setVariables
!#     -> Sets the conservative variables from an UCD import
!#
!# </purpose>
!##############################################################################

module euler_basic

  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use problem
  use statistics
  use thermodynamics
  use triangulation
  use ucd

  implicit none

  private
  public :: euler_getNVAR
  public :: euler_getVariable
  public :: euler_getVarInterleaveFormat
  public :: euler_getVarBlockFormat
  public :: euler_setVariables

!<constants>

!<constantblock description="Global type of mass matrix">

  ! no mass matrix, i.e., no time derivative
  integer, parameter, public :: MASS_ZERO       = 0
  
  ! consistent mass matrix
  integer, parameter, public :: MASS_CONSISTENT = 1

  ! lumped mass matrix
  integer, parameter, public :: MASS_LUMPED     = 2

!</constantblock>


!<constantblock description="Global type of initial solution">

  ! zero initial solution
  integer, parameter, public :: SOLUTION_ZERO     = 0

  ! analytical initial solution
  integer, parameter, public :: SOLUTION_ANALYTIC = 1

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

  ! Scalar dissipation adopting dimensional splitting
  integer, parameter, public :: DISSIPATION_SCALAR_DSPLIT = -DISSIPATION_SCALAR

  ! Tensorial dissipation adopting dimensional splitting
  integer, parameter, public :: DISSIPATION_TENSOR_DSPLIT = -DISSIPATION_TENSOR

  ! Rusanov flux adopting dimensional splitting
  integer, parameter, public :: DISSIPATION_RUSANOV_DSPLIT = -DISSIPATION_RUSANOV

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

contains

  !*****************************************************************************

!<function>

  pure function euler_getNVAR(rproblemLevel) result(NVAR)

!<description>
    ! This function returns the number of flow variables
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel
!</input>

!<result>
    ! number of variables
    integer :: NVAR
!</result>
!</function>

    NVAR = rproblemLevel%rtriangulation%ndim + 2
    
  end function euler_getNVAR

  !*****************************************************************************

!<subroutine>

  subroutine euler_getVariable(rvectorBlock, cvariable, rvectorScalar)

!<description>
    ! This subroutine extracts a single variable from the vector of
    ! conservative variables which is stored in interleave of block format
!</description>

!<input>
    ! Vector of conservative variables
    type(t_vectorBlock), intent(in) :: rvectorBlock

    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable
!</input>

!<inputoutput>
    ! Extracted single variable
    type(t_vectorScalar), intent(inout) :: rvectorScalar
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
      call euler_getVarInterleaveFormat(neq, nvar, cvariable,&
          p_Ddata, p_Dvalue)

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
      call euler_getVarBlockFormat(neq, nvar, cvariable,&
          p_Ddata, p_Dvalue)

    end if
    
    ! Attach spatial discretization from first subvector
    rvectorScalar%p_rspatialDiscr =>&
        rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
    
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
    integer, intent(in) :: neq

    ! Number of variables
    integer, intent(in) :: nvar
    
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable
    
    ! Vector of conservative variables
    real(DP), dimension(nvar,neq), intent(in) :: Ddata
!</input>

!<output>
    ! Extracted single variable
    real(DP), dimension(:), intent(out) :: Dvalue
!</output>    
!</subroutine>
    
    ! local variables
    real(DP) :: p
    integer :: ieq
    
    
    select case (trim(adjustl(sys_upcase(cvariable))))
    case ('DENSITY')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(1, ieq)
      end do
      !$omp end parallel do

    case ('VELOCITY_MAGNITUDE')
      select case(nvar)
      case (NVAR1D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = abs(Ddata(2, ieq))/Ddata(1, ieq)
        end do
        !$omp end parallel do

      case (NVAR2D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = sqrt(Ddata(2, ieq)**2 +&
                             Ddata(3, ieq)**2)/Ddata(1, ieq)
        end do
        !$omp end parallel do

      case (NVAR3D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = sqrt(Ddata(2, ieq)**2 +&
                             Ddata(3, ieq)**2 +&
                             Ddata(4, ieq)**2)/Ddata(1, ieq)
        end do
        !$omp end parallel do
      end select
      
    case ('VELOCITY_X')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(2, ieq)/Ddata(1, ieq)
      end do
      !$omp end parallel do
      
    case ('VELOCITY_Y')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(3, ieq)/Ddata(1, ieq)
      end do
      !$omp end parallel do
      
    case ('VELOCITY_Z')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(4, ieq)/Ddata(1, ieq)
      end do
      !$omp end parallel do

    case ('MOMENTUM_X')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(2, ieq)
      end do
      !$omp end parallel do
      
    case ('MOMENTUM_Y')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(3, ieq)
      end do
      !$omp end parallel do
      
    case ('MOMENTUM_Z')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(4, ieq)
      end do
      !$omp end parallel do
      
    case ('ENERGY')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(nvar, ieq)/Ddata(1, ieq)
      end do
      !$omp end parallel do

    case ('EFFECTIVE_ENERGY')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(nvar, ieq)
      end do
      !$omp end parallel do
      
    case ('PRESSURE')
      select case (nvar)
      case (NVAR1D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq))
        end do
        !$omp end parallel do
        
      case (NVAR2D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq))
        end do
        !$omp end parallel do
        
      case (NVAR3D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq),&
              Ddata(3, ieq)/Ddata(1, ieq),&
              Ddata(4, ieq)/Ddata(1, ieq))
        end do
        !$omp end parallel do
      end select
      
    case ('MACHNUMBER')
      select case (nvar)
      case (NVAR1D)
        !$omp parallel do private(p)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(nvar, ieq)/Ddata(1, ieq), Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(1, ieq),&
              Ddata(2, ieq)/Ddata(1, ieq))
        end do
        !$omp end parallel do
        
      case (NVAR2D)
        !$omp parallel do private(p)
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
        !$omp end parallel do
        
      case (NVAR3D)
        !$omp parallel do private(p)
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
        !$omp end parallel do
      end select

    case DEFAULT
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = 0.0_DP
      end do
      !$omp end parallel do      
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
    integer, intent(in) :: neq

    ! Number of variables
    integer, intent(in) :: nvar
    
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable
    
    ! Vector of conservative variables
    real(DP), dimension(neq,nvar), intent(in) :: Ddata
!</input>

!<output>
    ! Extracted single variable
    real(DP), dimension(:), intent(out) :: Dvalue
!</output>    
!</subroutine>
      
    ! local variables
    real(DP) :: p
    integer :: ieq
    
    
    select case (trim(adjustl(sys_upcase(cvariable))))
    case ('DENSITY')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 1)
      end do
      !$omp end parallel do

    case ('VELOCITY_MAGNITUDE')
      select case(nvar)
      case (NVAR1D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = abs(Ddata(ieq, 2))/Ddata(ieq, 1)
        end do
        !$omp end parallel do

      case (NVAR2D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = sqrt(Ddata(ieq, 2)**2 +&
                             Ddata(ieq, 3)**2)/Ddata(ieq, 1)
        end do
        !$omp end parallel do

      case (NVAR3D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = sqrt(Ddata(ieq, 2)**2 +&
                             Ddata(ieq, 3)**2 +&
                             Ddata(ieq, 4)**2)/Ddata(ieq, 1)
        end do
        !$omp end parallel do
      end select
      
    case ('VELOCITY_X')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 2)/Ddata(ieq, 1)
      end do
      !$omp end parallel do
      
    case ('VELOCITY_Y')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 3)/Ddata(ieq, 1)
      end do
      !$omp end parallel do
      
    case ('VELOCITY_Z')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 4)/Ddata(ieq, 1)
      end do
      !$omp end parallel do

    case ('MOMENTUM_X')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 2)
      end do
      !$omp end parallel do
      
    case ('MOMENTUM_Y')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 3)
      end do
      !$omp end parallel do
      
    case ('MOMENTUM_Z')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, 4)
      end do
      !$omp end parallel do
      
    case ('ENERGY')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, nvar)/Ddata(ieq, 1)
      end do
      !$omp end parallel do

    case ('EFFECTIVE_ENERGY')
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Ddata(ieq, nvar)
      end do
      !$omp end parallel do
      
    case ('PRESSURE')
      select case (nvar)
      case (NVAR1D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1))
        end do
        !$omp end parallel do
        
      case (NVAR2D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1))
        end do
        !$omp end parallel do
        
      case (NVAR3D)
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1),&
              Ddata(ieq, 3)/Ddata(ieq, 1),&
              Ddata(ieq, 4)/Ddata(ieq, 1))
        end do
        !$omp end parallel do
      end select
      
    case ('MACHNUMBER')
      select case (nvar)
      case (NVAR1D)
        !$omp parallel do private(p)
        do ieq = 1, neq
          p = thdyn_pressure(GAMMA,&
              Ddata(ieq, nvar)/Ddata(ieq, 1), Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1))
          
          Dvalue(ieq) = thdyn_Machnumber(GAMMA,&
              p, Ddata(ieq, 1),&
              Ddata(ieq, 2)/Ddata(ieq, 1))
        end do
        !$omp end parallel do
        
      case (NVAR2D)
        !$omp parallel do private(p)
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
        !$omp end parallel do
        
      case (NVAR3D)
        !$omp parallel do private(p)
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
        !$omp end parallel do
      end select

    case DEFAULT
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = 0.0_DP
      end do
      !$omp end parallel do

    end select
    
  end subroutine euler_getVarBlockformat

  !*****************************************************************************

!<subroutine>
  
  subroutine euler_setVariables(rexport, rvectorBlock)

!<description>
    ! This subroutine sets the conservative variables from a UCD export
!</description>

!<input>
    ! UCD export
    type(t_ucdExport), intent(in) :: rexport
!</input>

!<inputoutput>
    ! Solution vector
    type(t_vectorBlock), intent(inout) :: rvectorBlock
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata, p_Ddensity
    real(DP), dimension(:), pointer :: p_Dvelocity_x, p_Dvelocity_y
    real(DP), dimension(:), pointer :: p_Denergy
    integer :: neq, nvar, nlength

    ! Check if we are in interleave of block format
    if (rvectorBlock%nblocks .eq. 1) then

      ! Set dimensions
      neq  = rvectorBlock%RvectorBlock(1)%NEQ
      nvar = rvectorBlock%RvectorBlock(1)%NVAR

      ! Set pointer
      call lsysbl_getbase_double(rvectorBlock, p_Ddata)

      ! Get density values
      call ucd_getVariable(rexport, 'density', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Ddensity(neq))
      call ucd_getVariable(rexport, 'density', p_Ddensity)

      ! Set density variable
      call setVarInterleaveFormat(neq, nvar, 1,&
          p_Ddensity, p_Ddata)


      ! Get velocity in x-direction
      call ucd_getVariable(rexport, 'velocity_X', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Dvelocity_x(neq))
      call ucd_getVariable(rexport, 'velocity_X', p_Dvelocity_x)

      ! Generate momentum in x-direction
      p_Dvelocity_x = p_Dvelocity_x * p_Ddensity

      ! Set momentum variable in x-direction
      call setVarInterleaveFormat(neq, nvar, 2,&
          p_Dvelocity_x, p_Ddata)

      
      ! Get velocity in y-direction
      call ucd_getVariable(rexport, 'velocity_Y', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Dvelocity_y(neq))
      call ucd_getVariable(rexport, 'velocity_y', p_Dvelocity_y)

      ! Generate momentum in y-direction
      p_Dvelocity_y = p_Dvelocity_y * p_Ddensity

      ! Set momentum variable in y-direction
      call setVarInterleaveFormat(neq, nvar, 3,&
          p_Dvelocity_y, p_Ddata)


      ! Get total energy
      call ucd_getVariable(rexport, 'energy', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Denergy(neq))
      call ucd_getVariable(rexport, 'energy', p_Denergy)

      ! Generate total energy
      p_Denergy = p_Denergy * p_Ddensity

      ! Set energy variable
      call setVarInterleaveFormat(neq, nvar, 4,&
          p_Denergy, p_Ddata)

      ! Deallocate temporal memory
      deallocate(p_Ddensity, p_Dvelocity_x, p_Dvelocity_y, p_Denergy)
      
    else

      ! Set dimensions
      neq  = int(rvectorBlock%NEQ/rvectorBlock%nblocks)
      nvar = rvectorBlock%nblocks

      ! Set pointer
      call lsysbl_getbase_double(rvectorBlock, p_Ddata)

      ! Get density values
      call ucd_getVariable(rexport, 'density', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Ddensity(neq))
      call ucd_getVariable(rexport, 'density', p_Ddensity)

      ! Set density variable
      call setVarBlockFormat(neq, nvar, 1,&
          p_Ddensity, p_Ddata)


      ! Get velocity in x-direction
      call ucd_getVariable(rexport, 'velocity_X', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Dvelocity_x(neq))
      call ucd_getVariable(rexport, 'velocity_X', p_Dvelocity_x)

      ! Generate momentum in x-direction
      p_Dvelocity_x = p_Dvelocity_x * p_Ddensity

      ! Set momentum variable in x-direction
      call setVarBlockFormat(neq, nvar, 2,&
          p_Dvelocity_x, p_Ddata)

      
      ! Get velocity in y-direction
      call ucd_getVariable(rexport, 'velocity_Y', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Dvelocity_y(neq))
      call ucd_getVariable(rexport, 'velocity_y', p_Dvelocity_y)

      ! Generate momentum in y-direction
      p_Dvelocity_y = p_Dvelocity_y * p_Ddensity

      ! Set momentum variable in y-direction
      call setVarBlockFormat(neq, nvar, 3,&
          p_Dvelocity_y, p_Ddata)


      ! Get total energy
      call ucd_getVariable(rexport, 'energy', nlength=nlength)
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'euler_setVariables')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(p_Denergy(neq))
      call ucd_getVariable(rexport, 'energy', p_Denergy)

      ! Generate total energy
      p_Denergy = p_Denergy * p_Ddensity

      ! Set energy variable
      call setVarBlockFormat(neq, nvar, 4,&
          p_Denergy, p_Ddata)

      ! Deallocate temporal memory
      deallocate(p_Ddensity, p_Dvelocity_x, p_Dvelocity_y, p_Denergy)

    end if

  contains

    ! Here, the working routine follow
    
    !**************************************************************
    ! Set variable stored in interleave format

    pure subroutine setVarInterleaveformat(neq, nvar, ivar, Dvalue, Ddata)

      integer, intent(in) :: neq, nvar, ivar
      real(DP), dimension(:), intent(in) :: Dvalue
      
      real(DP), dimension(nvar,neq), intent(inout) :: Ddata

      ! local variables
      integer :: ieq

      !$omp parallel do
      do ieq = 1, neq
        Ddata(ivar, ieq) = Dvalue(ieq)
      end do
      !$omp end parallel do

    end subroutine setVarInterleaveformat

    !**************************************************************
    ! Set variable stored in block format

    pure subroutine setVarBlockformat(neq, nvar, ivar, Dvalue, Ddata)

      integer, intent(in) :: neq, nvar, ivar
      real(DP), dimension(:), intent(in) :: Dvalue
      
      real(DP), dimension(neq,nvar), intent(inout) :: Ddata

      ! local variables
      integer :: ieq

      !$omp parallel do
      do ieq = 1, neq
        Ddata(ieq, ivar) = Dvalue(ieq)
      end do
      !$omp end parallel do

    end subroutine setVarBlockformat
   
  end subroutine euler_setVariables
 
end module euler_basic
