!##############################################################################
!# ****************************************************************************
!# <name> transport_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global vairables
!# which are required to solve a conservation law for a scalar variable.
!#
!# The following routines are available:
!#
!# 1.) transp_hasVelocityVector
!#     -> Checks if the velocity vector is given as explicit vector
!#
!# 2.) transp_setVariable
!#     -> Sets the scalar variable from an UCD import
!#
!# </purpose>
!##############################################################################

module transport_basic

  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use statistics
  use ucd

  implicit none

  private

  public :: transp_hasVelocityVector
  public :: transp_setVariable

!<constants>

!<constantblock description="Global type of mass matrix">

  ! no mass matrix, i.e., no time derivative
  integer, parameter, public :: MASS_ZERO       = 0

  ! consistent mass matrix
  integer, parameter, public :: MASS_CONSISTENT = 1

  ! lumped mass matrix
  integer, parameter, public :: MASS_LUMPED     = 2

!</constantblock>


!<constantblock description="Global type of flow velocities">

  ! velocity profile is defined externally
  integer, parameter, public :: VELOCITY_EXTERNAL          =-1

  ! zero velocity profile v=0
  integer, parameter, public :: VELOCITY_ZERO              = 0

  ! constant velocity profile v=v(x)
  integer, parameter, public :: VELOCITY_CONSTANT          = 1

  ! linear time-dependent velocity profile v=v(x,t)
  integer, parameter, public :: VELOCITY_TIMEDEP           = 2

  ! nonlinear Burgers' equation in space-time
  integer, parameter, public :: VELOCITY_BURGERS_SPACETIME = 3

  ! nonlinear Buckley-Leverett equation in space-time
  integer, parameter, public :: VELOCITY_BUCKLEV_SPACETIME = 4

  ! nonlinear Burgers' equation in 1D
  integer, parameter, public :: VELOCITY_BURGERS1D         = 5

  ! nonlinear Burgers' equation in 2D
  integer, parameter, public :: VELOCITY_BURGERS2D         = 6

  ! nonlinear Burgers' equation in 3D
  integer, parameter, public :: VELOCITY_BURGERS3D         = 7

  ! nonlinear Buckley-Leverett equation in 1D
  integer, parameter, public :: VELOCITY_BUCKLEV1D         = 8

!</constantblock>


!<constantblock description="Global type of diffusion">

  ! zero diffusion
  integer, parameter, public :: DIFFUSION_ZERO        = 0

  ! isotropic diffusion
  integer, parameter, public :: DIFFUSION_ISOTROPIC   = 1

  ! anisotropic diffusion
  integer, parameter, public :: DIFFUSION_ANISOTROPIC = 2

  ! variable diffusion
  integer, parameter, public :: DIFFUSION_VARIABLE    = 3

!</constantblock>


!<constantblock description="Global type of reaction">

  ! zero reaction
  integer, parameter, public :: REACTION_ZERO = 0

!</constantblock>


!<constantblock description="Global type of right-hand side">

  ! zero right-hand side
  integer, parameter, public :: RHS_ZERO     = 0

  ! analytical right-hand side
  integer, parameter, public :: RHS_ANALYTIC = 1

!</constantblock>

!<constantblock description="Global type of initial solution">

  ! zero solution
  integer, parameter, public :: SOLUTION_ZERO                   = 0

  ! analytical solution: given by pointwise values
  integer, parameter, public :: SOLUTION_ANALYTIC_POINTVALUE    = 1

  ! graymap profile for solution
  integer, parameter, public :: SOLUTION_GRAYMAP                = 2

  ! analytical solution: given by consistent L2-projection
  integer, parameter, public :: SOLUTION_ANALYTIC_L2_CONSISTENT = 3

  ! analytical solution: given by lumped L2-projection
  integer, parameter, public :: SOLUTION_ANALYTIC_L2_LUMPED     = 4

!</constantblock>


!<constantblock description="Global type of target functional">

  ! zero target functional
  integer, parameter, public :: TFUNC_ZERO     = 0

  ! volume integral target functional
  integer, parameter, public :: TFUNC_VOLINTG  = 1

  ! surface integral target functional
  integer, parameter, public :: TFUNC_SURFINTG = 2

  ! mixed volume + surface integral target functional
  integer, parameter, public :: TFUNC_MIXINTG  = 3

!</constantblock>


!<constantblock description="Global type of recovery-based error estimation">

  ! L2-projection
  integer, parameter, public :: ERREST_L2PROJECTION = 1

  ! Superconvergent patch recovery (vertex-based)
  integer, parameter, public :: ERREST_SPR_VERTEX   = 2

  ! Superconvergent patch recovery (element-based)
  integer, parameter, public :: ERREST_SPR_ELEMENT  = 3

  ! Superconvergent patch recovery (face-based)
  integer, parameter, public :: ERREST_SPR_FACE     = 4

  ! Limited averaging gradient recovery
  integer, parameter, public :: ERREST_LIMAVR       = 5

  ! Second-difference indicator (by Loehner)
  integer, parameter, public :: ERREST_SECONDDIFF   = 6

!</constantblock>


!<constantblock description="Global constants for error redistribution">

  ! Use error 'as is'
  integer, parameter, public :: ERREST_ASIS        = 0

  ! Equidistribution of error
  integer, parameter, public :: ERREST_EQUIDIST    = 1

  ! Logarithmic equidistribution of error
  integer, parameter, public :: ERREST_LOGEQUIDIST = 2

  ! Automatic treshold based on RMS
  integer, parameter, public :: ERREST_AUTORMS     = 3

!</constantblock>


!<constantblock description="Global types of perturbation parameters">

  ! Perturbation parameter is chosen as in the NITSOL package
  integer, parameter, public :: PERTURB_NITSOL  = 1

  ! Perturbation parameter is chosen as SQRT(machine precision)
  integer, parameter, public :: PERTURB_SQRTEPS = 2

!</constantblock>


!<constantblock description="Global types of boundary conditions">

  ! Homogeneous Neumann boundary conditions
  integer, parameter, public :: BDRC_HOMNEUMANN   = 1

  ! Inhomogeneous Neumann boundary conditions
  integer, parameter, public :: BDRC_INHOMNEUMANN = 2

  ! Dirichlet boundary conditions
  integer, parameter, public :: BDRC_DIRICHLET    = 3

  ! Penalty parameter for Dirichlet boundary conditions
  real(DP), parameter, public :: BDRC_DIRICHLET_PENALTY = 1e6

  ! Robin boundary conditions
  integer, parameter, public :: BDRC_ROBIN        = 4

  ! Flux boundary condition
  integer, parameter, public :: BDRC_FLUX         = 5
!</constantblock>

!<constantblock description="Global types of coordinate systems">

  ! Cartesian coordinate system (x,y,z)
  integer, parameter, public :: COORDS_CARTESIAN           = 0

  ! Axi-symmetric coordinate system (r,z) with symmetry around the
  ! z-axis (2D approximation of a 3D flow)
  integer, parameter, public :: COORDS_AXIALSYMMETRY       = 1

  ! Cylindrically symmtric coordinate system (r)
  ! (1D approximation of a 2D flow)
  integer, parameter, public :: COORDS_CYLINDRICALSYMMETRY = 2

  ! Spherically symmetric coordinate system (r)
  ! (1D approximation of a 3D flow)
  integer, parameter, public :: COORDS_SPHERICALSYMMETRY   = 3

!</constants>

contains

  !*****************************************************************************

!<function>

  pure function transp_hasVelocityVector(ivelocityType) result(bvector)

!<description>
    ! This function returns .true. if the velocity vector is given explicitly.
    ! If there is no velocity vector, then it returns .false.
!</description>

!<input>
    ! Type of velocity
    integer, intent(in) :: ivelocityType
!</input>

!<result>
    ! .true. if the velocity is given as explicit vector
    logical :: bvector
!</result>
!</function>

    select case(abs(ivelocityType))

    case(VELOCITY_EXTERNAL, VELOCITY_CONSTANT, VELOCITY_TIMEDEP)
      bvector = .true.

    case default
      bvector = .false.
    end select

  end function transp_hasVelocityVector

  !*****************************************************************************

!<subroutine>

  subroutine transp_setVariable(rexport, cvariable, rvectorBlock, iblock)

!<description>
    ! This subroutine sets the scalar variable from a UCD export
!</description>

!<input>
    ! UCD export
    type(t_ucdExport), intent(in) :: rexport

    ! Name of the veariable in the UCD export
    character(LEN=*), intent(in) :: cvariable

    ! OPTIONAL: Number of the scalar subvector
    ! where to set the variable
    integer, intent(in), optional :: iblock
!</input>

!<inputoutput>
    ! Solution vector
    type(t_vectorBlock), intent(inout) :: rvectorBlock
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: neq, nlength, jblock

    if ((rvectorBlock%nblocks .ne. 1) .and.&
        .not.present(iblock)) then
      call output_line ('Optional block number must be given for &
          &a multi-component vector', &
          OU_CLASS_ERROR,OU_MODE_STD,'transp_setVariable')
      call sys_halt()
    end if

    jblock = 1
    if (present(iblock)) jblock = iblock

    if (jblock .gt. rvectorBlock%nblocks) then
      call output_line ('Invalid block number', &
          OU_CLASS_ERROR,OU_MODE_STD,'transp_setVariable')
      call sys_halt()
    end if

    ! Set dimensions
    neq  = rvectorBlock%RvectorBlock(jblock)%NEQ

    ! Set pointer
    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(jblock), p_Ddata)

    ! Get length of values
    call ucd_getVariable(rexport, cvariable, nlength=nlength)
    if (nlength .ne. neq) then
      call output_line ('Invalid size of data', &
          OU_CLASS_ERROR,OU_MODE_STD,'transp_setVariable')
      call sys_halt()
    end if

    ! Get variable values
    call ucd_getVariable(rexport, cvariable, p_Ddata)

  end subroutine transp_setVariable

end module transport_basic
