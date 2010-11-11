!##############################################################################
!# ****************************************************************************
!# <name> mhd_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global
!# variables which are required to solve the compressible MHD equations
!#
!# The following routines are available:
!#
!# 1.) mhd_getNVAR
!#     -> Returns the number of variables depending on the spatial dimension
!#
!# 2.) mhd_getNVARtransformed
!#     -> Returns the number of variables after transformation
!#
!# 3.) mhd_getVariable
!#     -> Extracts a single variable from the vector of conservative
!#        variables stored in interleave or block format
!#
!# 4.) mhd_getVarInterleaveFormat
!#     -> Extracts a single variable from the scalar vector of
!#        conservative variables stored in interleave format
!#
!# 5.) mhd_getVarBlockFormat
!#     -> Extracts a single variable from the vector of conservative
!#        variables stored in block format
!#
!# 6.) mhd_setVariables
!#     -> Sets the conservative variables from an UCD import
!#
!# </purpose>
!##############################################################################

module mhd_basic

#include "mhd.h"

  use basicgeometry
  use boundarycondaux
  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use problem
  use statistics
  use triangulation
  use ucd

  implicit none

  private
  public :: mhd_getNVAR
  public :: mhd_getNVARtransformed
  public :: mhd_getVariable
  public :: mhd_getVarInterleaveFormat
  public :: mhd_getVarBlockFormat
  public :: mhd_setVariables

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


!<constantblock description="Global type of system coupling approach">

  ! time-dependent flow
  integer, parameter, public :: SYSTEM_SEGREGATED = 0

  ! steady-state flow
  integer, parameter, public :: SYSTEM_ALLCOUPLED = 1

!</constantblock>


!<constantblock description="Global type of dissipation">

  ! Zero dissipation
  integer, parameter, public :: DISSIPATION_ZERO           = 0

  ! Scalar dissipation
  integer, parameter, public :: DISSIPATION_SCALAR         = 1

  ! Tensorial dissipation
  integer, parameter, public :: DISSIPATION_ROE            = 2

  ! Rusanov flux
  integer, parameter, public :: DISSIPATION_RUSANOV        = 3

  ! Scalar dissipation adopting dimensional splitting
  integer, parameter, public :: DISSIPATION_SCALAR_DSPLIT  = -DISSIPATION_SCALAR

  ! Tensorial dissipation adopting dimensional splitting
  integer, parameter, public :: DISSIPATION_ROE_DSPLIT     = -DISSIPATION_ROE

  ! Rusanov flux adopting dimensional splitting
  integer, parameter, public :: DISSIPATION_RUSANOV_DSPLIT = -DISSIPATION_RUSANOV

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


!<constantblock description="Global types of perturbation parameters">

  ! Perturbation parameter is chosen as in the NITSOL package
  integer, parameter, public :: PERTURB_NITSOL  = 1

  ! Perturbation parameter is chosen as SQRT(machine precision)
  integer, parameter, public :: PERTURB_SQRTEPS = 2

!</constantblock>


!<constantblock description="Global constants for number of variables">

  ! number of solution components in 1D
  integer, parameter, public :: NVAR1D = 7

  ! number of solution components in 2D
  integer, parameter, public :: NVAR2D = 8

  ! number of solution components in 3D
  integer, parameter, public :: NVAR3D = 8

!</constantblock>


!<constantblock description="Types of boundary conditions">

  integer, parameter, public :: BDRC_DUMMYALL  = 1
  integer, parameter, public :: BDRC_DUMMYNONE = 2

!</constantblock>

!</constants>

contains

  !*****************************************************************************

!<function>

  pure function mhd_getNVAR(rproblemLevel) result(NVAR)

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
    
    NVAR = merge(7, 8, rproblemLevel%rtriangulation%ndim .eq. NDIM1D)

  end function mhd_getNVAR
  
  !*****************************************************************************

!<function>

  pure function mhd_getNVARtransformed(rproblemLevel,&
      svariables) result(NVARtransformed)

!<description>
    ! This function returns the number of variables after transformation
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! list of variable names
    character(len=*), intent(in) :: svariables
!</input>

!<result>
    ! number of transformed variables
    integer :: NVARtransformed
!</result>
!</function>

    ! Check for supported lists of variables
    if (svariables .eq. 'density,energy,momentum') then

      NVARtransformed = 5

    elseif (svariables .eq. 'density,pressure,velocity') then

      NVARtransformed = 5

    elseif (svariables .eq. 'density,pressure') then

      NVARtransformed = 2

    elseif (svariables .eq. 'density,energy') then

      NVARtransformed = 2

    elseif (svariables .eq. 'momentum') then

      NVARtransformed = 3

    elseif (svariables .eq. 'velocity') then

      NVARtransformed = 3

    elseif (svariables .eq. 'magneticfield') then
      
      NVARtransformed = merge(2, 3, rproblemLevel%rtriangulation%ndim .eq. NDIM1D)

    else

      NVARtransformed = 1

    end if

  end function mhd_getNVARtransformed

  !*****************************************************************************

!<subroutine>

  subroutine mhd_getVariable(rvectorBlock, cvariable, rvectorScalar)

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
      call mhd_getVarInterleaveFormat(neq, nvar, cvariable,&
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
      call mhd_getVarBlockFormat(neq, nvar, cvariable,&
          p_Ddata, p_Dvalue)

    end if

    ! Attach spatial discretisation from first subvector
    rvectorScalar%p_rspatialDiscr =>&
        rvectorBlock%RvectorBlock(1)%p_rspatialDiscr

  end subroutine mhd_getVariable

  !*****************************************************************************

!<subroutine>

  subroutine mhd_getVarInterleaveFormat(neq, nvar, cvariable, Ddata, Dvalue)

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
    integer :: ieq


    if (trim(cvariable) .eq. 'density') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = DENSITY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'velocity_magnitude') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = VELOCITY_MAGNITUDE_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'magneticfield_magnitude') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MAGNETICFIELD_MAGNITUDE_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MAGNETICFIELD_MAGNITUDE_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'velocity_x') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = X_VELOCITY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'velocity_y') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Y_VELOCITY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'velocity_z') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Z_VELOCITY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'momentum_x') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = X_MOMENTUM_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'momentum_y') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Y_MOMENTUM_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'momentum_z') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Z_MOMENTUM_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'magneticfield_x') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = X_MAGNETICFIELD_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'magneticfield_y') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Y_MAGNETICFIELD_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if
      
    elseif (trim(cvariable) .eq. 'magneticfield_z') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Z_MAGNETICFIELD_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'energy') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = TOTAL_ENERGY_1T_FROM_CONSVAR(Ddata, nvar, ieq)/&
                      DENSITY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do
      
    elseif (trim(cvariable) .eq. 'total_energy') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = TOTAL_ENERGY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do
        
    elseif (trim(cvariable) .eq. 'internal_energy') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = INTERNAL_ENERGY_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = INTERNAL_ENERGY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = KINETIC_ENERGY_1T_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do
      
    elseif (trim(cvariable) .eq. 'total_pressure') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = TOTAL_PRESSURE_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = TOTAL_PRESSURE_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'pressure') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = PRESSURE_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = PRESSURE_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'machnumber') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MACH_NUMBER_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MACH_NUMBER_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'speedofsound') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = SPEED_OF_SOUND_1T_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = SPEED_OF_SOUND_1T_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    else
      
      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_getVarInterleaveFormat')
      call sys_halt()

    end if

  end subroutine mhd_getVarInterleaveFormat

  !*****************************************************************************

!<subroutine>

  subroutine mhd_getVarBlockformat(neq, nvar, cvariable, Ddata, Dvalue)

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
    integer :: ieq

    
    if (trim(cvariable) .eq. 'density') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = DENSITY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'velocity_magnitude') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = VELOCITY_MAGNITUDE_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do
      
    elseif (trim(cvariable) .eq. 'magneticfield_magnitude') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MAGNETICFIELD_MAGNITUDE_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MAGNETICFIELD_MAGNITUDE_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'velocity_x') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = X_VELOCITY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'velocity_y') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Y_VELOCITY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'velocity_z') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Z_VELOCITY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'momentum_x') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = X_MOMENTUM_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'momentum_y') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Y_MOMENTUM_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'momentum_z') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = Z_MOMENTUM_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'magneticfield_x') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = X_MAGNETICFIELD_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = X_MAGNETICFIELD_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'magneticfield_y') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Y_MAGNETICFIELD_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Y_MAGNETICFIELD_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'magneticfield_z') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Z_MAGNETICFIELD_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = Z_MAGNETICFIELD_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'energy') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = TOTAL_ENERGY_1L_FROM_CONSVAR(Ddata, nvar, ieq)/&
                      DENSITY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do
      
    elseif (trim(cvariable) .eq. 'total_energy') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = TOTAL_ENERGY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'internal_energy') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = INTERNAL_ENERGY_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = INTERNAL_ENERGY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      !$omp parallel do
      do ieq = 1, neq
        Dvalue(ieq) = KINETIC_ENERGY_1L_FROM_CONSVAR(Ddata, nvar, ieq)
      end do
      !$omp end parallel do

    elseif (trim(cvariable) .eq. 'total_pressure') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = TOTAL_PRESSURE_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = TOTAL_PRESSURE_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'pressure') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = PRESSURE_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = PRESSURE_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'machnumber') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MACH_NUMBER_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = MACH_NUMBER_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    elseif (trim(cvariable) .eq. 'speedofsound') then
      if (nvar .eq. NVAR1D) then
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = SPEED_OF_SOUND_1L_FROM_CONSVAR_1D(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        do ieq = 1, neq
          Dvalue(ieq) = SPEED_OF_SOUND_1L_FROM_CONSVAR(Ddata, nvar, ieq)
        end do
        !$omp end parallel do
      end if

    else

      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_getVarBlockformat')
      call sys_halt()
      
    end if

  end subroutine mhd_getVarBlockformat

  !*****************************************************************************

!<subroutine>

  subroutine mhd_setVariables(rexport, rvectorBlock)

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
    real(DP), dimension(:), pointer :: p_Dvelocity, p_Dmagneticfield, p_Denergy
    integer :: neq, nvar, nlength, idim, ieq
    logical :: bisInterleaveFormat

    ! Check if we are in interleave of block format
    if (rvectorBlock%nblocks .eq. 1) then
      ! Set dimensions
      neq  = rvectorBlock%RvectorBlock(1)%NEQ
      nvar = rvectorBlock%RvectorBlock(1)%NVAR
      bisInterleaveFormat = .true.
    else
      ! Set dimensions
      neq  = int(rvectorBlock%NEQ/rvectorBlock%nblocks)
      nvar = rvectorBlock%nblocks
      bisInterleaveFormat = .false.
    end if
    
    ! Set pointer
    call lsysbl_getbase_double(rvectorBlock, p_Ddata)
    
    ! Get density values
    call ucd_getVariable(rexport, 'density', nlength=nlength)

    if (nlength .ne. neq) then
      call output_line ('Invalid size of data', &
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_setVariables')
      call sys_halt()
    end if
    
    ! Allocate temporal memory
    allocate(p_Ddensity(neq))
    call ucd_getVariable(rexport, 'density', p_Ddensity)
    
    ! Set density variable
    if (bisInterleaveFormat) then
      call setVarInterleaveFormat(neq, nvar, 1, p_Ddensity, p_Ddata)
    else
      call setVarBlockFormat(neq, nvar, 1, p_Ddensity, p_Ddata)
    end if
    
    ! Get velocity values
    do idim = 1, (nvar-2)/2
      select case(idim)
      case (NDIM1D)
        call ucd_getVariable(rexport, 'momentum_X', nlength=nlength)
        if (nlength .eq. -1)&
            call ucd_getVariable(rexport, 'velocity_X', nlength=nlength)
      case (NDIM2D)
        call ucd_getVariable(rexport, 'momentum_Y', nlength=nlength)
        if (nlength .eq. -1)&
            call ucd_getVariable(rexport, 'velocity_Y', nlength=nlength)
      case (NDIM3D)
        call ucd_getVariable(rexport, 'momentum_Z', nlength=nlength)
        if (nlength .eq. -1)&
            call ucd_getVariable(rexport, 'velocity_Z', nlength=nlength)
      end select
      
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_setVariables')
        call sys_halt()
      end if
      
      ! Allocate temporal memory
      allocate(p_Dvelocity(neq))
      select case(idim)
      case (NDIM1D)
        call ucd_getVariable(rexport, 'momentum_X', p_Dvelocity, nlength)
        if (nlength .eq. -1)&
            call ucd_getVariable(rexport, 'velocity_X', p_Dvelocity)
      case (NDIM2D)
        call ucd_getVariable(rexport, 'momentum_Y', p_Dvelocity, nlength)
        if (nlength .eq. -1)&
            call ucd_getVariable(rexport, 'velocity_Y', p_Dvelocity)
      case (NDIM3D)
        call ucd_getVariable(rexport, 'momentum_Z', p_Dvelocity, nlength)
        if (nlength .eq. -1)&
            call ucd_getVariable(rexport, 'velocity_Z', p_Dvelocity)
      end select
      
      ! Generate momentum (if required)
      if (nlength .eq. -1) then
        do ieq = 1, neq
          p_Dvelocity(ieq) = p_Dvelocity(ieq) * p_Ddensity(ieq)
        end do
      end if
      
      ! Set momentum variable
      if (bisInterleaveFormat) then
        call setVarInterleaveFormat(neq, nvar, idim+1, p_Dvelocity, p_Ddata)
      else
        call setVarBlockFormat(neq, nvar, idim+1, p_Dvelocity, p_Ddata)
      end if
      
      ! Deallocate temporal memory
      deallocate(p_Dvelocity)
    end do

    ! Get magnetic field values
    do idim = 1, (nvar-2)/2
      select case(idim)
      case (NDIM1D)
        call ucd_getVariable(rexport, 'magneticfield_X', nlength=nlength)
      case (NDIM2D)
        call ucd_getVariable(rexport, 'magneticfield_Y', nlength=nlength)
      case (NDIM3D)
        call ucd_getVariable(rexport, 'magneticfield_Z', nlength=nlength)
      end select
      
      if (nlength .ne. neq) then
        call output_line ('Invalid size of data', &
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_setVariables')
        call sys_halt()
      end if
      
      ! Allocate temporal memory
      allocate(p_Dmagneticfield(neq))
      select case(idim)
      case (NDIM1D)
        call ucd_getVariable(rexport, 'magneticfield_X', p_Dmagneticfield)
      case (NDIM2D)
        call ucd_getVariable(rexport, 'magneticfield_Y', p_Dmagneticfield)
      case (NDIM3D)
        call ucd_getVariable(rexport, 'magneticfield_Z', p_Dmagneticfield)
      end select

      ! Set magnetic field variable
      if (bisInterleaveFormat) then
        call setVarInterleaveFormat(neq, nvar, nvar/2+idim, p_Dmagneticfield, p_Ddata)
      else
        call setVarBlockFormat(neq, nvar, nvar/2+idim, p_Dmagneticfield, p_Ddata)
      end if
      
      ! Deallocate temporal memory
      deallocate(p_Dmagneticfield)
    end do
    
    ! Get total energy
    call ucd_getVariable(rexport, 'total_energy', nlength=nlength)
    if (nlength .eq. -1)&
        call ucd_getVariable(rexport, 'energy', nlength=nlength)

    if (nlength .ne. neq) then
      call output_line ('Invalid size of data', &
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_setVariables')
      call sys_halt()
    end if
    
    ! Allocate temporal memory
    allocate(p_Denergy(neq))
    call ucd_getVariable(rexport, 'total_energy', p_Denergy, nlength)
    if (nlength .eq. -1)&
        call ucd_getVariable(rexport, 'energy', p_Denergy)
    
    ! Generate total energy
    if (nlength .eq. 1) then
      do ieq = 1, neq
        p_Denergy(ieq) = p_Denergy(ieq) * p_Ddensity(ieq)
      end do
    end if
    
    ! Set energy variable
    if (bisInterleaveFormat) then
      call setVarInterleaveFormat(neq, nvar, nvar, p_Denergy, p_Ddata)
    else
      call setVarBlockFormat(neq, nvar, nvar, p_Denergy, p_Ddata)
    end if
      
    ! Deallocate temporal memory
    deallocate(p_Ddensity, p_Denergy)

  contains

    ! Here, the working routine follow

    !**************************************************************
    ! Set variable stored in interleave format

#ifndef USE_OPENMP
    pure &
#endif

    subroutine setVarInterleaveformat(neq, nvar, ivar, Dvalue, Ddata)

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

#ifndef USE_OPENMP
    pure &
#endif

    subroutine setVarBlockformat(neq, nvar, ivar, Dvalue, Ddata)

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

  end subroutine mhd_setVariables

end module mhd_basic
