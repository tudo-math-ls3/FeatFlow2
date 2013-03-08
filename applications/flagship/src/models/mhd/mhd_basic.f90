!##############################################################################
!# ****************************************************************************
!# <name> mhd_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global
!# variables which are required to solve the compressible ideal MHD equations
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
!# 4.) mhd_setVariables
!#     -> Sets the conservative variables from an UCD import
!#
!# 5.) mhd_convertVariableNondim = mhd_convertVariableNondim1 /
!#                                 mhd_convertVariableNondimDP /
!#                                 mhd_convertVariableNondimSP /
!#     -> Non-dimensionalizes a physical variable using reference values
!#
!# 6.) mhd_convertVariableDim = mhd_convertVariableDim1 /
!#                              mhd_convertVariableDimDP /
!#                              mhd_convertVariableDimSP
!#     -> Convert a non-dimensionalized variable into its physical quantity
!#
!# </purpose>
!##############################################################################

module mhd_basic

#include "../../flagship.h"
#include "mhd.h"

!$use omp_lib
  use basicgeometry
  use boundarycondaux
  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use problem
  use statistics
  use storage
  use triangulation
  use ucd

  ! Modules from MHD model
  use mhd_basic1d
  use mhd_basic2d
  use mhd_basic3d

  implicit none

  private
  public :: mhd_getNVAR
  public :: mhd_getNVARtransformed
  public :: mhd_getVariable
  public :: mhd_setVariables
  public :: mhd_convertVariableNonDim
  public :: mhd_convertVariableDim

  interface mhd_convertVariableNonDim
    module procedure mhd_convertVariableNonDim1
    module procedure mhd_convertVariableNonDimDP
    module procedure mhd_convertVariableNonDimSP
  end interface mhd_convertVariableNonDim

  interface mhd_convertVariableDim
    module procedure mhd_convertVariableDim1
    module procedure mhd_convertVariableDimDP
    module procedure mhd_convertVariableDimSP
  end interface mhd_convertVariableDim

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


!<constantblock description="Types of boundary conditions">

  ! Supersonic outlet boundary condition
  ! No boundary conditions are prescribed at all
  
  integer, parameter, public :: BDRC_SUPEROUTLET  = 11

!</constantblock>


!<constantblock description="Global types of UCD export">

  ! Standard UCD export
  integer, parameter, public :: UCDEXPORT_STD             = 0

  ! UCD export on discontinuous P1 finite elements
  integer, parameter, public :: UCDEXPORT_P1DISCONTINUOUS = 1

  ! UCD export on continuous P1 finite elements
  integer, parameter, public :: UCDEXPORT_P1CONTINUOUS    = 2

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

  subroutine mhd_getVariable(rvectorBlock, cvariable, rvectorScalar, Imask)

!<description>
    ! This subroutine extracts a single variable from the vector of
    ! conservative variables which is stored in interleave of block format
!</description>

!<input>
    ! Vector of conservative variables
    type(t_vectorBlock), intent(in) :: rvectorBlock

    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! OPTIONAL: integer mask array
    ! If present only those entries of the destination vector are
    ! computed which are given by the integer mask.
    integer, dimension(:), intent(in), optional :: Imask
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
      select case(rvectorBlock%p_rblockDiscr%ndimension)
      case(NDIM1D)
        call mhd_getVarInterleaveFormat1d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM2D)
        call mhd_getVarInterleaveFormat2d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM3D)
        call mhd_getVarInterleaveFormat3d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      end select
      
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
      select case(rvectorBlock%p_rblockDiscr%ndimension)
      case(NDIM1D)
        call mhd_getVarBlockFormat1d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM2D)
        call mhd_getVarBlockFormat2d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM3D)
        call mhd_getVarBlockFormat3d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      end select
      
    end if

    ! Attach spatial discretisation from first subvector
    rvectorScalar%p_rspatialDiscr =>&
        rvectorBlock%RvectorBlock(1)%p_rspatialDiscr

  end subroutine mhd_getVariable

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
    if (nlength .eq. -1) then
      call ucd_getVariable(rexport, 'energy', p_Denergy)
    
      ! Generate total energy
      if (nlength .eq. 1) then
        do ieq = 1, neq
          p_Denergy(ieq) = p_Denergy(ieq) * p_Ddensity(ieq)
        end do
      end if
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

    pure subroutine setVarInterleaveformat(neq, nvar, ivar, Dvalue, Ddata)

      integer, intent(in) :: neq, nvar, ivar
      real(DP), dimension(:), intent(in) :: Dvalue

      real(DP), dimension(nvar,neq), intent(inout) :: Ddata

      ! local variables
      integer :: ieq

      do ieq = 1, neq
        Ddata(ivar, ieq) = Dvalue(ieq)
      end do

    end subroutine setVarInterleaveformat

    !**************************************************************
    ! Set variable stored in block format

    pure subroutine setVarBlockformat(neq, nvar, ivar, Dvalue, Ddata)

      integer, intent(in) :: neq, nvar, ivar
      real(DP), dimension(:), intent(in) :: Dvalue

      real(DP), dimension(neq,nvar), intent(inout) :: Ddata

      ! local variables
      integer :: ieq

      do ieq = 1, neq
        Ddata(ieq, ivar) = Dvalue(ieq)
      end do

    end subroutine setVarBlockformat

  end subroutine mhd_setVariables

  !*****************************************************************************

!<subroutine>

  subroutine mhd_convertVariableNonDim1(rvectorScalar, cvariable,&
                                          density_ref, velocity_ref, length_ref)

!<description>
    ! This subroutine non-dimensionalises the physical variable given
    ! in rvectorScalar based on the reference values density_ref,
    ! velocity_ref and length_ref and overwrites rvectorScalar.
!</description>

!<input>
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Reference values
    real(DP), intent(in) :: density_ref,velocity_ref,length_ref
!</input>

!<inputoutput>
    ! On input: scalar physical variable
    ! On output: non-dimensionalised variable
    type(t_vectorScalar), intent(inout) :: rvectorScalar
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    if (rvectorScalar%NEQ .eq. 0) then
      call output_line ('Cannot convert empty vector', &
          OU_CLASS_WARNING,OU_MODE_STD,'mhd_convertVariableNonDim1')
      return
    end if
    
    select case(rvectorScalar%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double(rvectorScalar, p_Ddata)
      call mhd_convertVariableNonDimDP(p_Ddata, cvariable,&
          density_ref, velocity_ref, length_ref)

    case (ST_SINGLE)
      call lsyssc_getbase_single(rvectorScalar, p_Fdata)
      call mhd_convertVariableNonDimSP(p_Fdata, cvariable,&
          real(density_ref,SP), real(velocity_ref,SP), real(length_ref,SP))

    case default
      call output_line ('Unspported data type', &
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_convertVariableNonDim1')
      call sys_halt()
    end select

  end subroutine mhd_convertVariableNonDim1

  !*****************************************************************************

!<subroutine>

  subroutine mhd_convertVariableNonDimDP(Ddata, cvariable,&
                                             density_ref, velocity_ref, length_ref)

!<description>
    ! This subroutine non-dimensionalises the physical variable given
    ! in the double-valued array Ddata based on the reference values
    ! density_ref, velocity_ref and length_ref and overwrites Ddata.
!</description>

!<input>
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Reference values
    real(DP), intent(in) :: density_ref,velocity_ref,length_ref
!</input>

!<inputoutput>
    ! On input: scalar physical variable
    ! On output: non-dimensionalised variable
    real(DP), dimension(:), intent(inout) :: Ddata
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ieq

    if (trim(cvariable) .eq. 'density') then
      ! density = density_phys / density_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/density_ref
      end do

    elseif (trim(cvariable) .eq. 'velocity') then
      ! velocity = velocity_phys / velocity_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/velocity_ref
      end do
      
    elseif (trim(cvariable) .eq. 'momentum') then
      ! momentum = momentum_phys / (density_ref*velocity_ref)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/(density_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'energy') then
      ! energy = energy_phys / (density_ref*velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'total_energy') then
      ! energy = energy_phys / (velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'internal_energy') then
      ! energy = energy_phys / (velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      ! energy = energy_phys / (velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/(velocity_ref*velocity_ref)
      end do
      
    elseif (trim(cvariable) .eq. 'pressure') then
      ! pressure = pressure_phys / (density_ref*velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'machnumber') then
      ! No conversion needed

    elseif (trim(cvariable) .eq. 'speedofsound') then
      ! soundspeed = soundspeed_phys / velocity_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/velocity_ref
      end do

    elseif (trim(cvariable) .eq. 'coordinate') then
      ! coordinate = coordinate_phys / length_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)/length_ref
      end do

    elseif (trim(cvariable) .eq. 'time') then
      ! time = velocity_ref * time_phys / length_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = velocity_ref*Ddata(ieq)/length_ref
      end do

    else
      
      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_convertVariableNonDimDP')
      call sys_halt()
      
    end if

  end subroutine mhd_convertVariableNonDimDP

  !*****************************************************************************

!<subroutine>

  subroutine mhd_convertVariableNonDimSP(Fdata, cvariable,&
                                             density_ref, velocity_ref, length_ref)

!<description>
    ! This subroutine non-dimensionalises the physical variable given
    ! in the single-valued array Fdata based on the reference values
    ! density_ref, velocity_ref and length_ref and overwrites Fdata.
!</description>

!<input>
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Reference values
    real(SP), intent(in) :: density_ref,velocity_ref,length_ref
!</input>

!<inputoutput>
    ! On input: scalar physical variable
    ! On output: non-dimensionalised variable
    real(SP), dimension(:), intent(inout) :: Fdata
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ieq

    if (trim(cvariable) .eq. 'density') then
      ! density = density_phys / density_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/density_ref
      end do

    elseif (trim(cvariable) .eq. 'velocity') then
      ! velocity = velocity_phys / velocity_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/velocity_ref
      end do
      
    elseif (trim(cvariable) .eq. 'momentum') then
      ! momentum = momentum_phys / (density_ref*velocity_ref)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/(density_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'energy') then
      ! energy = energy_phys / (density_ref*velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'total_energy') then
      ! energy = energy_phys / (velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'internal_energy') then
      ! energy = energy_phys / (velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      ! energy = energy_phys / (velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/(velocity_ref*velocity_ref)
      end do
      
    elseif (trim(cvariable) .eq. 'pressure') then
      ! pressure = pressure_phys / (density_ref*velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'machnumber') then
      ! No conversion needed

    elseif (trim(cvariable) .eq. 'speedofsound') then
      ! soundspeed = soundspeed_phys / velocity_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/velocity_ref
      end do

    elseif (trim(cvariable) .eq. 'coordinate') then
      ! coordinate = coordinate_phys / length_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)/length_ref
      end do

    elseif (trim(cvariable) .eq. 'time') then
      ! time = velocity_ref * time_phys / length_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = velocity_ref*Fdata(ieq)/length_ref
      end do

    else
      
      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_convertVariableNonDimSP')
      call sys_halt()
      
    end if

  end subroutine mhd_convertVariableNonDimSP

  !*****************************************************************************

!<subroutine>

  subroutine mhd_convertVariableDim1(rvectorScalar, cvariable,&
                                       density_ref, velocity_ref, length_ref)

!<description>
    ! This subroutine converts a non-dimensional variable given in
    ! rvectorScalar into physical quantities based on the reference
    ! values density_ref, velocity_ref and length_ref and
    ! overwrites rvectorScalar.
!</description>

!<input>
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Reference values
    real(DP), intent(in) :: density_ref,velocity_ref,length_ref
!</input>

!<inputoutput>
    ! On input: non-dimensionalised variable
    ! On output: scalar physical variable
    type(t_vectorScalar), intent(inout) :: rvectorScalar
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    if (rvectorScalar%NEQ .eq. 0) then
      call output_line ('Cannot convert empty vector', &
          OU_CLASS_WARNING,OU_MODE_STD,'mhd_convertVariableDim1')
      return
    end if

    select case(rvectorScalar%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double(rvectorScalar, p_Ddata)
      call mhd_convertVariableDimDP(p_Ddata, cvariable,&
          density_ref, velocity_ref, length_ref)

    case (ST_SINGLE)
      call lsyssc_getbase_single(rvectorScalar, p_Fdata)
      call mhd_convertVariableDimSP(p_Fdata, cvariable,&
          real(density_ref,SP), real(velocity_ref,SP), real(length_ref,SP))

    case default
      call output_line ('Unspported data type', &
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_convertVariableDim1')
      call sys_halt()
    end select

  end subroutine mhd_convertVariableDim1

  !*****************************************************************************

!<subroutine>

  subroutine mhd_convertVariableDimDP(Ddata, cvariable,&
                                          density_ref, velocity_ref, length_ref)

!<description>
    ! This subroutine converts a non-dimensional variable given in the
    ! double-valued array Ddata into physical quantities based on the
    ! reference values density_ref, velocity_ref and length_ref and
    ! overwrites Ddata.
!</description>

!<input>
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Reference values
    real(DP), intent(in) :: density_ref,velocity_ref,length_ref
!</input>

!<inputoutput>
    ! On input: non-dimensionalised variable
    ! On output: scalar physical variable
    real(DP), dimension(:), intent(inout) :: Ddata
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ieq

    if (trim(cvariable) .eq. 'density') then
      ! density_phys = density * density_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*density_ref
      end do

    elseif (trim(cvariable) .eq. 'velocity') then
      ! velocity_phys = velocity * velocity_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*velocity_ref
      end do
      
    elseif (trim(cvariable) .eq. 'momentum') then
      ! momentum_phys = momentum * (density_ref*velocity_ref)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*(density_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'energy') then
      ! energy_phys = energy * (density_ref*velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'total_energy') then
      ! energy_phys = energy * (velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'internal_energy') then
      ! energy_phys = energy * (velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      ! energy_phys = energy * (velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*(velocity_ref*velocity_ref)
      end do
      
    elseif (trim(cvariable) .eq. 'pressure') then
      ! pressure_phys = pressure * (density_ref*velocity_ref^2)
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'machnumber') then
      ! No conversion needed

    elseif (trim(cvariable) .eq. 'speedofsound') then
      ! soundspeed_phys = soundspeed * velocity_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*velocity_ref
      end do

    elseif (trim(cvariable) .eq. 'coordinate') then
      ! coordinate_phys = coordinate * length_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = Ddata(ieq)*length_ref
      end do

    elseif (trim(cvariable) .eq. 'time') then
      ! time_phys = length_ref * time / velocity_ref
      do ieq=1,size(Ddata)
        Ddata(ieq) = length_ref*Ddata(ieq)/velocity_ref
      end do

    else
      
      print *, trim(cvariable)

      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_convertVariableDimDP')
      call sys_halt()
      
    end if
      
  end subroutine mhd_convertVariableDimDP

  !*****************************************************************************

!<subroutine>

  subroutine mhd_convertVariableDimSP(Fdata, cvariable,&
                                          density_ref, velocity_ref, length_ref)

!<description>
    ! This subroutine converts a non-dimensional variable given in the
    ! single-valued array Fdata into physical quantities based on the
    ! reference values density_ref, velocity_ref and length_ref and
    ! overwrites Fdata.
!</description>

!<input>
    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Reference values
    real(SP), intent(in) :: density_ref,velocity_ref,length_ref
!</input>

!<inputoutput>
    ! On input: non-dimensionalised variable
    ! On output: scalar physical variable
    real(SP), dimension(:), intent(inout) :: Fdata
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ieq

    if (trim(cvariable) .eq. 'density') then
      ! density_phys = density * density_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*density_ref
      end do

    elseif (trim(cvariable) .eq. 'velocity') then
      ! velocity_phys = velocity * velocity_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*velocity_ref
      end do
      
    elseif (trim(cvariable) .eq. 'momentum') then
      ! momentum_phys = momentum * (density_ref*velocity_ref)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*(density_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'energy') then
      ! energy_phys = energy * (density_ref*velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'total_energy') then
      ! energy_phys = energy * (velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'internal_energy') then
      ! energy_phys = energy * (velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*(velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      ! energy_phys = energy * (velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*(velocity_ref*velocity_ref)
      end do
      
    elseif (trim(cvariable) .eq. 'pressure') then
      ! pressure_phys = pressure * (density_ref*velocity_ref^2)
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*(density_ref*velocity_ref*velocity_ref)
      end do

    elseif (trim(cvariable) .eq. 'machnumber') then
      ! No conversion needed

    elseif (trim(cvariable) .eq. 'speedofsound') then
      ! soundspeed_phys = soundspeed * velocity_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*velocity_ref
      end do

    elseif (trim(cvariable) .eq. 'coordinate') then
      ! coordinate_phys = coordinate * length_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = Fdata(ieq)*length_ref
      end do

    elseif (trim(cvariable) .eq. 'time') then
      ! time_phys = length_ref * time / velocity_ref
      do ieq=1,size(Fdata)
        Fdata(ieq) = length_ref*Fdata(ieq)/velocity_ref
      end do

    else
      
      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_convertVariableDimSP')
      call sys_halt()
      
    end if
      
  end subroutine mhd_convertVariableDimSP

end module mhd_basic
