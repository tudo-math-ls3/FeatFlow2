!##############################################################################
!# ****************************************************************************
!# <name> hydro_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global variables
!# which are required to solve the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!# 1.) hydro_getNVAR
!#     -> Returns the number of variables depending on the spatial dimension
!#
!# 2.) hydro_getNVARtransformed
!#     -> Returns the number of variables after transformation
!#
!# 3.) hydro_getVariable
!#     -> Extracts a single variable from the vector of conservative
!#        variables stored in interleave or block format
!#
!# 4.) hydro_setVariables
!#     -> Sets the conservative variables from an UCD import
!#
!# </purpose>
!##############################################################################

module hydro_basic

#include "hydro.h"

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
  use triangulation
  use ucd

  ! Modules from hydrodynamic model
  use hydro_basic1d
  use hydro_basic2d
  use hydro_basic3d
  
  implicit none

  private
  public :: hydro_getNVAR
  public :: hydro_getNVARtransformed
  public :: hydro_getVariable
  public :: hydro_setVariables

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

  ! Free-slip and symmetry plane boundary condition
  ! The normal component of the velocity vector is set to zero
  !
  ! V*n = 0

  integer, parameter, public :: BDRC_FREESLIP    = 1

  ! Relaxed free-slip and symmetry plane boundary condition
  ! The normal component of the velocity vector is "approching" zero
  !
  ! (V-c*Vold)*n = 0, where   0 < c <= 1

  integer, parameter, public :: BDRC_RLXFREESLIP = 2

  ! Viscous wall boundary condition
  ! The velocity vector is set to zero
  !
  ! V = 0

  integer, parameter, public :: BDRC_VISCOUSWALL  = 3

  ! Free stream boundary condition using characteristics
  ! These boundary conditions can be used for both subsonic and
  ! supersonic in- and outflow where the characteristics are either
  ! set from free stream quantities for ingoing characteristics or
  ! adopted from the interior values for outgoing characteristics.

  integer, parameter, public :: BDRC_FREESTREAM   = 4

  ! Subsonic inlet boundary condition
  ! At a subsonic inlet, the recommended boundary condition is to specify
  ! the total temperature and total pressure as well as the flow angle.

  integer, parameter, public :: BDRC_SUBINLET     = 5

  ! Subsonic outlet boundary condition
  ! At a subsonic outlet the recommended boundary condition is to specify
  ! the static pressure.

  integer, parameter, public :: BDRC_SUBOUTLET    = 6

  ! Massflow inlet boundary condition
  ! This boundary condition can be prescribed at a subsonic inflow boundary
  ! which requires a given value of mass flow. It is in principal identical
  ! to the total states inflow boundary condition, where the value of total
  ! pressure is a function of mass flow. First, the velocity normal to the
  ! boundary is found from the value of desired mass flow and density
  !
  ! (1) V_n = Massflow/rho
  !
  ! Next, the Lavel number is defined as
  !
  ! (2) Laval = SQRT(U_n^2+V_n^2)/c
  !
  ! Finally, the total pressure is calculated as
  !
  ! (3) p0 = p*(1-(gamma-1)/(gamma+1)Laval^2)^(-gamma)/(gamma-1)

  integer, parameter, public :: BDRC_MASSINLET    = 7

  ! Massflow outlet boundary condition
  ! This boundary condition can be prescribed at a subsonic outflow boundary
  ! which requires a given value of mass flow. It is in principal identical
  ! to the static pressure inflow boundary condition where the value of
  ! static pressure is a function of mass flow. For each value of total
  ! state and given mass flow, there is a corresponding value of static pressure.
  ! Mathematically, this can be achieved by solving the implicit equation
  ! for the velocity V_n:
  !
  ! (1) Massflow = p0*V_n*(1-(gamma-1)/(2*c0^2)*V_n^2)^1/(gamma-1)*S,
  !
  ! where the total states are known from the solution. Next, the Laval number
  ! is computed from the velocity value as follows
  !
  ! (2) Laval = V_n/(c0*SQRT(2/(gamma+1))),
  !
  ! where c0 is the total speed of sound. Finally, the static pressure is
  ! computed from the Laval number as follows
  !
  ! (3) p = p0*(1-(gamma-1)/(gamma+1)*Laval^2)^gamma/(gamma-1)

  integer, parameter, public :: BDRC_MASSOUTLET   = 8

  ! Mach outflow boundary condition
  ! This condition is similar to the mass flow outflow boundary condition. It is
  ! basically the static pressure outflow boundary condition, where the value of
  ! static pressure is expressed as function of desired Mach number M and known
  ! value of total pressure p0
  !
  ! p = p0*(1+(gamma-1)/2*M^2)^(-gamma/(gamma-1))

  integer, parameter, public :: BDRC_MACHOUTLET   = 9

  ! Supersonic inlet boundary condition
  ! All boundary conditions are prescribed by the free stream quantities

  integer, parameter, public :: BDRC_SUPERINLET   = 10

  ! Supersonic outlet boundary condition
  ! No boundary conditions are prescribed at all

  integer, parameter, public :: BDRC_SUPEROUTLET  = 11

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

  pure function hydro_getNVAR(rproblemLevel) result(NVAR)

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

  end function hydro_getNVAR
  
  !*****************************************************************************

!<function>

  pure function hydro_getNVARtransformed(rproblemLevel,&
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

      NVARtransformed = hydro_getNVAR(rproblemLevel)

    elseif (svariables .eq. 'density,pressure,velocity') then

      NVARtransformed = hydro_getNVAR(rproblemLevel)

    elseif (svariables .eq. 'density,pressure') then

      NVARtransformed = 2

    elseif (svariables .eq. 'density,energy') then

      NVARtransformed = 2

    elseif (svariables .eq. 'momentum') then

      NVARtransformed = rproblemLevel%rtriangulation%ndim

    elseif (svariables .eq. 'velocity') then

      NVARtransformed = rproblemLevel%rtriangulation%ndim

    else

      NVARtransformed = 1

    end if

  end function hydro_getNVARtransformed

  !*****************************************************************************

!<subroutine>

  subroutine hydro_getVariable(rvectorBlock, cvariable, rvectorScalar, Imask)

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
        call hydro_getVarInterleaveFormat1d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM2D)
        call hydro_getVarInterleaveFormat2d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM3D)
        call hydro_getVarInterleaveFormat3d(neq, nvar, cvariable,&
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
        call hydro_getVarBlockFormat1d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM2D)
        call hydro_getVarBlockFormat2d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      case(NDIM3D)
        call hydro_getVarBlockFormat3d(neq, nvar, cvariable,&
            p_Ddata, p_Dvalue, Imask)
      end select
      
    end if

    ! Attach spatial discretisation from first subvector
    rvectorScalar%p_rspatialDiscr =>&
        rvectorBlock%RvectorBlock(1)%p_rspatialDiscr

  end subroutine hydro_getVariable

  !*****************************************************************************

!<subroutine>

  subroutine hydro_setVariables(rexport, rvectorBlock)

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
    real(DP), dimension(:), pointer :: p_Dvelocity, p_Denergy
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
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_setVariables')
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
    do idim = 1, nvar-2
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
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_setVariables')
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
       
    ! Get total energy
    call ucd_getVariable(rexport, 'total_energy', nlength=nlength)
    if (nlength .eq. -1)&
        call ucd_getVariable(rexport, 'energy', nlength=nlength)

    if (nlength .ne. neq) then
      call output_line ('Invalid size of data', &
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_setVariables')
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

  end subroutine hydro_setVariables

end module hydro_basic
