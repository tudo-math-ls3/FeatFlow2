!##############################################################################
!# ****************************************************************************
!# <name> sse_base_sse </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for the
!# SSE problem.
!# </purpose>
!##############################################################################

module sse_base_sse

  use basicgeometry, only : NDIM2D
  use fsystem
  use genoutput
  use lineariser
  use linearsystemscalar
  use linearsystemblock
  use paramlist
  use spatialdiscretisation

  use sse_base
  
  implicit none

  private
  public :: sse_initParamSSE
  public :: sse_infoSSE
  
  public :: sse_bottomProfile
  public :: sse_bottomStress
  public :: sse_eddyViscosity
  public :: sse_calcVelocity

!<constants>

!<constantblock description="Constants for problem subtypes">

  ! SSE Alex benchmark
  integer, parameter, public :: SSE_ALEX    = 0

  ! SSE Marchi benchmark
  integer, parameter, public :: SSE_MARCHI  = 1

  ! SSE Walters benchmark
  integer, parameter, public :: SSE_WALTERS = 2

  ! SSE Winant benchmark
  integer, parameter, public :: SSE_WINANT  = 3
!</constantblock>

!<constantblock description="Constants for bottom profiles">

  ! SSE linear bottom profile
  integer, parameter, public :: SSE_BOTTOMPROFILE_LINEAR    = 0

  ! SSE constant bottom profile
  integer, parameter, public :: SSE_BOTTOMPROFILE_CONSTANT  = 1

  ! SSE parabolic bottom profile
  integer, parameter, public :: SSE_BOTTOMPROFILE_PARABOLIC = 2
!</constantblock>

!<constantblock description="Constants for stress">

  ! SSE variable stress
  integer, parameter, public :: SSE_STRESS_VARIABLE     = 0

  ! SSE constant stress
  integer, parameter, public :: SSE_STRESS_CONSTANT     = 1

  ! SSE stress proportional to bathymetry
  integer, parameter, public :: SSE_STRESS_PROPORTIONAL = 2
  !</constantblock>

!<constantblock description="Constants for eddy viscosity">

  ! SSE variable eddy viscosity
  integer, parameter, public :: SSE_VISCOSITY_VARIABLE     = 0

  ! SSE constant eddy viscosity
  integer, parameter, public :: SSE_VISCOSITY_CONSTANT     = 1

  ! SSE eddy viscosity proportional to bathymetry
  integer, parameter, public :: SSE_VISCOSITY_PROPORTIONAL = 2
!</constantblock>
  
!</constants>

!<publicvars>

  ! Length of the channel (=L)
  real(DP), public :: dlength             = 0.0_DP
  
  ! Convergent length of the channel (=Lb)
  real(DP), public :: dlengthB            = 0.0_DP

  ! Width of the entrance of the channel (=B)
  real(DP), public :: dwidth              = 0.0_DP

  ! Type of the width of the channel
  integer, public :: iwidthType           = 0

  ! Type of the bed profile
  integer, public :: ibathymetryType      = 0

  ! Mean depth of the channel (H and H0)
  real(DP), public :: dheight             = 0.0_DP
  real(DP), public :: dheight0            = 0.0_DP

  ! Depth at the end in the scaled domain (=a = dheight0/dheight)
  real(DP), public :: dheightRatio        = 0.0_DP

  ! Constant forcing at open boundary (=M2)
  real(DP), public :: dforcing            = 0.0_DP

  ! Coriolis acceleration (=f)
  real(DP), public :: dcoraccel           = 0.0_DP
  
  ! Frequency of the tidal constituent (=omega)
  real(DP), public :: dtidalfreq          = 0.0_DP

  ! Type of the bottom stress
  integer, public :: istressType          = 0

  ! Bottom stress (=s0)
  real(DP), public :: dstress             = 0.0_DP

  ! Type of the vertical eddy viscosity
  integer, public :: iviscosityType       = 0

  ! Vertical eddy viscosity (=Av0)
  real(DP), public :: dviscosity          = 0.0_DP

  ! Gravitational acceleration (=g)
  real(DP), public :: dgravaccel          = 0.0_DP

  ! Flag indicating the existence of an analytical solution
  logical, public :: bhasAnalyticSolution = .false.
  
!</publicvars>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamSSE(cproblemtype,rparlist)

!<description>
    ! This subroutine initialises the global parameters of the SSE problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
    
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sconfig,ssection

    ! Read config section
    ssection = sse_getSection(cproblemtype)
    call parlst_getvalue_string(rparlist, ssection, 'problemconfig', sconfig)
    
    ! Read parameters from parameter list (non-existing parameters are replaced by maximum value)
    call parlst_getvalue_double(rparlist,  sconfig, 'dcoraccel',      dcoraccel,   SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dforcing',       dforcing,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dgravaccel',     dgravaccel,  SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dheight',        dheight,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dheight0',       dheight0,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dheightRatio',   dheightRatio,SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dlength',        dlength,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dlengthB',       dlengthB,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dstress',        dstress,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dtidalfreq',     dtidalfreq,  SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dviscosity',     dviscosity,  SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sconfig, 'dwidth',         dwidth,      SYS_MAXREAL_DP)
    call parlst_getvalue_int(rparlist,     sconfig, 'ibathymetryType',ibathymetryType, SYS_MAXINT)
    call parlst_getvalue_int(rparlist,     sconfig, 'istressType',    istressType,     SYS_MAXINT)
    call parlst_getvalue_int(rparlist,     sconfig, 'iviscosityType', iviscosityType,  SYS_MAXINT)
    call parlst_getvalue_int(rparlist,     sconfig, 'iwidthType',     iwidthType,      SYS_MAXINT)
    call parlst_getvalue_logical(rparlist, sconfig, 'bhasAnalyticSolution', bhasAnalyticSolution, .false.)
        
  end subroutine sse_initParamSSE

  ! ***************************************************************************

!<subroutine>

  subroutine sse_infoSSE(cproblemtype,cproblemsubtype)

!<description>   
    ! This subroutine outputs information about the global parameters
    ! of the SSE problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
    
    ! Problem subtype
    integer, intent(in) :: cproblemsubtype
!</input>
!</subroutine>

    select case (cproblemtype)
    case (SSE_SCALAR)
      call output_line('PROBLEMTYPE........: SSE problem')
      call output_line('FORMULATION........: second-order equation')
    case (SSE_SYSTEM1)
      call output_line('PROBLEMTYPE........: SSE problem')
      call output_line('FORMULATION........: first-order formulation No.1')
    case (SSE_SYSTEM2)
      call output_line('PROBLEMTYPE........: SSE problem')
      call output_line('FORMULATION........: first-order formulation No.2')
    case default
      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_infoSSE")
      call sys_halt()
    end select
    
    select case (cproblemsubtype)
    case (SSE_ALEX)
      call output_line('PROBLEMSUBTYPE.....: Alex benchmark')
    case (SSE_MARCHI)
      call output_line('PROBLEMSUBTYPE.....: Marchi benchmark')
    case (SSE_WALTERS)
      call output_line('PROBLEMSUBTYPE.....: Walters benchmark')
    case (SSE_WINANT)
      call output_line('PROBLEMSUBTYPE.....: Walters benchmark')
    case default
      call output_line("Invalid problem subtype", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_infoSSE")
      call sys_halt()
    end select

    call output_line('Av0................: '//trim(adjustl(sys_sdE(dviscosity,5))))
    call output_line('B..................: '//trim(adjustl(sys_sdE(dwidth,5))))
    call output_line('H0/H...............: '//trim(adjustl(sys_sdE(dheightRatio,5))))
    call output_line('H0.................: '//trim(adjustl(sys_sdE(dheight0,5))))
    call output_line('H..................: '//trim(adjustl(sys_sdE(dheight,5))))
    call output_line('L..................: '//trim(adjustl(sys_sdE(dlength,5))))
    call output_line('Lb.................: '//trim(adjustl(sys_sdE(dlengthB,5))))
    call output_line('M2.................: '//trim(adjustl(sys_sdE(dforcing,5))))
    call output_line('f..................: '//trim(adjustl(sys_sdE(dcoraccel,5))))
    call output_line('g..................: '//trim(adjustl(sys_sdE(dgravaccel,5))))
    call output_line('omega..............: '//trim(adjustl(sys_sdE(dtidalfreq,5))))
    call output_line('s0.................: '//trim(adjustl(sys_sdE(dstress,5))))
    call output_line('bathymetryType.....: '//trim(adjustl(sys_siL(ibathymetryType,1))))
    call output_line('stressType.........: '//trim(adjustl(sys_siL(istressType,1))))
    call output_line('viscosityType......: '//trim(adjustl(sys_siL(iviscosityType,1))))
    call output_line('widthType..........: '//trim(adjustl(sys_siL(iwidthType,1))))
    call output_line('hasAnalyticSol.....: '//trim(merge('yes','no ',bhasAnalyticSolution)))

  end subroutine sse_infoSSE
  
  ! ***************************************************************************

!<function>

  elemental function sse_bottomProfile(dx,dy) result(dh)

!<description>
    ! This function computes the bathymetry as a function of the
    ! Cartisian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Height of the bottom profile
    real(DP) :: dh
!</result>

!</function>

    select case(ibathymetryType)
    case(SSE_BOTTOMPROFILE_LINEAR)
      
      ! linear bottom profile
      dh = dheight0 + (dheight-dheight0)*(1-dx/dlength)
      
    case(SSE_BOTTOMPROFILE_CONSTANT)
      
      ! constant bottom profile
      dh = dheight
      
    case (SSE_BOTTOMPROFILE_PARABOLIC)
      
      ! parabolic bottom profile
      dh = dheight*(dheightRatio+(1.0_DP-dheightRatio)*(1.0-(dy/dwidth)**2))

    case default
      dh = 0.0_DP
      
    end select
    
  end function sse_bottomProfile

  ! ***************************************************************************

!<function>

  elemental function sse_bottomStress(dx,dy) result(ds)

!<description>
    ! This function computes the bottom stress as a function of the
    ! Cartesian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Bottom stress
    real(DP) :: ds
!</result>

!</function>

    ! local parameters
    real(DP), parameter :: dx1 = -15000.0_DP
    real(DP), parameter :: dx2 = -10000.0_DP
    real(DP), parameter :: ds1 =  0.1_DP
    real(DP), parameter :: ds2 =  0.00049_DP

    real(DP), parameter :: da = (ds1-ds2)/(dx1-dx2)
    real(DP), parameter :: db = (ds2*dx1-ds1*dx2)/(dx1-dx2)

    select case(istressType)
      
    case(SSE_STRESS_VARIABLE)
      ! variable stress
      if (dx .le. dx1) then
        ds = ds1
      elseif (dx .ge. dx2) then
        ds = ds2
      elseif ((dx .gt. dx1) .and. (dx .lt. dx2)) then
        ds = da*dx+db
      else
        ds = 0.0_DP
      end if
          
    case(SSE_STRESS_CONSTANT)
      ! constant stress
      ds = dstress

    case(SSE_STRESS_PROPORTIONAL)
      ! stress proportional to bathymetry
      ds = dstress * sse_bottomProfile(dx,dy) / dheight

    case default
      ds = 0.0_DP
    end select

  end function sse_bottomStress

  ! ***************************************************************************

!<function>

  elemental function sse_eddyViscosity(dx,dy) result(dAv)

!<description>
    ! This function computes the vertical eddy viscosity as a function
    ! of the Cartesian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Vertical eddy viscosity
    real(DP) :: dAv
!</result>

!</function>

    ! local parameters
    real(DP), parameter :: dx1 = -15000.0_DP
    real(DP), parameter :: dx2 = -10000.0_DP
    real(DP), parameter :: dAv1 = 1.0_DP
    real(DP), parameter :: dAv2 = 0.012_DP

    real(DP), parameter :: da = (dAv1-dAv2)/(dx1-dx2)
    real(DP), parameter :: db = (dAv2*dx1-dAv1*dx2)/(dx1-dx2)

    select case(iviscosityType)
    case(SSE_VISCOSITY_VARIABLE)
      ! variable eddy viscosity
      if (dx .le. dx1) then
        dAv = dAv1
      elseif (dx .ge. dx2) then
        dAv = dAv2
      elseif ((dx .gt. dx1) .and. (dx .lt. dx2)) then
        dAv = da*dx+db
      else
        dAv = 0.0_DP
      end if
          
    case(SSE_VISCOSITY_CONSTANT)
      ! constant eddy viscosity
      dAv = dviscosity

    case(SSE_VISCOSITY_PROPORTIONAL)
      ! eddy viscosity proportional to bathymetry
      dAv = dviscosity * sse_bottomProfile(dx,dy) / dheight

    case default
      dAv = 0.0_DP
    end select

  end function sse_eddyViscosity

  ! ***************************************************************************

!<subroutine>

  subroutine sse_calcVelocity(rvelocity,dz,rvector,rvectorGrad_SSEre,&
      rvectorGrad_SSEim,rvectorHessX_SSEre,rvectorHessX_SSEim,&
      rvectorHessY_SSEre,rvectorHessY_SSEim)

!<description>
  ! Calculates the vertical and horizontal velocities (u,v,w) from the
  ! first and (if present) second derivatives of the sea surface
  ! elevation N provided via rvector. The discretisation structure is
  ! provided via the problem structure rproblem.
!</description>

!<input>
  ! Sea surface elevation
  type(t_vectorBlock), intent(in) :: rvector

  ! Z-coordinate, where the velocity should be calculated
  real(DP), intent(in) :: dz

  ! Gradients of the sea surface elevation
  type(t_vectorBlock), intent(in) :: rvectorGrad_SSEre
  type(t_vectorBlock), intent(in) :: rvectorGrad_SSEim

  ! OPTIONAL: second derivatives of the sea surface elevation
  type(t_vectorBlock), intent(in), optional :: rvectorHessX_SSEre
  type(t_vectorBlock), intent(in), optional :: rvectorHessX_SSEim
  type(t_vectorBlock), intent(in), optional :: rvectorHessY_SSEre
  type(t_vectorBlock), intent(in), optional :: rvectorHessY_SSEim
!</input>

!<output>
  ! Velocity vector
  type(t_vectorBlock), intent(out) :: rvelocity
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_vectorBlock) :: rcoordsDOF
    real(DP), dimension(:), pointer :: p_DcoordsX,p_DcoordsY
    real(DP), dimension(:), pointer :: p_DsseX_SSEre,p_DsseX_SSEim
    real(DP), dimension(:), pointer :: p_DsseY_SSEre,p_DsseY_SSEim
    real(DP), dimension(:), pointer :: p_DvelU_SSEre,p_DvelU_SSEim
    real(DP), dimension(:), pointer :: p_DvelV_SSEre,p_DvelV_SSEim
    complex(DP) :: calpha1,calpha2,cr1,cr2,cSSEx,cSSEY,cvelU,cvelV
    real(DP) :: dAv,dh,ds
    integer, dimension(NDIM2D) :: Isize
    integer :: ipoint,i

    ! Compute coordinates of DOFs
    Isize = rvector%RvectorBlock(1)%NEQ
    call lsysbl_createVector(rcoordsDOF, Isize, .false.)
    call lin_calcDofCoords(rvector%RvectorBlock(1)%p_rspatialDiscr, rcoordsDOF)
    call lsyssc_getbase_double(rcoordsDOF%RvectorBlock(1), p_DcoordsX)
    call lsyssc_getbase_double(rcoordsDOF%RvectorBlock(2), p_DcoordsY)

    ! Create block discretisation structure
    call spdiscr_initBlockDiscr(rblockDiscr,6,&
        rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation,&
        rvector%p_rblockDiscr%p_rboundary)
    do i=1,6
      call spdiscr_duplicateDiscrSc(rvector%RvectorBlock(1)%p_rspatialDiscr,&
          rblockDiscr%RspatialDiscr(i), .true.)
    end do
    
    ! Set pointer to horizontal velocity values
    call lsysbl_createVector(rblockDiscr, rvelocity, .true.)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(1), p_DvelU_SSEre)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(2), p_DvelU_SSEim)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(3), p_DvelV_SSEre)
    call lsyssc_getbase_double(rvelocity%RvectorBlock(4), p_DvelV_SSEim)

    ! Set pointers to gradient values
    call lsyssc_getbase_double(rvectorGrad_SSEre%RvectorBlock(1), p_DsseX_SSEre)
    call lsyssc_getbase_double(rvectorGrad_SSEre%RvectorBlock(2), p_DsseY_SSEre)
    call lsyssc_getbase_double(rvectorGrad_SSEim%RvectorBlock(1), p_DsseX_SSEim)
    call lsyssc_getbase_double(rvectorGrad_SSEim%RvectorBlock(2), p_DsseY_SSEim)

    ! Compute horizontal velocities (U,V) from analytical expressions
    do ipoint = 1, size(p_DcoordsX)

      ! Compute bottom profile
      dh = sse_bottomProfile(p_DcoordsX(ipoint),p_DcoordsY(ipoint))
      
      ! Compute bottom stress
      ds = sse_bottomStress(p_DcoordsX(ipoint),p_DcoordsY(ipoint))
      
      ! Compute vertical eddy viscosity
      dAv = sse_eddyViscosity(p_DcoordsX(ipoint),p_DcoordsY(ipoint))

      ! Compute coefficients calpha1 and calpha2
      calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
      calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)
      
      ! Compute complex sea surface elevation
      cSSEx = cmplx(p_DsseX_SSEre(ipoint),p_DsseX_SSEim(ipoint))
      cSSEy = cmplx(p_DsseY_SSEre(ipoint),p_DsseY_SSEim(ipoint))

      ! Compute coefficient cr1
      cr1 = dgravaccel/(dAv*(calpha1**2))*&
          ((ds*cosh(calpha1*dz))/(calpha1*dAv*sinh(calpha1*dh)+&
                                           ds*cosh(calpha1*dh))-1.0_DP)*&
                                           (cSSEx+cimg*cSSEy)

      ! Compute coefficient cr2
      cr2 = dgravaccel/(dAv*(calpha2**2))*&
          ((ds*cosh(calpha2*dz))/(calpha2*dAv*sinh(calpha2*dh)+&
                                           ds*cosh(calpha2*dh))-1.0_DP)*&
                                           (cSSEx-cimg*cSSEy)

      ! Compute complex velocities
      cvelU = 0.5_DP*(cr1+cr2)
      cvelV = 0.5_DP*(cr1-cr2)/cimg
      
      ! And separate them into real and imaginary parts
      p_DvelU_SSEre(ipoint) = real(cvelU)
      p_DvelV_SSEre(ipoint) = real(cvelV)
      p_DvelU_SSEim(ipoint) = aimag(cvelU)
      p_DvelV_SSEim(ipoint) = aimag(cvelV)
    end do

    ! Release DOF coordinates
    call lsysbl_releaseVector(rcoordsDOF)
    
  end subroutine sse_calcVelocity
    
end module sse_base_sse
