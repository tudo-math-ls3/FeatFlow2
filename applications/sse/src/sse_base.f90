!##############################################################################
!# ****************************************************************************
!# <name> sse_base </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for solving
!# the elliptic equation for sea surface elevation.
!# </purpose>
!##############################################################################

module sse_base

  use basicgeometry
  use fsystem
  use lineariser
  use linearsystemblock
  use linearsystemscalar
  use spatialdiscretisation
!  use storage

  implicit none

  private
  public :: sse_bottomProfile
  public :: sse_bottomStress
  public :: sse_eddyViscosity
  public :: sse_calcVelocity

#ifdef USE_COMPILER_INTEL
  public :: sinh
  public :: cosh
#endif

!<constants>

!<constantblock description="Constants for problem type">

  ! Comput standard Poisson problem
  integer, parameter, public :: POISSON_SCALAR = 0

  ! Comput standard Poisson problem as first-order system
  integer, parameter, public :: POISSON_SYSTEM = 1

  ! Compute SSE solution from scalar problem
  integer, parameter, public :: SSE_SCALAR     = 2

  ! Compute SSE solution from first-order system ($\sigma=A\nabla N$)
  integer, parameter, public :: SSE_SYSTEM1    = 3

  ! Compute SSE solution from first-order system ($\sigma=\nabla N$)
  integer, parameter, public :: SSE_SYSTEM2    = 4

!</constantblock>

!<constantblock description="Constants for complex numbers">

  ! Real part of complex number
  complex(DP), parameter, public :: creal = cmplx(1.0_DP,0.0_DP)

  ! Imaginary part of complex number
  complex(DP), parameter, public :: cimg = cmplx(0.0_DP,1.0_DP)

!</constantblock>

#if defined(CASE_POISSON_DIRICHLET)
!<constantblock description="Problem parameters for pure Dirichlet Poisson benchmark">

  ! Scaling parameter
  real(DP), parameter, public :: dpoisson       = 1.0_DP

  ! Dirichlet boundary value (=u_D)
  real(DP), parameter, public :: ddirichlet     = 3.0_DP

#elif defined(CASE_POISSON_NEUMANN)

!<constantblock description="Problem parameters for mixed Dirichlet-Neumann Poisson benchmark">

  ! Scaling parameter
  real(DP), parameter, public :: dpoisson       = 1.0_DP

  ! Dirichlet boundary value (=u_D)
  real(DP), parameter, public :: ddirichlet     = 3.0_DP

  ! Neumann boundary value (=g)
  real(DP), parameter, public :: dneumann       = 0.0_DP

#else
#error 'Test case is undefined.'
#endif


#if defined(CASE_SSE_ALEX)
!<constantblock description="Problem parameters for Alex benchmark">

  ! Length of the channel (=L)
  real(DP), parameter, public :: dlength        = 6.0E4_DP
  
  ! Convergent length of the channel (=Lb)
  real(DP), parameter, public :: dlengthB       = 1.0e20_DP

  ! Width of the entrance of the channel (=B)
  real(DP), parameter, public :: dwidth         = 15000.0_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType      = 2

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = 0

  ! Mean depth of the channel (H and H0)
  real(DP), parameter, public :: dheight        = 7.7_DP
  real(DP), parameter, public :: dheight0       = 7.7_DP

  ! Depth at the end in the scaled domain (=a = dheight0/dheight)
  real(DP), parameter, public :: dheightRatio   = SYS_MAXREAL_DP

  ! Constant forcing at open boundary (=M2)
  real(DP), parameter, public :: dforcing       = 1.0_DP

  ! Coriolis acceleration (=f)
  real(DP), parameter, public :: dcoraccel      = 0.0_DP
  
  ! Frequency of the tidal constituent (=omega)
  real(DP), parameter, public :: dtidalfreq     = 1.405634296908185E-4_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType     = 1

  ! Bottom stress (=s0)
  real(DP), parameter, public :: dstress        = 0.0000098_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType  = 1

  ! Vertical eddy viscosity (=Av0)
  real(DP), parameter, public :: dviscosity     = 0.019_DP

  ! Flag indicating the existence of an analytical solution
  logical, parameter, public :: bhasAnalyticSolution = .true.

!</constantblock>
#elif defined(CASE_SSE_MARCHI)
!<constantblock description="Problem parameters for Marchi benchmark">

  ! Length of the channel (=L)
  real(DP), parameter, public :: dlength        = 1.0_DP
  
  ! Convergent length of the channel (=Lb)
  real(DP), parameter, public :: dlengthB       = SYS_MAXREAL_DP

  ! Width of the entrance of the channel (=B)
  real(DP), parameter, public :: dwidth         = 1.0_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType      = 1

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = SYS_MAXINT

  ! Mean depth of the channel (H and H0)
  real(DP), parameter, public :: dheight        = SYS_MAXREAL_DP
  real(DP), parameter, public :: dheight0       = SYS_MAXREAL_DP

  ! Depth at the end in the scaled domain (=a = dheight0/dheight)
  real(DP), parameter, public :: dheightRatio   = SYS_MAXREAL_DP

  ! Constant forcing at open boundary (=M2)
  real(DP), parameter, public :: dforcing       = SYS_MAXREAL_DP

  ! Coriolis acceleration (=f)
  real(DP), parameter, public :: dcoraccel      = SYS_MAXREAL_DP
  
  ! Frequency of the tidal constituent (=omega)
  real(DP), parameter, public :: dtidalfreq     = SYS_MAXREAL_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType     = SYS_MAXINT

  ! Bottom stress (=s0)
  real(DP), parameter, public :: dstress        = SYS_MAXREAL_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType  = SYS_MAXINT

  ! Vertical eddy viscosity (=Av0)
  real(DP), parameter, public :: dviscosity     = SYS_MAXREAL_DP

  ! Flag indicating the existence of an analytical solution
  logical, parameter, public :: bhasAnalyticSolution = .false.

!</constantblock>
#elif defined(CASE_SSE_WALTERS)
!<constantblock description="Problem parameters for Walters benchmark">

  ! Length of the channel (=L)
  real(DP), parameter, public :: dlength        = 4.0E3_DP
  
  ! Convergent length of the channel (=Lb)
  real(DP), parameter, public :: dlengthB       = SYS_MAXREAL_DP

  ! Width of the entrance of the channel (=B)
  real(DP), parameter, public :: dwidth         = 2.0E3_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType      = 1

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = 1

  ! Mean depth of the channel (H and H0)
  real(DP), parameter, public :: dheight        = 10.0_DP
  real(DP), parameter, public :: dheight0       = SYS_MAXREAL_DP

  ! Depth at the end in the scaled domain (=a = dheight0/dheight)
  real(DP), parameter, public :: dheightRatio   = SYS_MAXREAL_DP

  ! Constant forcing at open boundary (=M2)
  real(DP), parameter, public :: dforcing       = 0.1_DP

  ! Coriolis acceleration (=f)
  real(DP), parameter, public :: dcoraccel      = 0.0_DP
  
  ! Frequency of the tidal constituent (=omega)
  real(DP), parameter, public :: dtidalfreq     = 0.00174532925199433_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType     = 1

  ! Bottom stress (=s0)
  real(DP), parameter, public :: dstress        = 0.0_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType  = 1

  ! Vertical eddy viscosity (=Av0)
  real(DP), parameter, public :: dviscosity     = 0.012_DP

  ! Flag indicating the existence of an analytical solution
  logical, parameter, public :: bhasAnalyticSolution = .false.

!</constantblock>
#elif defined(CASE_SSE_WINANT)
!<constantblock description="Problem parameters for Winant benchmark">

  ! Length of the channel (=L)
  real(DP), parameter, public :: dlength        = 6.0E4_DP
  
  ! Convergent length of the channel (=Lb)
  real(DP), parameter, public :: dlengthB       = SYS_MAXREAL_DP

  ! Width of the entrance of the channel (=B)
  real(DP), parameter, public :: dwidth         = 1000.0_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType      = 1

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = 2

  ! Mean depth of the channel (H and H0)
  real(DP), parameter, public :: dheight        = 10.0_DP
  real(DP), parameter, public :: dheight0       = SYS_MAXREAL_DP

  ! Depth at the end in the scaled domain (=a = dheight0/dheight)
  real(DP), parameter, public :: dheightRatio   = 0.01_DP

  ! Constant forcing at open boundary (=M2)
  real(DP), parameter, public :: dforcing       = 1.0_DP

  ! Coriolis acceleration (=f)
  real(DP), parameter, public :: dcoraccel      = 0.0_DP
  
  ! Frequency of the tidal constituent (=omega)
  real(DP), parameter, public :: dtidalfreq     = 1.405634296908185E-4_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType     = 1

  ! Bottom stress (=s0)
  real(DP), parameter, public :: dstress        = 1.0E20_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType  = 1

  ! Vertical eddy viscosity (=Av0)
  real(DP), parameter, public :: dviscosity     = 1.0E-3_DP

  ! Flag indicating the existence of an analytical solution
  logical, parameter, public :: bhasAnalyticSolution = .false.

!</constantblock>
#else
#error 'Test case is undefined.'
#endif

!<constantblock description="Problem parameters for all benchmarks">

  ! Gravitational acceleration (=g)
  real(DP), parameter, public :: dgravaccel     = 10.0_DP

!</constantblock>

!</constants>

contains

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

    if     (ibathymetryType .eq. 0) then

      ! linear bottom profile
      dh = dheight0 + (dheight-dheight0)*(1-dx/dlength)

    elseif (ibathymetryType .eq. 1) then

      ! constant bottom profile
      dh = dheight

    elseif (ibathymetryType .eq. 2) then

      ! parabolic bottom profile
      dh = dheight*(dheightRatio+(1.0_DP-dheightRatio)*(1.0-(dy/dwidth)**2))
    end if

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

    if     (istressType .eq. 0) then

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
          
    elseif (istressType .eq. 1) then
      
      ! constant stress
      ds = dstress

    elseif (istressType .eq. 2) then

      ! stress proportional to bathymetry
      ds = dstress * sse_bottomProfile(dx,dy) / dheight

    end if

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

    if     (iviscosityType .eq. 0) then

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
          
    elseif (iviscosityType .eq. 1) then
      
      ! constant eddy viscosity
      dAv = dviscosity

    elseif (iviscosityType .eq. 2) then

      ! eddy viscosity proportional to bathymetry
      dAv = dviscosity * sse_bottomProfile(dx,dy) / dheight

    end if

  end function sse_eddyViscosity

  ! ***************************************************************************
#ifdef USE_COMPILER_INTEL
!<function>

  elemental function sinh(cx)

!<description>
    ! Complex valued hyperbolic sine functions (available in Fortran 2008)
!</description>

!<input>
    complex(DP), intent(in) :: cx
!</input>

!<result>
    complex(DP) :: sinh
!</result>
!</function>

    sinh = -cmplx(0.0_DP,1.0_DP) * sin(cmplx(0.0_DP,1.0_DP)*cx) 

  end function
#endif
  ! ***************************************************************************
#ifdef USE_COMPILER_INTEL
!<function>

  elemental function cosh(cx)

!<description>
    ! Complex valued hyperbolic sine functions (available in Fortran 2008)
!</description>

!<input>
    complex(DP), intent(in) :: cx
!</input>

!<result>
    complex(DP) :: cosh
!</result>
!</function>

    cosh = cos(cmplx(0.0_DP,1.0_DP)*cx) 

  end function
#endif

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
    
  end subroutine

end module sse_base
