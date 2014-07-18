!#########################################################################
!# ***********************************************************************
!# <name> pprocintegrals </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating integrals.
!#
!# The following routines can be found in this module:
!#
!# 1.) ppint_lineIntegral
!#     -> Calculate a line integral
!#
!# </purpose>
!#########################################################################

module pprocintegrals

!$ use omp_lib
  use fsystem
  use genoutput
  use cubature
  use mprimitives
  use basicgeometry
  use triangulation
  use triasearch
  use collection

  implicit none

  private

  public :: ppint_lineIntegral

contains

  !****************************************************************************

!<subroutine>

  subroutine ppint_lineIntegral (dvalue,rtriangulation,Dstart,Dend,&
      ccubature,nlevels,ffunctionRefSimple,rcollection)

!<description>
  ! Calculates a line integration on a function using a summed cubature rule.
!</description>

!<input>
  ! Underlying triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Starting point of the line.
  real(dp), dimension(:), intent(in) :: Dstart

  ! Ending point of the line.
  real(dp), dimension(:), intent(in) :: Dend

  ! Basic 1D cubature rule tu use for line integration.
  integer(I32), intent(in) :: ccubature

  ! Refinement of the cubature rule. >= 0.
  ! The line Dstart..Dend will be divided into 2**idegree intervals and a
  ! summed cubature formula will be applied.
  integer, intent(in) :: nlevels

  ! Callback-Function to integrate.
  include 'intf_functionScSimple.inc'

  ! OPTIONAL: Collection structure to pass to the callback functino.
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
  ! Value of the integral.
  real(DP), intent(out) :: dvalue
!</output>

!</subroutine>

    ! local variables
    integer(I32) :: ccub
    integer :: ncubp, ncdim,i,j,idim
    integer :: ifirstUnknownPoint,ilastUnknownPoint
    real(DP), dimension(:,:), allocatable :: Dpoints,DpointsReal
    real(DP), dimension(:,:,:), allocatable :: DpointsReal2
    real(DP), dimension(:), allocatable :: Domega
    real(DP), dimension(:,:), allocatable :: Dvalues
    integer, dimension(:), allocatable :: Ielements
    real(DP) :: dcubweight

    ! Determine the cubature formula
    ccub = cub_getSummedCubType(ccubature,max(0,nlevels))

    ! Get number of points and coordinate dimension
    ncubp = cub_igetNumPts(ccub)
    ncdim = cub_igetCoordDim(ccub)

    if (ncdim .ne. 1) then
      call output_line ('Invalid cubature rule. Must be 1D!', &
          OU_CLASS_ERROR,OU_MODE_STD,'ppint_lineIntegral')
      call sys_halt()
    end if

    ! Allocate arrays
    allocate(Domega(ncubp))
    allocate(Dpoints(ncdim,ncubp))

    ! Get cubature points/weights
    call cub_getCubature(ccub, Dpoints, Domega)

    ! Recompute into real-world coordinates
    allocate(DpointsReal(rtriangulation%ndim,ncubp))
    allocate(DpointsReal2(rtriangulation%ndim,1,ncubp))

    do i=1,ncubp
      do idim = 1,rtriangulation%ndim
        call mprim_linearRescale(Dpoints(1,i),&
            -1.0_DP,1.0_DP,Dstart(idim),Dend(idim),DpointsReal(idim,i))
        DpointsReal2(idim,1,i) = DpointsReal(idim,i)
      end do
    end do

    ! Reference coordinates no more needed.
    deallocate(Dpoints)

    ! Get the corresponding elements
    allocate(Ielements(ncubp))

    call tsrch_getElementsByRaytrace (rtriangulation,DpointsReal,Ielements,&
        ifirstUnknownPoint,ilastUnknownPoint)

    ! Sometimes, some elements may not be found (e.g. if the point is outside
    ! of the domain). We replace the correspoding function value by zero
    ! in this case and throw a warning. This is done by deleting the points
    ! and elements from the list.
    j = 0
    do i=1,ncubp
      if (Ielements(i) .ne. 0) then
        j = j+1
        DpointsReal2(:,1,j) = DpointsReal2(:,1,i)
        Ielements(j) = Ielements(i)
      else
        select case (ubound(DpointsReal2,1))
        case (NDIM1D)
          call output_line ("Element not found. Point "//trim(sys_siL(i,10))//"=("//&
              trim(sys_sdL(DpointsReal2(1,1,i),5))//")", &
              OU_CLASS_WARNING,OU_MODE_STD,"ppint_lineIntegral")
        case (NDIM2D)
          call output_line ("Element not found. Point "//trim(sys_siL(i,10))//"=("//&
              trim(sys_sdL(DpointsReal2(1,1,i),5))//","//trim(sys_sdL(DpointsReal2(2,1,i),5))//")", &
              OU_CLASS_WARNING,OU_MODE_STD,"ppint_lineIntegral")
        case (NDIM3D)
          call output_line ("Element not found. Point "//trim(sys_siL(i,10))//"=("//&
              trim(sys_sdL(DpointsReal2(1,1,i),5))//","//trim(sys_sdL(DpointsReal2(2,1,i),5))//&
              trim(sys_sdL(DpointsReal2(3,1,i),5))//")", &
              OU_CLASS_WARNING,OU_MODE_STD,"ppint_lineIntegral")
        end select
      end if
    end do

    ! Calculate the values in the points. One point per element.
    allocate(Dvalues(1,ncubp))

    call ffunctionRefSimple (j,1,Ielements,DpointsReal2,Dvalues,rcollection)

    ! Sum up all values to an integral
    dvalue = 0.0_DP

    dcubweight = 0.5_DP*sqrt((Dend(1)-Dstart(1))**2+(Dend(2)-Dstart(2))**2)
    do i=1,j
      dvalue = dvalue + Dvalues(1,i) * dcubweight * Domega(i)
    end do

    ! Release data
    deallocate(Dvalues)
    deallocate(DpointsReal)
    deallocate(DpointsReal2)
    deallocate(Ielements)
    deallocate(Domega)

  end subroutine

end module

