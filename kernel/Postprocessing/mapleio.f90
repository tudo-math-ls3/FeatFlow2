!#########################################################################
!# ***********************************************************************
!# <name> mapleio </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains different routines for MAPLE export.
!# This feature can be used e.g. for debugging purposes to quickly analyse
!# situaltions inside of the program.
!#
!# The following routines can be found in this module:
!#
!# 1.) mapleio_writePointArray
!#     -> Writes a MAPLE command to an output channel for creating a
!#        point array
!#
!# 2.) mapleio_writePolygonPlot
!#     -> Writes a MAPLE PLOT command to an output channel for creating a
!#        polygon
!#
!# 3.) mapleio_writePointPlot
!#     -> Writes a MAPLE PLOT command to an output channel for creating a
!#        point set
!'
!# </purpose>
!#########################################################################

module mapleio

!$use omp_lib
  use fsystem
  use storage
  use io
  use basicgeometry

  implicit none

  private

  interface mapleio_writePointPlot
    module procedure mapleio_writePointPlotSingle
    module procedure mapleio_writePointPlotMult
  end interface

  public :: mapleio_writePointArray
  public :: mapleio_writePolygonPlot
  public :: mapleio_writePointPlot

contains

  ! ***************************************************************************

!<subroutine>
  subroutine mapleio_writePointArray (ichannel,Dpoints,sname,iindex)

  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE array containing the points in Dpoints. The MAPLE variable
    ! name will get the name 'sname'. If
    ! iindex is specified, the array gets the name 'sname(iindex)'.
  !</description>

  !<input>
    ! Output channel number where to write the MAPLE output to.
    integer, intent(in) :: ichannel

    ! Point set that should be written out. This is a list of (2D or 2D)
    ! point coordinates that form the polygon.
    real(DP), dimension(:,:), intent(in) :: Dpoints

    ! Name of the variable to be used by MAPLE.
    character(len=*), intent(in) :: sname

    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    integer, intent(in), optional :: iindex
  !</input>

!</subroutine>

    integer :: ipoint,idim

    ! Create the first line
    if (.not. present(iindex)) then
      write (ichannel,'(A)') sname//':=['
    else
      write (ichannel,'(A)') sname//'('//trim(sys_siL(iindex,10))//'):=['
    end if

    ! Write the points
    do ipoint=1,ubound(Dpoints,2)

      write (ichannel,'(A)',ADVANCE='NO') &
        '     ['//trim(sys_sdL(Dpoints(1,ipoint),10))

      do idim=2,ubound(Dpoints,1)
        write (ichannel,'(A)',ADVANCE='NO') &
          ','//trim(sys_sdL(Dpoints(idim,ipoint),10))
      end do

      if (ipoint .eq. ubound(Dpoints,2)) then
        write (ichannel,'(A)',ADVANCE='NO') ']'
      else
        write (ichannel,'(A)',ADVANCE='YES') '],'
      end if

    end do

    write (ichannel,'(A)') '];';

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine mapleio_writePolygonPlot (ichannel,Dpoints,sname,iindex)

  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE PLOT command for a polygon defined  by the points in Dpoints.
    ! The MAPLE variable name will get the name 'sname'. If
    ! iindex is specified, the variable gets the name 'sname(iindex)'.
  !</description>

  !<input>
    ! Output channel number where to write the MAPLE output to.
    integer, intent(in) :: ichannel

    ! Polygon that should be created. This is a list of (2D or 2D)
    ! point coordinates that form the polygon.
    real(DP), dimension(:,:), intent(in) :: Dpoints

    ! Name of the variable to be used by MAPLE.
    character(len=*), intent(in) :: sname

    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    integer, intent(in), optional :: iindex
  !</input>

!</subroutine>

    integer :: ipoint,idim

    ! Create the first line
    if (.not. present(iindex)) then
      write (ichannel,'(A)') sname//':=polygonplot(['
    else
      write (ichannel,'(A)') sname//'('//trim(sys_siL(iindex,10))//&
        '):=polygonplot(['
    end if

    ! Write the points
    do ipoint=1,ubound(Dpoints,2)

      write (ichannel,'(A)',ADVANCE='NO') &
        '     ['//trim(sys_sdL(Dpoints(1,ipoint),10))

      do idim=2,ubound(Dpoints,1)
        write (ichannel,'(A)',ADVANCE='NO') &
          ','//trim(sys_sdL(Dpoints(idim,ipoint),10))
      end do

      if (ipoint .eq. ubound(Dpoints,2)) then
        write (ichannel,'(A)',ADVANCE='NO') ']'
      else
        write (ichannel,'(A)',ADVANCE='YES') '],'
      end if

    end do

    write (ichannel,'(A)')']);';

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine mapleio_writePointPlotSingle (ichannel,Dpoint,sname,iindex)

  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE PLOT command for a point defined by Dpoint.
    ! The MAPLE variable name will get the name 'sname'. If
    ! iindex is specified, the variable gets the name 'sname(iindex)'.
  !</description>

  !<input>
    ! Output channel number where to write the MAPLE output to.
    integer, intent(in) :: ichannel

    ! Point that should be written out.
    real(DP), dimension(:), intent(in) :: Dpoint

    ! Name of the variable to be used by MAPLE.
    character(len=*), intent(in) :: sname

    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    integer, intent(in), optional :: iindex
  !</input>

!</subroutine>

    integer :: idim

    ! Create the first line
    if (.not. present(iindex)) then
      write (ichannel,'(A)') sname//':=pointplot('
    else
      write (ichannel,'(A)') sname//'('//trim(sys_siL(iindex,10))//&
        '):=polygonplot('
    end if

    ! Write the points
    write (ichannel,'(A)',ADVANCE='NO') &
      '     ['//trim(sys_sdL(Dpoint(1),10))

    do idim=2,size(Dpoint)
      write (ichannel,'(A)',ADVANCE='NO') &
        ','//trim(sys_sdL(Dpoint(idim),10))
    end do

    write (ichannel,'(A)')']);';

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine mapleio_writePointPlotMult (ichannel,Dpoints,sname,iindex)

  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE PLOT command for a point set defined  by the points in Dpoints.
    ! The MAPLE variable name will get the name 'sname'. If
    ! iindex is specified, the variable gets the name 'sname(iindex)'.
  !</description>

  !<input>
    ! Output channel number where to write the MAPLE output to.
    integer, intent(in) :: ichannel

    ! Polygon that should be created. This is a list of (2D or 2D)
    ! point coordinates that form the polygon.
    real(DP), dimension(:,:), intent(in) :: Dpoints

    ! Name of the variable to be used by MAPLE.
    character(len=*), intent(in) :: sname

    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    integer, intent(in), optional :: iindex
  !</input>

!</subroutine>

    integer :: ipoint,idim

    ! Create the first line
    if (.not. present(iindex)) then
      write (ichannel,'(A)') sname//':=pointplot('
    else
      write (ichannel,'(A)') sname//'('//trim(sys_siL(iindex,10))//&
        '):=polygonplot(['
    end if

    ! Write the point
    write (ichannel,'(A)',ADVANCE='NO') &
      '     ['//sys_sd(Dpoints(1,ipoint),10)

    do idim=2,ubound(Dpoints,1)
      write (ichannel,'(A)',ADVANCE='NO') &
        ','//sys_sd(Dpoints(idim,ipoint),10)
    end do

    if (ipoint .eq. ubound(Dpoints,2)) then
      write (ichannel,'(A)',ADVANCE='YES') '],'
    else
      write (ichannel,'(A)',ADVANCE='NO') ']'
    end if

    write (ichannel,'(A)')');';

  end subroutine

end module
