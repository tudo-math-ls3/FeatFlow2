!##############################################################################
!# ****************************************************************************
!# <name> geometryoutput </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module adds support for the output of basic geometry objects
!# to files. This can be e.g. the output of lines/points to a file that
!# can be plotted by GNUplot or similar.
!#
!# The following routines can be found here:
!#
!# 1.) geoout_writeGnuplotPoint
!#     -> Writes a point to a Gnuplot file.
!#
!# 2.) geoout_writeGnuplotTria2D
!#     -> Writes a 2D triangular element to a Gnuplot file
!#
!# 3.) geoout_writeGnuplotQuad2D
!#     -> Writes a 2D quad element to a Gnuplot file
!# </purpose>
!##############################################################################
 
module geometryoutput

!$use omp_lib
  use fsystem
  use genoutput
  use basicgeometry
  use io
  
  implicit none
  
  private

  public :: geoout_writeGnuplotPoint
  public :: geoout_writeGnuplotTria2D
  public :: geoout_writeGnuplotQuad2D
  
contains

!************************************************************************
!<subroutine>

  subroutine geoout_writeGnuplotPoint (Dpoint,ifile,sfilename)
  
!<description>
  ! Writes a point Dpoint to a Gnuplot file. The point may be specified
  ! in 1D, 2D or 3D.
!</description>

!<input>
  ! The point to write out.
  real(dp), dimension(:), intent(in) :: Dpoint
  
  ! Handle of the file where to write out. If 0 is specified here,
  ! a filename sfilename must be specified where the data is written to.
  ! The data will be appendet to the file.
  integer, intent(in) :: ifile
  
  ! OPTIONAL: Filename of the file where to write to.
  ! This must be specified if ifile=0.
  character(len=*), intent(in), optional :: sfilename
!</subroutine>

    integer :: ihandle
    
    ! If ifile=0, open a new file and append.
    if (ifile .ne. 0) then
      ihandle = ifile
    else
      if (.not. present (sfilename)) then
        call output_line('No filename specified!',&
            OU_CLASS_ERROR,OU_MODE_STD,'geoout_writeGnuplotPoint')
        call sys_halt()
      end if
      
      if (sfilename .eq. '') then
        call output_line('No filename specified!',&
            OU_CLASS_ERROR,OU_MODE_STD,'geoout_writeGnuplotPoint')
        call sys_halt()
      end if
      
      call io_openFileForWriting(sfilename, ihandle, SYS_APPEND, bformatted=.true.)
    end if
    
    ! Write out the data, an empty line in front.
    write (ihandle,'(A)') ''
    select case (size(Dpoint))
    case (NDIM1D)
      write (ihandle,'(ES16.8E3)') Dpoint(1)
    case (NDIM2D)
      write (ihandle,'(2ES16.8E3)') Dpoint(1:2)
    case (NDIM3D:)
      write (ihandle,'(3ES16.8E3)') Dpoint(1:3)
    end select
    
    ! Close the file again
    if (ifile .eq. 0) then
      close (ihandle)
    end if

  end subroutine

!************************************************************************

!<subroutine>

  subroutine geoout_writeGnuplotTria2D (Dpoints,ifile,sfilename)
  
!<description>
  ! Writes a 2D tria element to a Gnuplot file.
!</description>

!<input>
  ! The points that form the element, given as 2-tuples for the X- and Y-
  ! coordinates.
  real(dp), dimension(:,:), intent(in) :: Dpoints
  
  ! Handle of the file where to write out. If 0 is specified here,
  ! a filename sfilename must be specified where the data is written to.
  ! The data will be appendet to the file.
  integer, intent(in) :: ifile
  
  ! OPTIONAL: Filename of the file where to write to.
  ! This must be specified if ifile=0.
  character(len=*), intent(in), optional :: sfilename
!</subroutine>

    integer :: ihandle
    
    ! If ifile=0, open a new file and append.
    if (ifile .ne. 0) then
      ihandle = ifile
    else
      if (.not. present (sfilename)) then
        call output_line('No filename specified!',&
            OU_CLASS_ERROR,OU_MODE_STD,'geoout_writeGnuplotPoint')
        call sys_halt()
      end if
      
      if (sfilename .eq. '') then
        call output_line('No filename specified!',&
            OU_CLASS_ERROR,OU_MODE_STD,'geoout_writeGnuplotPoint')
        call sys_halt()
      end if
      
      call io_openFileForWriting(sfilename, ihandle, SYS_APPEND, bformatted=.true.)
    end if
    
    ! Write out the data, an empty line in front.
    write (ihandle,'(A)') ''
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,1)
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,2)
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,3)
    ! Repeat the 1st point to close the polygon
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,1)
    
    ! Close the file again
    if (ifile .eq. 0) then
      close (ihandle)
    end if

  end subroutine

!************************************************************************

!<subroutine>

  subroutine geoout_writeGnuplotQuad2D (Dpoints,ifile,sfilename)
  
!<description>
  ! Writes a 2D quad element to a Gnuplot file.
!</description>

!<input>
  ! The points that form the element, given as 2-tuples for the X- and Y-
  ! coordinates.
  real(dp), dimension(:,:), intent(in) :: Dpoints
  
  ! Handle of the file where to write out. If 0 is specified here,
  ! a filename sfilename must be specified where the data is written to.
  ! The data will be appendet to the file.
  integer, intent(in) :: ifile
  
  ! OPTIONAL: Filename of the file where to write to.
  ! This must be specified if ifile=0.
  character(len=*), intent(in), optional :: sfilename
!</subroutine>

    integer :: ihandle
    
    ! If ifile=0, open a new file and append.
    if (ifile .ne. 0) then
      ihandle = ifile
    else
      if (.not. present (sfilename)) then
        call output_line('No filename specified!',&
            OU_CLASS_ERROR,OU_MODE_STD,'geoout_writeGnuplotPoint')
        call sys_halt()
      end if
      
      if (sfilename .eq. '') then
        call output_line('No filename specified!',&
            OU_CLASS_ERROR,OU_MODE_STD,'geoout_writeGnuplotPoint')
        call sys_halt()
      end if
      
      call io_openFileForWriting(sfilename, ihandle, SYS_APPEND, bformatted=.true.)
    end if
    
    ! Write out the data, an empty line in front.
    write (ihandle,'(A)') ''
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,1)
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,2)
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,3)
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,4)
    ! Repeat the 1st point to close the polygon
    write (ihandle,'(2ES16.8E3)') Dpoints(1:2,1)
    
    ! Close the file again
    if (ifile .eq. 0) then
      close (ihandle)
    end if

  end subroutine

end module
