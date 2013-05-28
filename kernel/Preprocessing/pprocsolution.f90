!#########################################################################
!# ***********************************************************************
!# <name> pprocsolution </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines to preprocess the solution vector
!#
!# The following routines can be found in this module:
!#
!# 1.) ppsol_readPGM
!#     -> Reads a Portable Graymap from file.
!#
!# 2.) ppsol_releasePGM
!#     -> Releases a Portable Graymap.
!#
!# 3.) ppsol_initArrayPGMDP / ppsol_initArrayPGMSP
!#     -> Initialises a 2D double array from a Portable Graymap image
!# </purpose>
!#########################################################################
module pprocsolution

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use io
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock

  implicit none

  private

!<types>
!<typeblock>

  ! A portable Graymap

  type t_pgm

    ! Identifier of the PGM
    character(LEN=2) :: id = ''

    ! Width of the image
    integer :: width       = 0

    ! Height of the image
    integer :: height      = 0

    ! Maximum gray-value
    integer :: maxgray     = 0

    ! Handle to the image data
    integer :: h_Idata     = ST_NOHANDLE
  end type t_pgm

  public :: t_pgm

!</typeblock>
!</types>

  public :: ppsol_readPGM
  public :: ppsol_releasePGM
  public :: ppsol_initArrayPGMDP
  public :: ppsol_initArrayPGMSP

contains

  ! ***************************************************************************

!<subroutine>

  subroutine ppsol_readPGM(ifile, sfile, rpgm)

!<description>
    ! This subroutine reads a portable graymap image from file.
!</description>

!<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
!</input>

!<output>
    ! portable graymap
    type(t_pgm), intent(out) :: rpgm
!</output>
!</subroutine>

    ! local variables
    character(LEN=SYS_STRLEN) :: sdata
    integer, dimension(:,:), pointer :: p_Idata
    integer, dimension(2) :: Isize
    integer :: iunit,ilinelen,istatus,ix,iy

    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, iunit, .true.)
      if (iunit .eq. -1) then
        call output_line('Could not open file '//trim(sfile)//'!',&
                         OU_CLASS_ERROR, OU_MODE_STD, 'ppsol_readPGM')
        call sys_halt()
      end if
    else
      iunit = ifile
    end if
    
    ! Read format identifier
    call io_readlinefromfile(iunit, sdata, ilinelen, istatus)
    call sys_tolower(sdata)
    select case(trim(sdata))
    case ('p2')
      rpgm%id='p2'
      
      ! Skip comments (if required)
      comment: do
        call io_readlinefromfile(iunit, sdata, ilinelen, istatus)
        if (sdata(1:1) .ne. '#') exit comment
      end do comment
      
      ! Read width and height
      read(sdata, *) rpgm%width, rpgm%height

      ! Read number of graylevels
      call io_readlinefromfile(iunit, sdata, ilinelen, istatus)
      read(sdata, *) rpgm%maxgray

      ! Allocate memory for image data
      Isize=(/rpgm%width, rpgm%height/)
      call storage_new('ppsol_readPGM', 'h_Idata',&
          Isize, ST_INT, rpgm%h_Idata, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2D(rpgm%h_Idata, p_Idata)
      
      do iy = 1, rpgm%height
        do ix = 1, rpgm%width
          call io_readlinefromfile(iunit, sdata, ilinelen, istatus)
          read(sdata,*) p_Idata(ix,iy)
          if (p_Idata(ix,iy) .gt. rpgm%maxgray) then
            call output_line('Image data exceeds range!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppsol_readPGM')
            call sys_halt()
          end if
        end do
      end do
      
    case default
      call output_line('Invalid PGM format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppsol_readPGM')
      call sys_halt()
    end select
    
    ! Close the file if necessary
    if (ifile .eq. 0) close(iunit)
   
  end subroutine ppsol_readPGM

  ! ***************************************************************************

!<subroutine>

  subroutine ppsol_releasePGM(rpgm)

!<description>
    ! This subroutine releases a portable graymap image.
!</description>

!<inputoutput>
    ! portable graymap image
    type(t_pgm), intent(inout) :: rpgm
!</inputoutput>
!</subroutine>

    rpgm%id      = ''
    rpgm%width   = 0
    rpgm%height  = 0
    rpgm%maxgray = 0

    call storage_free(rpgm%h_Idata)
  end subroutine ppsol_releasePGM

  ! ***************************************************************************

!<subroutine>

  subroutine ppsol_initArrayPGMDP(rpgm, Dpoints, Ddata, Dbounds)

!<description>
    ! Initialises a 2D double array by a Portable Graymap image
!</description>

!<input>
    ! portable graymap image
    type(t_pgm), intent(in) :: rpgm

    ! coordinates of 2D array
    real(DP), dimension(:,:), intent(in) :: Dpoints

    ! OPTIONAL: coordinates of the bounding box
    ! If not present, then the bounding box is calculated internally
    ! DIMENSION(2,2) Dpoints(X:Y,MIN:MAX)
    real(DP), dimension(:,:), intent(in), optional :: Dbounds
!</input>

!<output>
    ! double data array
    real(DP), dimension(:), intent(out) :: Ddata
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_Idata
    real(DP) :: x,y,xmin,ymin,xmax,ymax
    integer :: ipoint,npoints,ix,iy

    ! Set pointer for image data
    call storage_getbase_int2D(rpgm%h_Idata, p_Idata)

    ! Set number of points
    npoints = size(Ddata)

    if (present(Dbounds)) then
      xmin = Dbounds(1,1)
      ymin = Dbounds(2,1)
      xmax = Dbounds(1,2)
      ymax = Dbounds(2,2)
    else
      ! Determine minimum/maximum values of array
      xmin = huge(1.0_DP)
      xmax = -huge(1.0_DP)
      ymin = huge(1.0_DP)
      ymax = -huge(1.0_DP)

      do ipoint = 1, npoints
        xmin = min(xmin, Dpoints(1,ipoint))
        xmax = max(xmax, Dpoints(1,ipoint))
        ymin = min(ymin, Dpoints(2,ipoint))
        ymax = max(ymax, Dpoints(2,ipoint))
      end do
    end if

    ! Clear array
    call lalg_clearVector(Ddata)

    ! Fill array with scaled image data
    do ipoint = 1, npoints
      x = Dpoints(1,ipoint)
      y = Dpoints(2,ipoint)

      ix = 1+(rpgm%width-1)*(x-xmin)/(xmax-xmin)
      if (ix .lt. 1 .or. ix .gt. rpgm%width) cycle

      iy = rpgm%height-(rpgm%height-1)*(y-ymin)/(ymax-ymin)
      if (iy .lt. 1 .or. iy .gt. rpgm%height) cycle

      Ddata(ipoint) = real(p_Idata(ix,iy),DP)/real(rpgm%maxgray,DP)
    end do
  end subroutine ppsol_initArrayPGMDP

  ! ***************************************************************************

!<subroutine>

  subroutine ppsol_initArrayPGMSP(rpgm, Dpoints, Fdata, Dbounds)

!<description>
    ! Initialises a 2D single array by a Portable Graymap image
!</description>

!<input>
    ! portable graymap image
    type(t_pgm), intent(in) :: rpgm

    ! coordinates of 2D array
    real(DP), dimension(:,:), intent(in) :: Dpoints

    ! OPTIONAL: coordinates of the bounding box
    ! If not present, then the bounding box is calculated internally
    ! DIMENSION(2,2) Dpoints(X:Y,MIN:MAX)
    real(DP), dimension(:,:), intent(in), optional :: Dbounds
!</input>

!<output>
    ! single data array
    real(SP), dimension(:), intent(out) :: Fdata
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_Idata
    real(DP) :: x,y,xmin,ymin,xmax,ymax
    integer :: ipoint,npoints,ix,iy

    ! Set pointer for image data
    call storage_getbase_int2D(rpgm%h_Idata, p_Idata)

    ! Set number of points
    npoints = size(Fdata)

    if (present(Dbounds)) then
      xmin = Dbounds(1,1)
      ymin = Dbounds(2,1)
      xmax = Dbounds(1,2)
      ymax = Dbounds(2,2)
    else
      ! Determine minimum/maximum values of array
      xmin = huge(1.0_DP)
      xmax = -huge(1.0_DP)
      ymin = huge(1.0_DP)
      ymax = -huge(1.0_DP)

      do ipoint = 1, npoints
        xmin = min(xmin, Dpoints(1,ipoint))
        xmax = max(xmax, Dpoints(1,ipoint))
        ymin = min(ymin, Dpoints(2,ipoint))
        ymax = max(ymax, Dpoints(2,ipoint))
      end do
    end if

    ! Clear array
    call lalg_clearVector(Fdata)

    ! Fill array with scaled image data
    do ipoint = 1, npoints
      x = Dpoints(1,ipoint)
      y = Dpoints(2,ipoint)

      ix = 1+(rpgm%width-1)*(x-xmin)/(xmax-xmin)
      if (ix .lt. 1 .or. ix .gt. rpgm%width) cycle

      iy = rpgm%height-(rpgm%height-1)*(y-ymin)/(ymax-ymin)
      if (iy .lt. 1 .or. iy .gt. rpgm%height) cycle

      Fdata(ipoint) = real(p_Idata(ix,iy),SP)/real(rpgm%maxgray,SP)
    end do
  end subroutine ppsol_initArrayPGMSP

end module pprocsolution
