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
!# 3.) ppsol_initArrayPGM_Dble
!#     -> Initialises a 2D double array from a Portable Graymap image
!# </purpose>
!#########################################################################
module pprocsolution

  use fsystem
  use genoutput
  use storage
  use io
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
  public :: ppsol_initArrayPGM_Dble
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine ppsol_readPGM(ifile,sfile,rpgm)

!<description>
    ! This subroutine reads a portable graymap image from file.
!</description>

!<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer(I32), intent(IN)     :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
!</input>

!<output>
    ! portable graymap
    type(t_pgm), intent(OUT)     :: rpgm
!</output>
!</subroutine>
    
    ! local variables
    character(LEN=80) :: cbuffer
    integer, dimension(:,:), pointer :: p_Idata
    integer, dimension(2) :: Isize
    integer :: cf,ix,iy
    
    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, cf, .true.)
      if (cf .eq. -1) then
        print *, 'ppsol_readPGM: Could not open file '//trim(sfile)
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Read image ID
    read (cf,*) rpgm%id


    select case(rpgm%id)
    case ('P2','p2')
      ! Read the header
      call getNextEntryASCII(cbuffer); read(cbuffer,*) rpgm%width
      call getNextEntryASCII(cbuffer); read(cbuffer,*) rpgm%height
      call getNextEntryASCII(cbuffer); read(cbuffer,*) rpgm%maxgray
    
      ! Allocate memory for image data
      Isize=(/rpgm%width, rpgm%height/)
      call storage_new('ppsol_readPGM', 'h_Idata',&
          Isize, ST_INT, rpgm%h_Idata, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2D(rpgm%h_Idata, p_Idata)
      
      do iy = 1, rpgm%height
        do ix = 1, rpgm%width
          call getNextEntryASCII(cbuffer); read(cbuffer,*) p_Idata(ix,iy)
          if (p_Idata(ix,iy) .gt. rpgm%maxgray) then
            call output_line('Image data exceeds range!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppsol_readPGM')
            call sys_halt()
          end if
        end do
      end do

    case DEFAULT
      call output_line('Invalid PGM format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppsol_readPGM')
      call sys_halt()
    end select
      
    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  contains

    ! Here, the real working routines follow
    
    !**************************************************************
    ! Reads an item from file into the buffer.
    ! Items are separated by spaces, comma, (tabs ...??)
    ! Anything from character '#' till the end of line is ignored
    ! assuming Fortran's eor=eol ... seems to work

    subroutine getNextEntryASCII(cbuffer)
      
      character(LEN=*), intent(INOUT) :: cbuffer

      character(LEN=1) :: c
      integer          :: ipos

      ipos = 1
      
      ! Read until somthing not a space, comma or comment
      do
10      read(cf, FMT='(a)', ADVANCE='NO', end=99, EOR=10) c
        if ((c .ne. ' ') .and.&
            (c .ne. ',') .and.&
            (c .ne.'\t') .and.&
            (c .ne. '#')) exit
        if (c .eq. '#') then
          ! Skip the comment till the end of line
          do
            read(cf, FMT='(A)', ADVANCE='NO', end=99, EOR=10) c
          end do
        end if
      end do
      
      ! We have some significant character, 
      ! read until next space, comma or comment
      do
        cbuffer(ipos:ipos) = c
        ipos = ipos+1
        read(cf, FMT='(a)', ADVANCE='NO', end=99, EOR=99) c
        if ((c .eq. ' ') .or.&
            (c .eq. ',') .or.&
            (c .eq.'\t')) exit
        if (c .eq. '#') then
          ! Skip the comment till the end of line
          do
            read(cf, FMT='(A)', ADVANCE='NO', end=99, EOR=99) c
          end do
        end if
      end do
      ! End the read characters by one space
99    cbuffer(ipos:ipos) = ' '
    end subroutine getNextEntryASCII
  end subroutine ppsol_readPGM

  ! ***************************************************************************

!<subroutine>

  subroutine ppsol_releasePGM(rpgm)

!<description>
    ! This subroutine releases a portable graymap image.
!</description>

!<inputoutput>
    ! portable graymap image
    type(t_pgm), intent(INOUT) :: rpgm
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

  subroutine ppsol_initArrayPGM_Dble(rpgm, DvertexCoords, Ddata)

!<description>
    ! Initialises a 2D double array by a Portable Graymap image
!</description>

!<input>
    ! portable graymap image
    type(t_pgm), intent(IN)              :: rpgm

    ! vertex coordinates of 2D array
    real(DP), dimension(:,:), intent(IN) :: DvertexCoords
!</input>

!<output>
    ! double data array
    real(DP), dimension(:), intent(OUT)  :: Ddata
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_Idata
    real(DP) :: x,y,xmin,ymin,xmax,ymax
    integer :: ivt,nvt,ix,iy

    ! Set pointer for image data
    call storage_getbase_int2D(rpgm%h_Idata, p_Idata)

    ! Determine minimum/maximum values of array
    xmin = huge(DP); xmax = -huge(DP)
    ymin = huge(DP); ymax = -huge(DP)

    nvt = size(Ddata)
    do ivt = 1, nvt
      xmin = min(xmin, DvertexCoords(1,ivt))
      xmax = max(xmax, DvertexCoords(1,ivt))
      ymin = min(ymin, DvertexCoords(2,ivt))
      ymax = max(ymax, DvertexCoords(2,ivt))
    end do

    ! Clear array
    call lalg_clearVectorDble(Ddata)

    ! Fill array with scaled image data
    do ivt = 1, nvt
      x = DvertexCoords(1,ivt)
      y = DvertexCoords(2,ivt)
      
      ix = 1+(rpgm%width-1)*(x-xmin)/(xmax-xmin)
      if (ix .lt. 1 .or. ix .gt. rpgm%width) cycle

      iy = rpgm%height-(rpgm%height-1)*(y-ymin)/(ymax-ymin)
      if (iy .lt. 1 .or. iy .gt. rpgm%height) cycle

      Ddata(ivt) = real(p_Idata(ix,iy),DP)/real(rpgm%maxgray,DP)
    end do
  end subroutine ppsol_initArrayPGM_Dble

  ! ***************************************************************************

!<subroutine>

  subroutine ppsol_initArrayPGM_Sngl(rpgm, DvertexCoords, Fdata)

!<description>
    ! Initialises a 2D single array by a Portable Graymap image
!</description>

!<input>
    ! portable graymap image
    type(t_pgm), intent(IN)              :: rpgm

    ! vertex coordinates of 2D array
    real(DP), dimension(:,:), intent(IN) :: DvertexCoords
!</input>

!<output>
    ! single data array
    real(SP), dimension(:), intent(OUT)  :: Fdata
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_Idata
    real(DP) :: x,y,xmin,ymin,xmax,ymax
    integer :: ivt,nvt,ix,iy

    ! Set pointer for image data
    call storage_getbase_int2D(rpgm%h_Idata, p_Idata)

    ! Determine minimum/maximum values of array
    xmin = huge(DP); xmax = -huge(DP)
    ymin = huge(DP); ymax = -huge(DP)

    nvt = size(Fdata)
    do ivt = 1, nvt
      xmin = min(xmin, DvertexCoords(1,ivt))
      xmax = max(xmax, DvertexCoords(1,ivt))
      ymin = min(ymin, DvertexCoords(2,ivt))
      ymax = max(ymax, DvertexCoords(2,ivt))
    end do

    ! Clear array
    call lalg_clearVectorSngl(Fdata)

    ! Fill array with scaled image data
    do ivt = 1, nvt
      x = DvertexCoords(1,ivt)
      y = DvertexCoords(2,ivt)
      
      ix = 1+(rpgm%width-1)*(x-xmin)/(xmax-xmin)
      if (ix .lt. 1 .or. ix .gt. rpgm%width) cycle

      iy = rpgm%height-(rpgm%height-1)*(y-ymin)/(ymax-ymin)
      if (iy .lt. 1 .or. iy .gt. rpgm%height) cycle

      Fdata(ivt) = real(p_Idata(ix,iy),SP)/real(rpgm%maxgray,SP)
    end do
  end subroutine ppsol_initArrayPGM_Sngl
end module pprocsolution
