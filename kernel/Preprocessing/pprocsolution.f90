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
MODULE pprocsolution

  USE fsystem
  USE storage
  USE io
  USE linearsystemscalar
  USE linearsystemblock
  
  IMPLICIT NONE

!<types>
!<typeblock>

  ! A portable Graymap

  TYPE t_pgm

    ! Identifier of the PGM
    CHARACTER(LEN=2) :: id = ''

    ! Width of the image
    INTEGER :: width       = 0

    ! Height of the image
    INTEGER :: height      = 0

    ! Maximum gray-value
    INTEGER :: maxgray     = 0

    ! Handle to the image data
    INTEGER :: h_Idata     = ST_NOHANDLE
  END TYPE t_pgm
  
CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ppsol_readPGM(ifile,sfile,rpgm)

!<description>
    ! This subroutine reads a portable graymap image from file.
!</description>

!<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN)     :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
!</input>

!<output>
    ! portable graymap
    TYPE(t_pgm), INTENT(OUT)     :: rpgm
!</output>
!</subroutine>
    
    ! local variables
    CHARACTER(LEN=80) :: cbuffer
    INTEGER,      DIMENSION(:,:), POINTER :: p_Idata
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER :: cf,ix,iy
    
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForReading(sfile, cf, .TRUE.)
      IF (cf .EQ. -1) THEN
        PRINT *, 'ppsol_readPGM: Could not open file '//trim(sfile)
        CALL sys_halt()
      END IF
    ELSE
      cf = ifile
    END IF

    ! Read image ID
    READ (cf,*) rpgm%id


    SELECT CASE(rpgm%id)
    CASE ('P2','p2')
      ! Read the header
      CALL getNextEntryASCII(cbuffer); READ(cbuffer,*) rpgm%width
      CALL getNextEntryASCII(cbuffer); READ(cbuffer,*) rpgm%height
      CALL getNextEntryASCII(cbuffer); READ(cbuffer,*) rpgm%maxgray
    
      ! Allocate memory for image data
      Isize=(/rpgm%width, rpgm%height/)
      CALL storage_new('ppsol_readPGM', 'h_Idata',&
          Isize, ST_INT, rpgm%h_Idata, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int2D(rpgm%h_Idata, p_Idata)
      
      DO iy = 1, rpgm%height
        DO ix = 1, rpgm%width
          CALL getNextEntryASCII(cbuffer); READ(cbuffer,*) p_Idata(ix,iy)
          IF (p_Idata(ix,iy) .GT. rpgm%maxgray) THEN
            CALL output_line('Image data exceeds range!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppsol_readPGM')
            CALL sys_halt()
          END IF
        END DO
      END DO

    CASE DEFAULT
      CALL output_line('Invalid PGM format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppsol_readPGM')
      CALL sys_halt()
    END SELECT
      
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)

  CONTAINS

    ! Here, the real working routines follow
    
    !**************************************************************
    ! Reads an item from file into the buffer.
    ! Items are separated by spaces, comma, (tabs ...??)
    ! Anything from character '#' till the end of line is ignored
    ! assuming Fortran's eor=eol ... seems to work

    SUBROUTINE getNextEntryASCII(cbuffer)
      
      CHARACTER(LEN=*), INTENT(INOUT) :: cbuffer

      CHARACTER(LEN=1) :: c
      INTEGER          :: ipos

      ipos = 1
      
      ! Read until somthing not a space, comma or comment
      DO
10      READ(cf, FMT='(a)', ADVANCE='NO', END=99, EOR=10) c
        IF ((c .NE. ' ') .AND.&
            (c .NE. ',') .AND.&
            (c .NE.'\t') .AND.&
            (c .NE. '#')) EXIT
        IF (c .EQ. '#') THEN
          ! Skip the comment till the end of line
          DO
            READ(cf, FMT='(A)', ADVANCE='NO', END=99, EOR=10) c
          END DO
        END IF
      END DO
      
      ! We have some significant character, 
      ! read until next space, comma or comment
      DO
        cbuffer(ipos:ipos) = c
        ipos = ipos+1
        READ(cf, FMT='(a)', ADVANCE='NO', END=99, EOR=99) c
        IF ((c .EQ. ' ') .OR.&
            (c .EQ. ',') .OR.&
            (c .EQ.'\t')) EXIT
        IF (c .EQ. '#') THEN
          ! Skip the comment till the end of line
          DO
            READ(cf, FMT='(A)', ADVANCE='NO', END=99, EOR=99) c
          END DO
        END IF
      END DO
      ! End the read characters by one space
99    cbuffer(ipos:ipos) = ' '
    END SUBROUTINE getNextEntryASCII
  END SUBROUTINE ppsol_readPGM

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ppsol_releasePGM(rpgm)

!<description>
    ! This subroutine releases a portable graymap image.
!</description>

!<inputoutput>
    ! portable graymap image
    TYPE(t_pgm), INTENT(INOUT) :: rpgm
!</inputoutput>
!</subroutine>

    rpgm%id      = ''
    rpgm%width   = 0
    rpgm%height  = 0
    rpgm%maxgray = 0

    CALL storage_free(rpgm%h_Idata)
  END SUBROUTINE ppsol_releasePGM

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ppsol_initArrayPGM_Dble(rpgm, DvertexCoords, Ddata)

!<description>
    ! Initialises a 2D double array by a Portable Graymap image
!</description>

!<input>
    ! portable graymap image
    TYPE(t_pgm), INTENT(IN)              :: rpgm

    ! vertex coordinates of 2D array
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: DvertexCoords
!</input>

!<output>
    ! double data array
    REAL(DP), DIMENSION(:), INTENT(OUT)  :: Ddata
!</output>
!</subroutine>

    ! local variables
    INTEGER, DIMENSION(:,:), POINTER :: p_Idata
    REAL(DP)     :: x,y,xmin,ymin,xmax,ymax
    INTEGER(I32) :: ivt,nvt,ix,iy

    ! Set pointer for image data
    CALL storage_getbase_int2D(rpgm%h_Idata, p_Idata)

    ! Determine minimum/maximum values of array
    xmin = HUGE(DP); xmax = -HUGE(DP)
    ymin = HUGE(DP); ymax = -HUGE(DP)

    nvt = SIZE(Ddata)
    DO ivt = 1, nvt
      xmin = MIN(xmin, DvertexCoords(1,ivt))
      xmax = MAX(xmax, DvertexCoords(1,ivt))
      ymin = MIN(ymin, DvertexCoords(2,ivt))
      ymax = MAX(ymax, DvertexCoords(2,ivt))
    END DO

    ! Clear array
    CALL lalg_clearVectorDble(Ddata)

    ! Fill array with scaled image data
    DO ivt = 1, nvt
      x = DvertexCoords(1,ivt)
      y = DvertexCoords(2,ivt)
      
      ix = 1+(rpgm%width-1)*(x-xmin)/(xmax-xmin)
      IF (ix .LT. 1 .OR. ix .GT. rpgm%width) CYCLE

      iy = rpgm%height-(rpgm%height-1)*(y-ymin)/(ymax-ymin)
      IF (iy .LT. 1 .OR. iy .GT. rpgm%height) CYCLE

      Ddata(ivt) = REAL(p_Idata(ix,iy),DP)/REAL(rpgm%maxgray,DP)
    END DO
  END SUBROUTINE ppsol_initArrayPGM_Dble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ppsol_initArrayPGM_Sngl(rpgm, DvertexCoords, Fdata)

!<description>
    ! Initialises a 2D single array by a Portable Graymap image
!</description>

!<input>
    ! portable graymap image
    TYPE(t_pgm), INTENT(IN)              :: rpgm

    ! vertex coordinates of 2D array
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: DvertexCoords
!</input>

!<output>
    ! single data array
    REAL(SP), DIMENSION(:), INTENT(OUT)  :: Fdata
!</output>
!</subroutine>

    ! local variables
    INTEGER, DIMENSION(:,:), POINTER :: p_Idata
    REAL(DP)     :: x,y,xmin,ymin,xmax,ymax
    INTEGER(I32) :: ivt,nvt,ix,iy

    ! Set pointer for image data
    CALL storage_getbase_int2D(rpgm%h_Idata, p_Idata)

    ! Determine minimum/maximum values of array
    xmin = HUGE(DP); xmax = -HUGE(DP)
    ymin = HUGE(DP); ymax = -HUGE(DP)

    nvt = SIZE(Fdata)
    DO ivt = 1, nvt
      xmin = MIN(xmin, DvertexCoords(1,ivt))
      xmax = MAX(xmax, DvertexCoords(1,ivt))
      ymin = MIN(ymin, DvertexCoords(2,ivt))
      ymax = MAX(ymax, DvertexCoords(2,ivt))
    END DO

    ! Clear array
    CALL lalg_clearVectorSngl(Fdata)

    ! Fill array with scaled image data
    DO ivt = 1, nvt
      x = DvertexCoords(1,ivt)
      y = DvertexCoords(2,ivt)
      
      ix = 1+(rpgm%width-1)*(x-xmin)/(xmax-xmin)
      IF (ix .LT. 1 .OR. ix .GT. rpgm%width) CYCLE

      iy = rpgm%height-(rpgm%height-1)*(y-ymin)/(ymax-ymin)
      IF (iy .LT. 1 .OR. iy .GT. rpgm%height) CYCLE

      Fdata(ivt) = REAL(p_Idata(ix,iy),SP)/REAL(rpgm%maxgray,SP)
    END DO
  END SUBROUTINE ppsol_initArrayPGM_Sngl
END MODULE pprocsolution
