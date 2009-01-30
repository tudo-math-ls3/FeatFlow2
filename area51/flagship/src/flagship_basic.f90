!##############################################################################
!# ****************************************************************************
!# <name> flagship_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines for the main program
!#
!# The following routines are available:
!#
!# 1.) flagship_readParserFromFile
!#     -> Read expressions from file and initialize the function parser
!#
!# 2.) flagship_outputVectorScalar
!#     -> Write scalar vector to file in UCD format
!#
!# 3.) flagship_outputVectorBlock
!#     -> Write block vector to file in UCD format
!#
!# </purpose>
!##############################################################################

module flagship_basic

  use fsystem
  use fparser
  use io
  use problem
  use ucd

  implicit none

  private
  public :: flagship_readParserFromFile
  public :: flagship_initUCDexport

contains

  ! ***************************************************************************

!<subroutine>

  subroutine flagship_readParserFromFile(sfilename, ssectionname,&
                                         cvariables, rparser)

!<description>
    ! This subroutine initializes a vector profile from the parameters
    ! specified in the parameter file by calling a function parser w.r.t.
    ! the specified variable names
!</description>
    
!<input>
    ! name of parameter file
    character(LEN=*), intent(IN) :: sfilename

    ! name of the parameter section
    character(LEN=*), intent(IN) :: ssectionname

    ! symbolic variable names
    character(LEN=*), dimension(:), intent(IN) :: cvariables
!</input>

!<output>
    ! function parser
    type(t_fparser), intent(OUT) :: rparser
!</output>
!</subroutine>
    
    ! local variables
    character(SYS_STRLEN) :: skeyword
    character(1024) :: sdata,sexpression
    integer :: iunit,ios,ipos,idatalen,icomp,ncomp
    
    ! Try to open the file
    call io_openFileForReading(sfilename, iunit, .true.)

    ! Oops...
    if (iunit .eq. -1) then
      call output_line('Unable to open input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if

    ! Read through the input file until the given keyword is reached.
    ios = 0
    do while(ios .eq. 0)
      
      ! Read next line in file
      call io_readlinefromfile(iunit, sdata, idatalen, ios)
      if (ios .ne. 0) then
        call output_line('Unable to read KEYWORD from input file!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'flagship_readParserFromFile')
        call sys_halt()
      end if

      ! Check for keyword
      call sys_tolower(sdata(1:idatalen), skeyword)
      if (trim(adjustl(skeyword)) .eq. trim(adjustl(ssectionname))) exit
    end do
    
    ! We found the keyword. String NCOMP must be the next line to read.
    call io_readlinefromfile(iunit, sdata, idatalen, ios)
    if (ios .ne. 0) then
      call output_line('Unable to read data from input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if

    ! Check for keyword NCOMP
    call sys_tolower(sdata(1:idatalen), skeyword)
    if (trim(adjustl(skeyword)) .ne. 'ncomp') then
      call output_line('Syntax error in input file!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if
    
    ! Read value for NCOMP
    read(iunit,*,IOSTAT=ios) ncomp
    if (ios .ne. 0) then
      call output_line('Unable to read value of NCOMP from input file!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if   

    ! Create parser structure
    call fparser_create(rparser, ncomp)
    
    ! Read expressions into parser, relations and concatinations
    do icomp = 1, ncomp
      
      ! Read until expression is finished
      ios  = 0
      ipos = 1
      sexpression = " "

      do while(ios .eq. 0)

        ! Read next line in file
        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_readParserFromFile')
          call sys_halt()
        end if

        ! Append line to expression
        sexpression(ipos:) = sdata(1:idatalen)
        ipos = len_trim(sexpression)
        
        ! Check if expression is continued in the following line
        if (sexpression(max(1,ipos-2):ipos) .eq. '...') then
          ipos = ipos-2
        else
          exit
        end if
      end do

      call fparser_parseFunction(rparser, icomp, sexpression(1:ipos), cvariables)
    end do

    ! We are done. Close the unit.
    close(iunit)

  end subroutine flagship_readParserFromFile

  !*****************************************************************************

!<subroutine>

  subroutine flagship_initUCDexport(rproblemLevel, sfilename, ioutputUCD,&
                                    rexport, ifilenumber)

!<description> 
    ! This subroutine initializes the UCD exporter structure. If the
    ! optional parameter ifilenumber is given, the outputfile is named
    ! 'sfilename'.<ifilenumber>.'ext' where 'ext' is the file
    ! extension that corresponds to the UCD format.
!</description>

!<input>
    ! multigrid structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! name of the output file
    character(LEN=*), intent(IN) :: sfilename

    ! type of UCD output
    integer, intent(IN) :: ioutputUCD

    ! OPTIONAL: number of the output file
    integer, intent(IN), optional :: ifilenumber
!</input>

!<output>
    ! UCD export structure
    type(t_ucdExport), intent(OUT) :: rexport
!</output>
!</subroutine>

    select case(ioutputUCD)
    case (UCD_FORMAT_GMV)
      if (present(ifilenumber)) then
        call ucd_startGMV(rexport, UCD_FLAG_STANDARD,&
                          rproblemLevel%rtriangulation,&
                          trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.gmv')
      else
        call ucd_startGMV(rexport, UCD_FLAG_STANDARD,&
                          rproblemLevel%rtriangulation,&
                          trim(adjustl(sfilename))//'.gmv')
      end if

    case (UCD_FORMAT_BGMV)
      if (present(ifilenumber)) then
        call ucd_startBGMV(rexport, UCD_FLAG_STANDARD,&
                           rproblemLevel%rtriangulation,&
                           trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.gmv')
      else
        call ucd_startBGMV(rexport, UCD_FLAG_STANDARD,&
                           rproblemLevel%rtriangulation,&
                           trim(adjustl(sfilename))//'.gmv')
      end if

    case (UCD_FORMAT_AVS)
      if (present(ifilenumber)) then
        call ucd_startAVS(rexport, UCD_FLAG_STANDARD,&
                          rproblemLevel%rtriangulation,&
                          trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.avs')
      else
        call ucd_startAVS(rexport, UCD_FLAG_STANDARD,&
                          rproblemLevel%rtriangulation,&
                          trim(adjustl(sfilename))//'.avs')
      end if

    case (UCD_FORMAT_VTK)
      if (present(ifilenumber)) then
        call ucd_startVTK(rexport, UCD_FLAG_STANDARD,&
                          rproblemLevel%rtriangulation,&
                          trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.vtu')
      else
        call ucd_startVTK(rexport, UCD_FLAG_STANDARD,&
                          rproblemLevel%rtriangulation,&
                          trim(adjustl(sfilename))//'.vtu')
      end if

    case DEFAULT
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_initUCDexport')
      call sys_halt()
    end select

  end subroutine flagship_initUCDexport

end module flagship_basic
