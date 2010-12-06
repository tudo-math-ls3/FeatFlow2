!##############################################################################
!# ****************************************************************************
!# <name> textstream </name>
!# ****************************************************************************
!#
!# <purpose>
!# Encapsules a text file in a stream. Allows to read in a text file into
!# a stream object and read it line by line. Allows to create a stream, write
!# lines to it and produce a text file output from the stream.
!#
!# Shortly said, a textstream is a line based memory buffer for text files.
!#
!# The following routines can be found here:
!#
!# 1.) tstream_initStream
!#     -> Initialises a text stream
!#
!# 2.) tstream_doneStream
!#     -> Releases a text stream
!#
!# 3.) tstream_appendLine
!#     -> Appends a line to a text stream
!#
!# 4.) tstream_readLine
!#     -> Reads a text line
!#
!# 5.) tstream_peekLine
!#     -> Returns the next line which would be read by tstream_readLine
!#        without actually reading it.
!#
!# 6.) tstream_eof
!#     -> Checks if the read position is at the end of the stream
!#
!# 7.) tstream_pos
!#     -> Returns the number of the current line in reading
!#
!# 8.) tstream_length
!#     -> Returns the total number of lines in the buffer
!#
!# 9.) tstream_readFromFile
!#     -> Read a text file from disc and save it in a text stream
!#
!# 10.) tstream_writeToFile
!#     -> Writes the content of a text stream to a file on the disc.
!#
!# Auxiliary routines (not accessible from outside):
!#
!# 1.) tstream_writeLine
!#     -> Write data to the buffer
!#
!# </purpose>
!##############################################################################

module textstream

  use fsystem
  use io
  use genoutput
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! A text stream object. Used to store text files in memory.
  type t_textstream
  
    private
  
    ! Current read position
    integer :: ipositionRead = 1

    ! Current write position
    integer :: ipositionWrite = 1
    
    ! Number of lines in the buffer
    integer :: ilineCount = 0
    
    ! Current line number
    integer :: icurrentLine = 0
    
    ! Buffer of characters. The lines are divided by CHAR(11) characters.
    character, dimension(:), pointer :: p_Sbuf => null()
    
  end type

!</typeblock>

  public :: t_textstream

!</types>

!<constants>

!<constantblock>
  ! Size of chunks in the buffer (in characters)
  integer, parameter :: TSTREAM_BUFSIZE = 8192
  
  ! Maximum length of a line
  integer, parameter :: TSTREAM_LENLINEBUF = 1024
!</constantblock>

!</constants>

  public :: tstream_initStream
  public :: tstream_doneStream
  public :: tstream_appendLine
  public :: tstream_readLine
  public :: tstream_eof
  public :: tstream_length
  public :: tstream_peekLine
  public :: tstream_readFromFile
  public :: tstream_writeToFile

contains

! ****************************************************************************************

!<subroutine>

  subroutine tstream_reallocBuffer (rtextstream, isize)

!<description>
  ! Reallocates the stream buffer to size isize.
!</description>

!<input>
  ! New size of the buffer.
  integer, intent(in) :: isize
!</input>

!<inputoutput>
  ! A text stream object.
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!</subroutine>

    ! local variables
    character, dimension(:), pointer :: p_Sbuf => null()
    integer :: ilen

    ! Allocate, copy, release, replace... as usual.
    if (size(rtextstream%p_Sbuf) .eq. isize) return
    
    ilen = min(size(rtextstream%p_Sbuf),isize)

    allocate(p_Sbuf(isize))
    p_Sbuf(1:ilen) = rtextstream%p_Sbuf(1:ilen)
    deallocate(rtextstream%p_Sbuf)
    rtextstream%p_Sbuf => p_Sbuf

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_initStream (rtextstream)

!<description>
  ! Creates a new text stream for reading and writing.
!</description>

!<output>
  ! A text stream object.
  type(t_textstream), intent(out) :: rtextstream
!</output>

!</subroutine>

    ! Allocate memory for the buffer The other variables are set 
    ! by default initialisation.
    allocate(rtextstream%p_Sbuf(TSTREAM_BUFSIZE))
    
  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_reset (rtextstream)

!<description>
  ! Resets the read pointer of a stream.
!</description>

!<output>
  ! A text stream object.
  type(t_textstream), intent(out) :: rtextstream
!</output>

!</subroutine>

    ! Reset the pointers.
    rtextstream%ipositionRead = 1
    rtextstream%icurrentLine = 1

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_clear (rtextstream)

!<description>
  ! Clears a stream, removes all content.
!</description>

!<output>
  ! A text stream object.
  type(t_textstream), intent(out) :: rtextstream
!</output>

!</subroutine>

    ! Reset the pointers.
    rtextstream%ipositionRead = 1
    rtextstream%ipositionWrite = 1
    rtextstream%ilineCount = 0
    rtextstream%icurrentLine = 1

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_doneStream (rtextstream)

!<description>
  ! Releases a text stream object.
!</description>

!<inputoutput>
  ! A text stream object to be released.
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!</subroutine>

    ! Dellocate memory of the buffer.
    deallocate(rtextstream%p_Sbuf)
    rtextstream%ipositionRead = 1
    rtextstream%ipositionWrite = 1
    rtextstream%icurrentLine = 1
    rtextstream%ilineCount = 0

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_appendLine (rtextstream, sstring)

!<description>
  ! Writes a text line to the end of the buffer of the text stream.
!</description>

!<input>
  ! String to write to the buffer.
  character(len=*), intent(in) :: sstring
!</input>

!<inputoutput>
  ! A text stream object where to write to.
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!</subroutine>

    ! Adjust the write position and write.
    ! rtextstream%ipositionWrite = rtextstream%ilength+1
    call tstream_writeLine (rtextstream, sstring)

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_writeLine (rtextstream, sstring, btrim)

!<description>
  ! Writes a text line to the current write position in the text stream.
  ! The buffer will be truncated after the written text.
!</description>

!<input>
  ! String to write to the buffer.
  character(len=*), intent(in) :: sstring
  
  ! OPTIONAL: Trim the text before writing it to the buffer.
  ! If not specified, TRUE is assumed.
  logical, intent(in), optional :: btrim
!</input>

!<inputoutput>
  ! A text stream object where to write to.
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilen,inewsize,i
    
    ! Do not use .OR. with PRESENT() !
    if (.not. present(btrim)) then
      ilen = len(trim(sstring))
    else
      if (btrim) then
        ilen = len(trim(sstring))
      else
        ilen = len(sstring)
      end if
    end if
    
    ! If the buffer is not large enough, resize.
    if ((rtextstream%ipositionWrite+ilen+1) .gt. size(rtextstream%p_Sbuf)) then
      inewsize = ((rtextstream%ipositionWrite+ilen+1+TSTREAM_BUFSIZE-1) / &
          TSTREAM_BUFSIZE)*TSTREAM_BUFSIZE
      call tstream_reallocBuffer(rtextstream,inewsize)
    end if
    
    ! Append the data -- and a CHAR(11) as EOL character.
    do i=1,ilen
      rtextstream%p_Sbuf(rtextstream%ipositionWrite) = sstring(i:i)
      rtextstream%ipositionWrite = rtextstream%ipositionWrite + 1
    end do

    rtextstream%p_Sbuf(rtextstream%ipositionWrite) = char(11)
    rtextstream%ipositionWrite = rtextstream%ipositionWrite + 1
    
    rtextstream%ilineCount = rtextstream%ilineCount + 1

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_readLine (rtextstream, sstring, ilength)

!<description>
  ! Reads a text line from the current read position in the text stream.
!</description>

!<inputoutput>
  ! A text stream object where to read from.
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!<output>
  ! String from the buffer.
  character(len=*), intent(out) :: sstring
  
  ! OPTIONAL: Length of the string.
  integer, intent(out), optional :: ilength
!</output>

!</subroutine>

    ! local variables
    integer :: ilen,inewsize,i
    character :: c
    
    ! Maximum length of the string
    ilen = len(sstring)
  
    sstring = ""
    
    ! Read text until we find a CHAR(11) EOL character or the end of the stream.
    ilength = 0
    do 

      ! Cancel if we are at the end of the stream.
      if (rtextstream%ipositionRead .ge. rtextstream%ipositionWrite) exit

      c = rtextstream%p_Sbuf(rtextstream%ipositionRead)
      rtextstream%ipositionRead = rtextstream%ipositionRead + 1

      if (c .eq. char(11)) then
        rtextstream%ipositionRead = rtextstream%ipositionRead + 1
        rtextstream%icurrentLine = rtextstream%icurrentLine + 1

        exit
      end if
      
      i = i + 1
      if (present(ilength)) ilength = i
      
      sstring(i:i) = c
      
    end do

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_peekLine (rtextstream, sstring)

!<description>
  ! Reads a text line from the current read position in the text stream, but does not
  ! change the read position, so the text can be read again.
!</description>

!<inputoutput>
  ! A text stream object where to read from.
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!<output>
  ! String from the buffer.
  character(len=*), intent(out) :: sstring
!</output>

!</subroutine>

    ! local variables
    integer :: i,j
    
    i = rtextstream%ipositionRead
    j = rtextstream%icurrentLine
    call tstream_readLine (rtextstream, sstring)
    rtextstream%ipositionRead = i
    rtextstream%icurrentLine = j

  end subroutine

! ****************************************************************************************

!<function>

  logical function tstream_eof (rtextstream)

!<description>
  ! Returns true, if the read position is at the end of the stream.
!</description>

!<input>
  ! A text stream object to test.
  type(t_textstream), intent(in) :: rtextstream
!</input>

!<result>
  ! true, if the read position is at the end of the stream.
!</result>

!</function>

    tstream_eof = rtextstream%ipositionRead .ge. rtextstream%ipositionWrite

  end function

! ****************************************************************************************

!<function>

  integer function tstream_length (rtextstream)

!<description>
  ! Returns the number of lines in the buffer.
!</description>

!<input>
  ! A text stream object to test.
  type(t_textstream), intent(in) :: rtextstream
!</input>

!<result>
  ! Number of lines in the buffer.
!</result>

!</function>

    tstream_length = rtextstream%ilineCount

  end function

! ****************************************************************************************

!<function>

  integer function tstream_pos (rtextstream)

!<description>
  ! Returns the current line number.
!</description>

!<input>
  ! A text stream object to test.
  type(t_textstream), intent(in) :: rtextstream
!</input>

!<result>
  ! Current line number.
!</result>

!</function>

    tstream_pos = rtextstream%icurrentLine

  end function

! ****************************************************************************************

!<subroutine>

  subroutine tstream_readFromFile (rtextstream, sfilename, ilinecount)

!<description>
  ! Reads a text file into a stream buffer.
!</description>

!<inputoutput>
  ! A text stream object where to write to. The previous content of the
  ! stream is discarded.
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!<input>
  ! Filename of the file to read.
  character(len=*), intent(in) :: sfilename
!</input>

!<output>
  ! OPTIONAL: Number of lines read from the file.
  ! -1=File could not be opened.
  integer, intent(out), optional :: ilinecount
!</output>

!</subroutine>

    ! local variables
    integer :: iunit,ios,isbuflen,ilinenum
    character(len=TSTREAM_LENLINEBUF) :: sdata

    ! Reset the stream.
    call tstream_clear(rtextstream)

    ! Try to open the file
    call io_openFileForReading(sfilename, iunit)
    
    ! Oops...
    if (iunit .eq. -1) then
      if (present(ilinecount)) ilinecount = -1
      return
    end if
  
    ! Read all lines from the file
    ios = 0
    ilinenum = 0
    do while (ios .eq. 0) 
      
      ! Read a line from the file into sbuf
      call preadlinefromfile (iunit, sdata, isbuflen, ios)
      ilinenum = ilinenum + 1
      
      if (isbuflen .ne. 0) then
        call tstream_appendLine (rtextstream, sdata(1:isbuflen))
      end if

    end do
    
    ! Close the file, finish.
    close (iunit)
    
    if (present(ilinecount)) ilinecount = ilinenum
    
  contains
  
    ! Internal subroutine: Read a line from a text file.
    
    subroutine readlinefromfile (iunit, sdata, ilinelen, ios)
      
      ! The unit where to read from; must be connected to a file.
      integer, intent(in) :: iunit
    
      ! The string where to write data to
      character(len=*), intent(out) :: sdata
      
      ! Length of the output
      integer, intent(out) :: ilinelen
      
      ! Status of the reading process. Set to a value <> 0 if the end
      ! of the file is reached.
      integer, intent(out) :: ios
      
      ! local variables
      character :: c
      
      sdata = ''
      ilinelen = 0
      
      ! Read the data - as long as the line/file does not end.
      do
        
        ! Read a character.
        ! Unfortunately, Fortran forces me to use this dirty GOTO
        ! to decide processor-independently whether the line or
        ! the record ends.
        read (unit=iunit, fmt='(A1)', iostat=ios, advance='NO',&
            end=10, eor=20) c
        
        ! Do not do anything in case of an error
        if (ios .eq. 0) then
          
          ilinelen = ilinelen + 1
          sdata (ilinelen:ilinelen) = c
          
        end if
        
        ! Proceed to next character
        cycle
        
        ! End of file. 
10      ios = -1
        exit
        
        ! End of record = END OF LINE.
20      ios = 0
        exit
        
      end do
      
    end subroutine

  end subroutine

! ****************************************************************************************

!<subroutine>

  subroutine tstream_writeToFile (rtextstream, sfilename)

!<description>
  ! Writes the content of the whole stream to a text file.
!</description>

!<inputoutput>
  ! A text stream object which data should be written out. 
  type(t_textstream), intent(inout) :: rtextstream
!</inputoutput>

!<input>
  ! Filename of the file to write to.
  character(len=*), intent(in) :: sfilename
!</input>

!</subroutine>

    ! local variables
    integer :: i,ioldpos,ioldpos2,iunit
    character(len=TSTREAM_LENLINEBUF) :: sdata

    ! Try to open the file
    call io_openFileForWriting(sfilename, iunit, SYS_REPLACE)
    
    ! Oops...
    if (iunit .eq. -1) then
      call output_line ('Cannot open file for writing: '//trim(sfilename), &
                        OU_CLASS_ERROR,OU_MODE_STD,'tstream_writeToFile')
      call sys_halt()
    end if
  
    ! Write all lines to the file
    ioldpos = rtextstream%ipositionRead
    ioldpos2 = rtextstream%icurrentLine
    rtextstream%ipositionRead = 1
    rtextstream%icurrentLine = 1
    do 
      if (tstream_eof(rtextstream)) exit
      call tstream_readLine(rtextstream,sdata)
      write (iunit,'(A)') trim(sdata)
    end do
    rtextstream%ipositionRead = ioldpos
    rtextstream%icurrentLine = ioldpos2
    
    ! Close the file, finish.
    close (iunit)
    
  end subroutine

! ****************************************************************************************

end module
