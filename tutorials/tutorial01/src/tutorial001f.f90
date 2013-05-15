!##############################################################################
!# Tutorial 001f: Command line arguments
!##############################################################################

module tutorial001f

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  implicit none
  private
  
  public :: start_tutorial001f

contains

  ! ***************************************************************************

  subroutine start_tutorial001f
  
    integer :: iarg,nargc
    character(len=SYS_STRLEN) :: sargv, soption, svalue
    integer :: iformat

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 001f.")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Get the number of command line
    ! parameters and print them.
    ! =================================
    nargc = sys_ncommandLineArgs()

    call output_line ("ARGC = " // trim(sys_siL(nargc,10)) )
    call output_line ("ARGV =", bnolinebreak = (nargc .gt. 0) )
    
    do iarg=1,nargc
      ! Note that there is also an extended version of "sys_getcommandLineArg"
      ! which automatically parses options line "--option=value".
      call sys_getcommandLineArg(iarg,sargv)

      call output_line (" " // trim (sargv), bnolinebreak = (iarg .ne. nargc) )
    end do
    
    ! =================================
    ! Parse some prescribed parameters
    ! =================================
    
    ! -----------------------
    ! Single option, svalue not given.
    call sys_parseCommandLineArg ("--file=myfile.txt",soption)
    
    call output_lbrk()
    call output_line ("Parse-Test: '--file=myfile.txt'")
    call output_line ("   Option = " // trim(soption) )

    ! -----------------------
    ! Single option
    call sys_parseCommandLineArg ("'This is a string.'",soption,svalue,iformat)
    
    call output_lbrk()
    call output_line ("Parse-Test: '""This is a string.""'")
    call output_line ("   Option = " // trim(soption) )
    call output_line ("   Value  = " // trim(svalue) )
    call output_line ("   Format = " // trim(sys_siL(iformat,10)) )

    ! -----------------------
    ! Single option
    call sys_parseCommandLineArg ("-flag1",soption,svalue,iformat)
    
    call output_lbrk()
    call output_line ("Parse-Test: '-flag1'")
    call output_line ("   Option = " // trim(soption) )
    call output_line ("   Value  = " // trim(svalue) )
    call output_line ("   Format = " // trim(sys_siL(iformat,10)) )

    ! -----------------------
    ! Single option
    call sys_parseCommandLineArg ("--flag2",soption,svalue,iformat)
    
    call output_lbrk()
    call output_line ("Parse-Test: '--flag2'")
    call output_line ("   Option = " // trim(soption) )
    call output_line ("   Value  = " // trim(svalue) )
    call output_line ("   Format = " // trim(sys_siL(iformat,10)) )

    ! -----------------------
    ! Option with value
    call sys_parseCommandLineArg ("-file=myfile.txt",soption,svalue,iformat)
    
    call output_lbrk()
    call output_line ("Parse-Test: '-file=myfile.txt'")
    call output_line ("   Option = " // trim(soption) )
    call output_line ("   Value  = " // trim(svalue) )
    call output_line ("   Format = " // trim(sys_siL(iformat,10)) )

    ! -----------------------
    ! Option with value
    call sys_parseCommandLineArg ("--file=myfile.txt",soption,svalue,iformat)
    
    call output_lbrk()
    call output_line ("Parse-Test: '--file=myfile.txt'")
    call output_line ("   Option = " // trim(soption) )
    call output_line ("   Value  = " // trim(svalue) )
    call output_line ("   Format = " // trim(sys_siL(iformat,10)) )
    
  end subroutine

end module
