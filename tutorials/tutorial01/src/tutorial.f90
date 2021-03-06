!####################################################################
!# Collection of tutorials.
!#
!# The main routine calls all tutorial applications.
!####################################################################

program tutorial

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  ! Include the tutorial modules
  use tutorial001a
  use tutorial001b
  use tutorial001c
  use tutorial001d
  use tutorial001e
  use tutorial001f
  
  use tutorial002a
  use tutorial002b
  
  use tutorial003a
  use tutorial003b
  use tutorial003c
  use tutorial003d
  use tutorial003e
  
  use tutorial004a
  use tutorial004b
  use tutorial004c
  use tutorial004d
  use tutorial004e
  use tutorial004f
  use tutorial004g
  use tutorial004h
  use tutorial004i
  use tutorial004j
  use tutorial004k
  use tutorial004l
  use tutorial004m
  
  use tutorial005a
  use tutorial005b
  use tutorial005c
  
  use tutorial006a
  use tutorial006b
  use tutorial006c
  use tutorial006d
  use tutorial006e
  use tutorial006f
  use tutorial006g
  use tutorial006h
  use tutorial006i
  use tutorial006j
  use tutorial006k
  use tutorial006l
  use tutorial006m
  
  use tutorial007a
  use tutorial007b
  use tutorial007c
  
  use tutorial008a
  use tutorial008b
  use tutorial008c
  use tutorial008d
  use tutorial008e

  use tutorial009a
  use tutorial009b
  use tutorial009c
  
  use tutorial010a
  use tutorial010b
  use tutorial010c
  use tutorial010d
  use tutorial010e
  use tutorial010f
  use tutorial010g
  
  use tutorial011a
  use tutorial011b
  use tutorial011c
  use tutorial011d
  
  use tutorial012a
  use tutorial012b
  use tutorial012c
  use tutorial012d
  
  use tutorial013a
  use tutorial013b
  use tutorial013c
  
  use tutorial014a
  use tutorial014b
  use tutorial014c
  
  use tutorial015a
  use tutorial015b
  use tutorial015c
  
  use tutorial016a
  use tutorial016b
  use tutorial016c
  use tutorial016d
  use tutorial016e
  use tutorial016f
  
  use tutorial017a
  use tutorial017b
  use tutorial017c
  use tutorial017d
  use tutorial017e
  
  use tutorial018a
  use tutorial018b
  use tutorial018c
  use tutorial018d
  use tutorial018e
  use tutorial018f
  
  use tutorial019a
  use tutorial019b

  use tutorial020a
  use tutorial020b
  
  use tutorial021a
  use tutorial021b

  use tutorial022a
  use tutorial022b
  use tutorial022c
  use tutorial022d

  use tutorial023a
  use tutorial023b
  use tutorial023c
  use tutorial023d
  use tutorial023e
  use tutorial023f
  use tutorial023g

  implicit none
  
  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile

  ! =================================================================
  ! Main program
  ! =================================================================
  
  ! -----------------------------------------------------------------
  ! Initialisation of the feat library, output system and 
  ! memory management
  ! -----------------------------------------------------------------
  call sys_init()

  ! Initialise the output system.
  !
  ! Normally, we write all the output to the screen and to a file
  ! "./log/output.txt".
  ! In the case that environment variables "$logdir"/"$resultsfile" exists,
  ! we write all the output to that file. This can be used e.g. in
  ! regression tests to compare results to reference results.
  if (sys_getenv_string("LOGDIR",slogdir) .and. &
      sys_getenv_string("RESULTFILE",slogfile)) then
    call output_init (trim(slogdir)//"/"//trim(slogfile))
  else
    call output_init ("./log/output.txt")
  end if
  
  call storage_init(999, 100)

  ! -----------------------------------------------------------------
  ! Call the tutorial modules.
  ! -----------------------------------------------------------------
  call start_tutorial001a
  call start_tutorial001b
  call start_tutorial001c
  call start_tutorial001d
  call start_tutorial001e
  call start_tutorial001f
  
  call start_tutorial002a
  call start_tutorial002b
  
  call start_tutorial003a
  call start_tutorial003b
  call start_tutorial003c
  call start_tutorial003d
  call start_tutorial003e
  
  call start_tutorial004a
  call start_tutorial004b
  call start_tutorial004c
  call start_tutorial004d
  call start_tutorial004e
  call start_tutorial004f
  call start_tutorial004g
  call start_tutorial004h
  call start_tutorial004i
  call start_tutorial004j
  call start_tutorial004k
  call start_tutorial004l
  call start_tutorial004m
  
  call start_tutorial005a
  call start_tutorial005b
  call start_tutorial005c
  
  call start_tutorial006a
  call start_tutorial006b
  call start_tutorial006c
  call start_tutorial006d
  call start_tutorial006e
  call start_tutorial006f
  call start_tutorial006g
  call start_tutorial006h
  call start_tutorial006i
  call start_tutorial006j
  call start_tutorial006k
  call start_tutorial006l
  call start_tutorial006m
  
  call start_tutorial007a
  call start_tutorial007b
  call start_tutorial007c
  
  call start_tutorial008a
  call start_tutorial008b
  call start_tutorial008c
  call start_tutorial008d
  call start_tutorial008e
  
  call start_tutorial009a
  call start_tutorial009b
  call start_tutorial009c

  call start_tutorial010a
  call start_tutorial010b
  call start_tutorial010c
  call start_tutorial010d
  call start_tutorial010e
  call start_tutorial010f
  call start_tutorial010g
  
  call start_tutorial011a
  call start_tutorial011b
  call start_tutorial011c
  call start_tutorial011d

  call start_tutorial012a
  call start_tutorial012b
  call start_tutorial012c
  call start_tutorial012d

  call start_tutorial013a
  call start_tutorial013b
  call start_tutorial013c

  call start_tutorial014a
  call start_tutorial014b
  call start_tutorial014c

  call start_tutorial015a
  call start_tutorial015b
  call start_tutorial015c

  call start_tutorial016a
  call start_tutorial016b
  call start_tutorial016c
  call start_tutorial016d
  call start_tutorial016e
  call start_tutorial016f
  
  call start_tutorial017a
  call start_tutorial017b
  call start_tutorial017c
  call start_tutorial017d
  call start_tutorial017e

  call start_tutorial018a
  call start_tutorial018b
  call start_tutorial018c
  call start_tutorial018d
  call start_tutorial018e
  call start_tutorial018f
  
  call start_tutorial019a
  call start_tutorial019b

  call start_tutorial020a
  call start_tutorial020b

  call start_tutorial021a
  call start_tutorial021b

  call start_tutorial022a
  call start_tutorial022b
  call start_tutorial022c
  call start_tutorial022d

  call start_tutorial023a
  call start_tutorial023b
  call start_tutorial023c
  call start_tutorial023d
  call start_tutorial023e
  call start_tutorial023f
  call start_tutorial023g

  ! -----------------------------------------------------------------
  ! Clean up
  ! -----------------------------------------------------------------
  call storage_done()
  call output_done()

end program
