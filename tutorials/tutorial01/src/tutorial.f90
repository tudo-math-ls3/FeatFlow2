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

  implicit none
  
  ! =================================================================
  ! Main program
  ! =================================================================
  
  ! -----------------------------------------------------------------
  ! Initialisation of the feat library, output system and 
  ! memory management
  ! -----------------------------------------------------------------
  call system_init()
  call output_init ("")
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

  ! -----------------------------------------------------------------
  ! Clean up
  ! -----------------------------------------------------------------
  call storage_done()
  call output_done()

end program
