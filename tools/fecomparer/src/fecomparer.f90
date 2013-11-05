!####################################################################
!# Comparison of FE functions.
!####################################################################

program fecomparer

  ! Include basic Feat-2 modules
  use fparser
  use fsystem
  use genoutput
  use storage
  
  use fecompbase

  implicit none
  
  ! Local variables
  type(t_problem), pointer :: p_rproblem
  
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
  call fparser_init()
  
  ! -----------------------------------------------------------------
  ! Basic initialisations
  ! -----------------------------------------------------------------
  ! Create a problem structure
  allocate(p_rproblem)

  ! Read parameters, basic initialisation
  call base_init (p_rproblem)
  
  ! Create the mesh hierarchy
  call base_createMeshHierarchy (p_rproblem)
  
  ! Allocate and read solution vectors
  call base_allocateSolutions (p_rproblem)
  call base_readSolutions (p_rproblem)

  ! -----------------------------------------------------------------
  ! Do the computation
  ! -----------------------------------------------------------------
  
  select case (p_rproblem%ccomputation)
  
  case (0)
    ! Copy solutions 1 to 0.Postprocessing with 0.
    !
    ! 1.) Copy 1->0, Type cast.
    call base_vectorcopy (p_rproblem,1,0)
  
  case (2)
    ! Copy solutions 3/4 to 1/2. Calculate the difference 
    ! between 1 and 2 to 0. Postprocessing with 0.
    !
    ! 1.) Copy 3->1, Type cast
    call base_vectorcopy (p_rproblem,3,1)
  
    ! 2.) Copy 4->2, Type cast
    call base_vectorcopy (p_rproblem,4,2)
    
    ! 3.) Linear combination
    call base_vectorLinearComb (p_rproblem,2,1,0,-1.0_DP,1.0_DP)

  case (3)
    ! Copy solutions 3/4 to 1/2. Calculate the difference 
    ! between 1 and 2 to 0. Divide by 1. Postprocessing with 0.
    !
    ! 1.) Copy 3->1, Type cast
    call base_vectorcopy (p_rproblem,3,1)
  
    ! 2.) Copy 4->2, Type cast
    call base_vectorcopy (p_rproblem,4,2)
    
    ! 3.) Linear combination
    call base_vectorLinearComb (p_rproblem,2,1,0,-1.0_DP,1.0_DP)

    ! 4.) Divide by 1
    call base_vectorDivide (p_rproblem,0,1)

  end select
  
  ! -----------------------------------------------------------------
  ! Compute errors
  ! -----------------------------------------------------------------

  call base_computeErrors (p_rproblem)

  ! -----------------------------------------------------------------
  ! Write the output
  ! -----------------------------------------------------------------

  call base_writeOutput (p_rproblem)

  ! -----------------------------------------------------------------
  ! Clean up
  ! -----------------------------------------------------------------
  
  call base_done (p_rproblem)

  ! -----------------------------------------------------------------
  ! Clean up library
  ! -----------------------------------------------------------------
  call storage_done()
  call output_done()

end program
