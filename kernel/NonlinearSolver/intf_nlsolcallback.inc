  ! This file describes the interfaces of the callback routines of the nonlinear
  ! defect-correction loop solver.

  ! ***************************************************************************  
  
  INTERFACE 
  
    SUBROUTINE fcb_getDefect (ite,rx,rb,rd,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd. The routine accepts a pointer to a collection structure 
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
    INTEGER, INTENT(IN)                           :: ite
  
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled by the callback routine
    ! with data.
    TYPE(t_vectorBlock), INTENT(INOUT), TARGET    :: rd
  !</inputoutput>
  
    END SUBROUTINE
    
  END INTERFACE
  
  ! ***************************************************************************

  INTERFACE 
  
    SUBROUTINE fcb_precondDefect (rd,rx,rb,domega,bsuccess,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd. The routine accepts a pointer to a collection structure 
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<inputoutput>
    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    TYPE(t_vectorBlock), INTENT(INOUT), TARGET    :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    REAL(DP), INTENT(INOUT)                       :: domega
    
    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    LOGICAL, INTENT(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Current right hand side of the nonlinear system
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb
  !</input>

    END SUBROUTINE

  END INTERFACE
  
  ! ***************************************************************************
  
  INTERFACE 
  
    SUBROUTINE fcb_resNormCheck (ite,rx,rb,rd,bconvergence,bdivergence,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! Residual norm calculation & printing routine.
    ! This routine is called each time the norm of the residuum was calculated.
    ! It has to check the current residuum for convergence and/or divergence
    ! and can print the residuum to screen.
  !</description>

  !<inputoutput>
    ! Number of current iteration. Is set to 0 when the callback routine
    ! is called the first time. In this situation, rd describes the initial
    ! defect vector.
    INTEGER, INTENT(IN)                           :: ite
  
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb

    ! Defect vector b-A(x)x.
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
  !</inputoutput>
  
  !<output>
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is within a desired tolerance, so that the solver should treat
    ! the iteration as 'converged'.
    LOGICAL, INTENT(OUT)                        :: bconvergence
  
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is out of a desired tolerance, so that the solver should treat
    ! the iteration as 'diverged'.
    LOGICAL, INTENT(OUT)                        :: bdivergence
  !</output>
  
    END SUBROUTINE

  END INTERFACE
