  ! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines of the nonlinear
  ! defect-correction loop solver.

  ! ***************************************************************************
  
  interface
  
    subroutine fcb_getDefect (ite,rx,rb,rd,p_rcollection)
  
    use fsystem
    use linearsystemblock
    use collection
    
  !<description>
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
    integer, intent(in)                           :: ite
  
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled by the callback routine
    ! with data.
    type(t_vectorBlock), intent(inout), target    :: rd
  !</inputoutput>
  
    end subroutine
    
  end interface
  
  ! ***************************************************************************

  interface
  
    subroutine fcb_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
  
    use fsystem
    use linearsystemblock
    use collection
    
  !<description>
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<inputoutput>
    ! Number of current iteration.
    integer, intent(in)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(inout), target    :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it is changed again.
    real(DP), intent(inout)                       :: domega
    
    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(inout)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Current right hand side of the nonlinear system
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>

    end subroutine

  end interface
  
  ! ***************************************************************************
  
  interface
  
    subroutine fcb_resNormCheck (ite,rx,rb,rd,bconvergence,bdivergence,p_rcollection)
  
    use linearsystemblock
    use collection
    
  !<description>
    ! Residual norm calculation & printing routine.
    ! This routine is called each time the norm of the residuum was calculated.
    ! It has to check the current residuum for convergence and/or divergence
    ! and can print the residuum to screen.
  !</description>

  !<inputoutput>
    ! Number of current iteration.
    integer, intent(in)                           :: ite
  
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(in), target       :: rb

    ! Defect vector b-A(x)x.
    type(t_vectorBlock), intent(in), target       :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection
  !</inputoutput>
  
  !<output>
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is within a desired tolerance, so that the solver should treat
    ! the iteration as 'converged'.
    logical, intent(out)                        :: bconvergence
  
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is out of a desired tolerance, so that the solver should treat
    ! the iteration as 'diverged'.
    logical, intent(out)                        :: bdivergence
  !</output>
  
    end subroutine

  end interface
