!##############################################################################
!# ****************************************************************************
!# <name> heatcond_timeloop </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the actual nonstationary solver for the heat conduction
!# problem, i.e. the timeloop and time stepping. The following routines can
!# be found here:
!#
!# 1.) hc5_initparameters
!#     -> Initialise the parameters for the time stepping.
!#
!# 2.) hc5_timestep
!#     -> Calculate the solution of the next timestep.
!#
!# 3.) hc5_postprocessing
!#     -> Perform postprocessing (write GMV's,...)
!#
!# 4.) hc5_timeloop
!#     Start the nonstationary solver, perform the time stepping.
!# </purpose>
!##############################################################################

MODULE heatcond_timeloop

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE sortstrategy
  USE coarsegridcorrection
  USE ucd
  USE timestepping
  USE genoutput
  
  USE collection
  USE paramlist
    
  USE heatcond_callback
  
  USE heatcond_basic
  USE heatcond_matvec
  USE heatcond_boundarycondition
  USE heatcond_partridiscr
  USE heatcond_solver
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_initparameters (rparams,rproblem)
  
!<description>
  ! Reads the DAT file from disc into the parameter list rparams and
  ! initialises basic variables (number of levels, time stepping technique)
  ! in rproblem according to these settings.
!</description>

!<inputoutput>
  ! A parameter list structure accepting the parameters from the DAT file.
  TYPE(t_parlist), INTENT(INOUT) :: rparams

  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: cscheme,niterations
  REAL(DP) :: dtheta,dtstep,dtimemin,dtimemax

    ! Read the parameters from disc and put a reference to it
    ! to the collection
    CALL parlst_readfromfile(rparams, 'data/heatcond.dat')

    ! We want to solve our Laplace problem on level...
    CALL parlst_getvalue_int (rparams, 'GENERAL', 'NLMIN', rproblem%ilvmin, 7)
    CALL parlst_getvalue_int (rparams, 'GENERAL', 'NLMAX', rproblem%ilvmax, 7)
    
    ! Get the parameters for the time stepping scheme from the parameter list
    CALL parlst_getvalue_int (rparams, 'TIMESTEPPING', 'CSCHEME', cscheme, 0)
    CALL parlst_getvalue_int (rparams, 'TIMESTEPPING', 'NITERATIONS', niterations, 1000)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTHETA', dtheta, 1.0_DP)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTSTEP', dtstep, 0.1_DP)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTIMEMIN', dtimemin, 0.0_DP)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTIMEMAX', dtimemax, 1.0_DP)
    
    ! Initialise the time stepping in the problem structure
    CALL timstp_init (rproblem%rtimedependence%rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)
                     
    rproblem%rtimedependence%niterations = niterations
    
    rproblem%rtimedependence%dtimemin = dtimemin
    rproblem%rtimedependence%dtime = dtimemin
    rproblem%rtimedependence%dtimemax = dtimemax

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_timestep (rproblem,rvector,rrhs)
  
!<description>
  ! Performs one time step: $t^n -> t^n+1$. 
  ! Assembles system matrix and RHS vector. 
  ! Solves the corresponding time-step equation and returns the solution vector
  ! at the end of the time step.
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>
!</subroutine>

  ! local variables
    INTEGER :: ilvmin,ilvmax
    INTEGER :: i

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs
    TYPE(t_vectorBlock), TARGET :: rtempBlock

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rsmoother
    TYPE(t_linsolNode), POINTER :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1:rproblem%ilvmax) :: Rmatrices
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
    
    ! We have an equation of the type
    !
    !   d/dt u(x,t)  +  N(u(x,t))  =  f(x,t)
    !
    ! Which is discretised in time with a Theta scheme, leading to
    !
    !   $$ u_{n+1} + w_1*N(u_n+1) 
    !      =   u_n + w_2*N(u_n)  +  w_3*f_{n+1}  +  w_4*f_n $$
    !
    ! with k=time step size, u_{n+1} = u(.,t_{n+1}),etc., c.f. timestepping.f90.
    !
    ! The RHS of that equation therefore contains parts of the solution
    ! u_n, of the old RHS f_n and the new RHS f_{n+1}. At first, we make
    ! a weighted copy of the current RHS f_n to the 'global' RHS vector
    ! according to the time stepping scheme.
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
    p_rrhs    => rproblem%rrhs 
    
    ! Create a temporary vector we need for some preparations.
    CALL lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .FALSE.)
    
    ! Set up w_2*N(u_n) + w_4*f_n.
    
    CALL lsysbl_vectorLinearComb(rrhs,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightOldRHS,0.0_DP)
    
    ! Synchronise the sorting of the vectors according to the system matrix.
    ! We use the first subvector of rtempBlock as temporary data; it's
    ! large enough, as we have only one block.
    CALL lsysbl_synchroniseSortMatVec (p_rmatrix,p_rrhs,rtempBlock%RvectorBlock(1))
    CALL lsysbl_synchroniseSortMatVec (p_rmatrix,rvector,rtempBlock%RvectorBlock(1))
    
    CALL lsysbl_blockMatVec(rproblem%RlevelInfo(ilvmax)%rmatrixStatic,&
         rvector,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightMatrixRHS,&
         rproblem%rtimedependence%rtimestepping%dweightOldRHS)
         
    ! Add u_n -- or, more precisely, M u_n (with M being the mass matrix),
    ! since the discretisation with finite elements requires that.

    CALL lsysbl_blockMatVec(rproblem%RlevelInfo(ilvmax)%rmatrixMass,&
         rvector,p_rrhs,1.0_DP,1.0_DP)
         
    ! Switch to the next point in time.
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtime + &
          rproblem%rtimedependence%rtimestepping%dtstep
          
    ! Generate f_n+1 into the rrhs overwriting the previous RHS.
    CALL hc5_calcRHS (rproblem,rrhs)
    
    ! Add w_3*f_{n+1} to the current RHS. If necessary, unsort p_rrhs back before.
    CALL lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.FALSE.)
    
    CALL lsysbl_vectorLinearComb(rrhs,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightNewRHS,1.0_DP)

    ! That's it for the RHS vector.
    !
    ! The LHS "u_{n+1} + w_1*N(u_n+1)" results in the system matrix
    ! "M + w_1 N(.)" for the next linear system to solve. Set up that system
    ! on every level of the discretisation.
    
    DO i = ilvmin,ilvmax
      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixMass%RmatrixBlock(1,1),&
                                   rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                   LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
      CALL lsyssc_matrixLinearComb (rproblem%RlevelInfo(i)%rmatrixStatic%RmatrixBlock(1,1),&
                                    rproblem%rtimedependence%rtimestepping%dweightMatrixLHS,&
                                    rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    1.0_DP,&
                                    rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    .FALSE.,.FALSE.,.TRUE.,.TRUE.)
    END DO
    
    ! Discretise the boundary conditions at the new time instant
    CALL hc5_initDiscreteBC (rproblem)
          
    ! Implement boundary conditions into the RHS vector, the solution vector
    ! and the current system matrices.
    CALL hc5_implementBC (rproblem,rvector,p_rrhs,1.0_DP)
          
    ! Preparation of the linear system completed!
    !
    ! Attach the system matrices to the solver.
    
    p_rsolverNode => rproblem%p_rsolverNode
    
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    CALL linsol_setMatrices(p_rsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Initialise data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initData (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! Synchronise p_rrhs with the matrix so it's compatible to the linear system.
    CALL lsysbl_synchroniseSortMatVec (p_rmatrix,p_rrhs,rtempBlock%RvectorBlock(1))

    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    CALL linsol_solveAdaptively (p_rsolverNode,rvector,p_rrhs,rtempBlock)
    
    ! rvector is now u_n+1.
    !
    ! Release solver data.
    CALL linsol_doneData (p_rsolverNode)
    
    ! Unsort the vectors again in case they were resorted before calling 
    ! the solver.
    ! We use the first subvector of rtempBlock as temporary data; it's
    ! large enough, as we only have one block.
    CALL lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.FALSE.)
    CALL lsysbl_sortVectorInSitu (rvector,rtempBlock%RvectorBlock(1),.FALSE.)
    
    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)
    
    ! Finally tell the time stepping scheme that we completed the time step.
    CALL timstp_nextSubstep (rproblem%rtimedependence%rtimestepping)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_postprocessing (rproblem,rvector,iiteration,dtime)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<input>
  ! Number of current iteration
  INTEGER, INTENT(IN) :: iiteration
  
  ! Current simulation time
  REAL(DP), INTENT(IN) :: dtime
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! The current solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       'gmv/u5.gmv.'//TRIM(sys_si0L(iiteration,5)))
    CALL ucd_setSimulationTime (rexport,dtime)
    
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_timeloop (rproblem,rvector,rrhs)
  
!<description>
  ! Starts the time discretisation. Proceeds in time until the final time
  ! or the maximum number of time steps is reached.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! The initial solution vector. Is replaced by the final solution vector.
  ! Must be unsorted.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
  
  ! The initial RHS vector. Is replaced by the final RHS vector.
  ! Must be unsorted and without any boundary conditions implemented.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  INTEGER :: iiteration
  REAL(DP) :: dtime
  
    ! Let's start the timeloop
  
    iiteration = 1
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimemin
    DO WHILE ((iiteration .LE. rproblem%rtimedependence%niterations) .AND. &
              (rproblem%rtimedependence%dtime .LT. rproblem%rtimedependence%dtimemax))
              
      rproblem%rtimedependence%iiteration = iiteration
      
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line ('Time step '//TRIM(sys_siL(iiteration,6))// &
                        '     Time '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)))
      CALL output_lbrk ()
              
      ! Proceed to the next time step
      CALL hc5_timestep (rproblem,rvector,rrhs)
      
      ! Postprocessing. Write out the solution.
      CALL hc5_postprocessing (rproblem,rvector,iiteration,&
           rproblem%rtimedependence%dtime)
           
      iiteration = iiteration+1
    END DO

  END SUBROUTINE

END MODULE
