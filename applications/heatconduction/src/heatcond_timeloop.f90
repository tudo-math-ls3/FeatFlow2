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
!#     -> Perform postprocessing (write VTK`s,...)
!#
!# 4.) hc5_timeloop
!#     Start the nonstationary solver, perform the time stepping.
!# </purpose>
!##############################################################################

module heatcond_timeloop

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use sortstrategy
  use coarsegridcorrection
  use ucd
  use timestepping
  use genoutput
  
  use collection
  use paramlist
    
  use heatcond_callback
  
  use heatcond_basic
  use heatcond_matvec
  use heatcond_boundarycondition
  use heatcond_partridiscr
  use heatcond_solver
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_timestep (rproblem,rvector,rrhs)
  
!<description>
  ! Performs one time step: $t^n -> t^n+1$.
  ! Assembles system matrix and RHS vector.
  ! Solves the corresponding time-step equation and returns the solution vector
  ! at the end of the time step.
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  type(t_vectorBlock), intent(inout) :: rvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>
!</subroutine>

  ! local variables
    integer :: ilvmin,ilvmax
    integer :: i

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs
    type(t_vectorBlock), target :: rtempBlock

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1:rproblem%ilvmax) :: Rmatrices
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    type(t_linsolMGLevelInfo), pointer :: p_rlevelInfo
    
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
    call lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .false.)
    
    ! Set up w_2*N(u_n) + w_4*f_n.
    
    call lsysbl_vectorLinearComb(rrhs,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightOldRHS,0.0_DP)
    
    ! Synchronise the sorting of the vectors according to the system matrix.
    ! We use the first subvector of rtempBlock as temporary data; it is
    ! large enough, as we have only one block.
    call lsysbl_synchroniseSortMatVec (p_rmatrix,p_rrhs,rtempBlock%RvectorBlock(1))
    call lsysbl_synchroniseSortMatVec (p_rmatrix,rvector,rtempBlock%RvectorBlock(1))
    
    call lsysbl_blockMatVec(rproblem%RlevelInfo(ilvmax)%rmatrixStatic,&
         rvector,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightMatrixRHS,&
         rproblem%rtimedependence%rtimestepping%dweightOldRHS)
         
    ! Add u_n -- or, more precisely, M u_n (with M being the mass matrix),
    ! since the discretisation with finite elements requires that.

    call lsysbl_blockMatVec(rproblem%RlevelInfo(ilvmax)%rmatrixMass,&
         rvector,p_rrhs,1.0_DP,1.0_DP)
         
    ! Switch to the next point in time.
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtime + &
          rproblem%rtimedependence%rtimestepping%dtstep
          
    ! Generate f_n+1 into the rrhs overwriting the previous RHS.
    call hc5_calcRHS (rproblem,rrhs)
    
    ! Add w_3*f_{n+1} to the current RHS. If necessary, unsort p_rrhs back before.
    call lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.false.)
    
    call lsysbl_vectorLinearComb(rrhs,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightNewRHS,1.0_DP)

    ! That is it for the RHS vector.
    !
    ! The LHS "u_{n+1} + w_1*N(u_n+1)" results in the system matrix
    ! "M + w_1 N(.)" for the next linear system to solve. Set up that system
    ! on every level of the discretisation.
    
    do i = ilvmin,ilvmax
      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixMass%RmatrixBlock(1,1),&
                                   rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                   LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
      call lsyssc_matrixLinearComb (rproblem%RlevelInfo(i)%rmatrixStatic%RmatrixBlock(1,1),&
                                    rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    rproblem%rtimedependence%rtimestepping%dweightMatrixLHS,&
                                    1.0_DP,.false.,.false.,.true.,.true.)
    end do
    
    ! Discretise the boundary conditions at the new time instant
    call hc5_initDiscreteBC (rproblem)
          
    ! Implement boundary conditions into the RHS vector, the solution vector
    ! and the current system matrices.
    call hc5_implementBC (rproblem,rvector,p_rrhs,1.0_DP)
          
    ! Preparation of the linear system completed!
    !
    ! Attach the system matrices to the solver.
    
    p_rsolverNode => rproblem%p_rsolverNode
    
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    call linsol_setMatrices(p_rsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Initialise data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initData (p_rsolverNode,ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    
    ! Synchronise p_rrhs with the matrix so it is compatible to the linear system.
    call lsysbl_synchroniseSortMatVec (p_rmatrix,p_rrhs,rtempBlock%RvectorBlock(1))

    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rvector,p_rrhs,rtempBlock)
    
    ! rvector is now u_n+1.
    !
    ! Release solver data.
    call linsol_doneData (p_rsolverNode)
    
    ! Unsort the vectors again in case they were resorted before calling
    ! the solver.
    ! We use the first subvector of rtempBlock as temporary data; it is
    ! large enough, as we only have one block.
    call lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.false.)
    call lsysbl_sortVectorInSitu (rvector,rtempBlock%RvectorBlock(1),.false.)
    
    ! Release the temporary vector
    call lsysbl_releaseVector (rtempBlock)
    
    ! Finally tell the time stepping scheme that we completed the time step.
    call timstp_nextSubstep (rproblem%rtimedependence%rtimestepping)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_postprocessing (rproblem,rvector,iiteration,dtime)
  
!<description>
  ! Writes the solution into a VTK file.
!</description>

!<input>
  ! Number of current iteration
  integer, intent(in) :: iiteration
  
  ! Current simulation time
  real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! The current solution vector.
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    type(t_triangulation), pointer :: p_rtriangulation

    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       trim(rproblem%sucddir)//'/u5.vtk.'//trim(sys_si0L(iiteration,5)))
    call ucd_setSimulationTime (rexport,dtime)
    
    call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_timeloop (rproblem,rvector,rrhs)
  
!<description>
  ! Starts the time discretisation. Proceeds in time until the final time
  ! or the maximum number of time steps is reached.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem

  ! The initial solution vector. Is replaced by the final solution vector.
  ! Must be unsorted.
  type(t_vectorBlock), intent(inout) :: rvector
  
  ! The initial RHS vector. Is replaced by the final RHS vector.
  ! Must be unsorted and without any boundary conditions implemented.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

  integer :: iiteration
  real(DP) :: dtime
  
    ! Let us start the timeloop
  
    iiteration = 1
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimemin
    do while ((iiteration .le. rproblem%rtimedependence%niterations) .and. &
              (rproblem%rtimedependence%dtime .lt. rproblem%rtimedependence%dtimemax))
              
      rproblem%rtimedependence%iiteration = iiteration
      
      call output_separator(OU_SEP_MINUS)
      call output_line ('Time step '//trim(sys_siL(iiteration,6))// &
                        '     Time '//trim(sys_sdL(rproblem%rtimedependence%dtime,5)))
      call output_lbrk ()
              
      ! Proceed to the next time step
      call hc5_timestep (rproblem,rvector,rrhs)
      
      ! Postprocessing. Write out the solution.
      call hc5_postprocessing (rproblem,rvector,iiteration,&
           rproblem%rtimedependence%dtime)
           
      iiteration = iiteration+1
    end do

  end subroutine

end module
