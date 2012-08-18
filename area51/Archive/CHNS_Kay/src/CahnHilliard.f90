!##############################################################################
!# ****************************************************************************
!# <name> CahnHilliard_timeloop </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the actual nonstationary solver for the heat conduction
!# problem, i.e. the timeloop and time stepping. The following routines can
!# be found here:
!#
!# 1.) CH_initparameters
!#     -> Initialise the parameters for the time stepping.
!#
!# 2.) CH_timestep
!#     -> Calculate the solution of the next timestep.
!#
!# 3.) CH_postprocessing
!#     -> Perform postprocessing (write GMV's,...)
!#
!# 4.) CH_timeloop
!#     Start the nonstationary solver, perform the time stepping.
!#
!# 5.) CH_initTimeSteppingScheme
!#     Initialize the time stepping scheme
!#
!# </purpose>
!##############################################################################

module CahnHilliard

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
  use nonlinearsolver
    
   use CahnHilliard_callback
   use CahnHilliard_nonlinearcore
   use CahnHilliard_basic
   use CahnHilliard_matvec
   use CahnHilliard_boundarycondition
   use CahnHilliard_partridiscr

!   ! use NS modules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   use ccbasic
   use cccallback
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  subroutine CH_timestep (rCHproblem,rCHvector,rCHrhs, rtimestepping, &
            rCHnonlinearIteration,rnlSolver, rNSproblem, rNSvector, rNSrhs)
!<description>
  ! Performs one time step: $t^n -> t^n+1$.
  ! Assembles system matrix and RHS vector.
  ! Solves the corresponding time-step equation and returns the solution vector
  ! at the end of the time step.
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem
  
  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  type(t_vectorBlock), intent(INOUT) :: rCHvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  type(t_vectorBlock), intent(INOUT) :: rCHrhs

! MCai, for dealing with nonlinear iteration.
! Structure for the nonlinear iteration for solving the core equation.
  type(t_CHnonlinearIteration), intent(INOUT) :: rCHnonlinearIteration

  type(t_explicitTimestepping), intent(IN) :: rtimestepping
  ! A configuration stucture for the nonlinear solver
  type(t_nlsolNode), intent(INOUT) :: rnlSolver

!~~~~~~~~~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(INOUT) :: rNSproblem
  type(t_vectorBlock), intent(INOUT) :: rNSvector
  type(t_vectorBlock), intent(INOUT) :: rNSrhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!</inputoutput>
!</subroutine>

    ! local variables

    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1), TARGET :: RfilterChain
    type(t_filterChain), dimension(:), POINTER :: p_RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), POINTER :: p_rmatrix

 ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), POINTER :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), POINTER :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1:rCHproblem%NLMAX) :: Rmatrices
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    type(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo

 
    type(t_CHnonlinearIteration) :: rnonlinearIterationTmp
    type(t_nonlinearCHMatrix) :: rnonlinearCHMatrix
    
    type(t_vectorBlock) :: rtempVectorRhs, rtempVector, rCHvectorTmp
    real(DP), dimension(:), POINTER :: p_Ddata,p_Ddata2

    real(DP), dimension(2) :: Dresiduals
    integer, dimension(2) :: Cnorms
  

	! The new RHS will be set up in rtempVectorRhs. Assign the discretisation/
    ! boundary conditions of rrhs to that vector so that rtempVectorRhs
    ! acts as a RHS vector.
    
    call lsysbl_createVecBlockIndirect (rCHrhs, rtempVectorRhs, .false.)
    call lsysbl_assignDiscrIndirect(rCHrhs,rtempVectorRhs)

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
    
    ! Set up w_4*f_n.
    call lsysbl_vectorLinearComb(rCHrhs,rtempVectorRhs,&
         rtimestepping%dweightOldRHS,0.0_DP)
    
! We delete all sort part.....
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! For setting up M(u_n) + w_2*N(u_n), switch the sign of w_2 and call the method
    ! to calculate the Convection/Diffusion part of the nonlinear defect. This builds
    ! rtempVectorRhs = rtempVectorRhs - (-Mass)*u - (-w_2) (nu*Laplace*u + grad(u)u).
	! w2=-(1-theta) delta_t

    ! Initialise the matrix assembly structure rnonlinearCHMatrix to describe the
    ! matrix we want to have.
    ! Note that only the first equation is time dependent. We set rnonlinearCHmatrix as follows.
    rnonlinearCHMatrix%dalpha = -1.0_DP    ! Q: take care: 1, 0, or -1 ?
    rnonlinearCHMatrix%dgamma = -rtimestepping%dweightMatrixRHS
    rnonlinearCHMatrix%deta =   -rtimestepping%dweightMatrixRHS
    rnonlinearCHMatrix%dbeta = 0.0_DP
    rnonlinearCHMatrix%dtau = 0.0_DP
    rnonlinearCHMatrix%dtheta = 0.0_DP

    rnonlinearCHMatrix%p_rdiscretisation => &
        rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rdiscretisation
    rnonlinearCHMatrix%p_rmatrixLaplace => &
        rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixLaplace
    rnonlinearCHMatrix%p_rmatrixMass => &
        rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixMass

!     rnonlinearCHMatrix%p_rmatrixA => &
!         rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixA
!     rnonlinearCHMatrix%p_rmatrixB => &
!         rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixB
!     rnonlinearCHMatrix%p_rmatrixC => &
!         rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixC
!     rnonlinearCHMatrix%p_rmatrixD => &
!         rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixD

    ! call CH_nonlinearMatMul (rnonlinearCHMatrix,rCHvector,rtempVectorRhs,-1.0_DP,1.0_DP)
!	call lsysbl_copyVector(rCHvector, rCHvectorTmp)
!	call lsyssc_scaleVector(rCHvectorTmp%RvectorBlock(2), 0.0_DP)

   ! rtempVectorRhs=rtempVectorRhs-nonCHMat * rCHvector
   ! Instead, we need to use the following code

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!MCai, take care: Concentration dependent mobility
!   call CH_generateVarLaplace (rCHproblem,rCHvector, &
!          rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixB)
!   rnonlinearCHMatrix%p_rmatrixB=>&
!          rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! MCai, how to generate nonlinear matrix?
!   call CH_generateNonlinearMat (rCHproblem,rCHvector, &
!         rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixD)
!   rnonlinearCHMatrix%p_rmatrixD=>&
!         rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixD

   call CH_nonlinearMatMul (rnonlinearCHMatrix,rCHvector,&
                         rtempVectorRhs, -1.0_DP,1.0_DP, rCHvector,&
				         rNSproblem, rNSvector)

       ! Debug
!        Cnorms(:) = LINALG_NORMMAX
!        Dresiduals(:) = lsysbl_vectorNormBlock (rtempVectorRhs,Cnorms)
!        print *, Dresiduals(1)
!        print *, Dresiduals(2)
!        stop

!   print *, 'finishing assemble nonlinear RHS, no problem'
    ! -------------------------------------------
    ! Switch to the next point in time.
    rCHproblem%rtimedependence%dtime = rCHproblem%rtimedependence%rtimestepping%dcurrenttime&
                            + rCHproblem%rtimedependence%rtimestepping%dtstep

!~~MCai, Need to take care~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    ! Discretise the boundary conditions at the new point in time --
!    ! if the boundary conditions are nonconstant in time!
!    if (rproblem%iboundary .ne. 0) then
!      call cc_updateDiscreteBC (rproblem)
! Instead, we may need
!      call CH_updateDiscreteBC (rCHproblem)
!    end if
! Suggestion, read the cc code to make sure whether we need it.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    print *, 'there is no collection'
    ! Generate f_n+1 into the rCHrhs overwriting the previous RHS.
    call CH_calcRHS (rCHproblem,rCHrhs, &
	            rNSproblem, rNSvector, rNSrhs)
!MCai, we need to assemble the convetive term in the right hand side: take care

! Debug,
!         Cnorms(:) = LINALG_NORMMAX
!         Dresiduals(:) = lsysbl_vectorNormBlock (rCHrhs,Cnorms)
!         print *, Dresiduals(1)
!         print *, Dresiduals(2)
!         stop

    ! Add w_3*f_{n+1} to the rtempVectorRhs
    call lsysbl_vectorLinearComb(rCHrhs,rtempVectorRhs,&
         rCHproblem%rtimedependence%rtimestepping%dweightNewRHS,1.0_DP)
    ! That's it for the RHS vector.
    !
    ! The LHS "u_{n+1} + w_1*N(u_n+1)" results in the system matrix
    ! "M + w_1 N(.)" for the next linear system to solve. Set up that system
    ! on every level of the discretisation.

!~~~~~~~~~~~~~Shall we implement BC in solver core equation?~~~~~~~~~~~~~~~~~~~
!    ! Implement boundary conditions into the RHS vector, the solution vector
!    ! and the current system matrices.
     call CH_implementBC (rCHproblem,rCHvector,rtempVectorRhs,1.0_DP)

    rnonlinearIterationTmp = rCHnonlinearIteration
    rnonlinearIterationTmp%dalpha = 1.0_DP
	! A = \alpha M + \gamma Conv
    rnonlinearIterationTmp%dgamma = rCHproblem%rtimedependence%rtimestepping%dweightMatrixLHS
    rnonlinearIterationTmp%deta   = rCHproblem%rtimedependence%rtimestepping%dweightMatrixLHS
    rnonlinearIterationTmp%dbeta  = -1.0_DP
	! D = \tau N(.) + \theta Lap
    rnonlinearIterationTmp%dtau   = 1.0_DP
    rnonlinearIterationTmp%dtheta = 1.0_DP
    
    ! Update the preconditioner for the case, something changed (e.g.
    ! the boundary conditions).
    ! Note: The bstructuralChange-parameter is set to FALSE here.
    ! In case the template matrices changed (e.g. during a mesh adaption),
    ! the routine must be called with bstructuralChange=true!
	! For Cahn-Hilliard-Navier-Stokes eqn, the preconditioner also changes with time

!MCai, double check how to update preconditioner? Because preconditioner may need
! rCHvector in previous step. How to update preconditioner? Or preconditioner is fixed?
! MCai, in cc2d, the preconditioner does not change with time, but here we need
! the preconditioner changes with time: we also need rCHvector to assemble global sys.

! We need nonlinear solver here.
!MCai, we need to generate system matrix each time step, (it is different from heat eqn)
! we use updatePreconditioner to do this, CH_assembleMatrix need velocity vector...
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    print *, 'before update preconditioner'
    call CH_updatePreconditioner (rCHproblem,rnonlinearIterationTmp,&
                            rCHvector,rtempVectorRhs,.false.,.true.,&
							rNSproblem, rNSvector)
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! call the solver of the core equation to solve it using a nonlinear
    ! iteration.
 !   print *, 'before solve core eqn'
    call CH_solveCoreEquation (rCHproblem,rnonlinearIterationTmp,rnlSolver,&
                          rCHvector,rtempVectorRhs,rtempVector, &
                          rNSproblem, rNSvector, rNSrhs)
!-------------------------------------------------------------------------------
    ! rCHvector is now u_n+1.
    call lsysbl_releaseVector (rtempVectorRhs)
    call lsysbl_releaseVector (rtempVector)
!    call lsysbl_releaseVector (rCHvectorTmp)
    
    ! Finally tell the time stepping scheme that we completed the time step.
    call timstp_nextSubstep (rCHproblem%rtimedependence%rtimestepping)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine CH_initTimeSteppingScheme (rparams,rtimeStepping)
  
!<description>
  ! Initialises the time stepping scheme according to the parameters in the DAT file.
!</description>

!<input>
  ! A parameter list from a DAT file containing parameters that configure the
  ! time dependence.
  type(t_parlist), intent(IN) :: rparams
!</input>

!<inputoutput>
  ! The time stepping scheme structure to be initialised.
  type(t_explicitTimeStepping), intent(INOUT) :: rtimeStepping
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: cscheme
    real(DP) :: dtheta,dtstep,dtimemin

    ! Get the parameters for the time stepping scheme from the parameter list
    call parlst_getvalue_int (rparams, &
        'TIME-DISCRETISATION', 'ITIMESTEPSCHEME', cscheme, 0)
    call parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMESTEPTHETA', dtheta, 1.0_DP)
    call parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMESTEP', dtstep, 0.1_DP)
    call parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMEINIT', dtimemin, 0.0_DP)
    
    ! Initialise the time stepping in the problem structure
    call timstp_init (rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine CH_initparameters (rCHparams,rCHproblem)
  
!<description>
  ! Reads the DAT file from disc into the parameter list rCHparams and
  ! initialises basic variables (number of levels, time stepping technique)
  ! in rCHproblem according to these settings.
!</description>

!<inputoutput>
  ! A parameter list structure accepting the parameters from the DAT file.
  type(t_parlist), intent(INOUT) :: rCHparams

  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: cscheme,niterations
  real(DP) :: dtheta,dtstep,dtimemin,dtimemax

    ! Read the parameters from disc and put a reference to it
    ! to the collection
    ! Get the parameters for the time stepping scheme from the parameter list
    call  parlst_readfromfile(rCHparams, 'data/timediscr.dat')
    call parlst_getvalue_int (rCHparams, 'TIME-DISCRETISATION', 'itimeStepScheme', cscheme, 0)
    call parlst_getvalue_int (rCHparams, 'TIME-DISCRETISATION ', 'NITERATIONS', niterations, 1000)
    call parlst_getvalue_double (rCHparams, 'TIME-DISCRETISATION', 'DTIMESTEPTHETA', dtheta, 1.0_DP)
    call parlst_getvalue_double (rCHparams, 'TIME-DISCRETISATION', 'dtimeStep', dtstep, 0.1_DP)
    call parlst_getvalue_double (rCHparams, 'TIME-DISCRETISATION', 'DTIMEINIT', dtimemin, 0.0_DP)
    call parlst_getvalue_double (rCHparams, 'TIME-DISCRETISATION', 'DTIMEMAX', dtimemax, 1.0_DP)

    call  parlst_readfromfile(rCHparams, 'data/discretisation.dat')
    ! We want to solve our Laplace problem on level...
    call parlst_getvalue_int (rCHparams, 'CC-DISCRETISATION', 'NLMIN', rCHproblem%NLMIN, 7)
    call parlst_getvalue_int (rCHparams, 'CC-DISCRETISATION', 'NLMAX', rCHproblem%NLMAX, 7)
    
    ! Initialise the time stepping in the problem structure
    call timstp_init (rCHproblem%rtimedependence%rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)

    rCHproblem%rtimedependence%niterations = niterations

    rCHproblem%rtimedependence%dtimemin = dtimemin
    rCHproblem%rtimedependence%dtime = dtimemin
    rCHproblem%rtimedependence%dtimemax = dtimemax

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine CH_postprocessing (rCHproblem,rCHvector,iiteration,dtime)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<input>
  ! Number of current iteration
  integer, intent(IN) :: iiteration
  
  ! Current simulation time
  real(DP), intent(IN) :: dtime
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem
  
  ! The current solution vector.
  type(t_vectorBlock), intent(INOUT) :: rCHvector
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    real(DP), dimension(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    type(t_triangulation), POINTER :: p_rtriangulation

!~~~~~~~~~~~~~~~For the first block: for \phi~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rCHvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       'gmv/phi5.gmv.'//TRIM(sys_si0L(iiteration,5)))
    call ucd_setSimulationTime (rexport,dtime)
    
    call lsyssc_getbase_double (rCHvector%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol_phi',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)

!~~~~~~~~~~~~~~~For the second block: for w~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rCHvector%RvectorBlock(2)%p_rspatialDiscr%p_rtriangulation
    
    ! start the postprocessing for w
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       'gmv/chem5.gmv.'//TRIM(sys_si0L(iiteration,5)))
    call ucd_setSimulationTime (rexport,dtime)
    
    call lsyssc_getbase_double (rCHvector%RvectorBlock(2),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol_chemP',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
  end subroutine

end module
