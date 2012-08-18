!##############################################################################
!# ****************************************************************************
!# <name> AllenCahn </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the actual nonstationary solver for the Allen Cahn
!# problem, i.e. the timeloop and time stepping. The following routines can
!# be found here:
!#
!# 1.) AC5_initparameters
!#     -> Initialise the parameters for the time stepping.
!#
!# 2.) AC5_timestep
!#     -> Calculate the solution of the next timestep.
!#
!# 3.) AC5_postprocessing
!#     -> Perform postprocessing (write GMV's,...)
!# </purpose>
!##############################################################################

module AllenCahn

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

  use AllenCahn_callback
  use AllenCahn_basic
  use AllenCahn_nonlinearcore

  use AllenCahn_matvec
  use AllenCahn_boundarycondition
  use AllenCahn_partridiscr
  use AllenCahn_solver


  IMPLICIT NONE
  
CONTAINS
  ! ***************************************************************************
!<subroutine>

  subroutine AC5_timestep (rACproblem,rACvector,rACrhs,rNSproblem,rNSvector,rNSrhs)

!<description>
  ! Performs one time step: $t^n -> t^n+1$.
  ! Assembles system matrix and RHS vector.
  ! Solves the corresponding time-step equation and returns the solution vector
  ! at the end of the time step.
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
  
  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  type(t_vectorBlock), INTENT(INOUT) :: rACvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  type(t_vectorBlock), INTENT(INOUT) :: rACrhs

!Mcai
!~~~~~~~~~~~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(inout), optional :: rNSproblem
  type(t_vectorBlock), intent(inout), optional :: rNSvector
  type(t_vectorBlock), intent(inout), optional :: rNSrhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!</inputoutput>
!</subroutine>

  ! local variables
    INTEGER :: NLMIN,NLMAX
    INTEGER :: i

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    type(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), POINTER :: p_rmatrix
    type(t_vectorBlock), POINTER :: p_rrhs
    type(t_vectorBlock)  :: rtempBlock

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), POINTER :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), POINTER :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), DIMENSION(1:rACproblem%NLMAX) :: Rmatrices
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    type(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo

!~~~~~~~~~~~~~~local variables for AC problem~~~~~~~~~~~~~~~~~~~~~~~~~~
    !local variables for AC problem
    type(t_ACNonlinearIteration) :: rACnonlinearIteration
	! A tmp rnonlinearIterarion
    
    type(t_ACNonlinearIteration) :: rnonlinearIterationTmp
    ! The nonlinear solver configuration
    type(t_nlsolNode) :: rACnlSol

    ! Debug
    real(DP), dimension(:), pointer ::  p_vectordata

! Mcai
	! A parameter to determine, whether we use implicit scheme or explicit scheme for
	! Allen-Cahn problem, if it is 0, we treat convective term explictly. Otherwise,implicitly
	integer :: Implicit_Convective=0

!~~~~~~~~~~~~~~If we only have Allen-Cahn problem, activate the following one~~~~~~~~
    !MCai, For testing the Allen-Cahn solver only, we set rNSvector = 0 vector
	call lsysbl_clearVector(rNSvector)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    NLMIN = rACproblem%NLMIN
    NLMAX = rACproblem%NLMAX

! MCai~~~~~~~Geneate Poly matrix and Conv matrix on all levels~~~~~~~~~~~~~~~
    ! We first generate the matrix based on nonlinear term
!    call AC_assembleMatPoly(rACproblem, rACvector)

!~~~~~~~~~MCai, implicitly treat convection term?~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~think about how to implement this	*)*)*)*)*)*)*)**)*)*)*)*)*)*)*)*)*)*)
! we may also need to update the system matrix of AC problem by using new obtained
! velocity field in NS problem.

    ! We can also generate convetive matrix here
    call AC_assembleMatConv(rACproblem, rACvector, rNSproblem, rNSvector)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~We need to specify rACnonlinearIteration, rACnlSolver~~~~~~~
    call AC_getNonlinearSolver (rACnlSol)

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! our callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    call AC_initNonlinearLoop (rACproblem,rACproblem%NLMIN,rACproblem%NLMAX,&
        rACvector,rACrhs,rACnonlinearIteration)

!~~~~~~~~~~~~~~~call nonlinear solver here~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Set up both RHS and preconditioner for core equation:

    ! First, p_rrhs is to save RHS for core equation
    p_rrhs => rACproblem%rrhs
    call lsysbl_assignDiscrIndirect(rACproblem%rrhs, p_rrhs)

    ! Create a temporary vector we need for some preparations.
    call lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .FALSE.)

    ! Set up w_2*N(u_n) + w_4*f_n.
	! First, set up w_4*f_n
    call lsysbl_vectorLinearComb(rACrhs,p_rrhs,&
         rACproblem%rtimedependence%rtimestepping%dweightOldRHS,0.0_DP)

    ! Set up w_2*N(u_n). Unlike cc2d, we do not need rnonlinearACmatrix...
    ! We pass those coefficient in front of matrices through rACnonlinearIteration
    ! dalpha is for Mass matrix
    rACnonlinearIteration%dalpha=-1.0_DP
	! dtheta is for Laplace term
    rACnonlinearIteration%dtheta=-rACproblem%rtimedependence%rtimestepping%dweightMatrixRHS
	! dconv is for the convetive term
    rACnonlinearIteration%dconv=-rACproblem%rtimedependence%rtimestepping%dweightMatrixRHS
	! dgamma is for nonlinear term
    rACnonlinearIteration%dgamma=-rACproblem%rtimedependence%rtimestepping%dweightMatrixRHS

    call AC_nonlinearMatMul (rACproblem,rACvector,rACnonlinearIteration, &
                                           p_rrhs, -1.0_DP,1.0_DP)
        
    ! Switch to the next point in time.
    rACproblem%rtimedependence%dtime = rACproblem%rtimedependence%dtime + &
          rACproblem%rtimedependence%rtimestepping%dtstep

!~~~~~~~Mcai, we need to rewrite AC5_calRHS, because RHS may depending on rACvector
    ! Generate f_n+1 into the rACrhs overwriting the previous RHS.
	call AC5_calcRHS (rACproblem, rACvector, rACrhs,&
                          rNSproblem, rNSvector, rNSrhs)
! Debugging, to make sure that the integration of coeff_RHS3 is correct.
!     call lsyssc_getbaseVector_double(p_rrhs%rvectorBlock(1), p_vectordata)
!     print *, p_vectordata(1)
!     print *, p_vectordata(10)
!     print *, p_vectordata(100)
!     stop
    ! Add w_3*f_{n+1} to the p_rrhs
    call lsysbl_vectorLinearComb(rACrhs,p_rrhs,&
         rACproblem%rtimedependence%rtimestepping%dweightOldRHS,1.0_DP)
    

!     call AC5_calcRHS (rACproblem, rACvector, rACrhs,&
!                           rNSproblem, rNSvector, rNSrhs)
!
!    !print *, rACproblem%rtimedependence%rtimestepping%dweightOldRHS
!
!     call lsysbl_vectorLinearComb(rACrhs,p_rrhs,&
!          rACproblem%rtimedependence%rtimestepping%dweightOldRHS,1.0_DP)

    ! Discretise the boundary conditions at the new time instant
    call AC5_initDiscreteBC (rACproblem)

    ! Implement boundary conditions into the RHS vector, the solution vector
    ! and the current system matrices.
    call AC5_implementBC (rACproblem,rACvector,p_rrhs,1.0_DP)
	
    ! That's it for the RHS vector. (saved in p_rrhs)
    !
    ! The LHS "u_{n+1} + w_1*N(u_n+1)" results in the system matrix
    ! "M + w_1 N(.)" for the next linear system to solve. Set up that system
    ! on every level of the discretisation.

    ! we can use rnonlinearIterationTmp for assembling the matrix of LHS, or we directly
	! assemble system matrix for all levels and save them in rACproblem%RlevelInfo(i)%rmatrix
!    rnonlinearIterationTmp = rACnonlinearIteration
!	! A = \alpha M + \theta Laplace + \conv Conv + \gamma Nonlinear
    rnonlinearIterationTmp=rACnonlinearIteration
    rnonlinearIterationTmp%dalpha = 1.0_DP
    rnonlinearIterationTmp%dtheta = rACproblem%rtimedependence%rtimestepping%dweightMatrixLHS
    rnonlinearIterationTmp%dconv = rACproblem%rtimedependence%rtimestepping%dweightMatrixLHS
    rnonlinearIterationTmp%dgamma = rACproblem%rtimedependence%rtimestepping%dweightMatrixLHS

! The convection matrix is passed through rACproblem%RlevelInfo(i)%rmatrixConv
	!Preparing the system matrix for all levels, since we are going to use MG
!MCai, do we need it here?
    call AC_updatePreconditioner (rACproblem,rnonlinearIterationTmp,&
                                  rACvector)

! MCai, make sure that rnonlinearIterationTmp%dtheta is implemented correctly

! We have different ways to define p_rsolverNode, one is using the following
!    p_rsolverNode => rACproblem%p_rsolverNode
! The other is to use
!    p_rsolverNode => rACnonlinearIteration%p_rsolverNode
! remember to call
!    call linsol_doneData (p_rsolverNode)

    !print *, 'before core eqn, we have assembled system matrices and p_rsolverNode, saved in rnonlinearIterationTmp'
    call AC_solveCoreEquation (rACproblem,rnonlinearIterationTmp,rACnlSol,&
      rACvector,p_rrhs,rtempBlock)

    ! We need to clear system matrix
	do i= NLMIN, NLMAX
           call lsysbl_clearMatrix(rACproblem%RlevelInfo(i)%rmatrix)
	end do

    ! Release the temporary vector
    call lsysbl_clearVector(p_rrhs)
    call lsysbl_releaseVector (rtempBlock)
    call lsyssc_releaseVector(rnonlinearIterationTmp%p_rtempVectorSc)
    call lsyssc_releaseVector(rnonlinearIterationTmp%p_rtempVectorSc2)
    call linsol_doneData (rnonlinearIterationTmp%p_rsolverNode)
    ! there are memory problem,how to solve this
!     call linsol_doneStructure (rnonlinearIterationTmp%p_rsolverNode)
!     call linsol_releaseSolver (rnonlinearIterationTmp%p_rsolverNode)

    ! Finally tell the time stepping scheme that we completed the time step.
    call timstp_nextSubstep (rACproblem%rtimedependence%rtimestepping)

  end subroutine


  ! ***************************************************************************
!<subroutine>
  subroutine AC5_initparameters (rparams,rACproblem)

!<description>
  ! Reads the DAT file from disc into the parameter list rparams and
  ! initialises basic variables (number of levels, time stepping technique)
  ! in rACproblem according to these settings.
!</description>

!<inputoutput>
  ! A parameter list structure accepting the parameters from the DAT file.
  type(t_parlist), INTENT(INOUT) :: rparams

  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: cscheme,niterations
  real(DP) :: dtheta,dtstep,dtimemin,dtimemax

    ! Read the parameters from disc and put a reference to it
    ! to the collection
!    call parlst_readfromfile(rparams, 'data/AllenCahn.dat')

    
    ! Get the parameters for the time stepping scheme from the parameter list
    call parlst_readfromfile(rparams, 'data/timediscr.dat')
    call parlst_getvalue_int (rparams, 'TIME-DISCRETISATION', 'itimeStepScheme', cscheme, 0)
    call parlst_getvalue_int (rparams, 'TIME-DISCRETISATION ', 'NITERATIONS', niterations, 1000)
    call parlst_getvalue_double (rparams, 'TIME-DISCRETISATION', 'DTIMESTEPTHETA', dtheta, 1.0_DP)
    call parlst_getvalue_double (rparams, 'TIME-DISCRETISATION', 'dtimeStep', dtstep, 0.1_DP)
    call parlst_getvalue_double (rparams, 'TIME-DISCRETISATION', 'DTIMEINIT', dtimemin, 0.0_DP)
    call parlst_getvalue_double (rparams, 'TIME-DISCRETISATION', 'DTIMEMAX', dtimemax, 1.0_DP)


    call  parlst_readfromfile(rparams, 'data/discretisation.dat')
    ! We want to solve our Laplace problem on level...
    call parlst_getvalue_int (rparams, 'CC-DISCRETISATION', 'NLMIN', rACproblem%NLMIN, 7)
    call parlst_getvalue_int (rparams, 'CC-DISCRETISATION', 'NLMAX', rACproblem%NLMAX, 7)
 
    ! Initialise the time stepping in the problem structure
    call timstp_init (rACproblem%rtimedependence%rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)
                     
    rACproblem%rtimedependence%niterations = niterations
    
    rACproblem%rtimedependence%dtimemin = dtimemin
    rACproblem%rtimedependence%dtime = dtimemin
    rACproblem%rtimedependence%dtimemax = dtimemax

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine AC5_postprocessing (rACproblem,rvector,iiteration,dtime)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<input>
  ! Number of current iteration
  INTEGER, INTENT(IN) :: iiteration
  
  ! Current simulation time
  real(DP), INTENT(IN) :: dtime
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
  
  ! The current solution vector.
  type(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
     ! file name for output
     character(SYS_STRLEN) :: sfile,sfilename

    ! We need some more variables for postprocessing
    real(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    type(t_triangulation), POINTER :: p_rtriangulation

! Debug
!    write(*,*)'time dependent'
!    print *, dtime
!    print *, iiteration

    ! Basic filename
    sfilename='AC5_phi.gmv'


!    SFILENAMEUCD='gmv/AC5_phi.gmv'
    sfile = trim(adjustl(sfilename))//'.'//sys_si0(iiteration,5)
!    sfile = trim(adjustl(sfilename))//'.'//sys_si0(iiteration,5)
!    sfile = trim(adjustl(sfilename))//'.'//sys_si0(iiteration,5)


    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation

    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to GMV file:

!    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
!                       'gmv/ACu5.gmv.'//TRIM(sys_si0L(iiteration,5)))

    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       'gmv/'//sfile//TRIM(sys_si0L(iiteration,5)))
     
!    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
         
    call ucd_setSimulationTime (rexport,dtime)
    
    call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)

!    if (present(dtime)) then
      ! Update time stamp of last written out GMV.
!      rpostprocessing%dnextTimeUCD = rpostprocessing%dnextTimeUCD+dtimeDifferenceUCD
!      rpostprocessing%inextFileSuffixUCD = rpostprocessing%inextFileSuffixUCD + 1
!       iiteration=iiteration+1
!    end if

  end subroutine

  ! ***************************************************************************


end module
