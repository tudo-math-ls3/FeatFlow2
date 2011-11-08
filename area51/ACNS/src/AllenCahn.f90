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
!# 1.) AC_initparameters
!#     -> Initialise the parameters for the time stepping.
!#
!# 2.) AC_timestep
!#     -> Calculate the solution of the next timestep.
!#
!# 3.) AC_postprocessing
!#     -> Perform postprocessing (write GMV's,...)
!# </purpose>
!##############################################################################

MODULE AllenCahn

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

! use NS modules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  use ccbasic
  use cccallback
!  use ccnonstationary
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  use AllenCahn_callback
  use AllenCahn_basic
  use AllenCahn_matvec
  use AllenCahn_boundarycondition
  use AllenCahn_partridiscr
  use AllenCahn_solver
!  use AllenCahn_timeloop
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  subroutine AC_initparameters (rparams,rACproblem)

!<description>
  ! Reads the DAT file from disc into the parameter list rparams and
  ! initialises basic variables (number of levels, time stepping technique)
  ! in rACproblem according to these settings.
!</description>

!<inputoutput>
  ! A parameter list structure accepting the parameters from the DAT file.
  TYPE(t_parlist), INTENT(INOUT) :: rparams

  ! A problem structure saving problem-dependent information.
  TYPE(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
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

  subroutine AC_timestep (rACproblem,rACvector,rACrhs,rNSproblem,rNSvector,rNSrhs)
  
!<description>
  ! Performs one time step: $t^n -> t^n+1$.
  ! Assembles system matrix and RHS vector.
  ! Solves the corresponding time-step equation and returns the solution vector
  ! at the end of the time step.
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
  
  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rACvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rACrhs

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
    TYPE(t_matrixBlock), DIMENSION(1:rACproblem%NLMAX) :: Rmatrices
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo


! Mcai
	! A parameter to determine, whether we use implicit scheme or explicit scheme for
	! Allen-Cahn problem, if it is 0, we treat convective term explictly. Otherwise,implicitly
	integer :: Implicit_Convective=0


!  added~~~~for initialize the solution of AC problem
    real(DP), dimension(:,:), pointer ::  p_DvertexCoords
    real(DP), dimension(:), pointer ::  p_vectordata
    real(DP), dimension(:), pointer ::  p_data

!~~~~~~~~~~~~~~If we only have Allen-Cahn problem, activate the following one~~~~~~~~
    !MCai, For testing the Allen-Cahn solver only, we set rNSvector = 0 vector
	!call lsysbl_clearVector(rNSvector)
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
    call AC_assembleMatPoly(rACproblem, rACvector)

!~~~~~~~~~MCai, implicitly treat convection term?~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! we may also need to update the system matrix of AC problem by using new obtained
! velocity field in NS problem.

!~~think about how to implement this	*)*)*)*)*)*)*)**)*)*)*)*)*)*)*)*)*)*)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! We can also generate convetive matrix here
    call AC_assembleMatConv(rACproblem, rACvector, rNSproblem, rNSvector)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rmatrix => rACproblem%RlevelInfo(NLMAX)%rmatrix
    p_rrhs    => rACproblem%rrhs
    
    ! Create a temporary vector we need for some preparations.
    call lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .FALSE.)
    
    ! Set up w_2*N(u_n) + w_4*f_n.
    call lsysbl_vectorLinearComb(rACrhs,p_rrhs,&
         rACproblem%rtimedependence%rtimestepping%dweightOldRHS,0.0_DP)
!First, Laplace term
	call lsysbl_blockMatVec(rACproblem%RlevelInfo(NLMAX)%rmatrixStatic,&
         rACvector,p_rrhs,&
         rACproblem%rtimedependence%rtimestepping%dweightMatrixRHS,&
         rACproblem%rtimedependence%rtimestepping%dweightOldRHS)
!Then, Convective term
	call lsysbl_blockMatVec(rACproblem%RlevelInfo(NLMAX)%rmatrixConv,&
         rACvector,p_rrhs,&
         rACproblem%rtimedependence%rtimestepping%dweightMatrixRHS,&
         rACproblem%rtimedependence%rtimestepping%dweightOldRHS)
!MCai,
!Finally, Polynomial term: here we treated it through ????????????????????
! We can also do it through the following way:
! 	call lsysbl_blockMatVec(rACproblem%RlevelInfo(NLMAX)%rmatrixPoly,&
!          rACvector,p_rrhs,&
!          rACproblem%rtimedependence%rtimestepping%dweightMatrixRHS,&
!          rACproblem%rtimedependence%rtimestepping%dweightOldRHS)

    ! Add u_n -- or, more precisely, M u_n (with M being the mass matrix),
    ! since the discretisation with finite elements requires that.
    call lsysbl_blockMatVec(rACproblem%RlevelInfo(NLMAX)%rmatrixMass,&
         rACvector,p_rrhs,1.0_DP,1.0_DP)

    ! Switch to the next point in time.
    rACproblem%rtimedependence%dtime = rACproblem%rtimedependence%dtime + &
          rACproblem%rtimedependence%rtimestepping%dtstep

!~~~~~~~Mcai, we need to rewrite AC_calRHS, because RHS may depending on rACvector
    ! Generate f_n+1 into the rrhs overwriting the previous RHS.
!    call AC_calcRHS (rACproblem, rACvector, rACrhs)


! ! Debug
!     call lsyssc_getbaseVector_double(rACvector%rvectorBlock(1), p_vectordata)
!     call storage_getbase_double2D(rACproblem%RlevelInfo(rACproblem%NLMAX)%rtriangulation%h_DvertexCoords,p_DvertexCoords)
! !
!       do i=1,rACproblem%RlevelInfo(rACproblem%NLMAX)%rtriangulation%NVT
!            call AllenCahn_inicon(p_DvertexCoords(1,i),p_DvertexCoords(2,i), p_vectordata(i))
!       end do

    call AC_calcRHS (rACproblem, rACvector, rACrhs,&
                          rNSproblem, rNSvector, rNSrhs)

    call lsysbl_vectorLinearComb(rACrhs,p_rrhs,&
         rACproblem%rtimedependence%rtimestepping%dweightNewRHS,1.0_DP)

    !print *, rACproblem%rtimedependence%rtimestepping%dweightNewRHS
    !stop
    ! if \theta =0.5, we should use
    !call lsysbl_vectorLinearComb(rACrhs,p_rrhs,&
    !     rACproblem%rtimedependence%rtimestepping%dweightNewRHS,1.0_DP)

    ! For mass conservation
    call AC_conservativeRHS(rACproblem, rACvector, rACrhs)
    ! MCai note that it is  \delta t * \int_{\Omega} f(\phi)
     call lsysbl_vectorLinearComb(rACrhs,p_rrhs,&
                   rACproblem%rtimedependence%rtimestepping%dtstep,1.0_DP)
    ! That's it for the RHS vector. (saved in p_rrhs)
    !
    ! The LHS "u_{n+1} + w_1*N(u_n+1)" results in the system matrix
    ! "M + w_1 N(.)" for the next linear system to solve. Set up that system
    ! on every level of the discretisation.


    do i = NLMIN,NLMAX
      ! Mass matrix
      call lsyssc_duplicateMatrix (rACproblem%RlevelInfo(i)%rmatrixMass%RmatrixBlock(1,1),&
                                   rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                   LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
      ! Laplace Matrix
      call lsyssc_matrixLinearComb (rACproblem%RlevelInfo(i)%rmatrixStatic%RmatrixBlock(1,1),&
                                    rACproblem%rtimedependence%rtimestepping%dweightMatrixLHS,&
                                    rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    1.0_DP,&
                                    rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    .FALSE.,.FALSE.,.TRUE.,.TRUE.)
      ! Polynomical Matrix
! Here, the system matrix is time independent.
! We first add mass matrix from polynomial term first
!      call lsyssc_matrixLinearComb (rACproblem%RlevelInfo(i)%rmatrixPoly%RmatrixBlock(1,1),&
!                                   rACproblem%rtimedependence%rtimestepping%dweightMatrixLHS,&
!                                   rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
!                                   1.0_DP,&
!                                   rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
!                                   .FALSE.,.FALSE.,.TRUE.,.TRUE.)
      ! Then convective matrix

!    if (Implicit_Convective .eq. 1) then
      call lsyssc_matrixLinearComb (rACproblem%RlevelInfo(i)%rmatrixConv%RmatrixBlock(1,1),&
                                    rACproblem%rtimedependence%rtimestepping%dweightMatrixLHS,&
                                    rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    1.0_DP,&
                                    rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    .FALSE.,.FALSE.,.TRUE.,.TRUE.)
!    end if

    end do
    
    ! Discretise the boundary conditions at the new time instant
    call AC_initDiscreteBC (rACproblem)

    ! Implement boundary conditions into the RHS vector, the solution vector
    ! and the current system matrices.
    call AC_implementBC (rACproblem,rACvector,p_rrhs,1.0_DP)
    ! Preparation of the linear system completed!
    !
    ! Attach the system matrices to the solver.
    p_rsolverNode => rACproblem%p_rsolverNode
    
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(NLMIN:NLMAX) = rACproblem%RlevelInfo(NLMIN:NLMAX)%rmatrix
    call linsol_setMatrices(p_rsolverNode,Rmatrices(NLMIN:NLMAX))
    
    ! Initialise data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initData (p_rsolverNode,ierror)
    if (ierror .NE. LINSOL_ERR_NOERROR) STOP

    ! Synchronise p_rrhs with the matrix so it's compatible to the linear system.
    !call lsysbl_synchroniseSortMatVec (p_rmatrix,p_rrhs,rtempBlock%RvectorBlock(1))

    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. if b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rACvector,p_rrhs,rtempBlock)
    
     ! rvector is now u_n+1.
    !
    ! Release solver data.
    call linsol_doneData (p_rsolverNode)
 
    ! Release the temporary vector
    call lsysbl_releaseVector (rtempBlock)
    
    ! Finally tell the time stepping scheme that we completed the time step.
    call timstp_nextSubstep (rACproblem%rtimedependence%rtimestepping)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine AC_postprocessing (rACproblem,rvector,iiteration,dtime)
  
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
  TYPE(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
  
  ! The current solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
     ! file name for output
     character(SYS_STRLEN) :: sfile,sfilename

    ! We need some more variables for postprocessing
    real(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_triangulation), POINTER :: p_rtriangulation

! Debug
!    write(*,*)'time dependent'
!    print *, dtime
!    print *, iiteration

    ! Basic filename
    sfilename='AC_phi.gmv'


!    SFILENAMEUCD='gmv/AC_phi.gmv'
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


end MODULE
