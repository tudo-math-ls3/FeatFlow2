!##############################################################################
!# ****************************************************************************
!# <name> nonlinearoneshotspacetimesolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module aims to solve the fully space-time coupled Stokes and
!# Navier-Stokes optimal control problem with a one-shot approach.
!# As this type of problem is of elliptic nature with a very special
!# structure, a nonlinear defect correction loop is executed for the solution
!# where a space-time coupled multigrid preconditioner does the preconditioning
!# of the nonlinear defect.
!#
!# The following subroutines can be found here:
!#
!# 1.) nlstslv_initStdDiscrData
!#     -> Initialisation
!#
!# 2.) nlstslv_solve
!#     -> Applies a nonlinear loop to solve the system.
!#
!# </purpose>
!##############################################################################

module nonlinearoneshotspacetimesolver

  use fsystem
  use storage
  use genoutput
  use paramlist
  use basicgeometry
  use boundary
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  
  use spatialdiscretisation
  use timediscretisation
  
  use cubature
  use filtersupport
  use collection
  use statistics
  
  use analyticsolution
  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimehierarchy
  
  use constantsoptc
  use structuresoptc
  use structuresspacetimelinsol
  use structuresoptcspacetimenlsol
  use structuresoptflow
  
  use spacediscretisation
  
  use spatialbcdef
  use initmatrices
  use postprocessing
  use user_callback
  
  use timescalehierarchy
  use spacematvecassembly
  use spacetimevectors
  use spacetimelinearsystem
  use forwardbackwardsimulation
  use spacetimeinterlevelprojection
  use spacetimelinearsolver
  use optcanalysis
  use timeboundaryconditions
  use spacetimeneumannbc
  use spacetimedirichletbcc
  use timerhsevaluation
    
  implicit none
  
  private
  
  public :: nlstslv_initStdDiscrData
  public :: nlstslv_solve

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine nlstslv_initStdDiscrData (rsettings,ilevel,rdiscrData)

!<description>
  ! Initialises rdiscrData with standard values from rsettings.
!</description>

!<input>
  ! Global settings structure.
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Space-time level to take.
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Space-assembly structure.
  type(t_spaceTimeMatrixDiscrData), intent(out) :: rdiscrData
!</output>

!</subroutine>

    ! local variables
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    integer :: ispaceLevel
    integer :: itimeLevel

    rdiscrData%p_rphysicsPrimal => rsettings%rphysicsPrimal
    rdiscrData%p_rstabilPrimal => rsettings%rstabilPrimal
    rdiscrData%p_rstabilDual => rsettings%rstabilDual
    rdiscrData%p_rconstraints => rsettings%rsettingsOptControl%rconstraints
    
    call sth_getLevel (rsettings%rspaceTimeHierPrimal,ilevel,&
        p_rfeSpaceLevel)

    rdiscrData%p_rdiscrPrimal => p_rfeSpaceLevel%p_rdiscretisation

    call sth_getLevel (rsettings%rspaceTimeHierPrimalDual,ilevel,&
        p_rfeSpaceLevel,p_rtimeDiscr,ispaceLevel,itimeLevel)

    rdiscrData%p_rdiscrPrimalDual => p_rfeSpaceLevel%p_rdiscretisation
    rdiscrData%p_rtimeDiscr => p_rtimeDiscr
    rdiscrData%p_rspaceDiscr => rdiscrData%p_rdiscrPrimalDual
    
    rdiscrData%p_rstaticSpaceAsmTempl => &
        rsettings%rspaceAsmHierarchy%p_RasmTemplList(ispaceLevel)
    rdiscrData%p_rstaticSpaceAsmTemplOptC => &
        rsettings%rspaceAsmHierarchyOptC%p_RasmTemplList(ispaceLevel)
    
    rdiscrData%p_rsettingsOptControl => rsettings%rsettingsOptControl
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine nlstslv_assembleNeumannBC (rsettings,RneumannBC)
  
!<description>
  ! Sets up the Neumann boundary condition structures for all levels.
!</description>

!<input>
  ! Settings structure with global parametes.
  type(t_settings_optflow), intent(inout) :: rsettings
!</input>

!<inputoutput>
  ! Array with Neumann boundary definitions for nonlinear boundary conditions,
  ! for all levels.
  type(t_sptiNeumannBoundary), dimension(:), intent(inout), target :: RneumannBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev, nlevels
    
    ! How many levels do we have in the hierarchy?
    nlevels = size(RneumannBC)
    
    ! Loop through all levels to set up the matrices.
    do ilev = 1,nlevels
      ! Create the Neumann boundary conditions
      call stnm_assembleNeumannBoundary (rsettings%roptcBDC,RneumannBC(ilev),rsettings%rglobalData)
    end do
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine nlstslv_assembleDirichletBCC (rsettings,RdirichletBCC)
  
!<description>
  ! Sets up the Dirichlet boundary control boundary condition structures for all levels.
!</description>

!<input>
  ! Settings structure with global parametes.
  type(t_settings_optflow), intent(inout) :: rsettings
!</input>

!<inputoutput>
  ! Array with Neumann boundary definitions for nonlinear boundary conditions,
  ! for all levels.
  type(t_sptiDirichletBCCBoundary), dimension(:), intent(inout), target :: RdirichletBCC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev, nlevels
    
    ! How many levels do we have in the hierarchy?
    nlevels = size(RdirichletBCC)
    
    ! Loop through all levels to set up the matrices.
    do ilev = 1,nlevels
      ! Create the Neumann boundary conditions
      call stdbcc_assembleDirichletBCCBd (rsettings%roptcBDC,&
          RdirichletBCC(ilev),rsettings%rglobalData)
    end do
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine nlstslv_initPrecMatrices (rsettings,ctypeNonlinearIteration,rmatrix,&
      RprecMatrices,Rsolutions,RneumannBC,RdirichletBCC)
  
!<description>
  ! Sets up the preconditioner matrices on all levels according to
  ! the settings in the structures and the fine grid space-time matrix rmatrix.
!</description>

!<input>
  ! Settings structure with global parametes.
  type(t_settings_optflow), intent(inout) :: rsettings

  ! Type of the nonlinear iteration.
  integer, intent(in) :: ctypeNonlinearIteration
  
  ! Space-time matrix on the finest space-time mesh.
  type(t_ccoptspaceTimeMatrix), intent(in) :: rmatrix
!</input>

!<inputoutput>
  ! Array with space-time vectors with evaluation points of the nonlinearities
  ! on all levels.
  type(t_spaceTimeVector), dimension(:), intent(inout), target :: Rsolutions
  
  ! Array with Neumann boundary definitions for nonlinear boundary conditions,
  ! for all levels.
  type(t_sptiNeumannBoundary), dimension(:), intent(in), target :: RneumannBC

  ! Array with Neumann boundary definitions for Dirichlet control boundary conditions,
  ! for all levels.
  type(t_sptiDirichletBCCBoundary), dimension(:), intent(in), target :: RdirichletBCC
!</inputoutput>

!<output>
  ! Hierarchy of space-time matrices on all levels, used in the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(out) :: RprecMatrices
!</output>

!</subroutine>

    ! local variables
    integer :: nlevels,cmatrixType,ilev
    type(t_spaceTimeMatrixDiscrData) :: rdiscrData
    type(t_vectorBlock), dimension(:), allocatable :: p_RsolVectors
    
    ! Which type on nonlinear solver do we have. Should we assemble the Newton part?
    cmatrixType = MATT_OPTCONTROL
    select case (ctypeNonlinearIteration)
    case (2:)
      ! Assemble Newton
      cmatrixType = MATT_LINOPTCONTROL
    end select
    
    ! How many levels do we have in the hierarchy?
    nlevels = size(Rsolutions)

    ! Copy the solution of the topmost level to the corresponding
    ! evaluation point.
    call sptivec_copyVector (rmatrix%p_rsolution,Rsolutions(nlevels))
    
    ! Allocate temp vectors in space
    allocate(p_RsolVectors(nlevels))
    
    ! Loop through all levels to set up the matrices.
    do ilev = nlevels,1,-1
      ! Get the discretisation data of the level
      call nlstslv_initStdDiscrData (rsettings,ilev,rdiscrData)

      ! Initialise the corresponding space-time matrix
      call stlin_initSpaceTimeMatrix (&
          RprecMatrices(ilev),cmatrixType,rdiscrData,&
          Rsolutions(ilev),RneumannBC(ilev),RdirichletBCC(ilev),&
          rsettings%rglobalData,rsettings%rdebugFlags)

      ! Create a temp vector here for the interpolation
      call lsysbl_createVectorBlock(rdiscrData%p_rdiscrPrimalDual,p_RsolVectors(ilev))

      ! Project the solution of the higher level to the current level.
      if (ilev .lt. nlevels) then
        call sptipr_performInterpolation (rsettings%rprjHierSpaceTimePrimalDual,&
            ilev+1,Rsolutions(ilev),Rsolutions(ilev+1),p_RsolVectors(ilev),&
            p_RsolVectors(ilev+1))
      end if
      
    end do
    
    ! Release the temp vectors
    do ilev = 1,nlevels
      call lsysbl_releaseVector(p_RsolVectors(ilev))
    end do
    deallocate(p_RsolVectors)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine nlstslv_donePrecMatrices (RprecMatrices)
  
!<description>
  ! Cleans up the preconditioner matrices on all levels.
!</description>

!<inputoutput>
  ! Hierarchy of space-time matrices on all levels, to be cleaned up.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(out) :: RprecMatrices
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilev

    ! Loop through all levels to set up the matrices.
    do ilev = 1,size(RprecMatrices)
      call stlin_releaseSpaceTimeMatrix (RprecMatrices(ilev))
    end do

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine nlstslv_solve (rsettings,rnlstsolver,rpostproc,rx,rb,rd)
  
!<description>
  ! Solves a nonstationary space-time problem with parameters in rsettings
  ! and rnlstsolver.
!</description>

!<input>
  ! Settings structure with global parametes.
  type(t_settings_optflow), intent(inout) :: rsettings

  ! Discrete version of the space-time RHS.
  type(t_spacetimeVector), intent(in), target :: rb

  ! Postprocessing data.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</input>

!<inputoutput>
  ! Nonlinear solve structure representing the solver.
  type(t_nlstsolver), intent(inout) :: rnlstsolver

  ! A space-time vector defining the initial solution. Is replaced by a new
  ! solution vector.
  type(t_spacetimeVector), intent(inout), target :: rx

  ! A temporary space-time vector.
  type(t_spacetimeVector), intent(inout), target :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: nlevels, ierror, ilev
    real(dp) :: dinitDefNorm,ddefNorm,dlastDefNorm,dtempdef,delapsedReal
    real(DP), dimension(4) :: DerrorU,DerrorP,DerrorLambda,DerrorXi,Derror
    logical :: bnewtonAllowed, breassembleStructure
    
    ! Some statistical data
    type(t_timer) :: rtimeFactorisationStep,rtimerMGstep,rtimerIterate,rtimerPostproc

    ! Nonlinear space-time matrix representing the discrete operator, nonlinearity,...
    type(t_spaceTimeMatrixDiscrData) :: rdiscrData
    type(t_ccoptSpaceTimeMatrix) :: rmatrix
    type(t_sptiNeumannBoundary) :: rsptiNeumannBC
    type(t_sptiDirichletBCCBoundary) :: rsptiDirichletBCC

    ! Reset the timings
    call stat_clearTimer (rnlstsolver%rtimeSmoothing)
    call stat_clearTimer (rnlstsolver%rtimeSmoothingFinest)
    call stat_clearTimer (rnlstsolver%rtimeCoarseGridSolver)
    call stat_clearTimer (rnlstsolver%rtimeLinearAlgebra)
    call stat_clearTimer (rnlstsolver%rtimeProlRest)
    call stat_clearTimer (rnlstsolver%rtimeSpacePrecond)
    call stat_clearTimer (rnlstsolver%rtimePostprocessing)
    
    call stat_clearTimer (rnlstsolver%rtimerPreconditioner)
    call stat_clearTimer (rnlstsolver%rtimerNonlinear)
    
    ! Total number of available levels. The solution is computed
    ! on the topmost level.
    nlevels = rsettings%rspaceTimeHierPrimalDual%nlevels
    
    ! Get the discretisation data of the level
    call nlstslv_initStdDiscrData (rsettings,nlevels,rdiscrData)
    
    ! Initialise the Neumann boundary conditions on the maximum level.
    ! Used for nonlinear boundary conditions.
    call stnm_createNeumannBoundary (rdiscrData%p_rspaceDiscr,rdiscrData%p_rtimeDiscr,&
        rdiscrData%p_rstaticSpaceAsmTempl,rsptiNeumannBC)
    call stnm_assembleNeumannBoundary (rsettings%roptcBDC,rsptiNeumannBC,rsettings%rglobalData)

    ! Initialise the Dirichlet boundary control boundary conditions on the maximum level.
    ! Used for nonlinear boundary conditions.
    call stdbcc_createDirichletBCCBd (rdiscrData%p_rspaceDiscr,rdiscrData%p_rtimeDiscr,&
        rdiscrData%p_rstaticSpaceAsmTempl,rsptiDirichletBCC)
    call stdbcc_assembleDirichletBCCBd (rsettings%roptcBDC,rsptiDirichletBCC,rsettings%rglobalData)
    
    ! Initialise a space-time matrix of the corresponding system.
    ! The solution vector is our nonlinearity.
    call stlin_initSpaceTimeMatrix (&
        rmatrix,MATT_OPTCONTROL,rdiscrData,rx,rsptiNeumannBC,rsptiDirichletBCC,&
        rsettings%rglobalData,rsettings%rdebugFlags)
    
    ! ---------------------------------------------------------------
    ! Create the initial defect
    call output_line ('NLST-Solver: Assembling the initial defect...')

    ! Create the actual RHS (including the BC's) in rd and assemble the defect.
    call sptivec_copyVector (rb,rd)
    call trhsevl_implementBDCRHS (rsettings%rglobalData, rd, rsettings%roptcBDC)
    
    ! DEBUG!!!
    call stlin_spaceTimeMatVec (rmatrix, rx, rd, &
        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,&
        ddefNorm,rsettings%roptcBDC,rnlstsolver%ioutputLevel .ge. 2)
    !ddefNorm = 1.0_DP

    dinitDefNorm = ddefNorm
    dlastDefNorm = 0.0_DP
    if (dinitDefNorm .eq. 0.0_DP) then
      ! Trick to avoid div/0.
      dinitDefNorm = 1.0_DP
    end if
    if (rnlstsolver%ioutputLevel .ge. 2) call output_separator (OU_SEP_MINUS)
    
    if (rnlstsolver%ioutputLevel .ge. 1) &
      call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))

    ! ---------------------------------------------------------------
    ! Check if the Newton is allowed. This is the case if the
    ! max. number of fixed point iterations is > 0 => 1st iteration is
    ! fixed point.
    
    bnewtonAllowed = rnlstsolver%nmaxFixedPointIterations .eq. 0

    ! ---------------------------------------------------------------
    ! Start the nonlinear iteration
    
    ! Initialise Neumann boundary conditions.
    call nlstslv_assembleNeumannBC (rsettings,rnlstsolver%p_rsptiNeumannBC)

    ! Initialise Dirichlet boundary control boundary conditions.
    call nlstslv_assembleDirichletBCC (rsettings,rnlstsolver%p_rsptiDirichletBCC)
    
    ! Initialise the nonlinear matrices on all levels for the first iteration
    if (bnewtonAllowed .or. (rnlstsolver%ctypeNonlinearIteration .eq. 1)) then
      call nlstslv_initPrecMatrices (rsettings,rnlstsolver%ctypeNonlinearIteration,&
          rmatrix,rnlstsolver%p_RprecMatrices,rnlstsolver%p_Rsolutions,&
          rnlstsolver%p_rsptiNeumannBC,rnlstsolver%p_rsptiDirichletBCC)
    else
      ! Newton not allowed. Select fixed point.
      !if (rnlstsolver%ioutputLevel .ge. 1) then
      !  call output_line ("Adaptive Newton: Selecting fixed point iteration.")
      !end if
      call nlstslv_initPrecMatrices (rsettings,1,&
          rmatrix,rnlstsolver%p_RprecMatrices,rnlstsolver%p_Rsolutions,&
          rnlstsolver%p_rsptiNeumannBC,rnlstsolver%p_rsptiDirichletBCC)
    end if
    
    ! Pass the matrices to the linear solver.
    call sptils_setMatrix(rnlstsolver%p_rspaceTimePrec,&
        rnlstsolver%p_RprecMatrices(size(rnlstsolver%p_RprecMatrices)))

    ! Pass matrices of lower levels to the MG solver if present
    if (associated(rnlstsolver%p_rmgSolver)) then
      do ilev = 1,size(rnlstsolver%p_RprecMatrices)-1
        call sptils_setMatrixMultigrid(rnlstsolver%p_rspaceTimePrec,&
            rnlstsolver%p_RprecMatrices(ilev),ilev)
      end do
    end if
    
    ! Initialise the structure of the linear subsolver(s).
    call stat_startTimer (rnlstsolver%rtimeFactorisation)
    call sptils_initStructure (rnlstsolver%p_rspaceTimePrec,ierror)
    call stat_stopTimer (rnlstsolver%rtimeFactorisation)
    if (ierror .gt. 0) then
      call output_line ('sptils_initStructure failed! Matrix structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'nlstslv_solve')
      call sys_halt()
    end if

    ! Here we go...

    call stat_startTimer (rnlstsolver%rtimerNonlinear)

    rnlstsolver%nnonlinearIterations = 0
    rnlstsolver%nlinearIterations = 0
    rnlstsolver%nlinearIterationsSpace = 0
    rnlstsolver%ncoarsegridIterations = 0
    
    do while ((rnlstsolver%nnonlinearIterations .lt. rnlstsolver%nminIterations) .or. &
              ((((ddefNorm .gt. rnlstsolver%depsRel*dinitDefNorm) .or. &
                 (ddefNorm .ge. rnlstsolver%depsAbs)) .and. &
                (abs(ddefNorm-dlastDefNorm) .ge. rnlstsolver%depsDiff*dlastDefNorm)) &
              .and. (rnlstsolver%nnonlinearIterations .lt. rnlstsolver%nmaxIterations)))

      ! Reset the timings    
      call stat_clearTimer (rtimerIterate)
      call stat_clearTimer (rtimerPostproc)

      ! Measure time for the current iterate.
      call stat_startTimer (rtimerIterate)
    
      if (rnlstsolver%cpostprocessIterates .ne. 1) then
        ! Postprocessing of the current iterate.
        call stat_startTimer (rtimerPostproc)
        if (rnlstsolver%ioutputLevel .ge. 1) then
          call output_separator (OU_SEP_MINUS)
          call output_line ("Postprocessing of the current iterate.")
        end if
        
        select case (rnlstsolver%cpostprocessIterates)
        case (1)
          call optcpp_postprocessSpaceTimeVec (rpostproc,rx,rb,&
              rsettings%rsettingsOptControl,rsettings)
        case (2)
          call optcpp_postprocessSpaceTimeVec (rpostproc,rx,rb,&
              rsettings%rsettingsOptControl,rsettings,&
              rnlstsolver%nnonlinearIterations)
        end select
        
        call stat_stopTimer (rtimerPostproc)
      end if
                  
      rnlstsolver%nnonlinearIterations = rnlstsolver%nnonlinearIterations+1
      
      if (rnlstsolver%ctypeNonlinearIteration .eq. CCNLS_INEXACTNEWTON) then
      
        ! Determine the stopping criterion for the inexact Newton.
        ! This is an adaptive stopping criterion depending on the current
        ! defect. In detail, we calculate:
        !
        !   |b-Ax_{i+1}|         ( |b-Ax_i| ) exp             ( |b-Ax_i| )
        !   ------------ = min { ( -------- )     , depsrel * ( -------- ) }
        !     |b-Ax_0|           ( |b-Ax_0| )                 ( |b-Ax_0| )
        !
        ! see e.g. [Michael Hinze, Habilitation, p. 51].
        ! We furthermore introduce the additional stopping criterion
        !
        !   |b-Ax_{i+1}|
        !   ------------ = depsrel(solver)*dinexactNewtonEpsRel
        !     |b-Ax_0|
        !
        ! which stops the iteration in order to prevent the linear
        ! solver from solving too much further than necessary for the
        ! nonlinear solver.
        !
        ! Switch off the relative stopping criterion in the linear solver:
        
        rnlstsolver%p_rspaceTimePrec%depsRel = 0.0_DP
        rnlstsolver%p_rspaceTimePrec%istoppingcriterion = SPTILS_STOP_STANDARD
        
        ! Calculate the new absolute stopping criterion:
        
        dtempdef = ddefNorm / dinitDefNorm
        
        rnlstsolver%p_rspaceTimePrec%depsAbs = &
            
            min(&
        
              ! Try to gain quadratic convergence, but get at most only a little but
              ! more digits than depsRel. If depsAbs is lower, take that as a bound.
              max(&
                  min(0.1_DP * rnlstsolver%depsRel * dinitDefNorm,&
                      0.1_DP * rnlstsolver%depsAbs),&
                  dtempDef**rnlstsolver%dinexactNewtonExponent * dinitDefNorm),&
                
              ! Get at least dinexactNewtonEpsRel.
              rnlstsolver%dinexactNewtonEpsRel*dtempdef * dinitDefNorm)
        
        if (associated(rnlstsolver%p_rcgrSolver)) then
          ! For the coarse grid solver, we choose the same stopping criterion.
          ! But just for safetyness, the coarse grid solver should gain at least
          ! one digit!
          rnlstsolver%p_rcgrSolver%istoppingcriterion = SPTILS_STOP_STANDARD
          rnlstsolver%p_rcgrSolver%depsRel = 1.0E-1_DP
          rnlstsolver%p_rcgrSolver%depsAbs = rnlstsolver%p_rspaceTimePrec%depsAbs
        end if
      
      end if

      ! Value of the functional
      call stat_startTimer (rtimerPostproc)
      call optcana_nonstatFunctional (rsettings%rglobalData,rsettings%rphysicsPrimal,&
          rsettings%rsettingsOptControl%rconstraints,&
          rx,rsettings%rsettingsOptControl%rtargetFunction,&
          rsettings%rsettingsOptControl%dalphaC,rsettings%rsettingsOptControl%dbetaC,&
          rsettings%rsettingsOptControl%dgammaC,&
          Derror)
      if (rnlstsolver%ioutputLevel .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(1),10)))
        call output_line ('||u||         = '//trim(sys_sdEL(Derror(2),10)))
        call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
        call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(4),10)))
      end if

      ! If error analysis has to be performed, we can calculate
      ! the real error.
      if (rpostproc%icalcError .eq. 1) then
        if (rnlstsolver%ioutputLevel .ge. 1) then
          call optcana_analyticalError (rsettings%rglobalData,&
              rsettings%rsettingsOptControl%rconstraints,&
              rx,rpostproc%ranalyticRefFunction,&
              DerrorU,DerrorP,DerrorLambda,DerrorXi,.true.)
              
          call output_line ('||y-y0||_[0,T]           = '//trim(sys_sdEL(DerrorU(1),10)))
          call output_line ('||y-y0||_[0,T)           = '//trim(sys_sdEL(DerrorU(2),10)))
          call output_line ('||y-y0||_(0,T]           = '//trim(sys_sdEL(DerrorU(3),10)))
          call output_line ('||y-y0||_(0,T)           = '//trim(sys_sdEL(DerrorU(4),10)))
          
          call output_line ('||p-p0||_[0,T]           = '//trim(sys_sdEL(DerrorP(1),10)))
          call output_line ('||p-p0||_[0,T)           = '//trim(sys_sdEL(DerrorP(2),10)))
          call output_line ('||p-p0||_(0,T]           = '//trim(sys_sdEL(DerrorP(3),10)))
          call output_line ('||p-p0||_(0,T)           = '//trim(sys_sdEL(DerrorP(4),10)))
          
          call output_line ('||lambda-lambda0||_[0,T] = '//trim(sys_sdEL(DerrorLambda(1),10)))
          call output_line ('||lambda-lambda0||_[0,T) = '//trim(sys_sdEL(DerrorLambda(2),10)))
          call output_line ('||lambda-lambda0||_(0,T] = '//trim(sys_sdEL(DerrorLambda(3),10)))
          call output_line ('||lambda-lambda0||_(0,T) = '//trim(sys_sdEL(DerrorLambda(4),10)))

          call output_line ('||xi-xi0||_[0,T]         = '//trim(sys_sdEL(DerrorXi(1),10)))
          call output_line ('||xi-xi0||_[0,T)         = '//trim(sys_sdEL(DerrorXi(2),10)))
          call output_line ('||xi-xi0||_(0,T]         = '//trim(sys_sdEL(DerrorXi(3),10)))
          call output_line ('||xi-xi0||_(0,T)         = '//trim(sys_sdEL(DerrorXi(4),10)))
        end if
      end if
      
      call stat_stopTimer (rtimerPostproc)

      if (rnlstsolver%ioutputLevel .ge. 1) then
        if (rnlstsolver%ctypeNonlinearIteration .eq. CCNLS_INEXACTNEWTON) then
          call output_lbrk ()
          call output_line ('Inexact Newton: Stopping criterion = '//&
              trim(sys_sdEL(rnlstsolver%p_rspaceTimePrec%depsAbs,10)))
        end if
        call output_separator (OU_SEP_EQUAL)
      end if
            
      ! Preconditioning of the defect: d=C^{-1}d
      !
      ! Re-initialise the preconditioner matrices on all levels.
!      ! for iteration >= 2.
!      if (rnlstsolver%nnonlinearIterations .gt. 1) then

      breassembleStructure = .false.

      select case (rnlstsolver%ctypeNonlinearIteration)
      case (1)
        ! Standard fixed point iteraion
        call nlstslv_donePrecMatrices (rnlstsolver%p_RprecMatrices)
        call nlstslv_initPrecMatrices (rsettings,rnlstsolver%ctypeNonlinearIteration,&
            rmatrix,rnlstsolver%p_RprecMatrices,rnlstsolver%p_Rsolutions,&
            rnlstsolver%p_rsptiNeumannBC,rnlstsolver%p_rsptiDirichletBCC)

        ! The structure of the matrices does not change,
        ! so we don't have to call sptils_initStructure again.

      case (2:)
        ! Newton iteration, probably fallback to fixed point iteration.
        ! Newton allowed -> initialisation as above
        
        ! 2.) Is Newton to be allowed?
        if ((rnlstsolver%nnonlinearIterations .gt. rnlstsolver%nmaxFixedPointIterations) .or. &
            (ddefNorm .le. rnlstsolver%depsRelFixedPoint*dinitDefNorm)) then
          
          ! Yep, has to be activated.

          call nlstslv_donePrecMatrices (rnlstsolver%p_RprecMatrices)
          call nlstslv_initPrecMatrices (rsettings,rnlstsolver%ctypeNonlinearIteration,&
              rmatrix,rnlstsolver%p_RprecMatrices,rnlstsolver%p_Rsolutions,&
              rnlstsolver%p_rsptiNeumannBC,rnlstsolver%p_rsptiDirichletBCC)

          ! Is Newton already active?
          if (.not. bnewtonAllowed) then
          
            ! Activate Newton, reasseble matrix structure.
            bnewtonAllowed = .true.
            breassembleStructure = .true.
          
          end if
          
        else
        
          ! No, select fixed point iteration
          if (rnlstsolver%ioutputLevel .ge. 1) then
            call output_line ("Adaptive Newton: Selecting fixed point iteration.")
          end if
          
          call nlstslv_donePrecMatrices (rnlstsolver%p_RprecMatrices)
          call nlstslv_initPrecMatrices (rsettings,1,&
              rmatrix,rnlstsolver%p_RprecMatrices,rnlstsolver%p_Rsolutions,&
              rnlstsolver%p_rsptiNeumannBC,rnlstsolver%p_rsptiDirichletBCC)

          ! Is Newton currently active?
          if (.not. bnewtonAllowed) then
          
            ! Dectivate Newton, reasseble matrix structure.
            bnewtonAllowed = .false.
            breassembleStructure = .true.
          
          end if
        
        end if
        
      end select

      if (breassembleStructure) then
      
        call stat_startTimer (rnlstsolver%rtimeFactorisation)
        call sptils_doneStructure (rnlstsolver%p_rspaceTimePrec)
        call sptils_initStructure (rnlstsolver%p_rspaceTimePrec,ierror)
        call stat_stopTimer (rnlstsolver%rtimeFactorisation)
        if (ierror .gt. 0) then
          call output_line ('sptils_initStructure failed! Matrix structure invalid!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'nlstslv_solve')
          call sys_halt()
        end if

      end if
        
!      end if

      call stat_clearTimer (rtimeFactorisationStep)
      call stat_startTimer (rtimeFactorisationStep)
      call sptils_initData (rnlstsolver%p_rspaceTimePrec,ierror)
      call stat_stopTimer (rtimeFactorisationStep)
      if (ierror .gt. 0) then
        call output_line ('sptils_initData failed! Matrix singular!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'nlstslv_solve')
        call sys_halt()
      end if

      call stat_clearTimer (rtimerMGStep)
      call stat_startTimer (rtimerMGStep)
      call sptils_precondDefect (rnlstsolver%p_rspaceTimePrec,rd)
      call stat_stopTimer (rtimerMGStep)
      
      call sptils_doneData (rnlstsolver%p_rspaceTimePrec)
      
      if (rnlstsolver%ioutputLevel .ge. 1) &
        call output_lbrk ()

      ! Count the number of linear iterations and the time for
      ! preconditioning
      rnlstsolver%nlinearIterations = &
          rnlstsolver%nlinearIterations + rnlstsolver%p_rspaceTimePrec%iiterations
      rnlstsolver%nlinearIterationsSpace = rnlstsolver%nlinearIterationsSpace + &
          rnlstsolver%p_rspaceTimePrec%niteLinSolveSpace

      call stat_addtimers (rtimerMGStep,rnlstsolver%rtimerPreconditioner)
      call stat_addtimers (rtimeFactorisationStep,rnlstsolver%rtimeFactorisation)
      call stat_addtimers (rnlstsolver%p_rspaceTimePrec%rtimeSpacePrecond,&
          rnlstsolver%rtimeSpacePrecond)
      call stat_addtimers (rtimerPostProc,rnlstsolver%rtimePostprocessing)

      ! Sum up statistical multigrid data for if we have it.
      if (associated(rnlstsolver%p_rmgSolver)) then
        call stat_addtimers (rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSmoothing,&
            rnlstsolver%rtimeSmoothing)
        call stat_addtimers (rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSmoothingFinest,&
            rnlstsolver%rtimeSmoothingFinest)
        call stat_addtimers (rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeCoarseGridSolver,&
            rnlstsolver%rtimeCoarseGridSolver)
        call stat_addtimers (rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeLinearAlgebra,&
            rnlstsolver%rtimeLinearAlgebra)
        call stat_addtimers (rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeProlRest,&
            rnlstsolver%rtimeProlRest)
            
        rnlstsolver%ncoarsegridIterations = rnlstsolver%ncoarsegridIterations + &
            rnlstsolver%p_rspaceTimePrec%p_rsubnodeMultigrid%niteLinSolveCoarse

        ! Print some statistical output
        if (rnlstSolver%ioutputLevel .ge. 1) then
          call output_line ("Time for non-c.grid solving : "//&
              sys_sdL(rnlstsolver%p_rmgSolver%rtimeTotal%delapsedReal-&
                      rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeCoarseGridSolver%delapsedReal,10))
          call output_line ("Time for smoothing          : "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSmoothing%delapsedReal,10))
          call output_line ("Time for smoothing finest g.: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSmoothingFinest%delapsedReal,10))
          call output_line ("Time for coarse grid solving: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeCoarseGridSolver%delapsedReal,10))
          call output_line ("Time for linear algebra     : "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeLinearAlgebra%delapsedReal,10))
          call output_line ("Time for prol/rest          : "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeProlRest%delapsedReal,10))
          call output_line ("Time for def. asm. coarse gr: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSpaceDefectAssemblyCoarse%delapsedReal,10))
          call output_line ("Time for mat. asm. coarse gr: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSpaceMatrixAssemblyCoarse%delapsedReal,10))
          call output_line ("Time for def. asm. fine grid: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSpaceDefectAssemblyFine%delapsedReal,10))
          call output_line ("Time for mat. asm. fine grid: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSpaceMatrixAssemblyFine%delapsedReal,10))
          call output_line ("Time for def. asm. finest g.: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSpaceDefectAssemblyFinest%delapsedReal,10))
          call output_line ("Time for mat. asm. finest g.: "//&
              sys_sdL(rnlstsolver%p_rmgSolver%p_rsubnodeMultigrid%rtimeSpaceMatrixAssemblyFinest%delapsedReal,10))
          call output_lbrk ()
        end if
      end if
      
      if (rnlstSolver%ioutputLevel .ge. 1) then
        call output_line ("Time for prec. in space     : "//&
            sys_sdL(rnlstsolver%p_rspaceTimePrec%rtimeSpacePrecond%delapsedReal,10))
        call output_line ("Time for def. asm. in space : "//&
            sys_sdL(rnlstsolver%p_rspaceTimePrec%rtimeSpaceDefectAssembly%delapsedReal,10))
        call output_line ("Time for mat. asm. in space : "//&
            sys_sdL(rnlstsolver%p_rspaceTimePrec%rtimeSpaceMatrixAssembly%delapsedReal,10))
      end if

      call sptivec_vectorLinearComb (rd,rx,rnlstsolver%domega,1.0_DP)
      
      ! Normalise the primal and dual pressure to integral mean value zero
      ! where no Neumann boundary is present.
      call tbc_pressureToL20 (rsettings%roptcBDC,rx,rsettings%rglobalData)
      
      ! Remember the last defect norm for the stopping criterion
      dlastDefNorm = ddefNorm
      
      call output_separator (OU_SEP_EQUAL)
      if (rnlstsolver%ioutputLevel .ge. 2) call output_line('Nonlinear defect:')

      ! Assemble the new defect: d=b-Ax
      call sptivec_copyVector (rb,rd)
      call stlin_spaceTimeMatVec (rmatrix, rx, rd, &
          -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,&
          ddefNorm,rsettings%roptcBDC,rnlstsolver%ioutputLevel .ge. 2)
          
      if (rnlstsolver%ioutputLevel .ge. 1) &
          call output_separator (OU_SEP_MINUS)
      
      call stat_stopTimer (rtimerIterate)

      call output_line ('Iteration: '//trim(sys_siL(rnlstsolver%nnonlinearIterations,10))//&
          '                          Defect of supersystem:  '//adjustl(sys_sdEP(ddefNorm,20,10)))
      call output_line ('Time for computation of this iterate:       '//&
          trim(sys_sdL(rtimerIterate%delapsedReal-rtimerPostproc%delapsedReal,10)))
      call output_line ('Time for postprocessing in this iterate:    '//&
          trim(sys_sdL(rtimerPostproc%delapsedReal,10)))
      call output_line ('Total time for this nonlinear step:         '//&
          trim(sys_sdL(rtimerIterate%delapsedReal,10)))
          
      ! Measure the current wallclock
      call stat_sampleTimer(rnlstsolver%rtimerNonlinear,delapsedReal)
      call output_line ('Total time for nonlinear iteration:         '//trim(sys_sdL(delapsedReal,10)))
          
    end do

    if (rnlstsolver%ioutputLevel .ge. 1) &
      call output_separator (OU_SEP_EQUAL)

    call stat_stopTimer (rnlstsolver%rtimerNonlinear)
    
    ! Release the space-time preconditioner
    call sptils_doneStructure(rnlstsolver%p_rspaceTimePrec)
    
    ! Release the preconditioner matrices.
    call nlstslv_donePrecMatrices (rnlstsolver%p_RprecMatrices)
    
    ! Release other stuff
    call stnm_releaseNeumannBoundary (rsptiNeumannBC)
    call stdbcc_releaseDirichletBCCBd (rsptiDirichletBCC)
    call stlin_releaseSpaceTimeMatrix(rmatrix)
    
    ! Decrease rnlstsolver%nnonlinearIterations if the DO-loop was completely processed;
    ! iglobIter = nmaxIterations+1 in that case!
    if (rnlstsolver%nnonlinearIterations .gt. rnlstsolver%nmaxIterations) &
      rnlstsolver%nnonlinearIterations = rnlstsolver%nmaxIterations
    
    if (rnlstsolver%ioutputLevel .ge. 1) &
      call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))

    ! Value of the functional
    if (rnlstsolver%ioutputLevel .ge. 1) &
      call output_separator (OU_SEP_MINUS)
    call stat_clearTimer (rtimerPostproc)
    call stat_startTimer (rtimerPostproc)
    call optcana_nonstatFunctional (rsettings%rglobalData,rsettings%rphysicsPrimal,&
        rsettings%rsettingsOptControl%rconstraints,&
        rx,rsettings%rsettingsOptControl%rtargetFunction,&
        rsettings%rsettingsOptControl%dalphaC,rsettings%rsettingsOptControl%dbetaC,&
        rsettings%rsettingsOptControl%dgammaC,&
        Derror)
    if (rnlstsolver%ioutputLevel .ge. 1) then
      call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(1),10)))
      call output_line ('||u||         = '//trim(sys_sdEL(Derror(2),10)))
      call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
      call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(4),10)))
    end if

    ! If error analysis has to be performed, we can calculate
    ! the real error.
    if (rpostproc%icalcError .eq. 1) then
      if (rnlstsolver%ioutputLevel .ge. 1) then
        call optcana_analyticalError (rsettings%rglobalData,&
            rsettings%rsettingsOptControl%rconstraints,&
            rx,rpostproc%ranalyticRefFunction,&
            DerrorU,DerrorP,DerrorLambda,DerrorXi,.true.)
            
        call output_line ('||y-y0||_[0,T]           = '//trim(sys_sdEL(DerrorU(1),10)))
        call output_line ('||y-y0||_[0,T)           = '//trim(sys_sdEL(DerrorU(2),10)))
        call output_line ('||y-y0||_(0,T]           = '//trim(sys_sdEL(DerrorU(3),10)))
        call output_line ('||y-y0||_(0,T)           = '//trim(sys_sdEL(DerrorU(4),10)))
        
        call output_line ('||p-p0||_[0,T]           = '//trim(sys_sdEL(DerrorP(1),10)))
        call output_line ('||p-p0||_[0,T)           = '//trim(sys_sdEL(DerrorP(2),10)))
        call output_line ('||p-p0||_(0,T]           = '//trim(sys_sdEL(DerrorP(3),10)))
        call output_line ('||p-p0||_(0,T)           = '//trim(sys_sdEL(DerrorP(4),10)))
        
        call output_line ('||lambda-lambda0||_[0,T] = '//trim(sys_sdEL(DerrorLambda(1),10)))
        call output_line ('||lambda-lambda0||_[0,T) = '//trim(sys_sdEL(DerrorLambda(2),10)))
        call output_line ('||lambda-lambda0||_(0,T] = '//trim(sys_sdEL(DerrorLambda(3),10)))
        call output_line ('||lambda-lambda0||_(0,T) = '//trim(sys_sdEL(DerrorLambda(4),10)))

        call output_line ('||xi-xi0||_[0,T]         = '//trim(sys_sdEL(DerrorXi(1),10)))
        call output_line ('||xi-xi0||_[0,T)         = '//trim(sys_sdEL(DerrorXi(2),10)))
        call output_line ('||xi-xi0||_(0,T]         = '//trim(sys_sdEL(DerrorXi(3),10)))
        call output_line ('||xi-xi0||_(0,T)         = '//trim(sys_sdEL(DerrorXi(4),10)))
      end if
    end if
    call stat_stopTimer (rtimerPostproc)
    call stat_addtimers (rtimerPostProc,rnlstsolver%rtimePostprocessing)

  end subroutine


!  ! ***************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE cc_solveSupersysDirect (rproblem, rspaceTimeDiscr, rx, rd, &
!      rtempvectorX, rtempvectorB, rtempvectorD)
!
!!<description>
!  ! This routine assembles and solves the time-space coupled supersystem:
!  ! $Ax=b$. The RHS vector is generated on-the-fly.
!  ! The routine generates the full matrix in memory and solves with UMFPACK,
!  ! so it should only be used for debugging!
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!
!  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
!  ! coupled space-time matrix.
!  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!!</input>
!
!!<inputoutput>
!  ! A space-time vector defining the current solution.
!  ! Is replaced by the new solution
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx
!
!  ! A temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorX
!
!  ! A second temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorB
!
!  ! A third temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVectorD
!
!  ! A space-time vector that receives the defect.
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    INTEGER :: isubstep,ilevel,ierror,i
!    TYPE(t_matrixBlock) :: rblockTemp
!    TYPE(t_matrixScalar) :: rmassLumped
!    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
!    TYPE(t_vectorBlock) :: rxGlobal, rbGlobal, rdGlobal
!    TYPE(t_vectorBlock) :: rxGlobalSolve, rbGlobalSolve, rdGlobalSolve
!    TYPE(t_matrixBlock) :: rglobalA
!    TYPE(t_linsolNode), POINTER :: rsolverNode
!    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
!    integer, DIMENSION(:), ALLOCATABLE :: Isize
!    integer, DIMENSION(6) :: Isize2
!
!    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd
!
!    ! If the following constant is set from 1.0 to 0.0, the primal system is
!    ! decoupled from the dual system!
!    REAL(DP), PARAMETER :: dprimalDualCoupling = 0.0 !1.0
!
!    ilevel = rspaceTimeDiscr%ilevel
!
!    ! Calculate the lumped mass matrix of the FE space -- we need it later!
!    CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass,&
!        rmassLumped, LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
!    CALL lsyssc_lumpMatrixScalar (rmassLumped,LSYSSC_LUMP_STD)
!
!    ! ----------------------------------------------------------------------
!    ! 1.) Generate the global RHS vector
!
!    CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
!    CALL lsysbl_getbase_double (rtempVectorB,p_Db)
!    CALL lsysbl_getbase_double (rtempVectorD,p_Dd)
!
!    DO isubstep = 0,rspaceTimeDiscr%niterations
!
!      ! Current point in time
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%dtstep
!
!      ! Generate the RHS of that point in time.
!      CALL cc_generateBasicRHS (rproblem,rtempVectorB)
!
!      ! Multiply the RHS of the dual velocity by dtstep according to the time
!      ! discretisation -- except for if we are in the last timestep!
!      !
!      ! In the last timestep, if GAMMA is =0, we have no terminal condition
!      ! and thus have to force the RHS to 0!
!      ! Otherwise, the terminat condition is multiplied with dgammaC.
!      IF (isubstep .NE. rspaceTimeDiscr%niterations) THEN
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(4),rspaceTimeDiscr%dtstep)
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(5),rspaceTimeDiscr%dtstep)
!      ELSE
!        ! Multiply -z by gamma, that's it.
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(4),rspaceTimeDiscr%dgammaC)
!        CALL lsyssc_scaleVector (rtempVectorB%RvectorBlock(5),rspaceTimeDiscr%dgammaC)
!      END IF
!
!      ! Initialise the collection for the assembly process with callback routines.
!      ! Basically, this stores the simulation time in the collection if the
!      ! simulation is nonstationary.
!      CALL user_initCollectForAssembly (rproblem,rproblem%rcollection)
!
!      ! Discretise the boundary conditions at the new point in time --
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL cc_updateDiscreteBC (rproblem)
!      END IF
!
!      ! Implement the boundary conditions into the RHS.
!      ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
!      ! as Dirichlet values mustn't be multiplied with GAMMA!
!      CALL vecfil_discreteBCsol (rtempVectorB)
!      CALL vecfil_discreteFBCsol (rtempVectorB)
!
!      CALL sptivec_setTimestepData(rd, isubstep, rtempVectorB)
!
!      ! Clean up the collection (as we are done with the assembly, that's it.
!      CALL user_doneCollectForAssembly (rproblem,rproblem%rcollection)
!
!    END DO
!
!    ! Release the mass matr, we don't need it anymore
!    CALL lsyssc_releaseMatrix (rmassLumped)
!
!    ! ----------------------------------------------------------------------
!    ! 2.) Generate the matrix A
!    !
!    ! Create a global matrix:
!    CALL lsysbl_createEmptyMatrix (rglobalA,6*(rspaceTimeDiscr%niterations+1))
!
!    ! Loop through the substeps
!
!    DO isubstep = 0,rspaceTimeDiscr%niterations
!
!      ! Current point in time
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%dtstep
!
!      ! -----
!      ! Discretise the boundary conditions at the new point in time --
!      ! if the boundary conditions are nonconstant in time!
!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!        CALL cc_updateDiscreteBC (rproblem)
!      END IF
!
!      ! The first and last substep is a little bit special concerning
!      ! the matrix!
!      IF (isubstep .EQ. 0) THEN
!
!        ! We are in the first substep
!
!        ! -----
!
!        ! Create a matrix that applies "-M" to the dual velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(4,4)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(5,5)%dscaleFactor = -1.0_DP
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!
!        ! Include "-M" in the global matrix at position (1,2).
!        CALL insertMatrix (rblockTemp,rglobalA,7,1)
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----
!
!        ! Now the hardest -- or longest -- part: The diagonal matrix.
!        !
!        ! Generate the basic system matrix level rspaceTimeDiscr%ilevel
!        ! Will be modified by cc_assembleLinearisedMatrices later.
!        CALL cc_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!
!        CALL cc_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%ilevel,rtempVectorX,rtempVectorB,&
!            rnonlinearIterationTmp,'CC2D-NONLINEAR')
!
!        ! Set up all the weights in the core equation according to the current timestep.
!        rnonlinearIterationTmp%diota1 = 1.0_DP
!        rnonlinearIterationTmp%diota2 = 0.0_DP
!
!        rnonlinearIterationTmp%dkappa1 = 1.0_DP
!        rnonlinearIterationTmp%dkappa2 = 0.0_DP
!
!        rnonlinearIterationTmp%dalpha1 = 0.0_DP
!        rnonlinearIterationTmp%dalpha2 = 1.0_DP
!
!        rnonlinearIterationTmp%dtheta1 = 0.0_DP
!        rnonlinearIterationTmp%dtheta2 = rspaceTimeDiscr%dtstep
!
!        rnonlinearIterationTmp%dgamma1 = 0.0_DP
!        rnonlinearIterationTmp%dgamma2 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%cequation,DP)
!
!        rnonlinearIterationTmp%deta1 = 0.0_DP
!        rnonlinearIterationTmp%deta2 = rspaceTimeDiscr%dtstep
!
!        rnonlinearIterationTmp%dtau1 = 0.0_DP
!        rnonlinearIterationTmp%dtau2 = 1.0_DP
!
!        rnonlinearIterationTmp%dmu1 = 0.0_DP
!        rnonlinearIterationTmp%dmu2 = -rspaceTimeDiscr%dtstep
!
!        ! Assemble the system matrix on level rspaceTimeDiscr%ilevel.
!        ! Include the boundary conditions into the matrices.
!        CALL cc_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
!            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
!
!        ! Insert the system matrix for the dual equation to our global matrix.
!        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
!            rglobalA,1,1)
!
!      ELSE IF (isubstep .EQ. rspaceTimeDiscr%niterations) THEN
!
!        ! We are in the last substep
!
!        ! -----
!
!        ! Create a matrix that applies "-M" to the primal velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(1,1)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(2,2)%dscaleFactor = -1.0_DP
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!
!        ! Include "-M" in the global matrix at position (2,1).
!        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----
!
!        ! Now the hardest -- or longest -- part: The diagonal matrix.
!        !
!        ! Generate the basic system matrix level rspaceTimeDiscr%ilevel
!        ! Will be modified by cc_assembleLinearisedMatrices later.
!        CALL cc_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!
!        CALL cc_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%ilevel,rtempVectorX,rtempVectorB,&
!            rnonlinearIterationTmp,'CC2D-NONLINEAR')
!
!        ! Set up all the weights in the core equation according to the current timestep.
!        rnonlinearIterationTmp%diota1 = 0.0_DP
!        rnonlinearIterationTmp%diota2 = 0.0_DP
!
!        rnonlinearIterationTmp%dkappa1 = 0.0_DP
!        rnonlinearIterationTmp%dkappa2 = 1.0_DP
!
!        rnonlinearIterationTmp%dalpha1 = 1.0_DP
!        rnonlinearIterationTmp%dalpha2 = 1.0_DP
!
!        rnonlinearIterationTmp%dtheta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%dtheta2 = 0.0_DP
!
!        rnonlinearIterationTmp%dgamma1 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%cequation,DP)
!        rnonlinearIterationTmp%dgamma2 = 0.0_DP
!
!        rnonlinearIterationTmp%deta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%deta2 = 0.0_DP
!
!        rnonlinearIterationTmp%dtau1 = 1.0_DP
!        rnonlinearIterationTmp%dtau2 = 0.0_DP
!
!        rnonlinearIterationTmp%dmu1 = dprimalDualCoupling * &
!            rspaceTimeDiscr%dtstep / rspaceTimeDiscr%dalphaC
!        rnonlinearIterationTmp%dmu2 = -rspaceTimeDiscr%dgammaC
!
!        ! Assemble the system matrix on level rspaceTimeDiscr%ilevel.
!        ! Include the boundary conditions into the matrices.
!        CALL cc_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
!            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
!
!        ! Insert the system matrix for the dual equation to our global matrix.
!        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
!            rglobalA,isubstep*6+1,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!      ELSE
!
!        ! We are sonewhere in the middle of the matrix. There is a substep
!        ! isubstep+1 and a substep isubstep-1!
!
!        ! -----
!
!        ! Create a matrix that applies "-M" to the dual velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(4,4)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        rblockTemp%RmatrixBlock(5,5)%dscaleFactor = -1.0_DP
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!
!        ! Include that in the global matrix below the diagonal
!        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+7,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----
!
!        ! Create a matrix that applies "-M" to the primal velocity and include it
!        ! to the global matrix.
!        CALL lsysbl_createMatBlockByDiscr (rtempVectorX%p_rblockDiscretisation,&
!            rblockTemp)
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(1,1)%dscaleFactor = -1.0_DP
!        CALL lsyssc_duplicateMatrix (rspaceTimeDiscr%p_rlevelInfo%rmatrixMass, &
!            rblockTemp%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        rblockTemp%RmatrixBlock(2,2)%dscaleFactor = -1.0_DP
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        CALL matfil_discreteBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteBC)
!        CALL matfil_discreteFBC (rblockTemp,rproblem%RlevelInfo(ilevel)%p_rdiscreteFBC)
!
!        ! Include that in the global matrix above the diagonal
!        CALL insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        CALL lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----
!
!        ! Now the hardest -- or longest -- part: The diagonal matrix.
!        !
!        ! Generate the basic system matrix level rspaceTimeDiscr%ilevel
!        ! Will be modified by cc_assembleLinearisedMatrices later.
!        CALL cc_generateStaticSystemMatrix (rproblem%RlevelInfo(ilevel), &
!            rproblem%RlevelInfo(ilevel)%rmatrix,.FALSE.)
!
!        ! Set up a core equation structure and assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!
!        CALL cc_initNonlinearLoop (&
!            rproblem,rproblem%NLMIN,rspaceTimeDiscr%ilevel,rtempVectorX,rtempVectorB,&
!            rnonlinearIterationTmp,'CC2D-NONLINEAR')
!
!        ! Set up all the weights in the core equation according to the current timestep.
!        rnonlinearIterationTmp%diota1 = 0.0_DP
!        rnonlinearIterationTmp%diota2 = 0.0_DP
!
!        rnonlinearIterationTmp%dkappa1 = 0.0_DP
!        rnonlinearIterationTmp%dkappa2 = 0.0_DP
!
!        rnonlinearIterationTmp%dalpha1 = 1.0_DP
!        rnonlinearIterationTmp%dalpha2 = 1.0_DP
!
!        rnonlinearIterationTmp%dtheta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%dtheta2 = rspaceTimeDiscr%dtstep
!
!        rnonlinearIterationTmp%dgamma1 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%cequation,DP)
!        rnonlinearIterationTmp%dgamma2 = rspaceTimeDiscr%dtstep * REAL(1-rproblem%cequation,DP)
!
!        rnonlinearIterationTmp%deta1 = rspaceTimeDiscr%dtstep
!        rnonlinearIterationTmp%deta2 = rspaceTimeDiscr%dtstep
!
!        rnonlinearIterationTmp%dtau1 = 1.0_DP
!        rnonlinearIterationTmp%dtau2 = 1.0_DP
!
!        rnonlinearIterationTmp%dmu1 = dprimalDualCoupling * &
!            rspaceTimeDiscr%dtstep / rspaceTimeDiscr%dalphaC
!        rnonlinearIterationTmp%dmu2 = -rspaceTimeDiscr%dtstep
!
!        ! Assemble the system matrix on level rspaceTimeDiscr%ilevel.
!        ! Include the boundary conditions into the matrices.
!        CALL cc_assembleLinearisedMatrices (rnonlinearIterationTmp,rproblem%rcollection,&
!            .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
!
!        ! Insert the system matrix for the dual equation to our global matrix.
!        CALL insertMatrix (rnonlinearIterationTmp%RcoreEquation(ilevel)%p_rmatrix,&
!            rglobalA,isubstep*6+1,isubstep*6+1)
!
!      END IF
!
!    END DO
!
!    ! Update structural information of the global matrix.
!    CALL lsysbl_updateMatStrucInfo (rglobalA)
!
!    ! Write the global matrix to a file.
!    !CALL matio_writeBlockMatrixHR(rglobalA,'MATRIX',.TRUE.,0,'matrix.txt','(E13.2)')
!
!    ! Get the global solution/rhs/temp vector.
!    CALL sptivec_convertSupervecToVector (rx, rxGlobal)
!    CALL sptivec_convertSupervecToVector (rd, rbGlobal)
!    CALL lsysbl_createVecBlockIndirect (rbGlobal,rdGlobal,.FALSE.)
!
!    ! Initialise the UMFPACK solver.
!    CALL linsol_initUMFPACK4 (rsolverNode)
!
!    ! Add the matrices
!    Rmatrices(1) = rglobalA
!    CALL linsol_setMatrices (rsolverNode,Rmatrices)
!
!    ! Init the solver
!    CALL linsol_initStructure (rsolverNode,ierror)
!    CALL linsol_initData (rsolverNode,ierror)
!
!    ! Reshape the x,b and d-vector to fit to our matrix.
!    CALL lsysbl_deriveSubvector(rxGlobal,rxGlobalSolve,bshare=.TRUE.)
!    CALL lsysbl_deriveSubvector(rbGlobal,rbGlobalSolve,bshare=.TRUE.)
!    CALL lsysbl_deriveSubvector(rdGlobal,rdGlobalSolve,bshare=.TRUE.)
!    ALLOCATE(Isize(rxGlobal%nblocks*6))
!    Isize2= (/rtempVectorB%RvectorBlock(1)%NEQ,&
!              rtempVectorB%RvectorBlock(2)%NEQ,&
!              rtempVectorB%RvectorBlock(3)%NEQ,&
!              rtempVectorB%RvectorBlock(4)%NEQ,&
!              rtempVectorB%RvectorBlock(5)%NEQ,&
!              rtempVectorB%RvectorBlock(6)%NEQ &
!             /)
!    DO i=0,rxGlobal%nblocks-1
!      Isize(i*6+1:i*6+6) = Isize2(1:6)
!    END DO
!    CALL lsysbl_enforceStructureDirect (Isize,rxGlobalSolve)
!    CALL lsysbl_enforceStructureDirect (Isize,rbGlobalSolve)
!    CALL lsysbl_enforceStructureDirect (Isize,rdGlobalSolve)
!
!    ! Solve
!    CALL lsysbl_getbase_double (rxGlobalSolve,p_Dx)
!    CALL lsysbl_getbase_double (rbGlobalSolve,p_Db)
!    CALL linsol_solveAdaptively (rsolverNode,rxGlobalSolve,rbGlobalSolve,rdGlobalSolve)
!
!    ! Release
!    CALL lsysbl_releaseVector (rxGlobalSolve)
!    CALL lsysbl_releaseVector (rbGlobalSolve)
!    CALL lsysbl_releaseVector (rdGlobalSolve)
!
!    CALL linsol_releaseSolver (rsolverNode)
!    CALL lsysbl_releaseVector (rdGlobal)
!    CALL lsysbl_releaseVector (rbGlobal)
!
!    ! Remember the solution
!    CALL sptivec_convertVectorToSupervec (rxGlobal, rx)
!
!    ! Release the global matrix
!    CALL lsysbl_releaseMatrix (rglobalA)
!
!    CALL lsysbl_releaseVector (rxGlobal)
!
!  CONTAINS
!
!    SUBROUTINE insertMatrix (rsource,rdest,ileft,itop)
!
!    ! Includes rsource into rdest at position ileft,itop
!    TYPE(t_matrixBlock), INTENT(IN) :: rsource
!    TYPE(t_matrixBlock), INTENT(INOUT) :: rdest
!    INTEGER, INTENT(IN) :: ileft
!    INTEGER, INTENT(IN) :: itop
!
!    INTEGER :: i,j
!
!    DO j=1,rsource%ndiagBlocks
!      DO i=1,rsource%ndiagBlocks
!        IF (lsysbl_isSubmatrixPresent (rsource,i,j)) THEN
!          IF (lsysbl_isSubmatrixPresent (rdest,i+itop-1,j+ileft-1)) THEN
!            CALL lsyssc_releaseMatrix (rdest%RmatrixBlock(j+ileft-1,i+itop-1))
!          END IF
!          CALL lsyssc_duplicateMatrix (rsource%RmatrixBlock(i,j),&
!              rdest%RmatrixBlock(i+itop-1,j+ileft-1),&
!              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        END IF
!      END DO
!    END DO
!
!    END SUBROUTINE
!
!  END SUBROUTINE
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_applyUmfpackToSupersystem (rproblem, rspaceTimeMatrix, rx, rd, &
!      rtempvectorX, rtempvectorB, rtempvectorD)
!
!!<description>
!  ! This routine assembles and solves the time-space coupled supersystem:
!  ! $Ax=b$. The RHS vector is generated on-the-fly.
!  ! The routine generates the full matrix in memory and solves with UMFPACK,
!  ! so it should only be used for debugging!
!  ! The routine assembles the global matrix with the time step scheme
!  ! specified in dtimeStepTheta in the problem structure.
!  !
!  ! rspaceTimeMatrix specifies the global space time matrix of the system.
!  ! In this routine, there is no processing of any nonlinearity, so only
!  ! linear problems can be solved (e.g. Stokes).
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  type(t_problem), intent(INOUT), target :: rproblem
!
!  ! The definition of the global space time matrix that should be used
!  ! for solving the system.
!  type(t_ccoptSpaceTimeMatrix), intent(IN),target :: rspaceTimeMatrix
!!</input>
!
!!<inputoutput>
!  ! A space-time vector defining the current solution.
!  ! Is replaced by the new solution
!  type(t_spacetimeVector), intent(INOUT) :: rx
!
!  ! A temporary vector in the size of a spatial vector.
!  type(t_vectorBlock), intent(INOUT) :: rtempVectorX
!
!  ! A second temporary vector in the size of a spatial vector.
!  type(t_vectorBlock), intent(INOUT) :: rtempVectorB
!
!  ! A third temporary vector in the size of a spatial vector.
!  type(t_vectorBlock), intent(INOUT) :: rtempVectorD
!
!  ! A space-time vector that receives the defect.
!  type(t_spacetimeVector), intent(INOUT) :: rd
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: isubstep,ilevel,ierror,i
!    type(t_matrixBlock) :: rblockTemp
!    type(t_vectorBlock) :: rxGlobal, rbGlobal, rdGlobal
!    type(t_vectorBlock) :: rxGlobalSolve, rbGlobalSolve, rdGlobalSolve
!    type(t_matrixBlock) :: rglobalA
!    type(t_linsolNode), pointer :: rsolverNode,p_rpreconditioner
!    type(t_matrixBlock), dimension(1) :: Rmatrices
!    integer, dimension(:), allocatable :: Isize
!    integer, dimension(6) :: Isize2
!    real(DP) :: dtheta,dtime
!    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
!    type(t_matrixBlock) :: rmatrix
!    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
!    type(t_vectorBlock), dimension(3) :: rtempVectorSol
!    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol
!    type(t_ccPreconditionerSpecials) :: rprecSpecials
!    type(t_discreteBC) :: rdiscreteBC
!    type(t_discreteFBC) :: rdiscreteFBC
!
!    real(DP), dimension(:),pointer :: p_Dx, p_Db, p_Dd
!
!    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscr
!
!    ilevel = p_rspaceTimeDiscr%ilevel
!
!    ! Theta-scheme identifier.
!    ! =1: impliciz Euler.
!    ! =0.5: Crank Nicolson
!    dtheta = rproblem%rtimedependence%dtimeStepTheta
!
!    ! Initialise the collection for the assembly.
!    call user_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)
!
!    ! Initialise the boundary conditions
!    call bcasm_initDiscreteBC(rdiscreteBC)
!    call bcasm_initDiscreteFBC(rdiscreteFBC)
!
!    ! Generate the RHS for the initial condition.
!    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rinitialCondRHS,.false.)
!    call sptivec_getTimestepData (rx, 1, rtempVectorX)
!    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
!    call stlin_generateInitCondRHS (rproblem,p_rspaceTimeDiscr,&
!        rtempVectorX,rinitialCondRHS)
!
!    ! Assemble the space-time RHS into rd.
!    call trhsevl_assembleRHS (rproblem, p_rspaceTimeDiscr, rd, .true.)
!
!    ! Implement the initial condition into the RHS/solution.
!    call tbc_implementInitCondRHS (rproblem, rd, rinitialCondRHS, rtempvectorD)
!    call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvectorD)
!
!    ! Release the rhs vector with the init. condition again.
!    call lsysbl_releaseVector (rinitialCondRHS)
!    call lsysbl_releaseVector (rinitialCondSol)
!
!    ! ----------------------------------------------------------------------
!    ! 2.) Generate the matrix A
!    !
!    ! Create a global matrix:
!    call lsysbl_createEmptyMatrix (rglobalA,6*(p_rspaceTimeDiscr%rtimeDiscr%nintervals+1))
!
!    ! Get a temporary system matrix. Use the default preconditioner specials, that's
!    ! enough for us.
!    call cc_allocPrecSystemMatrix (rproblem,rprecSpecials,&
!        p_rspaceTimeDiscr%ilevel,rproblem%nlmin,rproblem%nlmax,&
!        p_rspaceTimeDiscr%p_rlevelInfo,CCMASM_MTP_AUTOMATIC,&
!        rblockTemp)
!
!    ! Solution vectors -- for the nonlinearity (if there is one).
!    ! For the previous (1), current (2) and next (3) time step.
!    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVectorSol(1),.true.)
!    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVectorSol(2),.true.)
!    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVectorSol(3),.true.)
!
!    ! Load the last solution vectors for handling the nonlinearity.
!    ! rtempVectorSol(2) holds the data for the current timestep.
!    ! rtempVectorSol(1) holds the data from the previous and
!    ! rtempVectorSol(3) that of the next timestep.
!    if (associated(rspaceTimeMatrix%p_rsolution)) then
!      call sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, &
!          0, rtempVectorSol(2))
!    end if
!
!    ! Loop through the substeps
!
!    do isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
!
!      ! Current point in time
!      dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*p_rspaceTimeDiscr%rtimeDiscr%dtstep
!
!      ! -----
!      ! Discretise the boundary conditions at the new point in time.
!      call bcasm_clearDiscreteBC(rdiscreteBC)
!      call bcasm_clearDiscreteFBC(rdiscreteFBC)
!      call sbc_assembleBDconditions (rproblem,dtime,&
!          p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!          CCSPACE_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
!      call sbc_assembleFBDconditions (rproblem,dtime,&
!          p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!          CCSPACE_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
!
!      ! The first and last substep is a little bit special concerning
!      ! the matrix!
!      if (isubstep .eq. 0) then
!
!        ! We are in the first substep
!
!        ! -----
!
!        ! Get the evaluation point for the nonlinearity in the next timestep
!        ! into rtempVectorSol(3).
!        if (associated(rspaceTimeMatrix%p_rsolution)) then
!          call sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, &
!              1+isubstep+1, rtempVectorSol(3))
!        end if
!
!        ! The diagonal matrix.
!
!        ! Set up the matrix weights of that submatrix.
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,0,rnonlinearSpatialMatrix)
!
!        ! Assemble the matrix. No 'previous' solution vector.
!        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
!            rmatrix,rnonlinearSpatialMatrix,&
!            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
!
!        ! Assemble the system matrix on level p_rspaceTimeDiscr%ilevel.
!        ! Include the boundary conditions into the matrices.
!        !CALL cc_assembleLinearisedMatrices (&
!        !    rnonlinearIterationTmp,rproblem%rcollection,&
!        !    .FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
!
!        ! Insert the system matrix for the dual equation to our global matrix.
!        call insertMatrix (rmatrix,rglobalA,1,1)
!
!        ! -----
!
!        ! Create the matrix
!        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
!        ! and include that into the global matrix for the primal velocity.
!
!        ! Set up the matrix weights of that submatrix.
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,1,rnonlinearSpatialMatrix)
!
!        ! Assemble the matrix
!        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
!            rmatrix,rnonlinearSpatialMatrix,&
!            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
!
!        call lsysbl_duplicateMatrix (&
!            rmatrix,&
!            rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        call matfil_discreteBC (rblockTemp,rdiscreteBC)
!        call matfil_discreteFBC (rblockTemp,rdiscreteFBC)
!
!        ! We don't need submatrix (1,1) to (2,2).
!        rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP
!
!        ! Include that in the global matrix below the diagonal
!        call insertMatrix (rblockTemp,rglobalA,isubstep*6+1+6,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        call lsysbl_releaseMatrix (rblockTemp)
!
!      else if (isubstep .lt. p_rspaceTimeDiscr%NEQtime-1) then
!
!        ! We are sonewhere in the middle of the matrix. There is a substep
!        ! isubstep+1 and a substep isubstep-1!
!
!        ! -----
!
!        ! Create the matrix
!        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
!        ! and include that into the global matrix for the primal velocity.
!
!        ! Set up the matrix weights of that submatrix.
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,-1,rnonlinearSpatialMatrix)
!
!        ! Assemble the matrix
!        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
!            rmatrix,rnonlinearSpatialMatrix,&
!            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
!
!        call lsysbl_duplicateMatrix (&
!            rmatrix,rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!        ! We don't need submatrix (4,4) to (5,5).
!        rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        call matfil_discreteBC (rblockTemp,rdiscreteBC)
!        call matfil_discreteFBC (rblockTemp,rdiscreteFBC)
!
!        ! Include that in the global matrix below the diagonal
!        call insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        call lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----
!
!        ! Now the diagonal matrix.
!
!        ! Assemble the nonlinear defect.
!        ! We use explicit Euler, so the weights are easy.
!
!        ! Set up the matrix weights of that submatrix.
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,0,rnonlinearSpatialMatrix)
!
!        ! Assemble the matrix
!        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
!            rmatrix,rnonlinearSpatialMatrix,&
!            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
!
!        ! Insert the system matrix for the dual equation to our global matrix.
!        call insertMatrix (rmatrix,&
!            rglobalA,isubstep*6+1,isubstep*6+1)
!
!        ! -----
!
!        ! Create the matrix
!        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
!        ! and include that into the global matrix for the dual velocity.
!
!        ! Set up the matrix weights of that submatrix.
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,1,rnonlinearSpatialMatrix)
!
!        ! Assemble the matrix
!        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
!            rmatrix,rnonlinearSpatialMatrix,&
!            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
!
!        call lsysbl_duplicateMatrix (&
!            rmatrix,rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!        ! We don't need submatrix (1,1) to (2,2).
!        rblockTemp%RmatrixBlock(1:2,1:2)%dscaleFactor = 0.0_DP
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        call matfil_discreteBC (rblockTemp,rdiscreteBC)
!        call matfil_discreteFBC (rblockTemp,rdiscreteFBC)
!
!        ! Include that in the global matrix above the diagonal
!        call insertMatrix (rblockTemp,rglobalA,isubstep*6+1+6,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        call lsysbl_releaseMatrix (rblockTemp)
!
!      else
!
!        ! We are in the last substep
!
!        ! -----
!
!        ! Create the matrix
!        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
!        ! and include that into the global matrix for the dual velocity.
!
!        ! Set up the matrix weights of that submatrix.
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,-1,rnonlinearSpatialMatrix)
!
!        ! Assemble the matrix
!        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
!            rmatrix,rnonlinearSpatialMatrix,&
!            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
!
!        call lsysbl_duplicateMatrix (&
!            rmatrix,rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!        ! We don't need submatrix (4,4) to (5,5).
!        rblockTemp%RmatrixBlock(4:5,4:5)%dscaleFactor = 0.0_DP
!
!        ! Include the boundary conditions into that matrix.
!        ! Specify the matrix as 'off-diagonal' matrix because it's not on the
!        ! main diagonal of the supermatrix.
!        rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
!        call matfil_discreteBC (rblockTemp,rdiscreteBC)
!        call matfil_discreteFBC (rblockTemp,rdiscreteFBC)
!
!        ! Include that in the global matrix above the diagonal
!        call insertMatrix (rblockTemp,rglobalA,isubstep*6+1-6,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        call lsysbl_releaseMatrix (rblockTemp)
!
!        ! -----
!
!        ! The diagonal matrix.
!
!        ! Set up the matrix weights of that submatrix.
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,0,rnonlinearSpatialMatrix)
!
!        ! Assemble the matrix
!        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
!            rmatrix,rnonlinearSpatialMatrix,&
!            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
!
!        ! Insert the system matrix for the dual equation to our global matrix.
!        call insertMatrix (rmatrix,rglobalA,isubstep*6+1,isubstep*6+1)
!
!        ! Release the block mass matrix.
!        call lsysbl_releaseMatrix (rblockTemp)
!
!      end if
!
!      ! Shift the evaluation vectors: 1 <- 2 <- 3
!      if (associated(rspaceTimeMatrix%p_rsolution)) then
!        call lsysbl_copyVector (rtempVectorSol(2),rtempVectorSol(1))
!        call lsysbl_copyVector (rtempVectorSol(3),rtempVectorSol(2))
!      end if
!
!    end do
!
!    ! Release the temp matrix and vectors
!    call lsysbl_releaseMatrix (rmatrix)
!
!    ! Update structural information of the global matrix.
!    call lsysbl_updateMatStrucInfo (rglobalA)
!    call lsysbl_releaseVector (rtempVectorSol(3))
!    call lsysbl_releaseVector (rtempVectorSol(2))
!    call lsysbl_releaseVector (rtempVectorSol(1))
!
!    ! Write the global matrix to a file.
!    !CALL matio_writeBlockMatrixHR(rglobalA,'MATRIX',.TRUE.,0,'matrixcn.txt','(1X,E20.10)')
!    !'(E13.2)')
!
!    ! Get the global solution/rhs/temp vector.
!    call sptivec_convertSupervecToVector (rx, rxGlobal)
!    call sptivec_convertSupervecToVector (rd, rbGlobal)
!    call lsysbl_createVecBlockIndirect (rbGlobal,rdGlobal,.false.)
!
!    ! Initialise the UMFPACK solver.
!    !CALL linsol_initUMFPACK4 (rsolverNode)
!    call linsol_initUMFPACK4 (p_rpreconditioner)
!    call linsol_initDefCorr (rsolverNode,p_rpreconditioner)
!    rsolverNode%ioutputLevel = 2
!    rsolverNode%depsRel = 1E-10
!    rsolverNode%nmaxIterations = 10
!
!    ! Add the matrices
!    Rmatrices(1) = rglobalA
!    call linsol_setMatrices (rsolverNode,Rmatrices)
!
!    ! Init the solver
!    call linsol_initStructure (rsolverNode,ierror)
!    call linsol_initData (rsolverNode,ierror)
!
!    ! Reshape the x,b and d-vector to fit to our matrix.
!    call lsysbl_deriveSubvector(rxGlobal,rxGlobalSolve,bshare=.true.)
!    call lsysbl_deriveSubvector(rbGlobal,rbGlobalSolve,bshare=.true.)
!    call lsysbl_deriveSubvector(rdGlobal,rdGlobalSolve,bshare=.true.)
!    allocate(Isize(rxGlobal%nblocks*6))
!    Isize2= (/rtempVectorB%RvectorBlock(1)%NEQ,&
!              rtempVectorB%RvectorBlock(2)%NEQ,&
!              rtempVectorB%RvectorBlock(3)%NEQ,&
!              rtempVectorB%RvectorBlock(4)%NEQ,&
!              rtempVectorB%RvectorBlock(5)%NEQ,&
!              rtempVectorB%RvectorBlock(6)%NEQ &
!             /)
!    do i=0,rxGlobal%nblocks-1
!      Isize(i*6+1:i*6+6) = Isize2(1:6)
!    end do
!    call lsysbl_enforceStructureDirect (Isize,rxGlobalSolve)
!    call lsysbl_enforceStructureDirect (Isize,rbGlobalSolve)
!    call lsysbl_enforceStructureDirect (Isize,rdGlobalSolve)
!
!    ! Solve
!    call lsysbl_getbase_double (rxGlobalSolve,p_Dx)
!    call lsysbl_getbase_double (rbGlobalSolve,p_Db)
!    call linsol_solveAdaptively (rsolverNode,rxGlobalSolve,rbGlobalSolve,rdGlobalSolve)
!
!    ! DEBUG!!!
!    call lsysbl_blockMatVec (rglobalA,rxGlobalSolve,rbGlobalSolve,-1.0_DP,1.0_DP)
!    !CALL lsyssc_scalarMatVec (rglobalA%RmatrixBlock(16,13),&
!    !                          rxGlobalSolve%RvectorBlock(13),&
!    !                          rbGlobalSolve%RvectorBlock(16),1.0_DP,-1.0_DP)
!
!    ! Release
!    call lsysbl_releaseVector (rxGlobalSolve)
!    call lsysbl_releaseVector (rbGlobalSolve)
!    call lsysbl_releaseVector (rdGlobalSolve)
!
!    call linsol_releaseSolver (rsolverNode)
!    call lsysbl_releaseVector (rdGlobal)
!    call lsysbl_releaseVector (rbGlobal)
!
!    ! Remember the solution
!    call sptivec_convertVectorToSupervec (rxGlobal, rx)
!
!    ! Release the global matrix
!    call lsysbl_releaseMatrix (rglobalA)
!
!    call lsysbl_releaseVector (rxGlobal)
!
!    ! Release the BC's again.
!    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
!    call bcasm_releaseDiscreteBC(rdiscreteBC)
!
!  contains
!
!    subroutine insertMatrix (rsource,rdest,ileft,itop)
!
!    ! Includes rsource into rdest at position ileft,itop
!    type(t_matrixBlock), intent(IN) :: rsource
!    type(t_matrixBlock), intent(INOUT) :: rdest
!    integer, intent(IN) :: ileft
!    integer, intent(IN) :: itop
!
!    integer :: i,j
!
!    do j=1,rsource%nblocksPerRow
!      do i=1,rsource%nblocksPerCol
!        if (lsysbl_isSubmatrixPresent (rsource,i,j)) then
!          if (lsysbl_isSubmatrixPresent (rdest,i+itop-1,j+ileft-1)) then
!            call lsyssc_releaseMatrix (rdest%RmatrixBlock(j+ileft-1,i+itop-1))
!          end if
!          call lsyssc_duplicateMatrix (rsource%RmatrixBlock(i,j),&
!              rdest%RmatrixBlock(i+itop-1,j+ileft-1),&
!              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!        end if
!      end do
!    end do
!
!    end subroutine
!
!  end subroutine
!
!!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_solveSupersystemDirect (rproblem, rspaceTimeDiscr, rx, rb, rd)
!
!!<description>
!  ! This routine assembles and solves the time-space coupled supersystem:
!  ! $Ax=b$. The RHS vector is generated on-the-fly.
!  ! The routine generates the full matrix in memory and solves with UMFPACK,
!  ! so it should only be used for debugging!
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  type(t_problem), intent(INOUT) :: rproblem
!
!  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
!  ! coupled space-time matrix.
!  type(t_ccoptSpaceTimeDiscretisation), intent(IN),target :: rspaceTimeDiscr
!!</input>
!
!!<inputoutput>
!  ! A space-time vector defining the initial solution. Is replaced by a new
!  ! solution vector.
!  type(t_spacetimeVector), intent(INOUT) :: rx
!
!  ! A temporary space-time vector that receives the RHS during the calculation.
!  type(t_spacetimeVector), intent(INOUT) :: rb
!
!  ! A temporary space-time vector that receives the defect during the calculation.
!  type(t_spacetimeVector), intent(INOUT) :: rd
!!</inputoutput>
!
!!</subroutine>
!
!    ! The nonlinear solver configuration
!    integer :: isubstep,iglobIter
!    logical :: bneumann
!    character(LEN=SYS_STRLEN) :: sstring,slinearSolver
!    integer :: nminIterations,nmaxIterations
!
!    real(DP) :: ddefNorm,dinitDefNorm,depsRel,depsAbs
!
!    type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
!    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dx
!
!    ! A temporary vector in the size of a spatial vector.
!    type(t_vectorBlock) :: rtempVectorX
!
!    ! A second temporary vector in the size of a spatial vector.
!    type(t_vectorBlock) :: rtempVectorB
!
!    ! A third temporary vector for the nonlinear iteration
!    type(t_vectorBlock) :: rtempVector
!
!    ! Create temp vectors for X, B and D.
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVector,.true.)
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVectorX,.true.)
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVectorB,.true.)
!
!!    ! Implement the bondary conditions into all initial solution vectors
!!    DO isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals
!!
!!      ! Current point in time
!!      rproblem%rtimedependence%dtime = &
!!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%rtimeDiscr%dtstep
!!      rproblem%rtimedependence%itimestep = isubstep
!!
!!      ! -----
!!      ! Discretise the boundary conditions at the new point in time --
!!      ! if the boundary conditions are nonconstant in time!
!!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!!        CALL cc_updateDiscreteBC (rproblem)
!!      END IF
!!
!!      ! Implement the boundary conditions into the global solution vector.
!!      CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
!!
!!      ! DEBUG!!!
!!      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
!!
!!      CALL cc_implementBC (rproblem,rvector=rtempVectorX)
!!
!!      CALL sptivec_setTimestepData(rx, isubstep, rtempVectorX)
!!
!!    END DO
!
!    call tbc_implementBCsolution (rproblem,rspaceTimeDiscr,rx)
!
!    ! Generate the RHS for the initial condition.
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rinitialCondRHS,.false.)
!    call sptivec_getTimestepData (rx, 1, rtempVectorX)
!    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
!    call stlin_generateInitCondRHS (rproblem,rspaceTimeDiscr,&
!        rtempVectorX,rinitialCondRHS)
!
!    ddefNorm = 1.0_DP
!
!    ! ---------------------------------------------------------------
!    ! Set up the structure of the global space time matrix.
!    ! As we can only handle linear subproblems here, we don't have
!    ! to initialise the evaluation point of the nonlinearity.
!
!    rspaceTimeMatrix%p_rspaceTimeDiscr => rspaceTimeDiscr
!    rspaceTimeMatrix%cmatrixType = 0
!
!    ! ---------------------------------------------------------------
!    ! Solve the global space-time coupled system.
!    !
!    ! Get the initial defect: d=b-Ax
!    !CALL smva_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
!    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
!    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
!    !  rtempvectorX, rtempvectorB, rtempvector, .FALSE.)
!    call trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rb, .false.)
!
!    ! Implement the initial condition into the RHS.
!    call tbc_implementInitCondRHS (rproblem, rb, rinitialCondRHS, rtempvector)
!    call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)
!
!    ! Now work with rd, our 'defect' vector
!    call sptivec_copyVector (rb,rd)
!
!    ! Assemble the defect.
!    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
!    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!    dinitDefNorm = ddefNorm
!
!    call output_separator (OU_SEP_EQUAL)
!    call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
!    call output_separator (OU_SEP_EQUAL)
!
!    ! Call the routine to generate the global matrix and to solve the system.
!    call cc_applyUmfpackToSupersystem (rproblem, rspaceTimeMatrix, rx, rd, &
!      rtempvectorX, rtempvectorB, rtempVector)
!
!    ! Calculate the final defect
!    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rd, &
!    !  rtempvectorX, rtempvectorB, rtempvector,.FALSE.)
!    call trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rd, .false.)
!
!    ! Implement the initial condition into the RHS.
!    call tbc_implementInitCondRHS (rproblem, rd, rinitialCondRHS, rtempvector)
!    call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)
!
!    ! Assemble the defect
!    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
!    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!    call output_separator (OU_SEP_EQUAL)
!    call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
!    call output_separator (OU_SEP_EQUAL)
!
!    ! Normalise the primal and dual pressure to integral mean value zero
!    ! where no Neumann boundary is present.
!    call tbc_pressureToL20 (rproblem,rx,rtempVectorX)
!
!    call lsysbl_releaseVector (rtempVectorB)
!    call lsysbl_releaseVector (rtempVectorX)
!    call lsysbl_releaseVector (rtempVector)
!    call lsysbl_releaseVector (rinitialCondRHS)
!    call lsysbl_releaseVector (rinitialCondSol)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_precondDefectSupersystem (rproblem, rspaceTimeMatrix, rx, rb, rd, &
!      rtempvectorX,  rtempvectorD, rpreconditioner)
!
!!<description>
!  ! This routine performs preconditioning with the nonlinear super-defect
!  ! vector rd: $d = A^{-1} d$.
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  type(t_problem), intent(INOUT) :: rproblem
!
!  ! space time matrix structure defining the matrix A.
!  type(t_ccoptSpaceTimeMatrix), intent(INout) :: rspaceTimeMatrix
!
!  ! A space-time vector defining the current solution.
!  type(t_spacetimeVector), intent(IN) :: rx
!
!  ! A space-time vector defining the current RHS.
!  type(t_spacetimeVector), intent(IN) :: rb
!
!!</input>
!
!!<inputoutput>
!  ! A spatial preconditioner. This one is applied to each substep in the
!  ! global matrix.
!  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!
!  ! A temporary vector in the size of a spatial vector.
!  type(t_vectorBlock), intent(INOUT) :: rtempVectorX
!
!  ! A third temporary vector for the nonlinear iteration
!  type(t_vectorBlock), intent(INOUT) :: rtempVectorD
!
!  ! A space-time vector that receives the preconditioned defect.
!  type(t_spacetimeVector), intent(INOUT) :: rd
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: isubstep,ilevel,itejac
!    real(DP) :: dtheta,dtstep,ddefNorm,dinitDef,dtime
!    logical :: bsuccess
!    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
!    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
!    type(t_vectorBlock) :: rtempVectorX1,rtempVectorX3
!    type(t_spacetimeVector) :: rcurrentx,rcurrentd
!
!    type(t_sptilsNode), pointer         :: rsolverNode
!    type(t_ccoptSpaceTimeMatrix), dimension(1) :: RspaceTimeMatrices
!    integer :: ierror
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dx,p_Dd
!
!    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscr
!
!    rspaceTimeMatrix%cmatrixType=1
!
!    dtheta = rproblem%rtimedependence%dtimeStepTheta
!    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep
!
!    ! Level of the discretisation
!    ilevel = p_rspaceTimeDiscr%ilevel
!
!    ! Create two additional temp vectors
!    call lsysbl_createVecBlockIndirect (rtempVectorX,rtempVectorX1,.true.)
!    call lsysbl_createVecBlockIndirect (rtempVectorX,rtempVectorX3,.true.)
!
!    ! Create a temp vector for the current iterates of the preconditioner
!    call sptivec_initVectorDiscr (rcurrentx,rx%p_rtimeDiscretisation,rx%p_rblockDiscretisation)
!    call sptivec_initVectorDiscr (rcurrentd,rd%p_rtimeDiscretisation,rd%p_rblockDiscretisation)
!
!    ! Create the initial defect.
!    call sptivec_copyVector (rd,rcurrentd)
!    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rcurrentx, rcurrentd, &
!      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, dinitDef,rproblem%MT_outputLevel .ge. 2)
!    ddefNorm = dinitDef
!    call output_line("Block-JAC: Ite=0, ||res||="//adjustl(sys_sdEP(ddefNorm,20,10)))
!
!
!    !!! <!-- DEBUG
!    !call sptils_initUMFPACK4 (rproblem,rsolverNode)
!    !rsolverNode%p_rsubnodeUMFPACK4%cwriteMatrix = 1
!    !RspaceTimeMatrices(1) = rspaceTimeMatrix
!    !call sptils_setMatrices (rsolverNode,RspaceTimeMatrices)
!    !call sptils_initstructure (rsolverNode,ierror)
!    !call sptils_initdata (rsolverNode,ierror)
!    !CALL sptils_releaseSolver (rsolverNode)
!    !!! -->
!    !read *
!
!    ! ----------------------------------------------------------------------
!    ! We use a block-Jacobi scheme for preconditioning...
!    !
!    ! Loop until convergence.
!    do itejac = 1,100
!
!      ! Stopping criterion: 2 digits.
!      if (ddefNorm .le. dinitDef*1.0E-2_DP) exit
!
!      ! For this purpose, loop through the substeps.
!
!      do isubstep = 1,p_rspaceTimeDiscr%NEQtime
!
!        ! Current time step?
!        dtime = &
!            rproblem%rtimedependence%dtimeInit + (isubstep-1) * dtstep
!
!        call output_line ('Block-Jacobi preconditioning of timestep: '//&
!            trim(sys_siL(isubstep,10))//&
!            ' Time: '//trim(sys_sdL(dtime,10)))
!
!        ! -----
!        ! Discretise the boundary conditions at the new point in time --
!        ! if the boundary conditions are nonconstant in time!
!        call cc_updatePreconditionerBC (rproblem,rpreconditioner,dtime)
!
!        ! DEBUG!!!
!        call lsysbl_getbase_double (rtempVectorX,p_Dx)
!        call lsysbl_getbase_double (rtempVectorD,p_Dd)
!
!        ! Read in the RHS/solution/defect vector of the current timestep.
!        call sptivec_getTimestepData (rx, isubstep, rtempVectorX)
!        call sptivec_getTimestepData (rcurrentd, isubstep, rtempVectorD)
!
!        if (isubstep .gt. 1) then
!          call sptivec_getTimestepData (rx, isubstep-1, rtempVectorX1)
!        end if
!
!        if (isubstep .le. p_rspaceTimeDiscr%NEQtime-1) then
!          call sptivec_getTimestepData (rx, isubstep+1, rtempVectorX3)
!        end if
!
!        ! Set up the matrix weights for the diagonal matrix
!        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
!          isubstep,0,rnonlinearSpatialMatrix)
!
!        ! Perform preconditioning of the defect with the method provided by the
!        ! core equation module.
!        call cc_precondSpaceDefect (rpreconditioner,rnonlinearSpatialMatrix,rtempVectorD,&
!          rtempVectorX1,rtempVectorX,rtempVectorX3,&
!          bsuccess,rproblem%rcollection)
!
!        ! Save back the preconditioned defect.
!        call sptivec_setTimestepData (rcurrentd, isubstep, rtempVectorD)
!
!      end do
!
!      ! Add the correction to the current iterate
!      call sptivec_vectorLinearComb(rcurrentd,rcurrentx,0.7_DP,1.0_DP)
!
!      ! Create the new defect
!      call sptivec_copyVector (rd,rcurrentd)
!      call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rcurrentx, rcurrentd, &
!        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,.false.)
!      call output_line("Block-JAC: Ite="//trim(sys_siL(itejac,10))//&
!          ", ||res||="//adjustl(sys_sdEP(ddefNorm,20,10)))
!
!    end do ! itejac
!
!    ! Return the preconditioned defect
!    call sptivec_copyVector (rcurrentx,rd)
!
!    ! Release temp vectors
!    call lsysbl_releaseVector (rtempVectorX1)
!    call lsysbl_releaseVector (rtempVectorX3)
!    call sptivec_releaseVector (rcurrentx)
!    call sptivec_releaseVector (rcurrentd)
!
!    rspaceTimeMatrix%cmatrixType=0
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_solveSupersystemDefCorr (rproblem, rspaceTimeDiscr, rx, rb, rd,&
!      ctypePreconditioner)
!
!!<description>
!  ! This is a primitive space-time defect correction solver. It applies the
!  ! iteration
!  !    $$ rx = rx + C^{-1} ( rb - A rx ) $$
!  ! to a space time vector rx and a space time RHS vector rb.
!  ! Here, $C^{-1}$ is a spatial preconditioner (linear solver) that is applied
!  ! to each time step during the iteration. The configuration of this defect
!  ! correction loop is specified in the '[TIME-DEFCORR]' section in the DAT files.
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  type(t_problem), intent(INOUT) :: rproblem
!
!  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
!  ! coupled space-time matrix.
!  type(t_ccoptSpaceTimeDiscretisation), intent(IN),target :: rspaceTimeDiscr
!
!  ! Type of preconditioner to use for the space time system.
!  ! =1: Standard linear system.
!  ! =2: Newton iteration
!  integer, intent(IN) :: ctypePreconditioner
!!</input>
!
!!<inputoutput>
!  ! A space-time vector defining the initial solution. Is replaced by a new
!  ! solution vector.
!  type(t_spacetimeVector), intent(INOUT), target :: rx
!
!  ! A temporary space-time vector that receives the RHS during the calculation.
!  type(t_spacetimeVector), intent(INOUT) :: rb
!
!  ! A temporary space-time vector that receives the defect during the calculation.
!  type(t_spacetimeVector), intent(INOUT) :: rd
!!</inputoutput>
!
!!</subroutine>
!
!    ! The nonlinear solver configuration
!    integer :: isubstep,iglobIter
!    logical :: bneumann
!    character(LEN=SYS_STRLEN) :: sstring,slinearSolver
!    integer :: nminIterations,nmaxIterations
!    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol
!
!    real(DP) :: ddefNorm,dinitDefNorm,depsRel,depsAbs,dlastDefNorm,depsDiff
!
!    type(t_ccspatialPreconditioner) :: rpreconditioner
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dx
!
!    ! A temporary vector in the size of a spatial vector.
!    type(t_vectorBlock) :: rtempVectorX
!
!    ! A second temporary vector in the size of a spatial vector.
!    type(t_vectorBlock) :: rtempVectorB
!
!    ! A third temporary vector for the nonlinear iteration
!    type(t_vectorBlock) :: rtempVector
!
!    ! A structure for the global space-time matrix on the highest level.
!    type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix,rspaceTimePreconditioner
!
!    ! Create temp vectors for X, B and D.
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVector,.true.)
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVectorX,.true.)
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rtempVectorB,.true.)
!
!    ! Some preparations for the spatial preconditioner.
!    !
!    ! Get the name of the section containing the linear solver from the DAT file.
!    call parlst_getvalue_string (rproblem%rparamList, 'TIME-DEFCORR', &
!                                 'slinearSolver', sstring, '')
!    read(sstring,*) slinearSolver
!
!    ! Initialise the preconditioner for the preconditioning in every timestep.
!    ! Specify slinearSolver as the name of the section that configures the
!    ! spatial preconditioner. This is (up to now only) the name of a section
!    ! containing configuration data of a linear solver.
!    call cc_initSpacePreconditioner (rproblem,&
!        rproblem%NLMIN,rproblem%NLMAX,rpreconditioner)
!    call cc_configPreconditioner (rproblem,rpreconditioner,slinearSolver,&
!        ctypePreconditioner)
!
!!    ! Implement the bondary conditions into all initial solution vectors
!!    DO isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals
!!
!!      ! Current point in time
!!      rproblem%rtimedependence%dtime = &
!!          rproblem%rtimedependence%dtimeInit + isubstep*rspaceTimeDiscr%rtimeDiscr%dtstep
!!      rproblem%rtimedependence%itimestep = isubstep
!!
!!      ! -----
!!      ! Discretise the boundary conditions at the new point in time --
!!      ! if the boundary conditions are nonconstant in time!
!!      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!!        CALL cc_updateDiscreteBC (rproblem)
!!      END IF
!!
!!      ! Implement the boundary conditions into the global solution vector.
!!      CALL sptivec_getTimestepData(rx, isubstep, rtempVectorX)
!!
!!      ! DEBUG!!!
!!      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
!!
!!      CALL cc_implementBC (rproblem,rvector=rtempVectorX)
!!
!!      CALL sptivec_setTimestepData(rx, isubstep, rtempVectorX)
!!
!!    END DO
!    call tbc_implementBCsolution (rproblem,rspaceTimeDiscr,rx)
!
!    ! Generate the RHS for the initial condition.
!    call lsysbl_createVecBlockByDiscr (rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rinitialCondRHS,.false.)
!    call sptivec_getTimestepData (rx, 1, rtempVectorX)
!    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
!    call stlin_generateInitCondRHS (rproblem,rspaceTimeDiscr,&
!        rtempVectorX,rinitialCondRHS)
!
!    ddefNorm = 1.0_DP
!
!    ! ---------------------------------------------------------------
!    ! Set up the structure of the global space time matrix.
!    ! Set the evaluation point of the matrix to the current solution
!    ! vector.
!
!    rspaceTimeMatrix%p_rspaceTimeDiscr => rspaceTimeDiscr
!    rspaceTimeMatrix%cmatrixType = 0
!    rspaceTimeMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
!    rspaceTimeMatrix%p_rsolution => rx
!
!    ! Set up a structure for the matrix that serves as space-time
!
!    ! preconditioner. This is based on the space time matrix...
!
!    rspaceTimePreconditioner = rspaceTimeMatrix
!    if ((rproblem%rphysicsPrimal%cequation .eq. 0) .and. &
!        (ctypePreconditioner .eq. 1)) then
!      ! ...but may also be the Newton matrix!
!      rspaceTimePreconditioner%cmatrixType = 1
!    end if
!
!    ! ---------------------------------------------------------------
!    ! Solve the global space-time coupled system.
!    !
!    ! Get the initial defect: d=b-Ax
!    !CALL smva_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
!    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
!    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
!    !  rtempvectorX, rtempvectorB, rtempvector, .FALSE.)
!    call trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rb, .false.)
!
!    ! Implement the initial condition into the RHS.
!    call tbc_implementInitCondRHS (rproblem, rb, rinitialCondRHS, rtempvector)
!    call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)
!
!    ! Now work with rd, our 'defect' vector
!    call sptivec_copyVector (rb,rd)
!
!    ! Assemble the defect.
!    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
!    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!    dinitDefNorm = ddefNorm
!    dlastDefNorm = 0.0_DP
!
!    call output_separator (OU_SEP_EQUAL)
!    call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
!    call output_separator (OU_SEP_EQUAL)
!
!    ! Get some solver parameters for the iteration
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-DEFCORR', &
!                              'nminIterations', nminIterations, 0)
!
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-DEFCORR', &
!                              'nmaxIterations', nmaxIterations, 0)
!
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-DEFCORR', &
!                                 'depsRel', depsRel, 1.0E-5_DP)
!
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-DEFCORR', &
!                                 'depsDiff', depsDiff, 1.0E-5_DP)
!
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-DEFCORR', &
!                                 'depsAbs', depsAbs, 1.0E99_DP)
!
!    iglobIter = 0
!
!!    !CALL cc_solveSupersysDirect (rproblem, rspaceTimeDiscr, rx, rd, &
!!    !  rtempvectorX, rtempvectorB, rtempVector)
!!    CALL cc_solveSupersysDirectCN (rproblem, rspaceTimeDiscr, rx, rd, &
!!      rtempvectorX, rtempvectorB, rtempVector)
!
!    do while ((iglobIter .lt. nminIterations) .or. &
!              ((((ddefNorm .gt. depsRel*dinitDefNorm) .or. (ddefNorm .ge. depsAbs)) .and.&
!                (abs(ddefNorm-dlastDefNorm) .ge. depsDiff*dlastDefNorm)) .and. &
!               (iglobIter .lt. nmaxIterations)))
!
!      iglobIter = iglobIter+1
!
!      ! Preconditioning of the defect: d=C^{-1}d
!      call cc_precondDefectSupersystem (rproblem, rspaceTimePreconditioner, &
!          rx, rb, rd, rtempvectorX,  rtempVector, rpreconditioner)
!
!      ! Add the defect: x = x + omega*d
!      call sptivec_vectorLinearComb (rd,rx,1.0_DP,1.0_DP)
!
!      ! Assemble the new defect: d=b-Ax
!      call sptivec_copyVector (rb,rd)
!      !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
!      call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!      call output_separator (OU_SEP_EQUAL)
!      call output_line ('Iteration: '//trim(sys_siL(iglobIter,10))//&
!          '              Defect of supersystem: '//adjustl(sys_sdEP(ddefNorm,20,10)))
!      call output_separator (OU_SEP_EQUAL)
!
!    end do
!
!    ! ---------------------------------------------------------------
!    ! Release the preconditioner of the nonlinear iteration
!    call cc_doneSpacePreconditioner (rpreconditioner)
!
!    !CALL smva_assembleDefectSupersystem (rproblem, rspaceTimeDiscr, rx, rd, &
!    !    rtempvectorX, rtempvectorB, rtempVector, ddefNorm)
!    !CALL cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rd, &
!    !  rtempvectorX, rtempvectorB, rtempvector,.FALSE.)
!    call trhsevl_assembleRHS (rproblem, rspaceTimeDiscr, rd, .false.)
!
!    ! Implement the initial condition into the RHS.
!    call tbc_implementInitCondRHS (rproblem, rd, rinitialCondRHS, rtempvector)
!    call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)
!
!    ! Assemble the defect
!    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, rx, rd, ddefNorm)
!    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!    call output_separator (OU_SEP_EQUAL)
!    call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
!    call output_separator (OU_SEP_EQUAL)
!
!    ! Normalise the primal and dual pressure to integral mean value zero
!    ! where no Neumann boundary is present.
!    call tbc_pressureToL20 (rproblem,rx,rtempVectorX)
!
!    call lsysbl_releaseVector (rtempVectorB)
!    call lsysbl_releaseVector (rtempVectorX)
!    call lsysbl_releaseVector (rtempVector)
!    call lsysbl_releaseVector (rinitialCondRHS)
!    call lsysbl_releaseVector (rinitialCondSol)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_calcOptimalCorrection (rproblem,rspaceTimeDiscr,rx,rb,rd,rtemp,rtemp2,domega)
!
!!<description>
!  ! Performs a line search algorithm to calculate the optimal correction
!  ! factor for the Newton iteration iteration.
!!</description>
!
!!<input>
!
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  type(t_problem), intent(INOUT) :: rproblem
!
!  ! Space-time discretisation
!  type(t_ccoptSpaceTimeDiscretisation), intent(IN), target :: rspaceTimeDiscr
!
!  ! Current solution vector
!  type(t_spacetimeVector), intent(INOUT), target :: rx
!
!  ! Current RHS vector
!  type(t_spacetimeVector), intent(INOUT), target :: rb
!
!  ! Correction vector
!  type(t_spacetimeVector), intent(INOUT), target :: rd
!
!  ! Temporary vector
!  type(t_spacetimeVector), intent(INOUT), target :: rtemp
!
!  ! Temporary vector
!  type(t_spacetimeVector), intent(INOUT), target :: rtemp2
!
!!</input>
!
!!<output>
!
!  ! Optimal correction factor
!  real(DP), intent(out) :: domega
!
!!</output>
!
!!</subroutine>
!
!    real(DP) :: df0, df0p, df, dt, df2
!
!    ! We implement an Armijo step length control for the Newton iteration here.
!    ! Our functional is given by the function calcfunctional below.
!    !
!    ! The functional we want to minimise here is
!    !
!    !  f(t) = F(x+t*dx)dx
!    !
!    ! with F(u)=-Laplace(u)+...-rhs being the residual of the Navier--Stokes equation.
!    !
!    ! At first, get f(0) and f'(0). We note that here, we have
!    !  f(0) = F(x)dx and f0' = DF(x)dx*dx = -F(x)dx
!    ! because Newton computes dx as solution of DF(x)dx = -F(x) !
!    df0 = calcfunctional (rproblem,rspaceTimeDiscr,0.0_DP,rx,rb,rd,rtemp,rtemp2)
!    df0p = -df0
!    !df0p = calcderiv (rproblem,rspaceTimeDiscr,rx,rb,rd,rtemp,rtemp2)
!
!    ! Now search for f(x+td) with t=1,1/2,1/4,1/8,...
!    ! until f <= f0 + t f0'.
!    dt = 1.0_DP
!    df = calcfunctional (rproblem,rspaceTimeDiscr,dt,rx,rb,rd,rtemp,rtemp2)
!    df2 = calcfunctional (rproblem,rspaceTimeDiscr,0.001_DP,rx,rb,rd,rtemp,rtemp2)
!
!    if (df0p .lt. 0.0_DP) then
!      do while (df .gt. df0 + 0.1_DP*dt*df0p)
!        dt = 0.5_DP*dt
!        df = calcfunctional (rproblem,rspaceTimeDiscr,dt,rx,rb,rd,rtemp,rtemp2)
!      end do
!
!    else
!      print *,'************** Derivative positive *******************'
!    end if
!
!    domega = dt
!
!  contains
!
!    subroutine calcresidual (rproblem,rspaceTimeDiscr,dalpha,rx,rb,rd,rF,rtemp)
!
!    ! With u=rx+alpha rd, this calculates the residuaö
!    !    res(alpha) := -Laplace(u) + ... - RHS
!    ! of the Navier--Stokes equation.
!
!    ! A problem structure that provides information about matrices on all
!    ! levels as well as temporary vectors.
!    type(t_problem), intent(INOUT) :: rproblem
!
!    ! Space-time discretisation
!    type(t_ccoptSpaceTimeDiscretisation), intent(IN), target :: rspaceTimeDiscr
!
!    ! Relaxation parameter
!    real(DP), intent(in) :: dalpha
!
!    ! Current solution vector
!    type(t_spacetimeVector), intent(INOUT), target :: rx
!
!    ! Current RHS vector
!    type(t_spacetimeVector), intent(INOUT), target :: rb
!
!    ! Correction vector
!    type(t_spacetimeVector), intent(INOUT), target :: rd
!
!    ! Temporary vector
!    type(t_spacetimeVector), intent(INOUT), target :: rtemp
!
!    ! Vector obtaining the residual vector
!    type(t_spacetimeVector), intent(INOUT), target :: rF
!
!      ! Local variables
!      type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
!      real(DP) :: ddef
!
!      ! Create a space time matrix for the maximum level which serves as basis for
!      ! setting up the defect. This is the actual system matrix, thus it doesn't
!      ! contain any Newton parts!
!      rspaceTimeMatrix%p_rspaceTimeDiscr => rspaceTimeDiscr
!      rspaceTimeMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
!
!      ! Create the point where to evaluate
!      !if (dalpha .ne. 0.0_DP) then
!        rspaceTimeMatrix%p_rsolution => rtemp
!        call sptivec_copyVector(rx,rtemp)
!        call sptivec_vectorLinearComb (rd,rtemp,dalpha,1.0_DP)
!      !else
!      !  rspaceTimeMatrix%p_rsolution => rx
!      !end if
!
!      ! Evaluate to rF.
!      call sptivec_copyVector(rb,rF)
!      call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, &
!        rspaceTimeMatrix%p_rsolution, rF, &
!        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddef)
!      print *,'Def=',ddef
!
!    end subroutine
!
!    ! -------------------------------------------------------------------------
!
!    subroutine applyderiv (rproblem,rspaceTimeDiscr,rx,ry,rF)
!
!    ! Calculates rd = A'(rx)ry
!
!    ! A problem structure that provides information about matrices on all
!    ! levels as well as temporary vectors.
!    type(t_problem), intent(INOUT) :: rproblem
!
!    ! Space-time discretisation
!    type(t_ccoptSpaceTimeDiscretisation), intent(IN), target :: rspaceTimeDiscr
!
!    ! Current solution vector where to evaluate the matrix
!    type(t_spacetimeVector), intent(INOUT), target :: rx
!
!    ! Solution vector
!    type(t_spacetimeVector), intent(INOUT), target :: ry
!
!    ! Vector obtaining A'(x)y
!    type(t_spacetimeVector), intent(INOUT), target :: rF
!
!      ! Local variables
!      type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
!      real(DP) :: ddef
!
!      ! Create a space time matrix for the maximum level which serves as basis for
!      ! setting up the defect. This is the actual system matrix, thus it doesn't
!      ! contain any Newton parts!
!      rspaceTimeMatrix%p_rspaceTimeDiscr => rspaceTimeDiscr
!      rspaceTimeMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
!      rspaceTimeMatrix%cmatrixType = 1
!
!      ! Create the point where to evaluate
!      rspaceTimeMatrix%p_rsolution => rx
!
!      ! Evaluate to rF.
!      call sptivec_clearVector(rF)
!      call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, &
!        ry, rF, 1.0_DP, 0.0_DP, 0,ddef)
!      print *,'Def=',ddef
!
!    end subroutine
!
!    ! -------------------------------------------------------------------------
!
!    real(DP) function calcfunctional (rproblem,rspaceTimeDiscr,dalpha,rx,rb,rd,rtemp,rtemp2)
!
!    ! Calculates the functional
!    !    f(alpha) := F(rx+alpha rd)rd
!    ! with F(x) being the residual of the nonlinear (Navier-)Stokes equation.
!
!    ! A problem structure that provides information about matrices on all
!    ! levels as well as temporary vectors.
!    type(t_problem), intent(INOUT) :: rproblem
!
!    ! Space-time discretisation
!    type(t_ccoptSpaceTimeDiscretisation), intent(IN), target :: rspaceTimeDiscr
!
!    ! Relaxation parameter
!    real(DP), intent(in) :: dalpha
!
!    ! Current solution vector
!    type(t_spacetimeVector), intent(INOUT), target :: rx
!
!    ! Current RHS vector
!    type(t_spacetimeVector), intent(INOUT), target :: rb
!
!    ! Correction vector
!    type(t_spacetimeVector), intent(INOUT), target :: rd
!
!    ! Temporary vector
!    type(t_spacetimeVector), intent(INOUT), target :: rtemp
!
!    ! Vector obtaining the residual vector
!    type(t_spacetimeVector), intent(INOUT), target :: rtemp2
!
!      ! Calculate F(x+omega*d)d. Note that calcresidual calculates the
!      ! residual RHS-NavSt(...), where the functional needs the defect
!      ! NavSt(...)-RHS -- so we would have to spend an additional '-'-sign.
!      ! On the other hand, there is a '-' sign included in rd which cancels
!      ! against the '-' sign we would have to add -- so nothing is to do here.
!      call calcresidual (rproblem,rspaceTimeDiscr,dalpha,rx,rb,rd,rtemp,rtemp2)
!      calcfunctional = -sptivec_scalarProduct (rtemp, rd)
!
!    end function
!
!
!    ! -------------------------------------------------------------------------
!
!    real(DP) function calcderiv (rproblem,rspaceTimeDiscr,rx,rb,rd,rtemp,rtemp2)
!
!    ! Calculates the functional
!    !    f(alpha) := F(rx+alpha rd)rd
!    ! with F(x) being the residual of the nonlinear (Navier-)Stokes equation.
!
!    ! A problem structure that provides information about matrices on all
!    ! levels as well as temporary vectors.
!    type(t_problem), intent(INOUT) :: rproblem
!
!    ! Space-time discretisation
!    type(t_ccoptSpaceTimeDiscretisation), intent(IN), target :: rspaceTimeDiscr
!
!    ! Current solution vector
!    type(t_spacetimeVector), intent(INOUT), target :: rx
!
!    ! Current RHS vector
!    type(t_spacetimeVector), intent(INOUT), target :: rb
!
!    ! Correction vector
!    type(t_spacetimeVector), intent(INOUT), target :: rd
!
!    ! Temporary vector
!    type(t_spacetimeVector), intent(INOUT), target :: rtemp
!
!    ! Vector obtaining the residual vector
!    type(t_spacetimeVector), intent(INOUT), target :: rtemp2
!
!      call applyderiv (rproblem,rspaceTimeDiscr,rx,rd,rtemp)
!      calcderiv = -sptivec_scalarProduct (rtemp, rd)
!
!    end function
!
!  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_solveSupersystemMultigrid (rproblem, RspaceTimeDiscr, rx, rb, rd,&
!      ctypePreconditioner)
!
!!<description>
!  ! This subroutine solves the nonstationary space time coupled (Navier-)Stokes
!  ! optimal control problem. For this purpose, a nonlinear defect correction
!  ! loop is executed. Every defect is preconditioned by a space-time coupled
!  ! preconditioner like Block-Jacobi, block SOR or whatever.
!  !
!  ! The caller must provide problem information in rproblem and a set of
!  ! matrix configurations for all time levels which are allowed for the solver
!  ! to use. If multiple time-levels are provided, space-time coupled multigrid
!  ! is used for preconditioning. matrix configurations must be provided in
!  ! RspaceTimeDiscr. The maximum level in this array defines the level where
!  ! the system is solved. This level must correspond to the vectors rx, rb and
!  ! rd.
!  !
!  ! The caller can provide an initial solution in rx. However, rb is
!  ! overwritten by a newly generated space-time coupled RHS vector.
!  ! rd is used as temporary vector.
!!</description>
!
!!<input>
!  ! A problem structure that provides information about matrices on all
!  ! levels as well as temporary vectors.
!  type(t_problem), intent(INOUT) :: rproblem
!
!  ! An array of t_ccoptSpaceTimeDiscretisation structure defining all the
!  ! levels of the coupled space-time the discretisation.
!  ! The solution/rhs rx/rb must correspond to the maximum level in this
!  ! array.
!  type(t_ccoptSpaceTimeDiscretisation), &
!      dimension(:), intent(IN), target :: RspaceTimeDiscr
!
!  ! Type of preconditioner to use for the space time system.
!  ! =1: Standard linear system.
!  ! =2: Newton iteration
!  integer, intent(IN) :: ctypePreconditioner
!!</input>
!
!!<inputoutput>
!  ! A space-time vector defining the initial solution. Is replaced by a new
!  ! solution vector.
!  type(t_spacetimeVector), intent(INOUT), target :: rx
!
!  ! A temporary space-time vector that receives the RHS during the calculation.
!  type(t_spacetimeVector), intent(INOUT), target :: rb
!
!  ! A temporary space-time vector that receives the defect during the calculation.
!  type(t_spacetimeVector), intent(INOUT) :: rd
!!</inputoutput>
!
!!</subroutine>
!
!    ! The nonlinear solver configuration
!    integer :: isubstep,iglobIter,ierror,ilev,ispacelev,nsmSteps,cspaceTimeSmoother
!    integer :: ilowerSpaceLevel,itemp,ctypeCoarseGridSolver,i,nlmax
!    integer :: iorderTimeProlRest
!    logical :: bneumann
!    integer :: icalcError
!    character(LEN=SYS_STRLEN) :: ssolutionExpressionY1,ssolutionExpressionY2
!    character(LEN=SYS_STRLEN) :: ssolutionExpressionP
!    character(LEN=SYS_STRLEN) :: ssolutionExpressionLAMBDA1,ssolutionExpressionLAMBDA2
!    character(LEN=SYS_STRLEN) :: ssolutionExpressionXI
!    real(DP), dimension(4) :: Derror
!    integer(I32) :: nminIterations,nmaxIterations
!    real(DP) :: depsRel,depsAbs,domega,domegaPrecond,depsDiff
!    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
!    type(t_ccspatialPreconditioner), dimension(size(RspaceTimeDiscr)) :: RspatialPrecond
!    type(t_ccspatialPreconditioner), dimension(size(RspaceTimeDiscr)) :: RspatialPrecondPrimal
!    type(t_ccspatialPreconditioner), dimension(size(RspaceTimeDiscr)) :: RspatialPrecondDual
!    type(t_sptiProjection), dimension(size(RspaceTimeDiscr)) :: RinterlevelProjection
!    type(t_ccoptSpaceTimeMatrix), dimension(size(RspaceTimeDiscr)) :: RspaceTimePrecondMatrix
!    type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
!    type(t_vectorBlock) :: rtempVecCoarse,rtempVecFine
!    type(t_vectorBlock) :: rinitialCondRHS,rinitialCondSol
!    integer :: ifbSORPartialUpdate
!
!    ! A solver node that identifies our solver.
!    type(t_sptilsNode), pointer :: p_rprecond,p_rsolverNode
!    type(t_sptilsNode), pointer :: p_rmgSolver,p_rsmoother,p_rcgrSolver
!    type(t_sptilsNode), pointer :: p_rpresmoother,p_rpostsmoother
!
!    ! TYPE(t_spaceTimeVector) :: rtemp
!
!    real(DP) :: ddefNorm,dinitDefNorm,dlastDefNorm,dtempdef,depsrelLinSol
!    real(DP) :: dinexactNewtonExponent,dinexactNewtonEpsRel
!
!    character(LEN=SYS_STRLEN) :: slinearSolver,sstring
!
!    integer :: inpos, innextrec
!    character :: ctest
!
!    type(t_timer) :: rtimerMGStep,rtimerNonlinear,rtimerPreconditioner
!    integer :: ilinearIterations
!
!    real(DP) :: domega2
!    type(t_spacetimeVector) :: rtempVec,rtempVec2
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dx
!
!    ! A temporary vector in the size of a spatial vector.
!    type(t_vectorBlock) :: rtempVectorX
!
!    ! A second temporary vector in the size of a spatial vector.
!    type(t_vectorBlock) :: rtempVectorB
!
!    ! A third temporary vector for the nonlinear iteration
!    type(t_vectorBlock) :: rtempVector
!
!    ! STATISTICS: Total time needed for smoothing operations
!    type(t_timer) :: rtimeSmoothing
!
!    ! STATISTICS: Total time needed for the coarse grid solver
!    type(t_timer) :: rtimeCoarseGridSolver
!
!    ! STATISTICS: Time needed for linear algebra stuff (matrix-vector,
!    ! vector-copy, prolongation/restriction,...)
!    type(t_timer) :: rtimeLinearAlgebra
!
!    ! STATISTICS: Time needed for prolongation/restriction
!    type(t_timer) :: rtimeProlRest
!
!    ! STATISTICS: Time for solving problems in space.
!    type(t_timer) :: rtimeSpacePrecond
!
!    ! STATISTICS: Time for initialisation / factorisation of the space time system;
!    ! global and in one step
!    type(t_timer) :: rtimeFactorisation,rtimeFactorisationStep
!
!    ! Get a poiter to the discretisation structure on the maximum
!    ! space/time level. That's the level where we solve here.
!    p_rspaceTimeDiscr => RspaceTimeDiscr(size(RspaceTimeDiscr))
!
!    ! Create temp vectors for X, B and D.
!    call lsysbl_createVecBlockByDiscr (&
!        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVector,.true.)
!    call lsysbl_createVecBlockByDiscr (&
!        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVectorX,.true.)
!    call lsysbl_createVecBlockByDiscr (&
!        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,rtempVectorB,.true.)
!
!    call sptivec_initVectorDiscr (rtempVec,rx%p_rtimeDiscretisation,rx%p_rblockDiscretisation)
!    call sptivec_initVectorDiscr (rtempVec2,rx%p_rtimeDiscretisation,rx%p_rblockDiscretisation)
!
!    ! Implement the bondary conditions into all initial solution vectors
!    call tbc_implementBCsolution (rproblem,p_rspaceTimeDiscr,rx,rtempvectorX)
!
!    call output_line ('NLST-Solver: Generating initial RHS...')
!
!    ! Generate the RHS for the initial condition.
!    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
!        rinitialCondRHS,.false.)
!    call sptivec_getTimestepData (rx, 1, rtempVectorX)
!    call lsysbl_copyVector (rtempVectorX,rinitialCondSol)
!    call stlin_generateInitCondRHS (rproblem,p_rspaceTimeDiscr,&
!        rtempVectorX,rinitialCondRHS)
!
!    ! We set up a space-time preconditioner, e.g. in the following configuration:
!    ! Main Preconditioner: Multigrid
!    !     -> Presmoother:  Block Jacobi/GS
!    !     -> Postsmoother: Block Jacobi/GS
!    !     -> CGr-Solver:   Defect correction (or UMFPACK)
!    !        -> Preconditioner: Block Jacobi
!    !
!    ! So we start creating the main solver: Multigrid.
!    !
!    ! Create as many time levels as specified by the length of RspatialPrecond.
!    call sptils_initMultigrid (rproblem,1,size(RspatialPrecond),p_rmgSolver)
!
!    ! Type of smoother to use?
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
!        'cspaceTimeSmoother', cspaceTimeSmoother, 0)
!
!    ! Type of coarse grid solver?
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                              'ctypeCoarseGridSolver', ctypeCoarseGridSolver, 0)
!
!    ! Type of prolongation/restriction in time
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
!        'iorderTimeProlRest', iorderTimeProlRest, -1)
!
!    if (iorderTimeProlRest .eq. -1) then
!      ! Automatic mode. Check the time stepping theta to determine
!      ! if we have 1st or 2nd order in time.
!      iorderTimeProlRest = 1
!      if (rproblem%rtimedependence%dtimeStepTheta .eq. 0.5_DP) then
!        ! 2nd order Crank Nicolson
!        iorderTimeProlRest = 2
!      end if
!    end if
!
!    ! Loop over the time levels.
!    do ilev=1,size(RspatialPrecond)
!
!      call output_line ('NLST-Solver: Initialising MG solver on level '//sys_siL(ilev,10))
!
!      ! Get the refinement level in space that belongs to this space-time level.
!      ispacelev = RspaceTimeDiscr(ilev)%ilevel
!
!      if (ilev .eq. 1) then
!        i = ctypeCoarseGridSolver
!
!        ! Get the name of the section containing the linear solver from the DAT file --
!        ! in case we use a linear solver as spatial preconditioner.
!        call parlst_getvalue_string (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                    'slinearSolver', sstring, '')
!        read(sstring,*) slinearSolver
!
!      else
!        i = cspaceTimeSmoother
!
!        ! Get the name of the section containing the linear solver from the DAT file --
!        ! in case we use a linear solver as spatial preconditioner.
!        call parlst_getvalue_string (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                    'slinearSolver', sstring, '')
!        read(sstring,*) slinearSolver
!
!      end if
!
!      select case (i)
!      case (0:)
!        ! Initialise the spatial preconditioner for Block Jacobi
!        ! (note: this is slightly expensive in terms of memory!
!        ! Probably we could use the same preconditioner for all levels,
!        ! but this has still to be implemeted and is a little bit harder!)
!        call cc_initSpacePreconditioner (rproblem,rproblem%nlmin,&
!            ispacelev,RspatialPrecond(ilev))
!        ! Specify slinearSolver as the name of the section that configures the
!        ! spatial preconditioner. This is (up to now only) the name of a section
!        ! containing configuration data of a linear solver.
!        call cc_configPreconditioner (rproblem,RspatialPrecond(ilev),&
!            slinearSolver,ctypePreconditioner)
!
!      case DEFAULT
!        print *,'Unknown preconditioner/smoother: ',i
!        call sys_halt()
!      end select
!
!      ! Generate an interlevel projection structure for that level.
!      ! Note that space restriction/prolongation must be switched off if
!      ! we are on the spatial coarse mesh!
!      if (ilev .gt. 1) then
!        ilowerSpaceLevel = RspaceTimeDiscr(ilev-1)%ilevel
!      else
!        ilowerSpaceLevel = ispacelev
!      end if
!      call sptipr_initProjection (rinterlevelProjection(ilev),&
!          RspaceTimeDiscr(ilev)%p_rlevelInfo%rdiscretisation,&
!          ilowerSpaceLevel .ne. ispacelev,iorderTimeProlRest)
!
!      if (ilev .eq. 1) then
!        ! ..., on the minimum level, create a coarse grid solver, ...
!
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                    'domega', domega, 1.0_DP)
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                    'ifbSORPartialUpdate', ifbSORPartialUpdate,0)
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!                                    'domega', domegaPrecond, 1.0_DP)
!
!        nullify(p_rprecond)
!        select case (ctypeCoarseGridSolver)
!        case (0)
!          ! Block Jacobi preconditioner
!          call sptils_initBlockJacobi (rproblem,p_rprecond,domega,RspatialPrecond(ilev))
!
!          ! Defect correction solver
!          call sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
!
!        case (1)
!          ! Block SOR preconditioner
!          call sptils_initBlockFBSOR (rproblem,p_rprecond,domega,&
!              domegaPrecond,RspatialPrecond(ilev),ifbSORPartialUpdate .ne. 0)
!
!          ! Defect correction solver
!          call sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
!
!        case (2)
!          ! Forward backward Gauss Seidel
!          call sptils_initBlockFBGS (rproblem,p_rprecond,&
!            domega,domegaPrecond,RspatialPrecond(ilev))
!
!          ! Defect correction solver
!          call sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
!
!        case (3)
!          ! CG with Block Jacobi as preconditioner
!          call sptils_initBlockJacobi (rproblem,p_rprecond,domegaPrecond,RspatialPrecond(ilev))
!
!          ! CG solver
!          call sptils_initCG (rproblem,p_rcgrSolver,p_rprecond)
!
!        case (4)
!          ! UMFACK Gauss elimination
!          call sptils_initUMFPACK4 (rproblem,p_rcgrSolver)
!          p_rcgrSolver%domega = domega
!
!          ! DEBUG!!!
!          !p_rcgrSolver%p_rsubnodeUMFPACK4%cwriteMatrix = 1
!          !p_rcgrSolver%p_rsubnodeUMFPACK4%ixfirst = 1
!          !p_rcgrSolver%p_rsubnodeUMFPACK4%iyfirst = 1
!
!        case (5)
!          ! Forward backward SOR as preconditioner
!          call sptils_initBlockFBSOR (rproblem,p_rprecond,&
!            domegaPrecond,domegaPrecond,RspatialPrecond(ilev),ifbSORPartialUpdate .ne. 0)
!          !CALL sptils_initBlockJacobi (rproblem,p_rprecond,domegaPrecond,RspatialPrecond(ilev))
!
!          ! BiCGStab solver
!          call sptils_initBiCGStab (rproblem,p_rcgrSolver,p_rprecond)
!
!        case (6)
!          ! UMFPACK
!          call sptils_initUMFPACK4 (rproblem,p_rprecond)
!
!          ! Defect correction solver
!          call sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
!          p_rcgrSolver%domega = domega
!
!        case (7)
!          ! Block SOR solver
!          call sptils_initBlockFBSOR (rproblem,p_rcgrSolver,domega,&
!              domegaPrecond,RspatialPrecond(ilev),ifbSORPartialUpdate .ne. 0)
!
!        case (8)
!          ! Forward backward solver as preconditioner
!          call sptils_initFBsim (rproblem,p_rprecond,1.0_DP)
!
!          ! Defect correction solver
!          call sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
!        case default
!          print *,'Unknown solver: ',ctypeCoarseGridSolver
!          stop
!        end select
!
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!            'nminIterations', p_rcgrSolver%nminIterations, 1)
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!            'nmaxIterations', p_rcgrSolver%nmaxIterations, 100)
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                    'depsRel', p_rcgrSolver%depsRel, 1E-5_DP)
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                    'depsAbs', p_rcgrSolver%depsAbs, 1E-5_DP)
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                    'depsDiff', p_rcgrSolver%depsDiff, 0.0_DP)
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                    'ddivRel', p_rcgrSolver%ddivRel, 1.0_DP)
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!                                'istoppingCriterion', p_rcgrSolver%istoppingCriterion, 0)
!
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEGRIDSOLVER', &
!            'ioutputLevel', p_rcgrSolver%ioutputLevel, 100)
!
!        ! If there's a subsolver, configure it.
!        if (associated(p_rprecond)) then
!          call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!              'nminIterations', p_rprecond%nminIterations, 1)
!          call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!              'nmaxIterations', p_rprecond%nmaxIterations, 100)
!          call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!                                      'depsRel', p_rprecond%depsRel, 1E-5_DP)
!          call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!                                      'depsAbs', p_rprecond%depsAbs, 1E-5_DP)
!          call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!                                      'depsDiff', p_rprecond%depsDiff, 0.0_DP)
!          call parlst_getvalue_double (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!                                      'ddivRel', p_rprecond%ddivRel, 1.0_DP)
!          call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!                                  'istoppingCriterion', p_rprecond%istoppingCriterion, 0)
!
!          call parlst_getvalue_int (rproblem%rparamList, 'TIME-COARSEPRECOND', &
!              'ioutputLevel', p_rprecond%ioutputLevel, 100)
!        end if
!
!        ! ...; finally initialise the level with that
!        call sptils_setMultigridLevel (p_rmgSolver,ilev,&
!                      rinterlevelProjection(ilev),&
!                      null(),null(),p_rcgrSolver)
!
!      else
!        ! ... on higher levels, create a smoother, ...
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                    'domega', domega, 1.0_DP)
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                    'domegaPrecond', domegaPrecond, 1.0_DP)
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                    'ifbSORPartialUpdate', ifbSORPartialUpdate, 0)
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
!            'nsmSteps', nsmSteps, 1)
!
!        ! DEBUG!!!
!!        IF (ilev .LE. SIZE(RspatialPrecond)-1) THEN
!!          cspaceTimeSmoother = 7
!!          ! nsmSteps = nsmSteps/2
!!        END IF
!        ! DEBUG!!!
!
!        select case (cspaceTimeSmoother)
!        case (0)
!          ! Block Jacobi
!          call sptils_initBlockJacobi (rproblem,p_rsmoother,&
!            domega,RspatialPrecond(ilev))
!        case (1)
!          ! Block SOR
!          call sptils_initBlockFBSOR (rproblem,p_rsmoother,&
!            domega,domegaPrecond,RspatialPrecond(ilev),ifbSORPartialUpdate .ne. 0)
!        case (2)
!          ! Block Forward-Backward Gauss-Seidel
!          call sptils_initBlockFBGS (rproblem,p_rsmoother,&
!            domega,domegaPrecond,RspatialPrecond(ilev))
!        case (3)
!          ! CG with Block Jacobi as preconditioner
!          call sptils_initBlockJacobi (rproblem,p_rprecond,&
!            domegaPrecond,RspatialPrecond(ilev))
!          call sptils_initCG (rproblem,p_rsmoother,p_rprecond)
!        case (4)
!          ! CG with Block Jacobi as preconditioner
!          call sptils_initBlockFBGS (rproblem,p_rprecond,domegaPrecond,&
!            domegaPrecond,RspatialPrecond(ilev))
!          call sptils_initCG (rproblem,p_rsmoother,p_rprecond)
!        case (6)
!          ! DefCorr with UMFPACK as preconditioner
!          call sptils_initUMFPACK4 (rproblem,p_rprecond)
!          call sptils_initDefCorr (rproblem,p_rsmoother,p_rprecond)
!        case (7)
!          ! BiCGStab with FBGS as preconditioner
!          call sptils_initBlockFBGS (rproblem,p_rprecond,&
!            domega,domegaPrecond,RspatialPrecond(ilev))
!            p_rprecond%nminIterations = 1
!            p_rprecond%nmaxIterations = 1
!          call sptils_initBiCGStab (rproblem,p_rsmoother,p_rprecond)
!        case DEFAULT
!          print *,'Unknown smoother: ',cspaceTimeSmoother
!          stop
!        end select
!
!        call sptils_convertToSmoother (p_rsmoother,nsmSteps,domega)
!
!        ! Switch off smoothing is set that way in the DAT file
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                  'ioutputLevel', p_rsmoother%ioutputLevel, 1)
!
!        p_rpresmoother => p_rsmoother
!        p_rpostsmoother => p_rsmoother
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                  'ipresmoothing', itemp, 1)
!        if (itemp .eq. 0) nullify(p_rpresmoother)
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                  'ipostsmoothing', itemp, 1)
!        if (itemp .eq. 0) nullify(p_rpostsmoother)
!
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                    'depsRel', p_rsmoother%depsRel, 0.0_DP)
!
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                    'depsAbs', p_rsmoother%depsAbs, 0.0_DP)
!
!        call parlst_getvalue_double (rproblem%rparamList, 'TIME-SMOOTHER', &
!                                    'depsDiff', p_rsmoother%depsDiff, 0.0_DP)
!
!        call parlst_getvalue_int (rproblem%rparamList, 'TIME-SMOOTHER', &
!            'istoppingCriterion', p_rsmoother%istoppingCriterion, 1)
!
!        ! DEBUG!!!
!!        IF (ilev .EQ. SIZE(RspatialPrecond)-1) THEN
!!          p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels(&
!!            SIZE(p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels)-1)%depsRelCycle = 1E-10_DP
!!          p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels(&
!!            SIZE(p_rmgSolver%p_rsubnodeMultigrid%p_Rlevels)-1)%depsAbsCycle = 1E-14_DP
!!        END IF
!        ! DEBUG!!!
!
!        if ((.not. associated(p_rpresmoother)) .and. &
!            (.not. associated(p_rpostsmoother))) then
!          call sptils_releaseSolver(p_rsmoother)
!        end if
!
!        ! ...; finally initialise the level with that
!        call sptils_setMultigridLevel (p_rmgSolver,ilev,&
!                      rinterlevelProjection(ilev),&
!                      p_rpresmoother,p_rpostsmoother,null())
!
!      end if
!
!    end do
!
!    ! Our main solver is MG now.
!    p_rsolverNode => p_rmgSolver
!
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
!                             'nmaxIterations', p_rsolverNode%nmaxIterations, 1)
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
!                             'nminIterations', p_rsolverNode%nminIterations, 10)
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
!                             'ioutputLevel', p_rsolverNode%ioutputLevel, 0)
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
!                             'icycle', p_rmgSolver%p_rsubnodeMultigrid%icycle, 0)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
!                                'depsRel', p_rmgSolver%depsRel, 1E-5_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
!                                'depsAbs', p_rmgSolver%depsAbs, 1E-5_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
!                                'depsDiff', p_rmgSolver%depsDiff, 0.0_DP)
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-MULTIGRID', &
!                             'istoppingCriterion', p_rmgSolver%istoppingCriterion, 0)
!
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
!                                'dalphaMin', p_rmgSolver%p_rsubnodeMultigrid%dalphaMin,&
!                                1.0_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-MULTIGRID', &
!                                'dalphaMax', p_rmgSolver%p_rsubnodeMultigrid%dalphaMax,&
!                                1.0_DP)
!
!    ! Save the relative stopping criterion of the linear solver; we
!    ! probably need it for the adaptive Newton.
!    depsrelLinSol = p_rmgSolver%depsRel
!
!    ! For the optimal coarse grid correction, we multiply the residuum of the
!    ! 3rd and 6th equation by -1; gives better results (matrix symmetry?!?).
!    ! Note that there is currently no INIT/DONE-Routine for the weights, so we
!    ! do that manually...
!    allocate(p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(6))
!    p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(:) = 1.0_DP
!    p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(3) = -1.0_DP
!    p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights(6) = -1.0_DP
!
!    ! Initialise the basic parameters of the system matrices on all levels.
!    do ilev=1,size(RspatialPrecond)
!
!      ! Pointer to the corresponding space time discretisation structure
!      RspaceTimePrecondMatrix(ilev)%p_rspaceTimeDiscr => RspaceTimeDiscr(ilev)
!
!      ! Configure the matrix type; standard or Newton matrix
!      select case (ctypePreconditioner)
!      case (CCPREC_LINEARSOLVER)
!        ! Standard system matrix
!        RspaceTimePrecondMatrix(ilev)%cmatrixType = 0
!
!      case (CCPREC_NEWTON,CCPREC_INEXACTNEWTON)
!        ! Newton matrix
!        RspaceTimePrecondMatrix(ilev)%cmatrixType = 1
!        RspaceTimePrecondMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
!      end select
!
!    end do
!
!    ! Allocate space-time vectors on all lower levels that hold the solution vectors
!    ! for the evaluation of the nonlinearity (if we have a nonlinearity).
!    do ilev=1,size(RspatialPrecond)-1
!      allocate(RspaceTimePrecondMatrix(ilev)%p_rsolution)
!      call sptivec_initVector (RspaceTimePrecondMatrix(ilev)%p_rsolution,&
!        RspaceTimeDiscr(ilev)%NEQtime,&
!        RspaceTimeDiscr(ilev)%p_rlevelInfo%rdiscretisation)
!    end do
!
!    ! On the maximum level attach the solution vector.
!    nlmax = ubound(RspaceTimeDiscr,1)
!    RspaceTimePrecondMatrix(nlmax)%p_rsolution => rx
!
!    ! Create a space time matrix for the maximum level which serves as basis for
!    ! setting up the defect. This is the actual system matrix, thus it doesn't
!    ! contain any Newton parts!
!    rspaceTimeMatrix%p_rspaceTimeDiscr => RspaceTimeDiscr(ilev)
!    rspaceTimeMatrix%p_rsolution => rx
!
!    ! That matrix also has to apply projection operators when being applied
!    ! to a vector -- in case control constraints are active.
!    rspaceTimeMatrix%ccontrolConstraints = rproblem%roptcontrol%ccontrolConstraints
!
!    call output_line ('NLST-Solver: Preparing space-time preconditioner...')
!
!    ! Attach matrix information to the linear solver
!    call sptils_setMatrices (p_rsolverNode,RspaceTimePrecondMatrix)
!
!    ! Initialise the space-time preconditioner
!    call stat_clearTimer (rtimeFactorisation)
!    call stat_startTimer (rtimeFactorisation)
!    call sptils_initStructure (p_rsolverNode,ierror)
!    call stat_stopTimer (rtimeFactorisation)
!
!    ddefNorm = 1.0_DP
!
!    ! ---------------------------------------------------------------
!    ! Solve the global space-time coupled system.
!    !
!    ! Get the initial defect: d=b-Ax
!    !CALL cc_assembleSpaceTimeRHS (rproblem, &
!    !  RspaceTimePrecondMatrix(nlmax)%p_rspaceTimeDiscr, rb, &
!    !  rtempvectorX, rtempvectorB, rtempvector, .FALSE.)
!    call output_line ('NLST-Solver: Assembling the RHS vector...')
!    call trhsevl_assembleRHS (rproblem, &
!      RspaceTimePrecondMatrix(nlmax)%p_rspaceTimeDiscr, rb, .false.)
!
!    ! Implement the initial condition into the RHS.
!    call tbc_implementInitCondRHS (rproblem, rb, rinitialCondRHS, rtempvector)
!    call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)
!
!    ! DEBUG!!!
!    !CALL sptivec_saveToFileSequence (rb,&
!    !    '(''./debugdata/initrhs.txt.'',I5.5)',.TRUE.)
!
!    ! Now work with rd, our 'defect' vector
!    call sptivec_copyVector (rb,rd)
!
!    call output_line ('NLST-Solver: Assembling the initial defect...')
!
!    ! Assemble the defect.
!    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, &
!    !    rx, rd, ddefNorm)
!    call sptivec_copyVector (rb,rd)
!    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,&
!      ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!    ! DEBUG!!!
!    !CALL sptivec_saveToFileSequence (rb,&
!    !    '(''./debugdata/initdef.txt.'',I5.5)',.TRUE.)
!
!    dinitDefNorm = ddefNorm
!    dlastDefNorm = 0.0_DP
!    if (dinitDefNorm .eq. 0.0_DP) then
!      ! Trick to avoid div/0.
!      dinitDefNorm = 1.0_DP
!    end if
!    call output_separator (OU_SEP_EQUAL)
!    call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
!    ! Value of the functional
!    call optcana_nonstatFunctional (rproblem,&
!        RspaceTimePrecondMatrix(size(RspatialPrecond))%p_rsolution,&
!        rtempVector,RspaceTimeDiscr(size(RspatialPrecond))%dalphaC,&
!        RspaceTimeDiscr(size(RspatialPrecond))%dgammaC,&
!        Derror)
!    call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(1),10)))
!    call output_line ('||u||         = '//trim(sys_sdEL(Derror(2),10)))
!    call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
!    call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(4),10)))
!
!    ! If error analysis has to be performed, we can calculate
!    ! the real error.
!    call parlst_getvalue_int (rproblem%rparamList,'TIME-POSTPROCESSING',&
!                              'icalcError',icalcError,0)
!    if (cinitialSolution .eq. 1) then
!      ! Get the expressions defining the solution
!      call parlst_getvalue_string (rproblem%rparamList,'TIME-POSTPROCESSING',&
!                                  'ssolutionExpressionY1',sstring,'''''')
!      read(sstring,*) ssolutionExpressionY1
!
!      call parlst_getvalue_string (rproblem%rparamList,'TIME-POSTPROCESSING',&
!                                  'ssolutionExpressionY2',sstring,'''''')
!      read(sstring,*) ssolutionExpressionY2
!
!      call parlst_getvalue_string (rproblem%rparamList,'TIME-POSTPROCESSING',&
!                                  'ssolutionExpressionP',sstring,'''''')
!      read(sstring,*) ssolutionExpressionP
!
!      call parlst_getvalue_string (rproblem%rparamList,'TIME-POSTPROCESSING',&
!                                  'ssolutionExpressionLAMBDA1',sstring,'''''')
!      read(sstring,*) ssolutionExpressionLAMBDA1
!
!      call parlst_getvalue_string (rproblem%rparamList,'TIME-POSTPROCESSING',&
!                                  'ssolutionExpressionLAMBDA2',sstring,'''''')
!      read(sstring,*) ssolutionExpressionLAMBDA2
!
!      call parlst_getvalue_string (rproblem%rparamList,'TIME-POSTPROCESSING',&
!                                  'ssolutionExpressionXi',sstring,'''''')
!      read(sstring,*) ssolutionExpressionXi
!
!      call optcana_analyticalError (rproblem,&
!          RspaceTimePrecondMatrix(size(RspatialPrecond))%p_rsolution,rtempVector,&
!          RspaceTimeDiscr(size(RspatialPrecond))%dalphaC,&
!          RspaceTimeDiscr(size(RspatialPrecond))%dgammaC,&
!          ssolutionExpressionY1,ssolutionExpressionY2,ssolutionExpressionP,&
!          ssolutionExpressionLAMBDA1,ssolutionExpressionLAMBDA2,ssolutionExpressionXI,&
!          Derror(1),Derror(2),Derror(3),Derror(4))
!      call output_line ('||y-y0||           = '//trim(sys_sdEL(Derror(1),10)))
!      call output_line ('||p-p0||           = '//trim(sys_sdEL(Derror(2),10)))
!      call output_line ('||lambda-lambda0|| = '//trim(sys_sdEL(Derror(3),10)))
!      call output_line ('||xi-xi0||         = '//trim(sys_sdEL(Derror(4),10)))
!    end if
!
!    call output_separator (OU_SEP_EQUAL)
!
!    ! Get configuration parameters from the DAT file
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-SOLVER', &
!                              'nminIterations', nminIterations, 1)
!    call parlst_getvalue_int (rproblem%rparamList, 'TIME-SOLVER', &
!                             'nmaxIterations', nmaxIterations, 10)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
!                                 'depsRel', depsRel, 1E-5_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
!                                 'depsAbs', depsAbs, 1E-5_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
!                                 'depsDiff', depsDiff, 0.0_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
!                                 'domega', domega, 1.0_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
!                                 'domega', domega, 1.0_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
!                                 'dinexactNewtonEpsRel', dinexactNewtonEpsRel, 1.0E-2_DP)
!    call parlst_getvalue_double (rproblem%rparamList, 'TIME-SOLVER', &
!                                 'dinexactNewtonExponent', dinexactNewtonExponent, 2.0_DP)
!
!    ! Initalise statistic variables
!    call stat_clearTimer (rtimeSmoothing)
!    call stat_clearTimer (rtimeCoarseGridSolver)
!    call stat_clearTimer (rtimeLinearAlgebra)
!    call stat_clearTimer (rtimeProlRest)
!    call stat_clearTimer (rtimeSpacePrecond)
!
!    call stat_clearTimer (rtimerPreconditioner)
!    call stat_clearTimer (rtimerNonlinear)
!    call stat_startTimer (rtimerNonlinear)
!
!    iglobIter = 0
!    ilinearIterations = 0
!
!    do while ((iglobIter .lt. nminIterations) .or. &
!              ((((ddefNorm .gt. depsRel*dinitDefNorm) .or. (ddefNorm .ge. depsAbs)) .and.&
!                (abs(ddefNorm-dlastDefNorm) .ge. depsDiff*dlastDefNorm)) &
!              .and. (iglobIter .lt. nmaxIterations)))
!
!      iglobIter = iglobIter+1
!
!      ! Project the solution down to all levels, so the nonlinearity
!      ! on the lower levels can be calculated correctly.
!      ! Use the memory in rtempvectorX and rtempvectorB as temp memory;
!      ! it's large enough.
!      call lsysbl_duplicateVector (rtempvectorX,rtempVecFine,&
!          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      call lsysbl_duplicateVector (rtempvectorB,rtempVecCoarse,&
!          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!      do ilev=size(RspatialPrecond)-1,1,-1
!        call lsysbl_enforceStructureDiscr(&
!            RspaceTimeDiscr(ilev+1)%p_rlevelInfo%rdiscretisation,rtempVecFine)
!
!        call lsysbl_enforceStructureDiscr(&
!            RspaceTimeDiscr(ilev)%p_rlevelInfo%rdiscretisation,rtempVecCoarse)
!
!        ! Interpolate down
!        call sptipr_performInterpolation (RinterlevelProjection(ilev+1),&
!            RspaceTimePrecondMatrix(ilev)%p_rsolution, RspaceTimePrecondMatrix(ilev+1)%p_rsolution,&
!            rtempVecCoarse,rtempVecFine)
!
!        ! Set boundary conditions
!        call tbc_implementBCsolution (rproblem,RspaceTimeDiscr(ilev),&
!            RspaceTimePrecondMatrix(ilev)%p_rsolution,rtempVecCoarse)
!      end do
!
!      call lsysbl_releaseVector(rtempVecFine)
!      call lsysbl_releaseVector(rtempVecCoarse)
!
!!      IF (rproblem%MT_outputLevel .GE. 2) THEN
!!        CALL output_line ('Writing solution to file '//&
!!            './ns/tmpsolution'//TRIM(sys_siL(iglobIter-1,10)))
!!        CALL sptivec_saveToFileSequence(&
!!            myRspaceTimeDiscr(SIZE(RspatialPrecond))%p_rsolution,&
!!            '(''./ns/tmpsolution'//TRIM(sys_siL(iglobIter-1,10))//'.'',I5.5)',&
!!            .TRUE.,rtempVectorX)
!!
!!        DO ilev=1,SIZE(RspatialPrecond)
!!          CALL cc_postprocSpaceTimeGMV(rproblem,myRspaceTimeDiscr(ilev),&
!!              myRspaceTimeDiscr(ilev)%p_rsolution,&
!!              './gmv/iteration'//TRIM(sys_siL(iglobIter-1,10))//'level'//TRIM(sys_siL(ilev,10))//&
!!              '.gmv')
!!        END DO
!!      END IF
!
!      if (ctypePreconditioner .eq. CCPREC_INEXACTNEWTON) then
!
!        ! Determine the stopping criterion for the inexact Newton.
!        ! This is an adaptive stopping criterion depending on the current
!        ! defect. In detail, we calculate:
!        !
!        !   |b-Ax_{i+1}|         ( |b-Ax_i| ) exp             ( |b-Ax_i| )
!        !   ------------ = min { ( -------- )     , depsrel * ( -------- ) }
!        !     |b-Ax_0|           ( |b-Ax_0| )                 ( |b-Ax_0| )
!        !
!        ! see e.g. [Michael Hinze, Habilitation, p. 51].
!        ! We furthermore introduce the additional stopping criterion
!        !
!        !   |b-Ax_{i+1}|
!        !   ------------ = depsrel(solver)*dinexactNewtonEpsRel
!        !     |b-Ax_0|
!        !
!        ! which stops the iteration in order to prevent the linear
!        ! solver from solving too much further than necessary for the
!        ! nonlinear solver.
!        !
!        ! Switch off the relative stopping criterion in the linear solver:
!
!        p_rsolverNode%depsRel = 0.0_DP
!
!        ! Calculate the new absolute stopping criterion:
!
!        dtempdef = ddefNorm / dinitDefNorm
!
!        p_rsolverNode%depsAbs = min(&
!            max(dinexactNewtonEpsRel*depsRel,dtempDef**dinexactNewtonExponent),&
!                dinexactNewtonEpsRel*dtempdef) * dinitDefNorm
!
!        ! For the coarse grid solver, we choose the same stopping criterion.
!        ! But just for safetyness, the coarse grid solver should gain at least
!        ! one digit!
!        p_rcgrSolver%depsRel = 1.0E-1_DP
!        p_rcgrSolver%depsAbs = p_rsolverNode%depsAbs
!
!      end if
!
!      if (rproblem%MT_outputLevel .ge. 1) then
!        ! Value of the functional
!        call optcana_nonstatFunctional (rproblem,&
!            rspaceTimeMatrix%p_rsolution,&
!            rtempVector,RspaceTimeDiscr(size(RspatialPrecond))%dalphaC,&
!            RspaceTimeDiscr(size(RspatialPrecond))%dgammaC,&
!            Derror)
!        call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(1),10)))
!        call output_line ('||u||         = '//trim(sys_sdEL(Derror(2),10)))
!        call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
!        call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(4),10)))
!
!        ! If error analysis has to be performed, we can calculate
!        ! the real error.
!        if (cinitialSolution .eq. 1) then
!          call optcana_analyticalError (rproblem,&
!              RspaceTimePrecondMatrix(size(RspatialPrecond))%p_rsolution,rtempVector,&
!              RspaceTimeDiscr(size(RspatialPrecond))%dalphaC,&
!              RspaceTimeDiscr(size(RspatialPrecond))%dgammaC,&
!              ssolutionExpressionY1,ssolutionExpressionY2,ssolutionExpressionP,&
!              ssolutionExpressionLAMBDA1,ssolutionExpressionLAMBDA2,ssolutionExpressionXI,&
!              Derror(1),Derror(2),Derror(3),Derror(4))
!          call output_line ('||y-y0||           = '//trim(sys_sdEL(Derror(1),10)))
!          call output_line ('||p-p0||           = '//trim(sys_sdEL(Derror(2),10)))
!          call output_line ('||lambda-lambda0|| = '//trim(sys_sdEL(Derror(3),10)))
!          call output_line ('||xi-xi0||         = '//trim(sys_sdEL(Derror(4),10)))
!        end if
!
!        if (ctypePreconditioner .eq. CCPREC_INEXACTNEWTON) then
!          call output_lbrk ()
!          call output_line ('Inexact Newton: Stopping criterion = '//&
!              trim(sys_sdEL(p_rsolverNode%depsAbs,10)))
!        end if
!        call output_separator (OU_SEP_EQUAL)
!      end if
!
!      call stat_clearTimer (rtimerMGStep)
!
!      ! DEBUG!!!
!      !CALL sptivec_copyVector (rd,rtemp)
!
!      ! Preconditioning of the defect: d=C^{-1}d
!      if (associated(p_rsolverNode)) then
!
!        call stat_clearTimer (rtimeFactorisationStep)
!        call stat_startTimer (rtimeFactorisationStep)
!        call sptils_initData (p_rsolverNode,ierror)
!        call stat_stopTimer (rtimeFactorisationStep)
!
!        !call sptivec_printVector (rd)
!
!        call stat_clearTimer (rtimerMGStep)
!        call stat_startTimer (rtimerMGStep)
!        call sptils_precondDefect (p_rsolverNode,rd)
!        call sptils_doneData (p_rsolverNode)
!
!        call stat_stopTimer (rtimerMGStep)
!
!        !call sptivec_printVector (rd)
!
!        ! Sum up time data for statistics.
!        call stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeSmoothing,&
!            rtimeSmoothing)
!        call stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeCoarseGridSolver,&
!            rtimeCoarseGridSolver)
!        call stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeLinearAlgebra,&
!            rtimeLinearAlgebra)
!        call stat_addtimers (p_rsolverNode%p_rsubnodeMultigrid%rtimeProlRest,&
!            rtimeProlRest)
!        call stat_addtimers (p_rsolverNode%rtimeSpacePrecond,rtimeSpacePrecond)
!
!        call output_lbrk ()
!        call output_line ("Time for smoothing          : "//&
!            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeSmoothing%delapsedReal,10))
!        call output_line ("Time for coarse grid solving: "//&
!            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeCoarseGridSolver%delapsedReal,10))
!        call output_line ("Time for linear algebra     : "//&
!            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeLinearAlgebra%delapsedReal,10))
!        call output_line ("Time for prol/rest          : "//&
!            sys_sdL(p_rsolverNode%p_rsubnodeMultigrid%rtimeProlRest%delapsedReal,10))
!        call output_lbrk ()
!        call output_line ("Time for prec. in space     : "//&
!            sys_sdL(p_rsolverNode%rtimeSpacePrecond%delapsedReal,10))
!
!        ! Count the number of linear iterations and the time for
!        ! preconditioning
!        ilinearIterations = ilinearIterations + p_rsolverNode%iiterations
!        call stat_addtimers (rtimerMGStep,rtimerPreconditioner)
!        call stat_addtimers (rtimeFactorisationStep,rtimeFactorisation)
!
!      end if
!
!      !CALL output_line('Linear defect:')
!      !CALL stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rd, rtemp, &
!      !  -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,.TRUE.)
!      !CALL output_separator (OU_SEP_MINUS)
!
!      ! Filter the defect for boundary conditions in space and time.
!      ! Normally this is done before the preconditioning -- but by doing it
!      ! afterwards, the initial conditions can be seen more clearly!
!      !CALL tbc_implementInitCondDefect (p_rspaceTimeDiscr,rd,rtempVector)
!      !CALL tbc_implementBCdefect (rproblem,p_rspaceTimeDiscr,rd,rtempVector)
!
!      ! Add the defect: x = x + omega*d
!      !call sptivec_vectorLinearComb (rd,rx,domega,1.0_DP)
!
!      !call cc_calcOptimalCorrection (rproblem,p_rspaceTimeDiscr,&
!      !    rx,rb,rd,rtempVec,rtempVec2,domega2)
!      domega2 = domega
!      call sptivec_vectorLinearComb (rd,rx,domega2,1.0_DP)
!
!      ! Normalise the primal and dual pressure to integral mean value zero
!      ! where no Neumann boundary is present.
!      call tbc_pressureToL20 (rproblem,rx,rtempVectorX)
!
!      ! Implement the initial condition to the new solution vector
!      call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvectorX)
!
!      ! Are bounds to the control active? If yes, restrict the control
!      ! to the allowed range.
!!      if (rproblem%roptcontrol%ccontrolContraints .eq. 1) then
!!        call cc_projectControl (rproblem,rx)
!!      end if
!
!      ! Call a parser that parses the script file commandfile.txt.
!      ! This allows in-program modification of the problem sructure.
!      rproblem%rdataOneshot%iglobIter = iglobIter
!      rproblem%rdataOneshot%ddefNorm = ddefNorm
!      rproblem%rdataOneshot%p_rx => rx
!      rproblem%rdataOneshot%p_rb => rb
!
!      call scr_readScript ('commandfile.txt',0,rproblem)
!
!      iglobIter = rproblem%rdataOneshot%iglobIter
!      ddefNorm = rproblem%rdataOneshot%ddefNorm
!
!      ! Remember the last defect norm for the stopping criterion
!      dlastDefNorm = ddefNorm
!
!      ! Assemble the new defect: d=b-Ax
!      call sptivec_copyVector (rb,rd)
!      !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, &
!      !    rx, rd, ddefNorm)
!      if (rproblem%MT_outputLevel .ge. 2) call output_line('Nonlinear defect:')
!
!      call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!        -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,&
!        ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!      ! Filter the defect for boundary conditions in space and time.
!      !CALL tbc_implementInitCondDefect (p_rspaceTimeDiscr,rd,rtempVector)
!      !CALL tbc_implementBCdefect (rproblem,p_rspaceTimeDiscr,rd,rtempVector)
!
!      call output_separator (OU_SEP_EQUAL)
!      call output_line ('Iteration: '//trim(sys_siL(iglobIter,10))//&
!          '                          Defect of supersystem:  '//adjustl(sys_sdEP(ddefNorm,20,10)))
!      call output_line ('Time for computation of this iterate: '//&
!          trim(sys_sdL(rtimerMGStep%delapsedReal+rtimeFactorisationStep%delapsedReal,10)))
!      call output_separator (OU_SEP_EQUAL)
!
!    end do
!
!    call stat_stopTimer (rtimerNonlinear)
!
!    ! Decrease iglobIter if the DO-loop was completely processed;
!    ! iglobIter = nmaxIterations+1 in that case!
!    if (iglobIter .gt. nmaxIterations) iglobIter = nmaxIterations
!
!    !CALL cc_assembleSpaceTimeRHS (rproblem, p_rspaceTimeDiscr, rd, &
!    !  rtempvectorX, rtempvectorB, rtempvector,.FALSE.)
!    call trhsevl_assembleRHS (rproblem, p_rspaceTimeDiscr, rd, .true.)
!
!    ! Implement the initial condition into the RHS/Solution.
!    call tbc_implementInitCondRHS (rproblem, rd, rinitialCondRHS, rtempvector)
!    call tbc_implementInitCond (rproblem, rx, rinitialCondSol, rtempvector)
!
!    ! Assemble the defect
!    !CALL cc_assembleSpaceTimeDefect (rproblem, rspaceTimeMatrix, &
!    !    rx, rd, ddefNorm)
!    call stlin_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, &
!      -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT,ddefNorm,rproblem%MT_outputLevel .ge. 2)
!
!    call output_separator (OU_SEP_EQUAL)
!    call output_line ('Defect of supersystem: '//sys_sdEP(ddefNorm,20,10))
!    ! Value of the functional
!    call optcana_nonstatFunctional (rproblem,&
!        rspaceTimeMatrix%p_rsolution,&
!        rtempVector,RspaceTimeDiscr(size(RspatialPrecond))%dalphaC,&
!        RspaceTimeDiscr(size(RspatialPrecond))%dgammaC,&
!        Derror)
!    call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(1),10)))
!    call output_line ('||u||         = '//trim(sys_sdEL(Derror(2),10)))
!    call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
!    call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(4),10)))
!
!    ! If the solution is given as analytical expression, we can calculate
!    ! the real error.
!    if (cinitialSolution .eq. 3) then
!      call optcana_analyticalError (rproblem,&
!          RspaceTimePrecondMatrix(size(RspatialPrecond))%p_rsolution,rtempVector,&
!          RspaceTimeDiscr(size(RspatialPrecond))%dalphaC,&
!          RspaceTimeDiscr(size(RspatialPrecond))%dgammaC,&
!          ssolutionExpressionY1,ssolutionExpressionY2,ssolutionExpressionP,&
!          ssolutionExpressionLAMBDA1,ssolutionExpressionLAMBDA2,ssolutionExpressionXI,&
!          Derror(1),Derror(2),Derror(3),Derror(4))
!      call output_line ('||y-y0||           = '//trim(sys_sdEL(Derror(1),10)))
!      call output_line ('||p-p0||           = '//trim(sys_sdEL(Derror(2),10)))
!      call output_line ('||lambda-lambda0|| = '//trim(sys_sdEL(Derror(3),10)))
!      call output_line ('||xi-xi0||         = '//trim(sys_sdEL(Derror(4),10)))
!    end if
!
!    call output_separator (OU_SEP_EQUAL)
!    call output_line ('Total computation time             = '// &
!        trim(sys_sdL(rtimerNonlinear%delapsedReal,10)))
!    call output_line ('Total time for factorisation       = '// &
!        trim(sys_sdL(rtimeFactorisation%delapsedReal,10)))
!    call output_line ('Total time for preconditioning     = '// &
!        trim(sys_sdL(rtimerPreconditioner%delapsedReal,10)))
!    call output_line ('#nonlinear iterations              = '//&
!        trim(sys_siL(iglobIter,10)))
!    call output_line ('#iterations preconditioner         = '//&
!        trim(sys_siL(ilinearIterations,10)))
!
!    call output_separator (OU_SEP_EQUAL)
!
!    call output_line ("Preconditioner statistics:")
!    call output_line ("Total time for smoothing           = "//&
!        sys_sdL(rtimeSmoothing%delapsedReal,10))
!    call output_line ("Total time for coarse grid solving = "//&
!        sys_sdL(rtimeCoarseGridSolver%delapsedReal,10))
!    call output_line ("Total time for linear algebra      = "//&
!        sys_sdL(rtimeLinearAlgebra%delapsedReal,10))
!    call output_line ("Total time for prol/rest           = "//&
!        sys_sdL(rtimeProlRest%delapsedReal,10))
!    call output_lbrk ()
!    call output_line ("Total time for prec. in space      = "//&
!        sys_sdL(rtimeSpacePrecond%delapsedReal,10))
!
!    call output_separator (OU_SEP_EQUAL)
!
!    ! Normalise the primal and dual pressure to integral mean value zero
!    ! where no Neumann boundary is present.
!    call tbc_pressureToL20 (rproblem,rx,rtempVectorX)
!
!    ! Release the multiplication weights for the energy minimisation.
!    deallocate(p_rmgSolver%p_rsubnodeMultigrid%p_DequationWeights)
!
!    ! Release the space-time and spatial preconditioner.
!    ! We don't need them anymore.
!    call sptils_releaseSolver (p_rsolverNode)
!
!    ! Release temp memory
!    ! CALL sptivec_releaseVector (rtemp)
!
!    ! Release the spatial preconditioner and temp vector on every level
!    do ilev=1,size(RspatialPrecond)
!      call cc_doneSpacePreconditioner (RspatialPrecondPrimal(ilev))
!      call cc_doneSpacePreconditioner (RspatialPrecondDual(ilev))
!      call cc_doneSpacePreconditioner (RspatialPrecond(ilev))
!    end do
!
!    do ilev=1,size(RspatialPrecond)-1
!      call sptivec_releaseVector (RspaceTimePrecondMatrix(ilev)%p_rsolution)
!      deallocate(RspaceTimePrecondMatrix(ilev)%p_rsolution)
!    end do
!
!    call sptivec_releaseVector (rtempVec)
!    call sptivec_releaseVector (rtempVec2)
!
!    call lsysbl_releaseVector (rtempVectorB)
!    call lsysbl_releaseVector (rtempVectorX)
!    call lsysbl_releaseVector (rtempVector)
!    call lsysbl_releaseVector (rinitialCondRHS)
!    call lsysbl_releaseVector (rinitialCondSol)
!
!  end subroutine
  
end module
