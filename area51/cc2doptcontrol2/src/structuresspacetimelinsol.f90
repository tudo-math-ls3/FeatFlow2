!##############################################################################
!# ****************************************************************************
!# <name> structuresspacetimelinsol </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying parameter structures for linear space-time preconditioners.
!#
!# Routines in this module:
!#
!# 1.) stlsinit_getSolver
!#     -> Create a linear space-time preconditioner based on parameter
!#        structures.
!#
!# Auxiliary routines:
!#
!# 1.) stlsinit_getMultiGridSolver
!#     -> Create a linear space-time MG structure
!#
!# 2.) stlsinit_getSingleGridSolver
!#     -> Create a linear space-time single grid structure
!# </purpose>
!##############################################################################

module structuresspacetimelinsol

  use fsystem
  use storage
  use genoutput
  use boundary
  use triangulation
  use paramlist
  use discretebc
  use discretefbc
  use fparser
  use linearsystemscalar
  use statistics
  
  use collection
  
  use structuresoptflow
  use spacetimelinearsolver

  implicit none
  
  private
  
!<typeblock>

  ! Space-time smoother structure
  type t_nlstprec_mgrsmooth

    ! Type of space-time smoother.
    ! =0: Block Jacobi
    ! =1: Block Forward-Backward SOR (domega, Preconditioner-domega)
    ! =2: Block Forward-Backward Gauss-Seidel
    ! =3: CG (Jacobi)
    ! =4: CG (Forward-Backward Gauss-Seidel)
    ! =6: Defect correction with UMFPACK preconditioning
    ! =7: BiCGStab with Block Forward-Backward Gauss-Seidel preconditioning
    integer :: cspaceTimeSmoother = 1

    ! Damping parameter
    real(dp) :: domega = 0.9 
    
    ! Relaxation parameter (e.g. in SOR)
    real(dp) :: drelax = 1.0

    ! Number of presmoothing steps
    integer :: nsmpre = 0

    ! Number of postsmoothing steps
    integer :: nsmpost = 1

    ! Relative convergence criterion. If this is reached, the smoother
    ! will stop the smoothing prematurely. Standard=0.0=always perform
    ! nsmSteps smoothing steps.
    real(dp) :: depsRel = 0.0

    ! Absolute convergence criterion. If this is reached, the smoother
    ! will stop the smoothing prematurely. Standard=0.0=always perform
    ! nsmSteps smoothing steps.
    real(dp) :: depsAbs = 0.0

    ! Limit for differences in the residuals
    real(dp) :: depsDiff           = 1E-5

    ! Type of stopping criterion, if nsmSteps is not reached.
    ! =0: Relative AND absolute stopping criterion must hold.
    ! =1: Relative OR absolute stopping criterion must hold.
    integer :: istoppingCriterion = 1

    ! Output level of the smoother
    integer :: ioutputLevel = 2

    ! For FBSOR-Type smoothers: Do a partial update of the solution
    ! during the iteration (only primal on forward sweep, only dual
    ! on backward sweep).
    integer :: ifbSORPartialUpdate = 0
    
    ! Reinitialisation counter, e.g. for BiCGStab.
    ! Every niteReinit iterations, the algorithm reinitialises.
    ! =0: no reinitialisation.
    integer :: niteReinit = 0
    
    ! Name of the section in the DAT file specifying a linear subsolver
    ! for the smoother.
    character(len=SYS_STRLEN) :: slinearSpaceSolver = ""

    ! Name of the section in the DAT file specifying an alternative linear subsolver
    ! for the smoother. ="": Not specified.
    character(len=SYS_STRLEN) :: slinearSpaceSolverAlternative = ""
  
  end type

!</typeblock>

  public :: t_nlstprec_mgrsmooth

!<typeblock>

  ! Space-time coarse grid solver structure
  
  type t_nlstprec_sgrsolver
  
    ! Type of coarse grid solver.
    ! =0: Iterative Defect correction with Block Jacobi prec.
    ! =1: Iterative Defect correction with Block SOR(domega, [TIME-COARSEPRECOND].domega)
    ! =2: Iterative Defect correction with Block Forward-Backward Gauss-Seidel precond.
    ! =3: Iterative CG
    ! =4: direct (UMFPACK Gauss elimination)
    ! =5: Iterative BiCGStab with Block Forward-Backward SOR
    ! =6: Defect correction with UMFPACK preconditioning
    ! =7: Pure block FBSOR solver
    ! =8: Defect correction with UMFPACK preconditioning
    ! =9: Defect correction with simple forward-backward solver.
    integer :: ctypeSolver = 5 

    ! Minimum number of time-iterations on the coarsest time-mesh
    integer :: nminIterations = 1 

    ! Maximum number of time-iterations on the coarsest time-mesh
    integer :: nmaxIterations = 100 

    ! Damping parameter
    real(dp) :: domega = 1.0 

    ! Relaxation parameter (e.g. in SOR)
    real(dp) :: drelax = 1.0

    ! Damping of residuals, i.e. reduction of relative error 
    ! on finest grid; smaller -> more iterations
    ! Not used if inexact Newton is used as nonlinear solver.
    real(dp) :: depsRel = 1E-2 

    ! Limit for residuals, i.e. absolute error on finest grid; 
    ! The linear solver stops if both, absolute error < depsAbs and 
    ! rel. error < depsRel
    ! Not used if inexact Newton is used as nonlinear solver.
    real(dp) :: depsAbs = 1E-0 

    ! Limit for differences in the residuals
    real(dp) :: depsDiff = 1E-5

    ! Maximum relative residual before the iteration is treated as 'divergent'
    real(dp) :: ddivRel = 1E20

    ! Output level of the solver.
    integer :: ioutputLevel = 2 
     
    ! Type of stopping criterion.
    ! =0: Relative AND absolute stopping criterion must hold.
    ! =1: Relative OR absolute stopping criterion must hold.
    integer :: istoppingCriterion = 0

    ! Reinitialisation counter, e.g. for BiCGStab.
    ! Every niteReinit iterations, the algorithm reinitialises.
    ! =0: no reinitialisation.
    integer :: niteReinit = 0

  end type

!</typeblock>

  public :: t_nlstprec_sgrsolver

!<typeblock>

  ! Space-time coarse grid and smoother preconditioner
  type t_nlstprec_sgrprec
  
    ! Minimum number of time-iterations on the coarsest time-mesh
    integer :: nminIterations = 1 

    ! Maximum number of time-iterations on the coarsest time-mesh
    integer :: nmaxIterations = 1 

    ! Damping parameter
    real(dp) :: domega = 1.0 

    ! Relaxation parameter (e.g. in SOR)
    real(dp) :: drelax = 1.0

    ! Damping of residuals, i.e. reduction of relative error 
    ! on finest grid; smaller -> more iterations
    real(dp) :: depsRel = 1E-13

    ! Limit for residuals, i.e. absolute error on finest grid; 
    ! The linear solver stops if both, absolute error < depsAbs and 
    ! rel. error < depsRel
    real(dp) :: depsAbs = 1E-0

    ! Limit for differences in the residuals
    real(dp) :: depsDiff = 1E-5

    ! Maximum relative residual before the iteration is treated as 'divergent'
    real(dp) :: ddivRel = 1E20

    ! Output level of the solver.
    ! If a subsolver is embedded in the solver, it receives the output level
    ! ioutputLevel-2.
    integer :: ioutputLevel = 0 

    ! Type of stopping criterion.
    ! =0: Relative AND absolute stopping criterion must hold.
    ! =1: Relative OR absolute stopping criterion must hold.
    integer :: istoppingCriterion = 0
  
    ! For FBSOR-Type smoothers: Do a partial update of the solution
    ! during the iteration (only primal on forward sweep, only dual
    ! on backward sweep).
    integer :: ifbSORPartialUpdate = 0
    
    ! Reference to a parameter list containing parameters of a linear
    ! subsolver (if one is needed).
    type(t_parlist), pointer :: p_rparlist => null()
    
    ! Name of the section in the DAT file specifying an alternative linear subsolver
    ! for the solver. ="": Not specified.
    character(len=SYS_STRLEN) :: slinearSpaceSolver = ""

    ! Name of the section in the DAT file specifying a linear subsolver
    ! for the solver.
    character(len=SYS_STRLEN) :: slinearSpaceSolverAlternative = ""
  end type

!</typeblock>

  public :: t_nlstprec_sgrprec

!<typeblock>

  ! Space-time multigrid solver structure.
  type t_nlstprec_mgrsolver

    ! Minimum number of time-MG sweeps per nonlinear time-iteration
    integer :: nminIterations = 1

    ! Maximum number of time-MG sweeps per nonlinear time-iteration
    integer :: nmaxIterations = 10

    ! Damping of residuals, i.e. reduction of relative error 
    ! on finest grid; smaller -> more iterations
    ! Not used if inexact Newton is used as nonlinear solver.
    real(dp) :: depsRel = 1E-2

    ! Limit for residuals, i.e. absolute error on finest grid; 
    ! The linear solver stops if both, absolute error < depsAbs and 
    ! rel. error < depsRel
    ! Not used if inexact Newton is used as nonlinear solver.
    real(dp) :: depsAbs = 1E-0

    ! Minimum relative difference between two iterates
    real(dp) :: depsDiff = 1E-5

    ! Type of stopping criterion.
    ! =0: Relative AND absolute stopping criterion must hold.
    ! =1: Relative OR absolute stopping criterion must hold.
    integer :: istoppingCriterion = 0

    ! Output level of the solver
    integer :: ioutputLevel = 2 

    ! Cycle. 0=F-cycle, 1=V-cycle, 2=W-cycle
    integer :: icycle = 1

    ! Minimum value for the adaptive coarse grid correction.
    ! Note: dalphamin=dalphamax deactivates adaptive coarse grid correction!
    ! Standard = 1.0
    real(dp) :: dalphamin = 1.0 

    ! Minimum value for the adaptive coarse grid correction.
    ! Note: dalphamin=dalphamax deactivates adaptive coarse grid correction!
    ! Standard = 1.0
    real(dp) :: dalphamax = 1.0 
  
  end type
  
!</typeblock>

  public :: t_nlstprec_mgrsolver

!<typeblock>

  ! Space-time multigrid solver structure.
  type t_settings_nlstprec

    ! Type of preconditioner.
    ! =0: Direct solver (only for testing purposes)
    ! =1: Defect correction on one level with forward-backward preconditioning.
    ! =2: Space-time multigrid solver.
    integer :: ctypeSolver = 0
  
    ! <!-- Parameters for single-grid preconditioners and/or coarse grid solvers -->

    ! Parameters of the space-time coarse grid solver (in multigrid iterations)
    ! and single grid solver (in one-grid iterations)
    type(t_nlstprec_sgrsolver) :: ronegridSolver
  
    ! Parameters of the space-time coarse grid / single grid preconditioner.
    type(t_nlstprec_sgrprec) :: ronegridPrecond

    ! <!-- Parameters for multigrid preconditioners -->
  
    ! Parameters of the space-time multigrid solver.
    type(t_nlstprec_mgrsolver) :: rmgSolver

    ! Parameters of the space-time smoother.
    type(t_nlstprec_mgrsmooth) :: rsmoother

    ! Parameters for the preconditioner of the space-time smoother.
    type(t_nlstprec_sgrprec) :: rprecsmoother

  end type
  
!</typeblock>

  public :: t_settings_nlstprec

!!<typeblock>
!
!  ! Simple nonlinear one-level defect correction solver in space/time.
!  type t_settings_nlstprec_defcorr
!    
!    ! Minimum number of iterations
!    integer :: nminIterations = 1
!
!    ! Maximum number of iterations
!    integer :: nmaxIterations = 100
!
!    ! Damping of residuals, i.e. reduction of relative error 
!    ! on finest grid; smaller -> more iterations
!    real(dp) :: depsRel        = 1E-5
!
!    ! Limit for residuals, i.e. absolute error on finest grid; 
!    ! The solver stops if both, absolute error < depsAbs and 
!    ! rel. error < depsRel
!    real(dp) :: depsAbs        = 1E-5
!    
!    ! Limit for differences in the residuals
!    real(dp) :: depsDiff           = 1E-5
!
!    ! Name of the section in the DAT file specifying a linear subsolver
!    ! for the solver.
!    character(len=SYS_STRLEN) :: slinearSpaceSolver = ""
!
!  end type
!
!!</typeblock>
!
!  public :: t_settings_nlstprec_defcorr

!</types>

  public :: stlsinit_getSolver

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlsinit_getSolver (rsettings,ispaceTimeLevel,rsolversettings,p_rsolver,&
      p_rmgsolver, p_rcgrsolver)

!<description>
  ! Creates a linear solver node p_rsolver based on the settings in
  ! rsettings and rprecsettings
!</description>

!<input>
  ! Global settings structure.
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Absolute level of the solver (relative to the global hierarchies in
  ! rsettings).
  integer, intent(in) :: ispaceTimeLevel

  ! Structure with parameters for the solver.
  type(t_settings_nlstprec), intent(in), target :: rsolversettings
!</input>

!<output>
  ! Pointer to a solver node representing a nonlinear space-time solver.
  type(t_sptilsNode), pointer :: p_rsolver
  
  ! Points to a multigrid solver node if one is involved.
  type(t_sptilsNode), pointer :: p_rmgsolver
  
  ! Points to a coarse grid solver node if a multigrid solver is involved.
  type(t_sptilsNode), pointer :: p_rcgrsolver
!</output>

!</subroutine>
    ! Based on the parameters in rsolversettings, create p_rsolver.
    p_rsolver => null()
    p_rmgsolver => null()
    p_rcgrsolver => null()
    select case (rsolversettings%ctypeSolver)
    case (0)
      ! Single grid preconditioner
      call stlsinit_getSingleGridSolver (rsettings,ispaceTimeLevel,&
          rsolversettings%ronegridSolver,rsolversettings%ronegridPrecond,p_rsolver)
          
    case (1)
      ! Multigrid preconditioner
      call stlsinit_getMultiGridSolver (rsettings,&
        ispaceTimeLevel,1,rsolversettings%rmgSolver,&
        rsolversettings%rsmoother,rsolversettings%rprecsmoother,&
        rsolversettings%ronegridSolver,rsolversettings%ronegridPrecond,&
        p_rsolver,p_rmgsolver,p_rcgrsolver)
      
    end select

  end subroutine  

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlsinit_getMultiGridSolver (rsettings,&
      ispaceTimeLevel,ispaceTimeLevelCoarse,rsolversettings,&
      rsmootherSettings,rprecSmootherSettings,&
      rcgrSolverSettings,rcgrPrecSettings,p_rsolver,p_rmgsolver,p_rcgrsolver)

!<description>
  ! Creates a linear solver node p_rsolver based on the settings in
  ! rsettings, rsolversettings and rprecsettings. The solver is a multigrid
  ! solver.
!</description>

!<input>
  ! Global settings structure.
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Absolute level of the solver (relative to the global hierarchies in
  ! rsettings).
  integer, intent(in) :: ispaceTimeLevel

  ! Absolute level of the coarse grid solver (relative to the global 
  ! hierarchies in rsettings).
  integer, intent(in) :: ispaceTimeLevelCoarse
  
  ! Structure with parameters for the solver.
  type(t_nlstprec_mgrsolver), intent(in), target :: rsolversettings
  
  ! Structure with parameters of the smoother
  type(t_nlstprec_mgrsmooth), intent(in), target :: rsmootherSettings

  ! Structure with parameters for the preconditioner of the smoother
  type(t_nlstprec_sgrprec), intent(in), target :: rprecSmootherSettings

  ! Structure with parameters of the coarse grid solver
  type(t_nlstprec_sgrsolver), intent(in), target :: rcgrSolverSettings

  ! Structure with parameters for the preconditioner of the coarse grid solver
  type(t_nlstprec_sgrprec), intent(in), target :: rcgrPrecSettings
!</input>

!<output>
  ! Pointer to a solver node representing a nonlinear space-time solver.
  type(t_sptilsNode), pointer :: p_rsolver
  
  ! Points to a multigrid solver node if one is involved.
  type(t_sptilsNode), pointer :: p_rmgsolver
  
  ! Points to a coarse grid solver node if a multigrid solver is involved.
  type(t_sptilsNode), pointer :: p_rcgrsolver
!</output>

!</subroutine>
    ! local variables
    integer :: ilev
    type(t_sptilsNode), pointer :: p_rsingleGridSolver,p_rpresmoother,p_rpostsmoother
    
    ! Multigrid preconditioner. Initialise p_rsolver as space-time MG.
    call sptils_initMultigrid (rsettings,ispaceTimeLevelCoarse,ispaceTimeLevel,p_rmgsolver)
    
    ! Initialise the hierarchy
    call sptils_setHierarchyMultigrid (p_rmgsolver,&
        rsettings%rspaceTimeHierPrimalDual,rsettings%rprjHierSpaceTimePrimalDual)
    
    ! The solver ís MG.
    p_rsolver => p_rmgsolver
        
    ! Transfer all parameters
    p_rsolver%nminIterations     = rsolversettings%nminIterations
    p_rsolver%nmaxIterations     = rsolversettings%nmaxIterations
    p_rsolver%depsRel            = rsolversettings%depsRel
    p_rsolver%depsAbs            = rsolversettings%depsAbs
    p_rsolver%depsDiff           = rsolversettings%depsDiff
    p_rsolver%ioutputLevel       = rsolversettings%ioutputLevel
    p_rsolver%istoppingCriterion = rsolversettings%istoppingCriterion
    p_rsolver%p_rsubnodeMultigrid%icycle = rsolversettings%icycle
    p_rsolver%p_rsubnodeMultigrid%dalphamin = rsolversettings%dalphamin
    p_rsolver%p_rsubnodeMultigrid%dalphamax = rsolversettings%dalphamax
    
    ! Initialise all levels.
    do ilev = ispaceTimeLevelCoarse,ispaceTimeLevel
    
      nullify(p_rpresmoother)
      nullify(p_rpostsmoother)
      nullify(p_rsingleGridSolver)

      ! We are on level ilev. If ilev=1, we have to initialise the
      ! coarse grid solver. On the other levels, we have to initialise
      ! the pre- and postsmoother.
      if (ilev .eq. 1) then
      
        call stlsinit_getSingleGridSolver (rsettings,ilev,&
            rcgrSolverSettings,rcgrPrecSettings,p_rsingleGridSolver)
        p_rcgrsolver => p_rsingleGridSolver
        
      else
      
        ! Pre=postsmoother?
        
        ! Same #steps= if yes, take the presmoother as postsmoother.
        if (rsmootherSettings%nsmpre .eq. rsmootherSettings%nsmpost) then
          
          call stlsinit_getSingleGridSmoother (rsettings,ilev,&
              rsmootherSettings,rprecSmootherSettings,p_rpresmoother)
              
          call sptils_convertToSmoother (p_rpresmoother,rsmootherSettings%nsmpre)
              
          p_rpostsmoother => p_rpresmoother
          
        else 
        
          ! Create the pre- and/or postsmoother.
          if (rsmootherSettings%nsmpre .gt. 0) then
            ! Take the solver as presmoother and create a postsmoother

            call stlsinit_getSingleGridSmoother (rsettings,ilev,&
                rsmootherSettings,rprecSmootherSettings,p_rpresmoother)
                
            call sptils_convertToSmoother (p_rpresmoother,rsmootherSettings%nsmpre)

          end if

          if (rsmootherSettings%nsmpost .gt. 0) then
            ! Take the solver as presmoother and create a postsmoother

            call stlsinit_getSingleGridSmoother (rsettings,ilev,&
                rsmootherSettings,rprecSmootherSettings,p_rpostsmoother)
            
            call sptils_convertToSmoother (p_rpostsmoother,rsmootherSettings%nsmpost)

          end if

        end if
      end if
    
      ! Initialise the level.
      call sptils_setLevelSubsolvers (p_rsolver,ilev,&
          p_rpresmoother,p_rpostsmoother,p_rsingleGridSolver)
    
    end do ! ilev

  end subroutine  

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlsinit_getSingleGridSmoother (rsettings,ispaceTimeLevel,&
      rsolversettings,rprecsettings,p_rsolver)

!<description>
  ! Creates a linear solver node p_rsolver based on the settings in
  ! rsettings, rsolversettings and rprecsettings. The solver can be used
  ! as smoother in a multigrid iteration.
!</description>

!<input>
  ! Global settings structure.
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Absolute level of the solver (relative to the global hierarchies in
  ! rsettings).
  integer, intent(in) :: ispaceTimeLevel
  
  ! Structure with parameters for the solver.
  type(t_nlstprec_mgrsmooth), intent(in), target :: rsolversettings

  ! Structure with parameters for the preconditioner (if there is one)
  type(t_nlstprec_sgrprec), intent(in), target :: rprecsettings
!</input>

!<output>
  ! Pointer to a solver node representing the smoother.
  type(t_sptilsNode), pointer :: p_rsolver
!</output>

!</subroutine>
    ! local variables
    type(t_sptilsNode), pointer :: p_rprecond

    ! Based on the parameters in rsolversettings, create p_rsolver.
    select case (rsolversettings%cspaceTimeSmoother)
    case (0)
      ! Defect correction with Block Jacobi preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initBlockJacobi (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction smoother.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case (1)
    
      ! Block SOR preconditioning.
      !
      ! Create a SOR preconditioner
      call sptils_initBlockFBSOR (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,rprecsettings%drelax,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction smoother.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case (2)
      ! BiCGStab with Block Jacobi preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initBlockJacobi (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStab (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)
          
      p_rsolver%p_rsubnodeBiCGStab%niteReinit = rsolversettings%niteReinit

    case (3)
      ! BiCGStab with Block SOR preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initBlockFBSOR (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,rprecsettings%drelax,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStab (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)
          
      p_rsolver%p_rsubnodeBiCGStab%niteReinit = rsolversettings%niteReinit

    case (8)
      ! Simple defect correction with UMFPACK preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initUMFPACK4 (rsettings,ispaceTimeLevel,p_rprecond)
      
      ! Set up some debug settings.
      p_rprecond%p_rsubnodeUMFPACK4%cwriteMatrix = &
          rsettings%rdebugFlags%cumfpackWriteMatrix
          
      p_rprecond%p_rsubnodeUMFPACK4%sfilename = &
          rsettings%rdebugFlags%sumfpackMatrixFilename

      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction smoother.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case (9)
      ! forward-backward preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rsolver)
          
    case (10)
      ! Defect correction with forward-backward preconditioning.
      !
      ! Create a forward-backward preconditioner.
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rprecond)
          
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction smoother.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case (11)
      ! BiCGStab with Forward backward preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rprecond)
          
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStab (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

      p_rsolver%p_rsubnodeBiCGStab%niteReinit = rsolversettings%niteReinit

    case (12)
      ! BiCGStab with right Forward backward preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rprecond)
          
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStabRight (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case default

      call output_line('Invalid smoother!',&
          OU_CLASS_ERROR,OU_MODE_STD,'stlsinit_getSingleGridSmoother')
      call sys_halt()
          
    end select

    ! Transfer all parameters
    p_rsolver%nminIterations     = 0
    p_rsolver%nmaxIterations     = 0
    p_rsolver%domega             = rsolversettings%domega
    p_rsolver%depsRel            = rsolversettings%depsRel
    p_rsolver%depsAbs            = rsolversettings%depsAbs
    p_rsolver%depsDiff           = rsolversettings%depsDiff
    p_rsolver%ioutputLevel       = rsolversettings%ioutputLevel
    p_rsolver%istoppingCriterion = rsolversettings%istoppingCriterion

  end subroutine  

  ! ***************************************************************************
  
!<subroutine>

  subroutine stlsinit_getSingleGridSolver (rsettings,ispaceTimeLevel,&
      rsolversettings,rprecsettings,p_rsolver)

!<description>
  ! Creates a linear solver node p_rsolver based on the settings in
  ! rsettings, rsolversettings and rprecsettings. The solver is a single-grid
  ! solver.
!</description>

!<input>
  ! Global settings structure.
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Absolute level of the solver (relative to the global hierarchies in
  ! rsettings).
  integer, intent(in) :: ispaceTimeLevel
  
  ! Structure with parameters for the solver.
  type(t_nlstprec_sgrsolver), intent(in), target :: rsolversettings

  ! Structure with parameters for the preconditioner (if there is one)
  type(t_nlstprec_sgrprec), intent(in), target :: rprecsettings
!</input>

!<output>
  ! Pointer to a solver node representing a nonlinear space-time solver.
  type(t_sptilsNode), pointer :: p_rsolver
!</output>

!</subroutine>
    ! local variables
    type(t_sptilsNode), pointer :: p_rprecond

    ! Based on the parameters in rsolversettings, create p_rsolver.
    select case (rsolversettings%ctypeSolver)
    case (0)
      ! Defect correction with Block Jacobi preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initBlockJacobi (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)
      
    case (1)
    
      ! Defect correction with Block SOR preconditioning.
      !
      ! Create a Block SOR preconditioner. We set the 2nd parameter to
      ! 1.0 which leads to a standard SOR preconditioner without
      ! GS relaxation.
      call sptils_initBlockFBSOR (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,rprecsettings%drelax,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)
      
    case (2)
      ! BiCGStab with Block Jacobi preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initBlockJacobi (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStab (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

      p_rsolver%p_rsubnodeBiCGStab%niteReinit = rsolversettings%niteReinit

    case (3)
      ! BiCGStab with Block SOR preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initBlockFBSOR (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          p_rprecond,rprecsettings%drelax,1.0_DP,rprecsettings%slinearSpaceSolverAlternative)
          
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStab (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

      p_rsolver%p_rsubnodeBiCGStab%niteReinit = rsolversettings%niteReinit

    case (4)
      ! BiCGStab with Forward backward preconditioning.
      !
      ! Create a Block_Jacobi preconditioner
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rprecond)
          
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStab (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)
      !call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
      !    p_rsolver,p_rprecond)

      p_rsolver%p_rsubnodeBiCGStab%niteReinit = rsolversettings%niteReinit

    case (6)
      ! Pure UMFPACK preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initUMFPACK4 (rsettings,ispaceTimeLevel,p_rsolver)
      
      ! Set up some debug settings.
      p_rsolver%p_rsubnodeUMFPACK4%cwriteMatrix = &
          rsettings%rdebugFlags%cumfpackWriteMatrix
      p_rsolver%p_rsubnodeUMFPACK4%sfilename = &
          rsettings%rdebugFlags%sumfpackMatrixFilename
      
    case (8)
      ! Simple defect correction with UMFPACK preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initUMFPACK4 (rsettings,ispaceTimeLevel,p_rprecond)
      
      ! Set up some debug settings.
      p_rprecond%p_rsubnodeUMFPACK4%cwriteMatrix = &
          rsettings%rdebugFlags%cumfpackWriteMatrix
      p_rprecond%p_rsubnodeUMFPACK4%sfilename = &
          rsettings%rdebugFlags%sumfpackMatrixFilename
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case (10)
      ! Simple defect correction with forward-backward preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rprecond)
          
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case (11)
      ! BiCGStab with forward-backward preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rprecond)
          
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStab (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)
          
      p_rsolver%p_rsubnodeBiCGStab%niteReinit = rsolversettings%niteReinit

    case (12)
      ! BiCGStab with right forward-backward preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initFBsim (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,p_rprecond)
          
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initBiCGStabRight (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)
          
    case (13)
      ! Simple defect correction with time-VANKA and forward-backward preconditioning.
      !
      ! Create a forward-backward solver.
      call sptils_initTimeVanka (rsettings,ispaceTimeLevel,&
          rprecsettings%p_rparlist,rprecsettings%slinearSpaceSolver,&
          rprecsettings%drelax,4,8,p_rprecond)
          
      p_rprecond%nminIterations = rprecsettings%nminIterations
      p_rprecond%nmaxIterations = rprecsettings%nmaxIterations
      p_rprecond%ioutputLevel = rprecsettings%ioutputLevel
      p_rprecond%domega = rprecsettings%domega
      
      ! Initialise the defect correction solver.
      call sptils_initDefCorr (rsettings,ispaceTimeLevel,&
          p_rsolver,p_rprecond)

    case default

      call output_line('Invalid single grid solver!',&
          OU_CLASS_ERROR,OU_MODE_STD,'stlsinit_getSingleGridSolver')
      call sys_halt()

    end select

    ! Transfer all parameters
    p_rsolver%nminIterations     = rsolversettings%nminIterations
    p_rsolver%nmaxIterations     = rsolversettings%nmaxIterations
    p_rsolver%domega             = rsolversettings%domega
    p_rsolver%depsRel            = rsolversettings%depsRel
    p_rsolver%depsAbs            = rsolversettings%depsAbs
    p_rsolver%depsDiff           = rsolversettings%depsDiff
    p_rsolver%ddivRel            = rsolversettings%ddivRel
    p_rsolver%ioutputLevel       = rsolversettings%ioutputLevel
    p_rsolver%istoppingCriterion = rsolversettings%istoppingCriterion

  end subroutine  

end module
