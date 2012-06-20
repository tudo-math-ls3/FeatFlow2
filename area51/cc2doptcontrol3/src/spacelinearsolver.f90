!##############################################################################
!# ****************************************************************************
!# <name> spacelinearsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains a set of routines to set up and maintain a linear solver in space.
!# Such solvers are used for preconditioning defects in space during
!# the loop over the solutions in time. They can be applied on every level
!# of a space discretisation.
!# </purpose>
!##############################################################################

module spacelinearsolver

  use fsystem
  use storage
  use genoutput
  use paramlist
  use collection
  use linearsolver
  use filtersupport
  use multilevelprojection
  use linearsystemscalar
  use linearsystemblock
  use coarsegridcorrection
  use linearsolverautoinitialise

  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  
  use structuresgeneral
  use constantsdiscretisation
  use structuresboundaryconditions
  
  implicit none

!<constants>

!<constantblock description = "Supported equations">

  ! General linear equation
  integer, parameter :: LSS_EQN_GENERAL = 0

  ! 2D Stokes / Navier-Stokes equations
  integer, parameter :: LSS_EQN_STNAVST2D = 1

!</constantblock>

!<constants>

!<types>

!<typeblock>

  ! Encapsules a linear solver including additional parameters
  type t_linsolSpace
  
    ! Type of solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Multigrid solver
    integer :: isolverType = 0
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of smoother.
    ! =0: general VANKA (slow, but independent of the discretisation and of the problem)
    ! =1: general VANKA; 'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 0, but slightly faster)
    ! =2: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =3: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 8, but faster)
    ! =4: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =5: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 10, but faster)
    integer :: ismootherType = 0
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of coarse grid solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Defect correction with diagonal VANKA preconditioning.
    ! =2: BiCGStab with diagonal VANKA preconditioning
    integer :: icoarseGridSolverType = 0

    ! Pointer to the solver node
    type (t_linsolNode), pointer :: p_rsolverNode => null()
    
    ! For Multigrid based solvers, pointer to the coarse grid solver
    type (t_linsolNode), pointer :: p_rcoarseGridSolver => null()

    ! A filter chain that is used for implementing boundary conditions or other
    ! things when invoking the linear solver.
    type(t_filterChain), dimension(6) :: RfilterChain
    
    ! Pointer to the system matrix.
    type(t_matrixBlock), pointer :: p_rmatrix => null()
  end type

!</typeblock>

!<typeblock>

  ! This structure encapsules a hierarchy of linear solvers for all
  ! levels in a space discrtisation hierarchy.
  type t_linsolHierarchySpace
  
    ! Minimum level
    integer :: nlmin = 0
    
    ! Maximum level
    integer :: nlmax = 0
  
    ! Underlying equation, the solvers in this structure support.
    integer :: cequation = LSS_EQN_GENERAL

    ! Hierarchy of FEM spaces associated with the linear solvers.
    type(t_feHierarchy), pointer :: p_feSpaceHierarchy => null()
    
    ! Pointer to the debug flags
    type(t_optcDebugFlags), pointer :: p_rdebugFlags => null()
    
    ! Array of linear solver structures. Element i in this array
    ! describes a linear solver which can be applied on level i
    ! of the above FE space hierarchy. The linear solver may be
    ! multigrid based; in this case, level 1 of the above FE space
    ! hierarchy corresponds to the coarse grid.
    type(t_linsolSpace), dimension(:), pointer :: p_RlinearSolvers => null()
  
    ! Pointer to interlevel projection hierarchy, for MG based solvers
    type(t_interlevelProjectionHier), pointer :: p_rprjHierarchy => null()

  end type

!</typeblock>

  public :: t_linsolHierarchySpace

!</types>

  ! Initialises a linear solver according to the settings in a parameter list.
  public :: lssh_initSolver

  ! Cleans up a linear solver according to the settings in a parameter list.
  public :: lssh_doneSolver

  ! Based on a FE space hierarchy and a parameter list, this subroutine
  ! creates a linear solver hierarchy.
  public :: lssh_createLinsolHierarchy

  ! Releases a given linear solver hierarchy.
  public :: lssh_releaseLinsolHierarchy

  ! Defines the system matrix for level ilevel.
  public :: lssh_setMatrix
  
  ! Initialises structural data for the solver at level ilevel.
  public :: lssh_initStructure
  
  ! Initialises calculation data for the solver at level ilevel.
  public :: lssh_initData
  
  ! Cleans up calculation data for the solver at level ilevel.
  public :: lssh_doneData
  
  ! Cleans up structural data for the solver at level ilevel.
  public :: lssh_doneStructure
  
  ! Applies preconditioning of rd with the solver on level ilevel.
  public :: lssh_precondDefect
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_initSolver (rsolver,cequation,nlevels,rprjHierarchy,&
      rparList,ssection,rdebugFlags)
  
!<description>
  ! Initialises a linear solver according to the settings in a parameter list.
!</description>

!<input>
  ! Underlying equation. One of the LSS_EQN_xxxx constants.
  integer, intent(in) :: cequation
  
  ! Number of available levels in the underlying hierarchy.
  integer, intent(in) :: nlevels

  ! Projection hierarchy, for MG based solvers
  type(t_interlevelProjectionHier), intent(in) :: rprjHierarchy

  ! Parameter list with solver parameters
  type(t_parlist), intent(in) :: rparList
  
  ! Section configuring the solver
  character(len=*), intent(in) :: ssection
  
  ! Debug flags
  type(t_optcDebugFlags), intent(in) :: rdebugFlags
!</input>

!<output>
  ! Solver node to initialise
  type(t_linsolSpace), intent(out) :: rsolver
!</output>
  
!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: ilev,nsm

    integer :: isolverType,ismootherType,icoarseGridSolverType
    character(LEN=SYS_STRLEN) :: sstring,ssolverSection,ssmootherSection
    character(LEN=SYS_STRLEN) :: scoarseGridSolverSection,spreconditionerSection
    type(t_linsolNode), pointer :: p_rpreconditioner, p_rsmoother
    type(t_linsolNode), pointer :: p_rsolverNode
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparList, ssection, p_rsection)
    
    if (.not. associated(p_rsection)) then
      call output_line ("Cannot create linear solver; no section '"//trim(ssection)//&
                        "'!", OU_CLASS_ERROR,OU_MODE_STD,"lssh_initSolver")
      call sys_halt()
    end if

    ! Get the parameters that configure the solver type
    
    call parlst_getvalue_int (p_rsection, "isolverType", isolverType, 1)
    call parlst_getvalue_int (p_rsection, "ismootherType", ismootherType, 3)
    call parlst_getvalue_int (p_rsection, "icoarseGridSolverType", &
        icoarseGridSolverType, 1)
        
    rsolver%isolverType = isolverType
    rsolver%ismootherType = ismootherType
    rsolver%icoarseGridSolverType = icoarseGridSolverType

    call parlst_getvalue_string (p_rsection, "ssolverSection", ssolverSection,&
        "",bdequote=.true.)
    call parlst_getvalue_string (p_rsection, "ssmootherSection", ssmootherSection,&
        "",bdequote=.true.)
    call parlst_getvalue_string (p_rsection, "scoarseGridSolverSection", &
        scoarseGridSolverSection,"",bdequote=.true.)
    
    ! Which type of solver do we have?
    select case (isolverType)
    
    ! ---------------------------------------------------------------
    ! UMFPACK - Gauss elimination
    ! ---------------------------------------------------------------
    case (0)
    
      ! This is the UMFPACK solver. Very easy to initialise. No parameters at all.
      ! This solver works for all types of equations.
      call linsol_initUMFPACK4 (p_rsolverNode)
    
    ! ---------------------------------------------------------------
    ! Multigrid solver
    ! ---------------------------------------------------------------
    case (1)
    
      ! At first, initialise the solver.
      call linsol_initMultigrid2 (p_rsolverNode,nlevels,&
          rsolver%RfilterChain)
      
      ! Init standard solver parameters and extended multigrid parameters
      ! from the DAT file.
      call linsolinit_initParams (p_rsolverNode,rparList,ssolverSection,&
          LINSOL_ALG_UNDEFINED)
      call linsolinit_initParams (p_rsolverNode,rparList,ssolverSection,&
          LINSOL_ALG_MULTIGRID2)
      
      select case (cequation)
      
      ! ---------------------------------------------------
      ! General equation
      ! ---------------------------------------------------
      case (-1)
        
        call output_line ("General equation not supported by multigrid.", &
            OU_CLASS_ERROR,OU_MODE_STD,"lssh_initSolver")
        call sys_halt()
      
      ! ---------------------------------------------------
      ! 2D Stokes / Navier-Stokes
      ! ---------------------------------------------------
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

        ! ---------------------------------------------------
        ! Coarse grid correction fine-tuning
        ! ---------------------------------------------------
        ! For the Navier-Stokes equations,
        ! manually trim the coarse grid correction in Multigrid to multiply the
        ! pressure equation with -1. This (un)symmetrises the operator and gives
        ! much better convergence rates.
        call cgcor_release(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)

        call cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,3)
        p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(3) &
            = -1.0_DP
        
        ! ---------------------------------------------------
        ! Coarse grid solver
        ! ---------------------------------------------------
        ! Ok, now we have to initialise all levels. First, we create a coarse
        ! grid solver and configure it.
        call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
        
        select case (icoarseGridSolverType)
        case (0)
          ! UMFPACK coarse grid solver. Easy.
          call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
          p_rlevelInfo%p_rcoarseGridSolver%p_rsubnodeUmfpack4%imatrixDebugOutput = &
              rdebugFlags%cwriteUmfpackMatrix
          p_rlevelInfo%p_rcoarseGridSolver%p_rsubnodeUmfpack4%smatrixName = "matrix.txt"
          
        case (1)
          ! Defect correction with diagonal VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_NAVST2D_DIAG)
        
          call parlst_getvalue_string (rparList, scoarseGridSolverSection, &
              "spreconditionerSection", sstring, "",bdequote=.true.)
          read (sstring,*) spreconditionerSection
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,p_rpreconditioner%calgorithm)
          
          ! Create the defect correction solver, attach VANKA as preconditioner.
          call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
              rsolver%RfilterChain)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,p_rpreconditioner%calgorithm)
          
        case (2)
          ! Defect correction with full VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_NAVST2D_FULL)
          
          call parlst_getvalue_string (rparList, scoarseGridSolverSection, &
              "spreconditionerSection", sstring, "", bdequote=.true.)
          read (sstring,*) spreconditionerSection
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,p_rpreconditioner%calgorithm)
          
          ! Create the defect correction solver, attach VANKA as preconditioner.
          call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
              rsolver%RfilterChain)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,p_rpreconditioner%calgorithm)
          
        case (3)
          ! BiCGStab with diagonal VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_NAVST2D_DIAG)
          
          call parlst_getvalue_string (rparList, scoarseGridSolverSection, &
            "spreconditionerSection", sstring, "", bdequote=.true.)
          read (sstring,*) spreconditionerSection
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,p_rpreconditioner%calgorithm)
          
          ! Create the defect correction solver, attach VANKA as preconditioner.
          call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
              rsolver%RfilterChain)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)
          
        case (4)
          ! BiCGStab with full VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_NAVST2D_FULL)
          
          call parlst_getvalue_string (rparList, scoarseGridSolverSection, &
            "spreconditionerSection", sstring, "", bdequote=.true.)
          read (sstring,*) spreconditionerSection
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,p_rpreconditioner%calgorithm)
          
          ! Create the defect correction solver, attach VANKA as preconditioner.
          call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
              rsolver%RfilterChain)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

        case (5)
          ! Defect correction with general VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
          
          call parlst_getvalue_string (rparList, scoarseGridSolverSection, &
            "spreconditionerSection", sstring, "", bdequote=.true.)
          read (sstring,*) spreconditionerSection
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,p_rpreconditioner%calgorithm)
          
          ! Create the defect correction solver, attach VANKA as preconditioner.
          call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
              rsolver%RfilterChain)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

        case (6)
          ! BiCGStab with general VANKA preconditioning.
          !
          ! Create VANKA and initialise it with the parameters from the DAT file.
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
          
          call parlst_getvalue_string (rparList, scoarseGridSolverSection, &
            "spreconditionerSection", sstring, "", bdequote=.true.)
          read (sstring,*) spreconditionerSection
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rpreconditioner,rparList,&
              spreconditionerSection,p_rpreconditioner%calgorithm)
          
          ! Create the defect correction solver, attach VANKA as preconditioner.
          call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
              rsolver%RfilterChain)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,rparList,&
              scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

        case default
        
          call output_line ("Unknown coarse grid solver.", &
              OU_CLASS_ERROR,OU_MODE_STD,"lssh_initSolver")
          call sys_halt()
            
        end select
        
        ! Save the reference to the coarse grid solver.
        rsolver%p_rcoarseGridSolver => p_rlevelInfo%p_rcoarseGridSolver
        
        ! ---------------------------------------------------
        ! Smoother
        ! ---------------------------------------------------
        ! Now after the coarse grid solver is done, we turn to the smoothers
        ! on all levels. Their initialisation is similar to the coarse grid
        ! solver. Note that we use the same smoother on all levels, for
        ! presmoothing as well as for postsmoothing.
        
        do ilev = 2,nlevels

          ! Initialise the smoothers.
          select case (ismootherType)
          
          case (0:10)

            nullify(p_rsmoother)
          
            ! This is some kind of VANKA smoother. Initialise the correct one.
            select case (ismootherType)
            case (0)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
            case (1)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
            case (2)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_NAVST2D_DIAG)
            case (3)
              call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_NAVST2D_FULL)
            case (4)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_NAVST2D_DIAG)
              call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                  rsolver%RfilterChain)
            case (5)
              call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_NAVST2D_FULL)
              call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                  rsolver%RfilterChain)
            end select
            
            ! Initialise the parameters -- if there are any.
            call linsolinit_initParams (p_rsmoother,rparList,&
                ssmootherSection,LINSOL_ALG_UNDEFINED)
            call linsolinit_initParams (p_rsmoother,rparList,&
                ssmootherSection,p_rsmoother%calgorithm)
            
            ! Convert to a smoother with a defined number of smoothing steps.
            call parlst_getvalue_int (rparList, ssmootherSection, &
                      "nsmoothingSteps", nsm, 4)
            call linsol_convertToSmoother (p_rsmoother,nsm)
            
            ! Put the smoother into the level info structure as presmoother
            ! and postsmoother
            call linsol_getMultigrid2Level (p_rsolverNode,ilev,p_rlevelInfo)
            p_rlevelInfo%p_rpresmoother => p_rsmoother
            p_rlevelInfo%p_rpostsmoother => p_rsmoother
            
            ! Set up the interlevel projection structure for the projection from/to
            ! the lower level.
            call linsol_initProjMultigrid2Level(p_rlevelInfo,&
                rprjHierarchy%p_Rprojection(ilev))
            
          case default
          
            call output_line ("Unknown smoother.", &
                OU_CLASS_ERROR,OU_MODE_STD,"lssh_initSolver")
            call sys_halt()
            
          end select
      
        end do

      ! ---------------------------------------------------
      ! Heat equation
      ! ---------------------------------------------------
      case (CCEQ_HEAT2D)

        ! ---------------------------------------------------
        ! Coarse grid solver
        ! ---------------------------------------------------
        ! Ok, now we have to initialise all levels. First, we create a coarse
        ! grid solver and configure it.
        call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
        
        select case (icoarseGridSolverType)
        case (0)
          ! UMFPACK coarse grid solver. Easy.
          call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
          p_rlevelInfo%p_rcoarseGridSolver%p_rsubnodeUmfpack4%imatrixDebugOutput = &
              rdebugFlags%cwriteUmfpackMatrix
          p_rlevelInfo%p_rcoarseGridSolver%p_rsubnodeUmfpack4%smatrixName = "matrix.txt"
          
        case default
        
          call output_line ("Unknown coarse grid solver.", &
              OU_CLASS_ERROR,OU_MODE_STD,"lssh_initSolver")
          call sys_halt()
            
        end select
        
        ! Save the reference to the coarse grid solver.
        rsolver%p_rcoarseGridSolver => p_rlevelInfo%p_rcoarseGridSolver
        
        ! ---------------------------------------------------
        ! Smoother
        ! ---------------------------------------------------
        ! Now after the coarse grid solver is done, we turn to the smoothers
        ! on all levels. Their initialisation is similar to the coarse grid
        ! solver. Note that we use the same smoother on all levels, for
        ! presmoothing as well as for postsmoothing.
        
        do ilev = 2,nlevels

          ! Initialise the smoothers.
          select case (ismootherType)
          
          case (0:10)

            nullify(p_rsmoother)
          
            ! This is some kind of VANKA smoother. Initialise the correct one.
            select case (ismootherType)
            case (0)
              call linsol_initJacobi (p_rsmoother)
            case (1)
              call linsol_initSOR (p_rsmoother)
            end select
            
            ! Initialise the parameters -- if there are any.
            call linsolinit_initParams (p_rsmoother,rparList,&
                ssmootherSection,LINSOL_ALG_UNDEFINED)
            call linsolinit_initParams (p_rsmoother,rparList,&
                ssmootherSection,p_rsmoother%calgorithm)
            
            ! Convert to a smoother with a defined number of smoothing steps.
            call parlst_getvalue_int (rparList, ssmootherSection, &
                      "nsmoothingSteps", nsm, 4)
            call linsol_convertToSmoother (p_rsmoother,nsm)
            
            ! Put the smoother into the level info structure as presmoother
            ! and postsmoother
            call linsol_getMultigrid2Level (p_rsolverNode,ilev,p_rlevelInfo)
            p_rlevelInfo%p_rpresmoother => p_rsmoother
            p_rlevelInfo%p_rpostsmoother => p_rsmoother
            
            ! Set up the interlevel projection structure for the projection from/to
            ! the lower level.
            call linsol_initProjMultigrid2Level(p_rlevelInfo,&
                rprjHierarchy%p_Rprojection(ilev))
            
          case default
          
            call output_line ("Unknown smoother.", &
                OU_CLASS_ERROR,OU_MODE_STD,"lssh_initSolver")
            call sys_halt()
            
          end select
      
        end do
      
      end select

    end select

    ! Put the final solver node to the solver structure.
    rsolver%p_rsolverNode => p_rsolverNode

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_doneSolver (rsolver)
  
!<description>
  ! Cleans up a linear solver according to the settings in a parameter list.
!</description>

!<inputoutput>
  ! Solver node to initialise
  type(t_linsolSpace), intent(inout) :: rsolver
!</inputoutput>
  
!</subroutine>

    type(t_linsolSpace) :: rtemplate

    ! Release the solver
    call linsol_releaseSolver(rsolver%p_rsolverNode)
    
    ! Overwrite with default settings
    rsolver = rtemplate

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_createLinsolHierarchy (&
      rlsshierarchy,rfespaceHierarchy,rprjHierarchy,&
      nlmin,nlmax,cequation,rparList,ssection,rdebugFlags)
  
!<description>
  ! Based on a FE space hierarchy and a parameter list, this subroutine
  ! creates a linear solver hierarchy.
!</description>

!<input>
  ! Underlying hierarchy of FE spaces.
  type(t_feHierarchy), target :: rfeSpaceHierarchy
  
  ! An interlevel projection hierarchy, for MG based solvers
  type(t_interlevelProjectionHier), intent(in), target :: rprjHierarchy

  ! Minimum level in the hierarchy rfeSpaceHierarchy.
  ! Standard = 1. =0: nlmin=nlmax
  integer, intent(in) :: nlmin

  ! Maximum level in the hierarchy. <=0: use level MAX+nlmax
  integer, intent(in) :: nlmax

  ! Underlying equation. One of the LSS_EQN_xxxx constants.
  integer, intent(in) :: cequation

  ! Parameter list with solver parameters
  type(t_parlist), intent(in) :: rparList
  
  ! Section configuring the solver
  character(len=*), intent(in) :: ssection

  ! Debug flags
  type(t_optcDebugFlags), intent(in), target :: rdebugFlags
!</input>

!<output>
  ! The hierarchy to be created
  type(t_linsolHierarchySpace), intent(out) :: rlssHierarchy
!</output>
  
!</subroutine>

    ! local variables
    integer :: i
  
    ! Save some structures
    rlssHierarchy%p_rdebugFlags => rdebugFlags
    rlssHierarchy%p_feSpaceHierarchy => rfeSpaceHierarchy
    rlssHierarchy%nlmin = nlmin
    rlssHierarchy%nlmax = nlmax
    
    if (rlssHierarchy%nlmax .le. 0) then
      rlssHierarchy%nlmax = rfeSpaceHierarchy%nlevels + rlssHierarchy%nlmax
    end if

    if (rlssHierarchy%nlmin .le. 0) then
      rlssHierarchy%nlmin = rfeSpaceHierarchy%nlevels + rlssHierarchy%nlmin
    end if
    
    ! Create an array of solvers
    allocate(rlssHierarchy%p_RlinearSolvers(rlssHierarchy%nlmin:rlssHierarchy%nlmax))
    
    ! Supported equation
    rlsshierarchy%cequation = cequation
    
    ! Projection hierarchy
    rlsshierarchy%p_rprjHierarchy => rprjHierarchy
    
    ! -----------------------------------------------------
    ! Solver specific initialisation
    ! -----------------------------------------------------
    ! Initialise the solvers
    do i=rlssHierarchy%nlmin,rlssHierarchy%nlmax
      call lssh_initSolver (&
          rlssHierarchy%p_RlinearSolvers(i),cequation,i,&
          rprjHierarchy,rparList,ssection,rdebugFlags)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_releaseLinsolHierarchy (rlsshierarchy)

!<description>
  ! Releases a given linear solver hierarchy.
!</description>
  
!<inputoutput>
  ! The hierarchy to release.
  type(t_linsolHierarchySpace), intent(inout) :: rlssHierarchy
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i

    ! Release the local linear solvers
    do i=rlssHierarchy%nlmin,rlssHierarchy%nlmax
      call lssh_doneSolver(rlssHierarchy%p_RlinearSolvers(i))
    end do

    ! DE-associate pointers, clenup
    deallocate(rlssHierarchy%p_RlinearSolvers)
    nullify(rlssHierarchy%p_feSpaceHierarchy)
    nullify(rlssHierarchy%p_rdebugFlags)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_setMatrix (rlsshierarchy,ilevel,rmatrix)

!<description>
  ! Defines the system matrix for level ilevel.
!</description>
  
!<input>
  ! Level of the hierarchy
  integer, intent(in) :: ilevel
  
  ! Matrix to be used for solving at that level.
  type(t_matrixBlock), intent(in), target :: rmatrix
!</input>
  
!<inputoutput>
  ! The underlying hierarchy.
  type(t_linsolHierarchySpace), intent(inout) :: rlssHierarchy
!</inputoutput>
  
!</subroutine>

    ! Remember that matrix
    rlssHierarchy%p_RlinearSolvers(ilevel)%p_rmatrix => rmatrix
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_initStructure (rlsshierarchy,ilevel,ierror)

!<description>
  ! Initialises structural data for the solver at level ilevel.
!</description>
  
!<input>
  ! Level of the hierarchy
  integer, intent(in) :: ilevel
!</input>
  
!<inputoutput>
  ! The underlying hierarchy.
  type(t_linsolHierarchySpace), intent(inout), target :: rlssHierarchy
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!</subroutine>
    
    ! local variables
    integer :: i
    type(t_matrixBlock), dimension(:), allocatable :: Rmatrices

    ! Pass the matrices to the solver.
    allocate(Rmatrices(1:ilevel))
    do i=1,ilevel
      call lsysbl_duplicateMatrix (rlsshierarchy%p_RlinearSolvers(i)%p_rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices (rlsshierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode,Rmatrices)

    do i=1,ilevel
      call lsysbl_releaseMatrix(Rmatrices(i))
    end do
    deallocate(Rmatrices)

    ! Initialise the solver node
    call linsol_initStructure (rlsshierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode,ierror)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_initData (rlsshierarchy,roptcBDCSpaceHierarchy,ilevel,ierror)

!<description>
  ! Initialises calculation data for the solver at level ilevel.
!</description>
  
!<input>
  ! Level of the hierarchy
  integer, intent(in) :: ilevel

  ! Boundary condition hierarchy with precomputed boundary conditions.
  type(t_optcBDCSpaceHierarchy), intent(in), target :: roptcBDCSpaceHierarchy
!</input>
  
!<inputoutput>
  ! The underlying hierarchy.
  type(t_linsolHierarchySpace), intent(inout) :: rlssHierarchy
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>
  
!</subroutine>

    ! lcoal variables
    integer :: ifilter
    type(t_linsolSpace), pointer :: p_rlinsolSpace
    type(t_optcBDCSpace), pointer :: p_roptcBDCSpace
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ! Get the solver structure of that level
    p_rlinsolSpace => rlssHierarchy%p_RlinearSolvers(ilevel)
    
    ! Get the corresponding boundary condition structure
    p_roptcBDCSpace => roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilevel)
    
    ! Initialise the filter chain according to cflags.
    ! This is partially equation dependent.
    
    ifilter = 0
    p_rlinsolSpace%RfilterChain(:)%ifilterType = FILTER_DONOTHING
    
    ! Filter for Dirichlet boundary conditions
    ifilter = ifilter + 1
    p_rlinsolSpace%RfilterChain(ifilter)%ifilterType = FILTER_DISCBCDEFREAL
    
    select case (rlsshierarchy%cequation)

    ! ---------------------------------------------------
    ! 2D Stokes / Navier-Stokes
    ! ---------------------------------------------------
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

      ! Pressure filter for iterative solvers
      select case (p_rlinsolSpace%isolverType)
      
      ! ---------------------------------------------------
      ! Multigrid
      ! ---------------------------------------------------
      case (1)
        
        ! Smoother or coarse grid solver?
        select case (ilevel)
        
        ! -------------------------------------------------
        ! Coarse grid solver
        ! -------------------------------------------------
        case (1)

          ! -----------------------------------------------
          ! Integral-mean-value-zero filter for pressure
          ! -----------------------------------------------
          if (p_roptcBDCSpace%rneumannBoundary%nregions .eq. 0) then
          
            ! Active if there is no Neumann boundary
          
            select case (p_rlinsolSpace%icoarseGridSolverType)
            
            ! -----------------------------------
            ! UMFPACK
            ! -----------------------------------
            case (0)
              call output_line (&
                  "UMFPACK coarse grid solver does not support filtering!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lssh_initData")
              call sys_halt()

            ! -----------------------------------
            ! Iterative solver
            ! -----------------------------------
            case default
              ifilter = ifilter + 1
              p_rlinsolSpace%RfilterChain(ifilter)%ifilterType = FILTER_TOL20
              p_rlinsolSpace%RfilterChain(ifilter)%itoL20component = 3
            end select
            
          end if
        
        ! -------------------------------------------------
        ! Smoother
        ! -------------------------------------------------
        case (2:)

          ! -----------------------------------------------
          ! Integral-mean-value-zero filter for pressure
          ! -----------------------------------------------
          if (p_roptcBDCSpace%rneumannBoundary%nregions .eq. 0) then
          
            ! Active if there is no Neumann boundary

            ifilter = ifilter + 1
            p_rlinsolSpace%RfilterChain(ifilter)%ifilterType = FILTER_TOL20
            p_rlinsolSpace%RfilterChain(ifilter)%itoL20component = 3
 
          end if
        
        end select ! ilevel
        
      end select ! Outer solver
      
      ! -----------------------------------
      ! DEBUG Flags
      ! -----------------------------------
      select case (p_rlinsolSpace%icoarseGridSolverType)
      
      ! -----------------------------------
      ! UMFPACK
      ! -----------------------------------
      case (0)
      
        ! Matrix name for debug output: Matrix to text file.
        call linsol_getMultigrid2Level (&
            rlssHierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode,1,p_rlevelInfo)
        p_rlevelInfo%p_rcoarseGridSolver%p_rsubnodeUmfpack4%smatrixName = &
            "matrix"//trim(rlssHierarchy%p_rdebugFlags%sstringTag)
        
      end select
    
    ! ---------------------------------------------------
    ! Heat equation
    ! ---------------------------------------------------
    case (CCEQ_HEAT2D)

      ! Pressure filter for iterative solvers
      select case (p_rlinsolSpace%isolverType)
      
      ! ---------------------------------------------------
      ! Multigrid
      ! ---------------------------------------------------
      case (1)
        
        ! Smoother or coarse grid solver?
        select case (ilevel)
        
        ! -------------------------------------------------
        ! Coarse grid solver
        ! -------------------------------------------------
        case (1)

          ! -----------------------------------------------
          ! Integral-mean-value-zero filter for solution
          ! -----------------------------------------------
          if (p_roptcBDCSpace%rdirichletBoundary%nregions .eq. 0) then
          
            ! Active if there is no Dirichlet boundary
          
            select case (p_rlinsolSpace%icoarseGridSolverType)
            
            ! -----------------------------------
            ! UMFPACK
            ! -----------------------------------
            case (0)
              call output_line (&
                  "UMFPACK coarse grid solver does not support filtering!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lssh_initData")
              call sys_halt()

            ! -----------------------------------
            ! Iterative solver
            ! -----------------------------------
            case default
              ifilter = ifilter + 1
              p_rlinsolSpace%RfilterChain(ifilter)%ifilterType = FILTER_TOL20
              p_rlinsolSpace%RfilterChain(ifilter)%itoL20component = 1
            end select
            
          end if
        
        ! -------------------------------------------------
        ! Smoother
        ! -------------------------------------------------
        case (2:)

          ! -----------------------------------------------
          ! Integral-mean-value-zero filter for the solution
          ! -----------------------------------------------
          if (p_roptcBDCSpace%rdirichletBoundary%nregions .eq. 0) then
          
            ! Active if there is no Dirichlet boundary

            ifilter = ifilter + 1
            p_rlinsolSpace%RfilterChain(ifilter)%ifilterType = FILTER_TOL20
            p_rlinsolSpace%RfilterChain(ifilter)%itoL20component = 1
 
          end if
        
        end select ! ilevel
        
      end select ! Outer solver
                
      ! -----------------------------------
      ! DEBUG Flags
      ! -----------------------------------
      select case (p_rlinsolSpace%icoarseGridSolverType)
      
      ! -----------------------------------
      ! UMFPACK
      ! -----------------------------------
      case (0)
      
        ! Matrix name for debug output: Matrix to text file.
        call linsol_getMultigrid2Level (&
            rlssHierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode,1,p_rlevelInfo)
        p_rlevelInfo%p_rcoarseGridSolver%p_rsubnodeUmfpack4%smatrixName = &
            "matrix"//trim(rlssHierarchy%p_rdebugFlags%sstringTag)
        
      end select
    
    end select ! equation

    ! Initialise the solver node
    call linsol_initData (&
        rlssHierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode,ierror)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_doneData (rlsshierarchy,ilevel)

!<description>
  ! Cleans up calculation data for the solver at level ilevel.
!</description>
  
!<input>
  ! Level of the hierarchy
  integer, intent(in) :: ilevel
!</input>
  
!<inputoutput>
  ! The underlying hierarchy.
  type(t_linsolHierarchySpace), intent(inout) :: rlssHierarchy
!</inputoutput>
  
!</subroutine>

    ! Initialise the solver node
    call linsol_doneData (&
        rlssHierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_doneStructure (rlsshierarchy,ilevel)

!<description>
  ! Cleans up structural data for the solver at level ilevel.
!</description>
  
!<input>
  ! Level of the hierarchy
  integer, intent(in) :: ilevel
!</input>
  
!<inputoutput>
  ! The underlying hierarchy.
  type(t_linsolHierarchySpace), intent(inout) :: rlssHierarchy
!</inputoutput>
  
!</subroutine>

    ! Initialise the solver node
    call linsol_doneStructure (&
        rlssHierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lssh_precondDefect (rlsshierarchy,ilevel,rd,p_rsolverNode)

!<description>
  ! Applies preconditioning of rd with the solver on level ilevel.
!</description>
  
!<input>
  ! Level of the hierarchy
  integer, intent(in) :: ilevel
!</input>
  
!<inputoutput>
  ! The underlying hierarchy.
  type(t_linsolHierarchySpace), intent(inout) :: rlssHierarchy
  
  ! The defect to apply preconditioning to.
  type(t_vectorBlock), intent(inout) :: rd
  
  ! OPTIONAL; If present, this is set to the solver node of the linear
  ! solver used to solve the system
  type(t_linsolNode), pointer, optional :: p_rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Initialise the solver node
    call linsol_precondDefect (&
        rlssHierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode,rd)
        
    if (present(p_rsolverNode)) then
      p_rsolverNode => rlssHierarchy%p_RlinearSolvers(ilevel)%p_rsolverNode
    end if
    
  end subroutine

end module
