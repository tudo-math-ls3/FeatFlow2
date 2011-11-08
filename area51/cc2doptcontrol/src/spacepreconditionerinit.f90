!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2nonlinearcoreinit </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains initialisation routines for the core equation
!# (see also cc2dminim2nonlinearcore).
!# These routines connect the "problem" structure with the "core equation"
!# structure. In detail, we have routines that initialise the preconditioner
!# and all the information structures that are used during the nonlinear
!# iteration.
!#
!# The following routines can be found here:
!#
!# 1.) cc_allocPrecSystemMatrix
!#     -> Allocates memory for the system matrix representing the
!#        core equation.
!#
!# 2.) cc_initSpacePreconditioner
!#     -> Initialises a spatial preconditioner structure with parameters from
!#        the DAT file. This is needed for preconditioning in space.
!#     -> Extension to cc_createSpacePreconditioner.
!#
!# 3.) cc_doneSpacePreconditioner
!#     -> Cleans up a spatial preconditioner structure initialised by
!#        cc_initSpacePreconditioner.
!#     -> Extension to cc_releaseSpacePreconditioner.
!#
!# 4.) cc_configPreconditioner
!#     -> Configures a preconditioner according to the actual situation
!#
!# 5.) cc_updatePreconditioner
!#     -> Updates the preconditioner if there was a change in the system matrices
!#
!# 7.) cc_getProlRest
!#     -> Auxiliary routine: Set up interlevel projection structure
!#        with information from INI/DAT files
!#
!# 8.) cc_updatePreconditionerBC
!#     -> Updates the boundary conditions in a preconditioner
!#
!# The module works in tight relationship to cc2dmediumm2nonlinearcore.
!# cc2dmediumm2nonlinearcodeinit provides the routines to initialise
!# preconditioner and important structures using the problem related
!# structure. This module cc2dmediumm2nonlinearcore on the other hand
!# contains only the 'main' worker routines that do the work of the
!# nonlinear iteration -- independent of the problem structure!
!#
!# Note that this module and the "nonlinearcore" module are the only modules
!# that 'know' the actual structure of the system matrix and how to link
!# it to the main problem! For the actual assembly of the matrices and defect
!# vectors, routines of the module spacematvecassembly are used.
!#
!# </purpose>
!##############################################################################

module spacepreconditionerinit

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use filtersupport
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use linearsolver
  
  use collection
  use convection
    
  use basicstructures
  use spacepreconditioner
  use spacetimelinearsystem
  
  implicit none
  

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_allocPrecSystemMatrix (rproblem,rprecSpecials,&
      ilev,nlmin,nlmax,rlevelInfo,cmatrixType,rmatrix)
  
!<description>
  ! Allocates memory for the system matrix in a preconditioner. rlevelInfo
  ! provides information about the level where the system matrix should be created.
  !
  ! Before this routine is called, the structure of all matrices in
  ! rlevelInfo must have been set up!
  !
  ! Memory for A33 is allocated, but the submatrix is switched off by
  ! the multiplication factor.
  ! Memory for A12 and A21 is only allocated if necessary.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem

  ! Current assembly level.
  integer, intent(in) :: ilev
  
  ! Minimum assembly level.
  integer, intent(in) :: nlmin
  
  ! Maximum assembly level.
  integer, intent(in) :: nlmax

  ! A level-info structure specifying the matrices of the problem.
  type(t_problem_lvl), intent(in), target :: rlevelInfo
  
  ! Type of matrix.
  ! =CCMASM_MTP_AUTOMATIC: standard matrix, A11=A22
  ! =CCMASM_MTP_DECOUPLED: Decoupled velocity matrices A11 and A22
  ! =CCMASM_MTP_FULLTENSOR: Full-tensor matrix with A11,A12,A21,A22 independent
  integer, intent(in) :: cmatrixType

  ! A t_ccPreconditionerSpecials structure that defines special parameters
  ! of the preconditioner. The choice of the preconditioner may influence
  ! the way the matrix must be set up.
  type(t_ccPreconditionerSpecials), intent(in) :: rprecSpecials
!</input>

!<output>
  ! A block matrix that receives the basic system matrix.
  type(t_matrixBlock), intent(out) :: rmatrix
!</output>

    ! local variables
  
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A t_nonlinearSpatialMatrix used for defining the matrix structure
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%rdiscretisation
    
    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
      
    ! Let's consider the global system in detail. It has roughly
    ! the following shape:
    !
    !    ( A11       B1  ) = ( A11  A12  A13 )
    !    (      A22  B2  )   ( A21  A22  A23 )
    !    ( B1^T B2^T .   )   ( A31  A32  A33 )
    !
    ! All matices may have multiplication factors in their front.
    !
    ! The structure of the matrices A11 and A22 of the global system matrix
    ! is governed by the template FEM matrix.
    ! Initialise them with the same structure, i.e. A11, A22 share (!) their
    ! structure (not the entries) with that of the template matrix.
    !
    ! We allocate the system matrix using the smva_assembleMatrix routine.
    ! For this purpose, we have to initialise a t_nonlinearSpatialMatrix structure
    ! which defines the shape of the matrix. We simply set the parameters
    ! of those terms which wshould appear in the matrix to a value <> 0,
    ! that's enough for the memory allocation.

    call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
        rlevelInfo%rdiscretisation,rlevelInfo%rstaticInfo)
    
    call cc_preparePrecondMatrixAssembly (rnonlinearSpatialMatrix,&
        ilev,nlmin,nlmax,rprecSpecials)
    
    ! Get a dummy structure for a full matrix.
    call cc_getFullMatrixDummy(rproblem%rphysicsPrimal,rnonlinearSpatialMatrix)
    
    ! As matrix flag we specify 0 here. This allocates a basic matrix which is
    ! modified later for our needs.
    call smva_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_AUTOMATIC,&
        rmatrix,rnonlinearSpatialMatrix)
                                  
    ! That's it, all submatrices are set up.
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initLinearSolver (rproblem,rpreconditioner,ssection)
  
!<description>
  ! This routine initialises a linear solver structure to be used
  ! for preconditioning of the nonlinear defect.
!</description>

!<input>
  ! Name of the section in the parameter list that contains the parameters
  ! of the linear solver.
  character(LEN=*), intent(IN) :: ssection
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! A preconditioner structure where to write data about the linear
  ! solver to.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    type(t_parlist), pointer :: p_rparamList
    integer :: nlevels, ilev, nsm
    
    integer :: isolverType,ismootherType,icoarseGridSolverType
    character(LEN=SYS_STRLEN) :: sstring,ssolverSection,ssmootherSection
    character(LEN=SYS_STRLEN) :: scoarseGridSolverSection,spreconditionerSection
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    type(t_linsolNode), pointer :: p_rpreconditioner, p_rsmoother
    type(t_linsolNode), pointer :: p_rsolverNode

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    p_rparamList => rproblem%rparamList
        
    call parlst_querysection(p_rparamList, ssection, p_rsection)
    
    if (.not. associated(p_rsection)) then
      call output_line ('Cannot create linear solver; no section '''//trim(ssection)//&
                        '''!', OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
      call sys_halt()
    end if
    
    ! Get the parameters that configure the solver type
    
    call parlst_getvalue_int (p_rsection, 'isolverType', isolverType, 1)
    call parlst_getvalue_int (p_rsection, 'ismootherType', ismootherType, 3)
    call parlst_getvalue_int (p_rsection, 'icoarseGridSolverType', &
        icoarseGridSolverType, 1)
        
    rpreconditioner%rprecSpecials%isolverType = isolverType
    rpreconditioner%rprecSpecials%ismootherType = ismootherType
    rpreconditioner%rprecSpecials%icoarseGridSolverType = icoarseGridSolverType

    call parlst_getvalue_string (p_rsection, 'ssolverSection', sstring,'')
    read (sstring,*) ssolverSection
    call parlst_getvalue_string (p_rsection, 'ssmootherSection', sstring,'')
    read (sstring,*) ssmootherSection
    call parlst_getvalue_string (p_rsection, 'scoarseGridSolverSection', sstring,'')
    read (sstring,*) scoarseGridSolverSection
    
    ! Which type of solver do we have?
    
    select case (isolverType)
    
    case (0)
    
      ! This is the UMFPACK solver. Very easy to initialise. No parameters at all.
      call linsol_initUMFPACK4 (p_rsolverNode)
    
    case (1)
    
      ! Multigrid solver. This is a little bit harder.
      !
      ! In a first step, initialise the main solver node for all our levels.
      nlevels = rpreconditioner%NLMAX - rpreconditioner%NLMIN + 1
      
      call linsol_initMultigrid2 (p_rsolverNode,nlevels,&
          rpreconditioner%RfilterChain)
      
      ! Manually trim the coarse grid correction in Multigrid to multiply the
      ! pressure equation with -1. This (un)symmetrises the operator and gives
      ! much better convergence rates.
      call cgcor_release(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
      call cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,6)
      p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(3) &
          = -1.0_DP
      p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(6) &
          = -1.0_DP

      ! Init standard solver parameters and extended multigrid parameters
      ! from the DAT file.
      call linsolinit_initParams (p_rsolverNode,p_rparamList,ssolverSection,&
          LINSOL_ALG_UNDEFINED)
      call linsolinit_initParams (p_rsolverNode,p_rparamList,ssolverSection,&
          LINSOL_ALG_MULTIGRID2)
      
      ! Ok, now we have to initialise all levels. First, we create a coarse
      ! grid solver and configure it.
      call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
      
      select case (icoarseGridSolverType)
      case (0)
        ! UMFPACK coarse grid solver. Easy.
        call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
        
      case (1)
        ! Defect correction with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
        
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%rprecSpecials%bneedVirtTransposedDonCoarse = .true.

      case (2)
        ! Defect correction with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
        
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%rprecSpecials%bneedVirtTransposedDonCoarse = .true.

      case (3)
        ! BiCGStab with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)
        
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%rprecSpecials%bneedVirtTransposedDonCoarse = .true.

      case (4)
        ! BiCGStab with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rpreconditioner%rprecSpecials%bneedVirtTransposedDonCoarse = .true.

      case (5)
        ! BiCGStab with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

      case (6)
        ! BiCGStab with diagonal VANKA preconditioning, new implementation
        ! for general elements
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

      case default
      
        call output_line ('Unknown coarse grid solver.', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
        call sys_halt()
          
      end select
      
      ! Now after the coarse grid solver is done, we turn to the smoothers
      ! on all levels. Their initialisation is similar to the coarse grid
      ! solver. Note that we use the same smoother on all levels, for
      ! presmoothing as well as for postsmoothing.
      
      do ilev = 2,nlevels

        ! Initialise the smoothers.
        select case (ismootherType)
        
        case (0:9)

          nullify(p_rsmoother)
        
          ! This is some kind of VANKA smoother. Initialise the correct one.
          select case (ismootherType)
          case (0)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
          case (1)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
          case (2)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

          case (3)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIRECT)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

          case (4)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

          case (5)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAGDIR)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

          case (6)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

          case (7)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

          case (8)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%RfilterChain)

          case (9)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)

          end select
          
          ! Initialise the parameters -- if there are any.
          call linsolinit_initParams (p_rsmoother,p_rparamList,&
              ssmootherSection,LINSOL_ALG_UNDEFINED)
          call linsolinit_initParams (p_rsmoother,p_rparamList,&
              ssmootherSection,p_rsmoother%calgorithm)
          
          ! Convert to a smoother with a defined number of smoothing steps.
          call parlst_getvalue_int (p_rparamList, ssmootherSection, &
                    'nsmoothingSteps', nsm, 4)
          call linsol_convertToSmoother (p_rsmoother,nsm)
          
          ! Put the smoother into the level info structure as presmoother
          ! and postsmoother
          call linsol_getMultigrid2Level (p_rsolverNode,ilev,p_rlevelInfo)
          p_rlevelInfo%p_rpresmoother => p_rsmoother
          p_rlevelInfo%p_rpostsmoother => p_rsmoother
          
          ! Set up the interlevel projection structure for the projection from/to
          ! the lower level.
          call linsol_initProjMultigrid2Level(p_rlevelInfo,&
              rpreconditioner%p_Rprojection(rpreconditioner%NLMIN-1+ilev))
          
        case default
        
          call output_line ('Unknown smoother.', &
              OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
          call sys_halt()
          
        end select
      
      end do

      ! Get information about adaptive matrix generation from INI/DAT files
      call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
          'iAdaptiveMatrix', rpreconditioner%rprecSpecials%iadaptiveMatrices, 0)
                                
      call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
          'dAdMatThreshold', rpreconditioner%rprecSpecials%dAdMatThreshold, 20.0_DP)

    case (2:3)
    
      ! VANKA smoother: 1 step defect correction with nmaxIterations steps VANKA.
      ! ismootherType defines the type of smoother to use.
      select case (ismootherType)
      
      case (0:9)

        nullify(p_rsmoother)
      
        ! This is some kind of VANKA smoother. Initialise the correct one.
        select case (ismootherType)
        case (0)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
        case (1)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
        case (2)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

        case (3)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIRECT)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

        case (4)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

        case (5)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAGDIR)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

        case (6)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%RfilterChain)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

        case (7)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%RfilterChain)

          ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
          rpreconditioner%rprecSpecials%bneedVirtTransposedD = .true.

        case (8)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%RfilterChain)

        case (9)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG2)

        end select
        
        ! Initialise the parameters -- if there are any.
        call linsolinit_initParams (p_rsmoother,p_rparamList,&
            ssmootherSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rsmoother,p_rparamList,&
            ssmootherSection,p_rsmoother%calgorithm)
        
        ! Convert to a smoother with a defined number of smoothing steps.
        call parlst_getvalue_int (p_rparamList, ssmootherSection, &
                  'nsmoothingSteps', nsm, 4)
        call linsol_convertToSmoother (p_rsmoother,nsm)
        
      case DEFAULT
      
        call output_line ('Unknown smoother.', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
        call sys_halt()
        
      end select

      ! Init defect correction, 1-step... with that smoother
      if (isolverType .eq. 2) then
        call linsol_initDefCorr (p_rsolverNode,p_rsmoother,rpreconditioner%RfilterChain)
      else
        call linsol_initBiCGStab (p_rsolverNode,p_rsmoother,rpreconditioner%RfilterChain)
      end if
      call linsolinit_initParams (p_rsolverNode,p_rparamList,&
          ssolverSection,LINSOL_ALG_UNDEFINED)
      call linsolinit_initParams (p_rsolverNode,p_rparamList,&
          ssolverSection,p_rsolverNode%calgorithm)
    
    end select

    ! Put the final solver node to the preconditioner structure.
    rpreconditioner%p_rsolverNode => p_rsolverNode

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneLinearSolver (rpreconditioner)
  
!<description>
  ! Releases information of the linear solver preconditioner from the
  ! structure rnonlinearIteration. Cleans up what was configured
  ! in cc_initLinearSolver.
!</description>

!<inputoutput>
  ! A preconditioner structure where to remove data about the linear
  ! solver from.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! Release the solver and all subsolvers.
    call linsol_releaseSolver(rpreconditioner%p_rsolverNode)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initSpacePreconditioner (rproblem,nlmin,nlmax,rpreconditioner)
  
!<description>
  ! Initialises the given spatial preconditioner structure rpreconditioner.
  ! Creates the structure with cc_createSpacePreconditioner and saves all
  ! problem dependent parameters and matrices in it.
  !
  ! This routine initialises only the basic structure. However, it does
  ! not set/initialise the type of preconditioner (Defect corection,
  ! Newton,...).
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! Minimum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners.
  integer, intent(IN) :: nlmin
  
  ! Maximum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioner. This is the level where the preconditioner is to be
  ! applied!
  integer, intent(IN) :: nlmax
!</input>

!<output>
  ! A spatial preconditioner structure to be initialised.
  type(t_ccspatialPreconditioner), intent(OUT) :: rpreconditioner
!</output>

!</subroutine>

    ! local variables
    integer :: ilevel

    ! Basic initialisation of the nonlinenar iteration structure.
    call cc_createSpacePreconditioner (rpreconditioner,nlmin,nlmax)
    
    ! Assign the matrix pointers in the nonlinear iteration structure to
    ! all our matrices that we want to use.
    do ilevel = nlmin,nlmax
    
      rpreconditioner%RcoreEquation(ilevel)%p_rstaticInfo => &
        rproblem%RlevelInfo(ilevel)%rstaticInfo

      rpreconditioner%RcoreEquation(ilevel)%p_rdiscrBlock => &
        rproblem%RlevelInfo(ilevel)%rdiscretisation

      ! Create some temp vectors.
      call lsysbl_createVectorBlock(&
        rpreconditioner%RcoreEquation(ilevel)%p_rdiscrBlock,&
        rpreconditioner%RcoreEquation(ilevel)%rtempVector1,.true.)

      call lsysbl_createVectorBlock(&
        rpreconditioner%RcoreEquation(ilevel)%p_rdiscrBlock,&
        rpreconditioner%RcoreEquation(ilevel)%rtempVector2,.true.)

      call lsysbl_createVectorBlock(&
        rpreconditioner%RcoreEquation(ilevel)%p_rdiscrBlock,&
        rpreconditioner%RcoreEquation(ilevel)%rtempVector3,.true.)

      ! Initialise basic boundary conditions on all levels
      call bcasm_initDiscreteBC(rpreconditioner%RcoreEquation(ilevel)%rdiscreteBC)
      call bcasm_initDiscreteFBC(rpreconditioner%RcoreEquation(ilevel)%rdiscreteFBC)
      
    end do
      
    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions.
    !
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter for defect vectors:
    rpreconditioner%RfilterChain(1)%ifilterType = &
        FILTER_DISCBCDEFREAL

    ! The second filter filters for boundary conditions of fictitious boundary
    ! components
    rpreconditioner%RfilterChain(2)%ifilterType = &
        FILTER_DISCBCDEFFICT

    ! By default, no filter at position 3/4. Can be Neumann boundary filter.
    rpreconditioner%RfilterChain(3)%ifilterType = FILTER_DONOTHING
    rpreconditioner%RfilterChain(4)%ifilterType = FILTER_DONOTHING
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneSpacePreconditioner (rpreconditioner)
  
!<description>
  ! Releases memory allocated in cc_initPrecoditioner..
  ! The routine automatically calls cc_releaseSpacePreconditioner to release
  ! internal parameters.
!</description>

!<inputoutput>
  ! A spatial preconditioner structure to be cleaned up.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>
    
    ! local variables
    integer :: i

    ! Which preconditioner do we have?
    select case (rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioning
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_INEXACTNEWTON)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the preconditioner matrix on every level
      do i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
        call lsysbl_releaseMatrix ( &
          rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner)
        
        call lsysbl_releaseVector (rpreconditioner%RcoreEquation(i)%rtempVector3)
        call lsysbl_releaseVector (rpreconditioner%RcoreEquation(i)%rtempVector2)
        call lsysbl_releaseVector (rpreconditioner%RcoreEquation(i)%rtempVector1)
        
        ! Release boundary conditions on all levels
        call bcasm_releaseDiscreteBC(rpreconditioner%RcoreEquation(i)%rdiscreteBC)
        call bcasm_releaseDiscreteFBC(rpreconditioner%RcoreEquation(i)%rdiscreteFBC)
        
      end do
      
      ! Release the temporary vector(s)
      call lsyssc_releaseVector (rpreconditioner%rtempVectorSc)
      
      ! Clean up data about the projection etc.
      do i=rpreconditioner%NLMAX,rpreconditioner%NLMIN+1,-1
        call mlprj_doneProjection(rpreconditioner%p_Rprojection(i))
      end do
      deallocate(rpreconditioner%p_Rprojection)

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      call linsol_doneStructure (rpreconditioner%p_rsolverNode)
      call cc_doneLinearSolver (rpreconditioner)
      
    case default
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

    call cc_releaseSpacePreconditioner (rpreconditioner)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_updatePreconditionerBC (rproblem,rpreconditioner,dtime)
  
!<description>
  ! Updates the boundary conditions in the preconditioner coprresponding to time
  ! dtime.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! Simulation time.
  real(dp), intent(in) :: dtime
!</input>

!<inputoutput>
  ! A spatial preconditioner structure to be initialised.
  type(t_ccspatialPreconditioner), intent(inout) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilevel

    ! Initialise the collection for the assembly.
    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

    ! Update the BC's on every level.
    do ilevel = rpreconditioner%nlmin,rpreconditioner%nlmax
    
      call bcasm_clearDiscreteBC(rpreconditioner%RcoreEquation(ilevel)%rdiscreteBC)
      call bcasm_clearDiscreteFBC(rpreconditioner%RcoreEquation(ilevel)%rdiscreteFBC)
      
      call cc_assembleBDconditions (rproblem,dtime,&
          rpreconditioner%RcoreEquation(ilevel)%p_rdiscrBlock,&
          CCDISCBC_PRIMALDUAL,rpreconditioner%RcoreEquation(ilevel)%rdiscreteBC,&
          rproblem%rcollection,rpreconditioner%RcoreEquation(ilevel)%bhasNeumann)
      call cc_assembleFBDconditions (rproblem,dtime,&
          rpreconditioner%RcoreEquation(ilevel)%p_rdiscrBlock,&
          CCDISCBC_PRIMALDUAL,rpreconditioner%RcoreEquation(ilevel)%rdiscreteFBC,&
          rproblem%rcollection)
      
    end do
      
    ! Adjust preconditioner specials (Neumann filter etc.)
    call cc_adjustPrecSpecials (rproblem,rpreconditioner)
      
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_configPreconditioner (rproblem,rpreconditioner,ssection,&
      ctypePreconditioner)
  
!<description>
  ! This routine prepares the preconditioner by means of the parameters
  ! in the DAT files and on the information in rproblem.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! Name of the section containing the configuration of the preconditioner.
  ! If ctypePreconditioner=CCPREC_LINEARSOLVER or =CCPREC_NEWTON, this must
  ! specify the name of the section that configures a linear solver.
  character(LEN=*), intent(IN) :: ssection
  
  ! Type of the preconditioner.
  ! =1: standard linear equation.
  ! =2: Newton matrix
  integer, intent(IN) :: ctypePreconditioner
!</input>

!<inputoutput>
  ! A spatial preconditioner structure to be initialised. Must have been
  ! created previously with initPreconditioner.
  ! This is configured according to the preconditioner as specified in
  ! the DAT files.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
    integer :: imaxmem
    character(LEN=PARLST_MLDATA) :: ssolverName
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Note in the structure which preconditioner we use and which section in
    ! the DAT file contains its parameters.
    rpreconditioner%ctypePreconditioning = ctypePreconditioner
    rpreconditioner%spreconditionerSection = ssection

    select case (rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioner
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_INEXACTNEWTON)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rpreconditioner%NLMIN
      NLMAX = rpreconditioner%NLMAX
      
      ! Get a pointer to the discretsation structure on the level
      ! where the preconditioner should act
      p_rdiscretisation => rproblem%RlevelInfo(NLMAX)%rdiscretisation
      
      ! The preconditioner is a linear solver, so ssection is the name of the section
      ! configuring the linear solver.
      ssolverName = ssection
      
      ! Initialise a standard interlevel projection structure for every level
      allocate(rpreconditioner%p_rprojection(NLMIN+1:NLMAX))
      do i=NLMIN+1,NLMAX
        call mlprj_initProjectionDiscr (rpreconditioner%p_Rprojection(i),&
            rproblem%RlevelInfo(i)%rdiscretisation)
      
        ! Initialise the projection structure with data from the INI/DAT
        ! files. This allows to configure prolongation/restriction.
        call cc_getProlRest (rpreconditioner%p_Rprojection(i), &
            rproblem%rparamList,  'CC-PROLREST')
      end do
      
      ! Initialise the linear solver as configured in the DAT file.
      call cc_initLinearSolver (rproblem,rpreconditioner,ssolverName)

      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...

      imaxmem = 0
      do i=NLMIN+1,NLMAX
        ! Pass the system metrices on the coarse/fine grid to
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = max(imaxmem,mlprj_getTempMemoryDirect (&
            rpreconditioner%p_Rprojection(i),&
            rproblem%RlevelInfo(i-1)% &
              rdiscretisation%RspatialDiscr(1:p_rdiscretisation%ncomponents),&
            rproblem%RlevelInfo(i)% &
              rdiscretisation%RspatialDiscr(1:p_rdiscretisation%ncomponents)))
      end do
      
      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      call lsyssc_createVector (rpreconditioner%rtempVectorSc,&
        max(imaxmem,dof_igetNDofGlobBlock(p_rdiscretisation)),.false.)
      
      ! Initialise the matrices.
      call cc_updatePreconditioner (rproblem,rpreconditioner,.true.,.true.)

    case DEFAULT
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_updatePreconditioner (rproblem,rpreconditioner,&
      binit,bstructuralUpdate)
  
!<description>
  ! This routine has to be called whenever the system matrices change.
  ! It initialises (depending on the system matrices) the matrices of the
  ! preconditioner or performs an update of them.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! First initialisation.
  ! Has to be set to TRUE on the first call. Initialises the preconditioner
  ! with the structure of the matrices.
  logical, intent(IN) :: binit

  ! Whether the structure of the system matrices is new.
  ! This variable has to be set to TRUE whenever there was a structure in
  ! the system matrices. This reinitialises the linear solver.
  logical, intent(IN) :: bstructuralUpdate
!</input>

!<inputoutput>
  ! A spatial preconditioner structure to be úpdated. Must have been
  ! created previously with initPreconditioner and configPreconditioner.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
    logical :: bphystranspose

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A pointer to the matrix of the preconditioner
    type(t_matrixBlock), pointer :: p_rmatrixPreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(rpreconditioner%NLMAX) :: Rmatrices
    
    ! Pointer to the template FEM matrix
    type(t_matrixScalar), pointer :: p_rmatrixTempateFEM
    
    select case (rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioner
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_INEXACTNEWTON)

      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Adjust the special preconditioner parameters to the current situation.
      call cc_adjustPrecSpecials (rproblem,rpreconditioner)

      ! Which levels have we to take care of during the solution process?
      NLMIN = rpreconditioner%NLMIN
      NLMAX = rpreconditioner%NLMAX
      
      ! Initialise the preconditioner matrices on all levels.
      do i=NLMIN,NLMAX
      
        if (binit .or. bstructuralUpdate) then
      
          ! Prepare the preconditioner matrices level i. This is
          ! basically the system matrix.
          ! Clean up the structure if necessary.
          if (rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner%NEQ .ne. 0) then
            call lsysbl_releaseMatrix (rpreconditioner%RcoreEquation(i)%&
                rmatrixPreconditioner)
          end if
        
          ! Allocate memory for the basic submatrices.
          call cc_allocPrecSystemMatrix (rproblem,rpreconditioner%rprecSpecials,&
              i,nlmin,nlmax,rproblem%RlevelInfo(i),CCMASM_MTP_AUTOMATIC,&
              rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner)
              
          ! Attach boundary conditions
          call lsysbl_assignDiscreteBC (&
              rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner,&
              rpreconditioner%RcoreEquation(i)%rdiscreteBC)
          
          call lsysbl_assignDiscreteFBC (&
              rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner,&
              rpreconditioner%RcoreEquation(i)%rdiscreteFBC)
              
        end if
        
        ! On the current level, set up a global preconditioner matrix.
      
        p_rmatrixPreconditioner => &
            rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner
        
        if (binit .or. bstructuralUpdate) then

          ! ----------------------------------------------------
          ! Should the linear solver use the Newton matrix?
          if ((rpreconditioner%ctypePreconditioning .eq. CCPREC_NEWTON) .or. &
              (rpreconditioner%ctypePreconditioning .eq. CCPREC_INEXACTNEWTON)) then
            ! That means, our preconditioner matrix must look like
            !
            !  A11  A12  B1
            !  A21  A22  B2
            !  B1^T B2^T 0
            !
            ! With A12, A21, A11, A22 independent of each other!
            ! Do we have that case? If not, we have to allocate memory
            ! for these matrices.
            p_rmatrixTempateFEM => rproblem%RlevelInfo(i)%rstaticInfo%rmatrixTemplateFEM
            
            ! If we have a Stokes problem, A12 and A21 don't exist, so we have nothing to
            ! do. In the Navier-Stokes case however, we have A12 and A21, so create them!
            ! Furthermore, there is to be a A51 and A42 submatrix allocated for the
            ! reactive coupling mass matrix R!
            if (rproblem%rphysicsPrimal%iequation .eq. 0) then
            
              if (p_rmatrixPreconditioner%RmatrixBlock(1,2)%cmatrixFormat &
                  .eq. LSYSSC_MATRIXUNDEFINED) then
                  
                call lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                  p_rmatrixPreconditioner%RmatrixBlock(1,2), &
                  LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
                ! Allocate memory for the entries; don't initialise the memory.
                call lsyssc_allocEmptyMatrix (&
                    p_rmatrixPreconditioner%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
                  
              end if

              if (p_rmatrixPreconditioner%RmatrixBlock(2,1)%cmatrixFormat &
                  .eq. LSYSSC_MATRIXUNDEFINED) then
                  
                ! Create a new matrix A21 in memory. create a new matrix
                ! using the template FEM matrix...
                call lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                  p_rmatrixPreconditioner%RmatrixBlock(2,1), &
                  LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
                ! Allocate memory for the entries; don't initialise the memory.
                call lsyssc_allocEmptyMatrix (&
                    p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)

              end if

              if (p_rmatrixPreconditioner%RmatrixBlock(5,1)%cmatrixFormat &
                  .eq. LSYSSC_MATRIXUNDEFINED) then
                  
                call lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                  p_rmatrixPreconditioner%RmatrixBlock(5,1), &
                  LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
                ! Allocate memory for the entries; don't initialise the memory.
                call lsyssc_allocEmptyMatrix (&
                    p_rmatrixPreconditioner%RmatrixBlock(5,1),LSYSSC_SETM_UNDEFINED)
                  
              end if

              if (p_rmatrixPreconditioner%RmatrixBlock(4,2)%cmatrixFormat &
                  .eq. LSYSSC_MATRIXUNDEFINED) then
                  
                ! Create a new matrix A21 in memory. create a new matrix
                ! using the template FEM matrix...
                call lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                  p_rmatrixPreconditioner%RmatrixBlock(4,2), &
                  LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
                ! Allocate memory for the entries; don't initialise the memory.
                call lsyssc_allocEmptyMatrix (&
                    p_rmatrixPreconditioner%RmatrixBlock(4,2),LSYSSC_SETM_UNDEFINED)

              end if

            end if

            ! A22 may share its entries with A11. If that's the case,
            ! allocate additional memory for A22, as the Newton matrix
            ! requires a separate A22!
            if (lsyssc_isMatrixContentShared( &
                p_rmatrixPreconditioner%RmatrixBlock(1,1), &
                p_rmatrixPreconditioner%RmatrixBlock(2,2)) ) then
              ! Release the matrix structure. As the matrix is a copy
              ! of another one, this will clean up the structure but
              ! not release any memory.
              call lsyssc_releaseMatrix ( &
                  p_rmatrixPreconditioner%RmatrixBlock(2,2))

              ! Create a new matrix A21 in memory. create a new matrix
              ! using the template FEM matrix...
              call lsyssc_duplicateMatrix ( &
                p_rmatrixPreconditioner%RmatrixBlock(1,1), &
                p_rmatrixPreconditioner%RmatrixBlock(2,2),&
                LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
              ! ... then allocate memory for the entries;
              ! don't initialise the memory.
              call lsyssc_allocEmptyMatrix (&
                  p_rmatrixPreconditioner%RmatrixBlock(2,2),&
                  LSYSSC_SETM_UNDEFINED)
            end if
            
          end if

          ! We add a zero diagonal matrix to the pressure block. This matrix
          ! is only used under rare circumstances, e.g. if we have a pure
          ! Dirichlet problem on that level and the solver does not support
          ! filtering. In this case, this matrix is used to ensure definiteness.
          !
          ! We ignore the scaling factor here as we only want to ensure that there's
          ! space available for the matrix.
          if (.not. lsysbl_isSubmatrixPresent(p_rmatrixPreconditioner,3,3,.true.)) then
            call lsyssc_createDiagMatrixStruc (p_rmatrixPreconditioner%RmatrixBlock(3,3),&
                p_rmatrixPreconditioner%RmatrixBlock(1,3)%NCOLS,LSYSSC_MATRIX9)
            call lsyssc_allocEmptyMatrix(&
                p_rmatrixPreconditioner%RmatrixBlock(3,3),LSYSSC_SETM_ZERO)
            p_rmatrixPreconditioner%RmatrixBlock(3,3)%dscaleFactor = 0.0_DP
          end if
          
        end if
        
        ! Under certain circumstances, the linear solver needs B^T-matrices.
        ! This is the case if
        ! - a direct solver (UMFPACK) is used on a level or
        ! - if the general VANKA preconditioner is used.
        ! In these cases, we create a separate copy of B1 and B2 and transpose them.
        ! Note that we do this only in that case when there is a 'structural update'
        ! (which means that the structure of the matrices have changed). Otherwise,
        ! B1^T and B2^T stay unchanged!
        !
        ! In case we have a pure-dirichlet problem, we activate the 3,3-submatrix
        ! It's probably needed for the preconditioner to make the pressure definite.
        if (binit .or. bstructuralUpdate) then
        
          bphystranspose = .false.
          select case (rpreconditioner%rprecSpecials%isolverType)
          case (0)
            ! UMFPACK solver.
            if (i .eq. NLMAX) then
              ! UMFPACK solver. Tweak the matrix on the max. level.
              ! Ignore the other levels.
              bphystranspose = .true.
              
              if (rpreconditioner%rprecSpecials%bneedPressureDiagonalBlock) then
                ! Activate the 3,3-block, UMFPACK needs it in case the pressure
                ! is globally indefinite.
                p_rmatrixPreconditioner%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
              end if
            end if

          case (1)
            ! Multigrid solver. Treat the matrix at the coarse level if there's
            ! UMFPACK chosen as coarse grid solver.
            if (i .eq. NLMIN) then
            
              if (rpreconditioner%rprecSpecials%icoarseGridSolverType .eq. 0) then
                bphystranspose = .true.
                
                if (rpreconditioner%rprecSpecials%bneedPressureDiagonalBlock) then
                  ! Activate the 3,3-block, UMFPACK needs it in case the pressure
                  ! is globally indefinite.
                  p_rmatrixPreconditioner%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
                end if
                
              else
                
                ! Tweak the matrix if the preconditioner needs transposed matrices.
                !
                ! Currently not implemented, as there is no configuration where a
                ! preconditioner needs transposed matrices...
                
              end if
              
            else
            
              ! On the other levels, tweak the matrix if the general VANKA is
              ! chosen as smoother; it needs transposed matrices.
              if ((rpreconditioner%rprecSpecials%ismootherType .eq. 0) .or. &
                  (rpreconditioner%rprecSpecials%ismootherType .eq. 1)) then
                bphystranspose = .true.
              end if
              
            end if

          end select
          
        end if
      
      end do
        
      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
      do i=NLMIN,NLMAX
        call lsysbl_duplicateMatrix ( &
          rpreconditioner%RcoreEquation(i)%rmatrixPreconditioner, &
          Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end do
      
      call linsol_setMatrices(&
          rpreconditioner%p_rsolverNode,Rmatrices(NLMIN:NLMAX))
          
      ! The solver got the matrices; clean up Rmatrices, it was only of temporary
      ! nature...
      do i=NLMIN,NLMAX
        call lsysbl_releaseMatrix (Rmatrices(i))
      end do
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      if (binit) then
        call linsol_initStructure (rpreconditioner%p_rsolverNode,ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
      else if (bstructuralUpdate) then
        call linsol_updateStructure (rpreconditioner%p_rsolverNode,ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
      end if
      
    case DEFAULT
      
      ! Unknown preconditioner
      call output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_updatePreconditioner')
      call sys_halt()
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getProlRest (rprojection, rparamList, sname)
  
!<description>
  ! Initialises an existing interlevel projection structure rprojection
  ! with parameters from the INI/DAT files. sname is the section in the
  ! parameter list containing parameters about prolongation restriction.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(IN) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the prolongation/restriction.
  character(LEN=*), intent(IN) :: sname
!</input>

!<output>
  ! An interlevel projection block structure containing an initial
  ! configuration of prolongation/restriction. The structure is modified
  ! according to the parameters in the INI/DAT file(s).
  type(t_interlevelProjectionBlock), intent(INOUT) :: rprojection
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: i1
    real(DP) :: d1

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, sname, p_rsection)

    if (.not. associated(p_rsection)) then
      ! We use the default configuration; stop here.
      return
    end if
    
    ! Now take a look which parameters appear in that section.

    ! Prolongation/restriction order for velocity components
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
      
      rprojection%RscalarProjection(:,4:5)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,4:5)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,4:5)%iinterpolationOrder = i1
    end if

    ! Prolongation/restriction order for pressure
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for pressure components
      rprojection%RscalarProjection(:,NDIM2D+1)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%iinterpolationOrder = i1
      
      rprojection%RscalarProjection(:,6)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,6)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,6)%iinterpolationOrder = i1
    end if
    
    ! Prolongation/restriction variant for velocity components
    ! in case of Q1~ discretisation
    call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
    
    if (i1 .ne. -1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1

      rprojection%RscalarProjection(:,4:5)%iprolVariant  = i1
      rprojection%RscalarProjection(:,4:5)%irestVariant  = i1
    end if
    
    ! Aspect-ratio indicator in case of Q1~ discretisation
    ! with extended prolongation/restriction
    call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
    
    if (i1 .ne. 1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
      
      rprojection%RscalarProjection(:,4:5)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,4:5)%irestARIndicatorEX3Y  = i1
    end if

    ! Aspect-ratio bound for switching to constant prolongation/restriction
    ! in case of Q1~ discretisation with extended prolongation/restriction
    call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
    
    if (d1 .ne. 20.0_DP) then
      rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
      
      rprojection%RscalarProjection(:,4:5)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,4:5)%drestARboundEX3Y  = d1
    end if

  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine cc_adjustPrecSpecials (rproblem,rpreconditioner)
  
!<description>
  ! This routine adjusts parameters in the rprecSpecials according to the current
  ! problem. This includes e.g. special 'tweaks' that must be done if the
  ! problem is a pure Dirichlet problem.
  !
  ! The routine must always be called if the situation changes during a
  ! simulation (e.g. if a nonstationary simulation proceeds to a new timestep
  ! and changes boundary conditions). It is usually called in
  ! cc_updatePreconditioner.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!<inputoutput>
  ! Nonlinear iteration structure.
  ! The t_ccPreconditionerSpecials substructure that receives information how to
  ! finally assembly the matrices such that everything in the callback routines
  ! will work.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! Do we have Neumann boundary?
    ! The bneumann flag on the max. level decides upon that.
    !
    ! The bhasNeumannBoundary flag of the higher level decides about that...
    if (rpreconditioner%RcoreEquation(rpreconditioner%nlmax)%bhasNeumann) then
      rpreconditioner%RfilterChain(3)%ifilterType = FILTER_DONOTHING
      rpreconditioner%RfilterChain(4)%ifilterType = FILTER_DONOTHING
      
      rpreconditioner%rprecSpecials%bneedPressureDiagonalBlock = .false.
    else
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rpreconditioner%RfilterChain(3)%ifilterType = FILTER_TOL20
      rpreconditioner%RfilterChain(3)%itoL20component = NDIM2D+1

      rpreconditioner%RfilterChain(4)%ifilterType = FILTER_TOL20
      rpreconditioner%RfilterChain(4)%itoL20component = 2*(NDIM2D+1)
      
      ! Matrices may have to be changed, depending on the solver.
      rpreconditioner%rprecSpecials%bneedPressureDiagonalBlock = .true.
    end if

  end subroutine

end module
