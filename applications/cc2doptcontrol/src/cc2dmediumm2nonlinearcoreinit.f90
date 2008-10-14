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
!# 1.) cc_allocSystemMatrix
!#     -> Allocates memory for the system matrix representing the 
!#        core equation.
!#
!# 2.) cc_initPreconditioner
!#     -> Initialises a spatial preconditioner structure with parameters from
!#        the DAT file. This is needed for preconditioning in space.
!#     -> Extension to cc_createPreconditioner.
!#
!# 3.) cc_donePreconditioner
!#     -> Cleans up a spatial preconditioner structure initialised by
!#        cc_initPreconditioner. 
!#     -> Extension to cc_releasePreconditioner.
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
!# Auxiliary routines, not to be called from outside:
!#
!# 1.) cc_checkAssembly
!#     -> Checks if the system matrices are compatible to the preconditioner.
!#        Set up some situation dependent 'tweak' flags for the assembly
!#        of the matrices.
!#
!# 2.) cc_finaliseMatrices
!#     -> Makes the system matrices compatible to the preconditioner if
!#        necessary.
!#
!# 3.) cc_unfinaliseMatrices 
!#     -> Reverts the changes of cc_finaliseMatrices and brings matrices
!#        into their original form.
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
!# vectors, routines of the module cc2dmediumm2matvecassembly are used.
!#
!# </purpose>
!##############################################################################

module cc2dmediumm2nonlinearcoreinit

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
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  
  use collection
  use convection
    
  use cc2dmediumm2basic
  use cc2dmediumm2nonlinearcore
  
  implicit none
  

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_allocSystemMatrix (rproblem,rlevelInfo,rmatrix)
  
!<description>
  ! Allocates memory for the system matrix. rlevelInfo provides information
  ! about the level where the system matrix should be created.
  !
  ! Before this routine is called, the structure of all matrices in
  ! rlevelInfo must have been set up!
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(IN) :: rproblem

  ! A level-info structure specifying the matrices of the problem.
  type(t_problem_lvl), intent(IN), target :: rlevelInfo
!</input>

!<output>
  ! A block matrix that receives the basic system matrix.
  type(t_matrixBlock), intent(OUT) :: rmatrix
!</output>

!</subroutine>

    ! local variables
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A t_ccmatrixComponents used for defining the matrix structure
    type(t_ccmatrixComponents) :: rmatrixAssembly
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%rdiscretisation
    
    ! Get a pointer to the template FEM matrix.
    p_rmatrixTemplateFEM => rlevelInfo%rmatrixTemplateFEM

    ! In the global system, there are two gradient matrices B1 and B2.
    ! Get a pointer to the template structure for these.
    p_rmatrixTemplateGradient => rlevelInfo%rmatrixTemplateGradient

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
    ! We allocate the system matrix using the cc_assembleMatrix routine.
    ! For this purpose, we have to initialise a t_ccmatrixComponents structure
    ! which defines the shape of the matrix. We simply set the parameters
    ! of those terms which wshould appear in the matrix to a value <> 0,
    ! that's enough for the memory allocation.
    
    rmatrixAssembly%dtheta1 = 1.0_DP   ! A velocity block
    rmatrixAssembly%deta1 = 1.0_DP     ! A gradient block
    rmatrixAssembly%dtau1 = 1.0_DP     ! A divergence block
    rmatrixAssembly%dkappa1 = 1.0_DP   ! Pressure block
    ! A Newton block, if we have Navier-Stokes. For the case
    ! that we use Newton
    rmatrixAssembly%dnewton1 = real(1-rproblem%iequation,DP)
    
    rmatrixAssembly%dtheta2 = 1.0_DP   ! A velocity block
    rmatrixAssembly%deta2 = 1.0_DP     ! A gradient block
    rmatrixAssembly%dtau2 = 1.0_DP     ! A divergence block
    ! A Newton block, if we have Navier-Stokes
    rmatrixAssembly%dnewton2 = real(1-rproblem%iequation,DP)
    rmatrixAssembly%dkappa2 = 1.0_DP   ! Pressure block
    
    rmatrixAssembly%p_rdiscretisation => rlevelInfo%rdiscretisation
    rmatrixAssembly%p_rmatrixTemplateFEM => rlevelInfo%rmatrixTemplateFEM
    rmatrixAssembly%p_rmatrixTemplateGradient => rlevelInfo%rmatrixTemplateGradient
    rmatrixAssembly%p_rmatrixStokes => rlevelInfo%rmatrixStokes
    rmatrixAssembly%p_rmatrixB1 => rlevelInfo%rmatrixB1
    rmatrixAssembly%p_rmatrixB2 => rlevelInfo%rmatrixB2
    rmatrixAssembly%p_rmatrixMass => rlevelInfo%rmatrixMass
    rmatrixAssembly%p_rmatrixIdentityPressure => rlevelInfo%rmatrixIdentityPressure

    call cc_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_AUTOMATIC,&
      rmatrix,rmatrixAssembly)
                                  
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
    type(t_linsolMGLevelInfo2), pointer :: p_rlevelInfo
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
          rpreconditioner%p_RfilterChain)
      
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
      call linsol_getMultigridLevel2 (p_rsolverNode,1,p_rlevelInfo)
      
      select case (icoarseGridSolverType)
      case (0)
        ! UMFPACK coarse grid solver. Easy.
        call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
        
      case (1)
        ! Defect correction with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
        
        call parlst_getvalue_string_direct (p_rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
        
      case (2)
        ! Defect correction with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
        
        call parlst_getvalue_string_direct (p_rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rpreconditioner%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
        
      case (3)
        ! BiCGSTab with diagonal VANKA preconditioning.
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
            rpreconditioner%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)
        
      case (4)
        ! BiCGSTab with full VANKA preconditioning.
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
            rpreconditioner%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

      case DEFAULT
      
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
        
        case (0:7)

          nullify(p_rsmoother)
        
          ! This is some kind of VANKA smoother. Initialise the correct one.
          select case (ismootherType)
          case (0)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
          case (1)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
          case (2)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
          case (3)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIRECT)
          case (4)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
          case (5)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAGDIR)
          case (6)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%p_RfilterChain)
          case (7)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rpreconditioner%p_RfilterChain)
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
          call linsol_getMultigridLevel2 (p_rsolverNode,ilev,p_rlevelInfo)
          p_rlevelInfo%p_rpresmoother => p_rsmoother
          p_rlevelInfo%p_rpostsmoother => p_rsmoother
          
          ! Set up the interlevel projection structure on all levels
          p_rlevelInfo%rinterlevelProjection = &
            rpreconditioner%p_rprojection
          
        case DEFAULT
        
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
      
      case (0:7)

        nullify(p_rsmoother)
      
        ! This is some kind of VANKA smoother. Initialise the correct one.
        select case (ismootherType)
        case (0)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
        case (1)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
        case (2)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
        case (3)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIRECT)
        case (4)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
        case (5)
          call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAGDIR)
        case (6)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOC)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%p_RfilterChain)
        case (7)
          call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVSTOCDIAG)
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
              rpreconditioner%p_RfilterChain)
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
        call linsol_initDefCorr (p_rsolverNode,p_rsmoother,rpreconditioner%p_RfilterChain)
      else
        call linsol_initBiCGStab (p_rsolverNode,p_rsmoother,rpreconditioner%p_RfilterChain)
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

  subroutine cc_initPreconditioner (rproblem,nlmin,nlmax,rpreconditioner)
  
!<description>
  ! Initialises the given spatial preconditioner structure rpreconditioner.
  ! Creates the structure with cc_createPreconditioner and saves all
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
    logical :: bneumann
    type(t_parlstSection), pointer :: p_rsection

    ! Basic initialisation of the nonlinenar iteration structure.
    call cc_createPreconditioner (rpreconditioner,nlmin,nlmax)    
    
    ! Assign the matrix pointers in the nonlinear iteration structure to
    ! all our matrices that we want to use.
    do ilevel = nlmin,nlmax
    
      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixStokes => &
        rproblem%RlevelInfo(ilevel)%rmatrixStokes
        
      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixB1 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB1
        
      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixB2 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB2

      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixMass => &
        rproblem%RlevelInfo(ilevel)%rmatrixMass

      rpreconditioner%RcoreEquation(ilevel)%p_rtempVector1 => &
        rproblem%RlevelInfo(ilevel)%rtempVector1

      rpreconditioner%RcoreEquation(ilevel)%p_rtempVector2 => &
        rproblem%RlevelInfo(ilevel)%rtempVector2

      rpreconditioner%RcoreEquation(ilevel)%p_rtempVector3 => &
        rproblem%RlevelInfo(ilevel)%rtempVector3

      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixIdentityPressure => &
        rproblem%RlevelInfo(ilevel)%rmatrixIdentityPressure

    end do
      
    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions.
    allocate(rpreconditioner%p_RfilterChain(4))
    
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter for defect vectors:
    rpreconditioner%p_RfilterChain(1)%ifilterType = &
        FILTER_DISCBCDEFREAL

    ! The second filter filters for boundary conditions of fictitious boundary
    ! components
    rpreconditioner%p_RfilterChain(2)%ifilterType = &
        FILTER_DISCBCDEFFICT
    
    ! Do we have Neumann boundary?
    !
    ! The bhasNeumannBoundary flag of the higher level decides about that...
    bneumann = rproblem%RlevelInfo(rproblem%NLMAX)%bhasNeumannBoundary
    rpreconditioner%p_RfilterChain(3)%ifilterType = FILTER_DONOTHING
    if (.not. bneumann) then
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rpreconditioner%p_RfilterChain(3)%ifilterType = FILTER_TOL20
      rpreconditioner%p_RfilterChain(3)%itoL20component = NDIM2D+1

      rpreconditioner%p_RfilterChain(4)%ifilterType = FILTER_TOL20
      rpreconditioner%p_RfilterChain(4)%itoL20component = 2*(NDIM2D+1)
    end if
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_donePreconditioner (rpreconditioner)
  
!<description>
  ! Releases memory allocated in cc_initPrecoditioner..
  ! The routine automatically calls cc_releasePreconditioner to release
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
          rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner)
        deallocate(rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner)

        ! Also release B1^T and B2^T everywhere where they exist.
        if (associated(rpreconditioner%RcoreEquation(i)%p_rmatrixB1T)) then
          call lsyssc_releaseMatrix (rpreconditioner%RcoreEquation(i)%p_rmatrixB1T)
          deallocate(rpreconditioner%RcoreEquation(i)%p_rmatrixB1T)
        end if

        if (associated(rpreconditioner%RcoreEquation(i)%p_rmatrixB2T)) then
          call lsyssc_releaseMatrix (rpreconditioner%RcoreEquation(i)%p_rmatrixB2T)
          deallocate(rpreconditioner%RcoreEquation(i)%p_rmatrixB2T)
        end if
      end do
      
      ! Release the temporary vector(s)
      call lsyssc_releaseVector (rpreconditioner%p_rtempVectorSc)
      deallocate(rpreconditioner%p_rtempVectorSc)
      
      ! Clean up data about the projection etc.
      call mlprj_doneProjection(rpreconditioner%p_rprojection)
      deallocate(rpreconditioner%p_rprojection)

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      call linsol_doneStructure (rpreconditioner%p_rsolverNode)
      call cc_doneLinearSolver (rpreconditioner)
      
    case DEFAULT
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

    ! Release the filter chain for the defect vectors.
    if (associated(rpreconditioner%p_RfilterChain)) &
      deallocate(rpreconditioner%p_RfilterChain)
      
    call cc_releasePreconditioner (rpreconditioner)

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
    integer :: i,icomponents
    integer(PREC_VECIDX) :: imaxmem
    character(LEN=PARLST_MLDATA) :: ssolverName,sstring
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
      
      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      allocate(rpreconditioner%p_rprojection)
      call mlprj_initProjectionDiscr (rpreconditioner%p_rprojection,p_rdiscretisation)
      
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      call cc_getProlRest (rpreconditioner%p_rprojection, &
          rproblem%rparamList,  'CC-PROLREST')
      
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
            rpreconditioner%p_rprojection,&
            rproblem%RlevelInfo(i-1)% &
              rdiscretisation%RspatialDiscr(1:p_rdiscretisation%ncomponents),&
            rproblem%RlevelInfo(i)% &
              rdiscretisation%RspatialDiscr(1:p_rdiscretisation%ncomponents)))
      end do
      
      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      allocate(rpreconditioner%p_rtempVectorSc)
      call lsyssc_createVector (rpreconditioner%p_rtempVectorSc,&
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
    logical :: btranspose
    character(LEN=PARLST_MLDATA) :: sstring,snewton

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
          ! Create an empty matrix structure or clean up the old one.
          if (.not. associated(rpreconditioner%RcoreEquation(i)%&
              p_rmatrixPreconditioner)) then
            allocate(rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner)
          else
            call lsysbl_releaseMatrix (rpreconditioner%RcoreEquation(i)%&
                p_rmatrixPreconditioner)
          end if
        
          ! Allocate memory for the basic submatrices.
          ! The B1^T and B2^T matrices are saved as 'virtual transpose' of B1 and
          ! B2 by default, we must change them later if necessary.
          call cc_allocSystemMatrix (rproblem,rproblem%RlevelInfo(i),&
              rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner)
        
          ! Attach boundary conditions
          rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner%p_rdiscreteBC &
              => rproblem%RlevelInfo(i)%p_rdiscreteBC
          
          rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner%p_rdiscreteBCfict &
              => rproblem%RlevelInfo(i)%p_rdiscreteFBC
        end if
        
        ! On the current level, set up a global preconditioner matrix.
      
        p_rmatrixPreconditioner => &
            rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner
        
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
            p_rmatrixTempateFEM => rproblem%RlevelInfo(i)%rmatrixTemplateFEM
            
            ! If we have a Stokes problem, A12 and A21 don't exist, so we have nothing to
            ! do. In the Navier-Stokes case however, we have A12 and A21, so create them!
            ! Furthermore, there is to be a A51 and A42 submatrix allocated for the
            ! reactive coupling mass matrix R!
            if (rproblem%iequation .eq. 0) then
            
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
        
          btranspose = .false.
          select case (rpreconditioner%rprecSpecials%isolverType)
          case (0)
            ! UMFPACK solver.
            if (i .eq. NLMAX) then
              ! UMFPACK solver. Tweak the matrix on the max. level.
              ! Ignore the other levels.
              btranspose = .true.
              
              if (rpreconditioner%rprecSpecials%bpressureGloballyIndefinite) then
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
                btranspose = .true.
                
                if (rpreconditioner%rprecSpecials%bpressureGloballyIndefinite) then
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
                btranspose = .true.
              end if              
              
            end if

          end select
          
          if (btranspose) then
            ! Allocate memory for B1^T and B2^T and transpose B1 and B2.
            if (.not. associated(rpreconditioner%RcoreEquation(i)%p_rmatrixB1T)) then
              allocate(rpreconditioner%RcoreEquation(i)%p_rmatrixB1T)
            else
              call lsyssc_releaseMatrix(rpreconditioner%RcoreEquation(i)%p_rmatrixB1T)
            end if
            call lsyssc_transposeMatrix (&
                rpreconditioner%RcoreEquation(i)%p_rmatrixB1,&
                rpreconditioner%RcoreEquation(i)%p_rmatrixB1T,LSYSSC_TR_ALL)

            if (.not. associated(rpreconditioner%RcoreEquation(i)%p_rmatrixB2T)) then
              allocate(rpreconditioner%RcoreEquation(i)%p_rmatrixB2T)
            else
              call lsyssc_releaseMatrix(rpreconditioner%RcoreEquation(i)%p_rmatrixB2T)
            end if
            call lsyssc_transposeMatrix (&
                rpreconditioner%RcoreEquation(i)%p_rmatrixB2,&
                rpreconditioner%RcoreEquation(i)%p_rmatrixB2T,LSYSSC_TR_ALL)
                
            ! B1^T and B2^T have the same structure, release unnecessary memory.
            call lsyssc_duplicateMatrix (&
                rpreconditioner%RcoreEquation(i)%p_rmatrixB1T,&
                rpreconditioner%RcoreEquation(i)%p_rmatrixB2T,&
                LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
          end if

        end if
      
      end do
        
      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
      do i=NLMIN,NLMAX
        call lsysbl_duplicateMatrix ( &
          rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner, &
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

    ! Check if we have Neumann boundary components. If not, the matrices
    ! may have to be changed, depending on the solver.
    rpreconditioner%rprecSpecials%bpressureGloballyIndefinite = &
      .not. rproblem%RlevelInfo(rproblem%NLMAX)%bhasNeumannBoundary

  end subroutine

end module
