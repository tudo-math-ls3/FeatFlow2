!##############################################################################
!# ****************************************************************************
!# <name> ccnonlinearcoreinit </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains initialisation routines for the core equation
!# (see also ccnonlinearcore).
!# These routines connect the "problem" structure with the "core equation"
!# structure. In detail, we have routines that initialise the preconditioner
!# and all the information structures that are used during the nonlinear
!# iteration.
!#
!# The following routines can be found here:
!#
!# 1.) cc_allocPrecSystemMatrix
!#     -> Allocates memory for the basic system matrix representing the
!#        core equation for the use in a preconditioner.
!#
!# 2.) cc_initNonlinearLoop
!#     -> Initialises a 'nonlinear iteration structure' with parameters from
!#        the DAT file. This is needed for solving the core equation.
!#     -> Extension to cc_createNonlinearLoop.
!#
!# 3.) cc_doneNonlinearLoop
!#     -> Cleans up a 'nonlinear iteration structure' initialised by
!#        cc_initNonlinearLoop.
!#     -> Extension to cc_releaseNonlinearLoop.
!#
!# 4.) cc_initPreconditioner
!#     -> Prepare preconditioner of nonlinear iteration
!#
!# 5.) cc_updatePreconditioner
!#     -> Updates the preconditioner if there was a change in the system matrices
!#
!# 6.) cc_releasePreconditioner
!#     -> Clean up preconditioner of nonlinear iteration
!#
!# 7.) cc_getProlRest
!#     -> Auxiliary routine: Set up interlevel projection structure
!#        with information from INI/DAT files
!#
!# Auxiliary routines, not to be called from outside:
!#
!# 1.) cc_adjustPrecSpecials
!#     -> Initialise some preconditioner 'specials' which must be
!#        respected when assembling matrices for the preconditioner.
!#
!# 2.) cc_initLinearSolver
!#     -> Initialise a linear solver for preconditioning
!#
!# 3.) cc_doneLinearSolver
!#     -> Release data allocated in cc_initLinearSolver
!#
!# The module works in tight relationship to ccnonlinearcore.
!# ccnlinearcodeinit provides the routines to initialise
!# preconditioner and important structures using the problem related
!# structure. This module ccnonlinearcore on the other hand
!# contains only the 'main' worker routines that do the work of the
!# nonlinear iteration -- independent of the problem structure!
!#
!# Note that this module and the "nonlinearcore" module are the only modules
!# that 'know' the actual structure of the system matrix and how to link
!# it to the main problem! For the actual assembly of the matrices and defect
!# vectors, routines of the module ccmatvecassembly are used.
!#
!# </purpose>
!##############################################################################

module ccnonlinearcoreinit

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use filtersupport
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use multilevelprojection
  use matrixrestriction
  use statistics
  
  use collection
  use convection
    
  use ccbasic
  use ccnonlinearcore
  
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
  ! =CCMASM_MTP_AUTOMATIC: standard matrix, A55=A66
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

!</subroutine>

    ! local variables
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A t_nonlinearCCMatrix used for defining the matrix structure
    type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
    
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%rdiscretisation
    
    ! Get a pointer to the template FEM matrix.
    p_rmatrixTemplateFEM => rlevelInfo%rasmTempl%rmatrixTemplateFEM

    ! In the global system, there are two gradient matrices B1 and B2.
    ! Get a pointer to the template structure for these.
    p_rmatrixTemplateGradient => rlevelInfo%rasmTempl%rmatrixTemplateGradient

    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
      
    ! Let us consider the global system in detail. It has roughly
    ! the following shape:
    !
      !    ( A11   A12   A13   0     A15   0      BS1  )   ( A11   A12   A13   A14   A15   A16    A17  )
      !    ( A21   A22   0     A24   0     A26    BS2  )   ( A21   A22   A23   A24   A25   A26    A27  )
      !    ( A31   0     A33   0     0     0      0    )   ( A31   A32   A33   A34   A35   A36    A37  )
      !    ( 0     A42   0     A44   0     0      0    )   ( A41   A42   A43   A44   A45   A46    A47  )
      !    ( 0     0     A53   0     A55   A56    BF1  ) = ( A51   A52   A53   A54   A55   A56    A57  )
      !    ( 0     0     0     A64   A65   A66    BF2  )   ( A61   A62   A63   A64   A65   A66    A67  )
      !    ( 0     0     DS1   DS2   DF1   DF2    0    )   ( A71   A72   A73   A74   A75   A76    A77  )
    !
    ! Indeed, the 2nd shape may occur if the matrix is used as Newton
    ! matrix using UMFPACK on a pure Dirichlet problem.
    ! (Newton=A56/A65 present, UMFPACK on pure Dirichlet=A33 present.)
    ! So we have to reserve memory for all submatrices.
    ! All matices may have multiplication factors in their front.
    !
    ! The structure of the matrices A55 and A66 of the global system matrix
    ! is governed by the template FEM matrix.
    ! Initialise them with the same structure, i.e. A55, A66 share (!) their
    ! structure (not the entries) with that of the template matrix.
    !
    ! We allocate the system matrix using the cc_assembleMatrix routine.
    ! For this purpose, we have to initialise a t_nonlinearCCMatrix structure
    ! which defines the shape of the matrix. We simply set the parameters
    ! of those terms which wshould appear in the matrix to a value <> 0,
    ! that is enough for the memory allocation.
    
    call cc_initNonlinMatrix (rnonlinearCCMatrix,rproblem,&
        rlevelInfo%rdiscretisation,rlevelInfo%rasmTempl,&
        rlevelInfo%rdynamicInfo)
    
    call cc_prepareNonlinMatrixAssembly (rnonlinearCCMatrix,&
        ilev,nlmin,nlmax,rprecSpecials)

    rnonlinearCCMatrix%dtheta = 1.0_DP   ! A velocity block
    rnonlinearCCMatrix%deta = 1.0_DP     ! A gradient block
    rnonlinearCCMatrix%dtau = 1.0_DP     ! A divergence block

    call cc_assembleMatrix (CCMASM_ALLOCMEM,cmatrixType,&
        rmatrix,rnonlinearCCMatrix,rproblem)
                                  
    ! That is it, all submatrices are set up.
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initNonlinearLoop (rproblem,nlmin,nlmax,rvector,rrhs,&
      rnonlinearIteration,sname)
  
!<description>
  ! Initialises the given nonlinear iteration structure rnonlinearIteration.
  ! Creates the structure with cc_createNonlinearLoop and saves all
  ! problem dependent parameters and matrices in it.
  ! The routine automatically calls cc_createNonlinearLoop to initialise
  ! the structure with internal parameters.
  ! Note: This does not initialise the preconditioner in that structure!
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! Minimum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners.
  integer, intent(in) :: nlmin
  
  ! Maximum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners. This level must correspond to rvector and rrhs.
  integer, intent(in) :: nlmax

  ! The solution vector which is modified later during the nonlinear iteration.
  type(t_vectorBlock), intent(in), target :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(in), target :: rrhs

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: sname
!</input>

!<output>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Is filled with data.
  type(t_ccnonlinearIteration), intent(out) :: rnonlinearIteration
!</output>

!</subroutine>

    ! local variables
    integer :: ilevel
    logical :: bneumann
    type(t_parlstSection), pointer :: p_rsection

    rnonlinearIteration%NLMIN = NLMIN
    rnonlinearIteration%NLMAX = NLMAX

    ! Initialise the matrix pointers on all levels that we have to maintain.
    allocate(rnonlinearIteration%RcoreEquation(NLMIN:NLMAX))

    rnonlinearIteration%MT_OutputLevel = rproblem%MT_OutputLevel
    
    ! Get the minimum/maximum damping parameter from the parameter list, save
    ! them to the nonlinear iteration structure (which is now initialised).
    call parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMin', rnonlinearIteration%domegaMin, 0.0_DP)
                              
    call parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMax', rnonlinearIteration%domegaMax, 2.0_DP)
    
    ! Save pointers to the RHS and solution vector
    rnonlinearIteration%p_rsolution => rvector
    rnonlinearIteration%p_rrhs => rrhs
    
    ! Set the preconditioner to 'nothing'
    rnonlinearIteration%rpreconditioner%ctypePreconditioning = -1
    
    ! Deactivate any 'tweak' flags in the final-assembly structure
    rnonlinearIteration%rprecSpecials%iadaptiveMatrices = 0
    
    ! Assign the matrix pointers in the nonlinear iteration structure to
    ! all our matrices that we want to use.
    do ilevel = nlmin,nlmax
      rnonlinearIteration%RcoreEquation(ilevel)%p_rdiscretisation => &
        rproblem%RlevelInfo(ilevel)%rdiscretisation

      rnonlinearIteration%RcoreEquation(ilevel)%p_rasmTempl => &
        rproblem%RlevelInfo(ilevel)%rasmTempl
      rnonlinearIteration%RcoreEquation(ilevel)%p_rdynamicInfo => &
        rproblem%RlevelInfo(ilevel)%rdynamicInfo
      rnonlinearIteration%RcoreEquation(ilevel)%p_rtempVector => &
        rproblem%RlevelInfo(ilevel)%rtempVector
    end do
      
    ! Clear auxiliary variables for the nonlinear iteration
    rnonlinearIteration%DresidualInit = 0.0_DP
    rnonlinearIteration%DresidualOld  = 0.0_DP
    
    call parlst_querysection(rproblem%rparamList, sname, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line ('Cannot create nonlinear solver; no section '''//&
          trim(sname)//'''!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initNonlinearLoop')
      call sys_halt()
    end if

    ! Get stopping criteria of the nonlinear iteration
    call parlst_getvalue_double (p_rsection, 'depsD_u', &
                                 rnonlinearIteration%DepsNL(1), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsD_vS', &
                                 rnonlinearIteration%DepsNL(2), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsD_vF', &
                                 rnonlinearIteration%DepsNL(3), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnonlinearIteration%DepsNL(4), 0.1_DP)
!     ------------------- ------------------------ ----------------------
    call parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnonlinearIteration%DepsNL(5), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsVSR', &
                                 rnonlinearIteration%DepsNL(6), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsVFR', &
                                 rnonlinearIteration%DepsNL(7), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnonlinearIteration%DepsNL(8), 0.1_DP)

    call parlst_getvalue_double (p_rsection, 'dDampingD', &
                                 rnonlinearIteration%DepsNL(9), 0.1_DP)
      
    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions. This filter chain is applied to each
    ! defect vector during the linear and nonlinear iteration.
    allocate(rnonlinearIteration%p_RfilterChain(3))
    
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter for defect vectors:
    rnonlinearIteration%p_RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! The second filter filters for boundary conditions of fictitious boundary
    ! components
    rnonlinearIteration%p_RfilterChain(2)%ifilterType = FILTER_DISCBCDEFFICT
    
    ! Do we have Neumann boundary?
    !
    ! The bhasNeumannBoundary flag of the higher level decides about that...
    bneumann = rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo%bhasNeumannBoundary
    rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_DONOTHING
    if (.not. bneumann) then
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_TOL20
      rnonlinearIteration%p_RfilterChain(3)%itoL20component = 3*NDIM2D+1
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneNonlinearLoop (rnonlinearIteration)
  
!<description>
  ! Releases memory allocated in cc_initNonlinearLoop.
  ! The routine automatically calls cc_releaseNonlinearLoop to release
  ! internal parameters.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that should be cleaned up.
  type(t_ccNonlinearIteration), intent(inout) :: rnonlinearIteration
!</inputoutput>

!</subroutine>
    
    ! Release the filter chain for the defect vectors.
    if (associated(rnonlinearIteration%p_RfilterChain)) &
      deallocate(rnonlinearIteration%p_RfilterChain)
      
    if (associated(rnonlinearIteration%RcoreEquation)) &
      deallocate(rnonlinearIteration%RcoreEquation)

    rnonlinearIteration%NLMIN = 0
    rnonlinearIteration%NLMAX = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initLinearSolver (rproblem,rnonlinearIteration,ssection)
  
!<description>
  ! This routine initialises a linear solver structure to be used
  ! for preconditioning of the nonlinear defect.
!</description>

!<input>
  ! Name of the section in the parameter list that contains the parameters
  ! of the linear solver.
  character(LEN=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! The nonlinear iteration structure to which a preconditioner should
  ! be initialised. The parameters for the linear solver are written to
  ! the rpreconditioner substructure.
  ! The rprecSpecials structure is configured according to the linear
  ! solver.
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration
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
        
    rnonlinearIteration%rprecSpecials%isolverType = isolverType
    rnonlinearIteration%rprecSpecials%ismootherType = ismootherType
    rnonlinearIteration%rprecSpecials%icoarseGridSolverType = icoarseGridSolverType

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
      nlevels = rnonlinearIteration%NLMAX - rnonlinearIteration%NLMIN + 1
      
      call linsol_initMultigrid2 (p_rsolverNode,nlevels,&
          rnonlinearIteration%p_RfilterChain)
      
      ! Manually trim the coarse grid correction in Multigrid to multiply the
      ! pressure equation with -1. This (un)symmetrises the operator and gives
      ! much better convergence rates.
      call cgcor_release(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
      call cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,3*NDIM2D+1)
      p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(7) &
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
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rnonlinearIteration%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
            
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rnonlinearIteration%rprecSpecials%bneedVirtTransposedDonCoarse = .true.
        
      case (2)
        ! BiCGStab with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rnonlinearIteration%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)
        
        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rnonlinearIteration%rprecSpecials%bneedVirtTransposedDonCoarse = .true.

      case (3)
        ! BiCGStab with full VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
        
        call parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        read (sstring,*) spreconditionerSection
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rnonlinearIteration%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

        ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
        rnonlinearIteration%rprecSpecials%bneedVirtTransposedDonCoarse = .true.

      case (4)
        ! BiCGStab with general VANKA preconditioning.
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
            rnonlinearIteration%p_RfilterChain)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        call linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)

      case default
      
        call output_line ('Unknown coarse grid solver!', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
        call sys_halt()

      end select
      
      ! Put the coarse grid solver node to the preconditioner structure.
      rnonlinearIteration%rpreconditioner%p_rcgrSolver => p_rlevelInfo%p_rcoarseGridSolver
      
      ! Now after the coarse grid solver is done, we turn to the smoothers
      ! on all levels. Their initialisation is similar to the coarse grid
      ! solver. Note that we use the same smoother on all levels, for
      ! presmoothing as well as for postsmoothing.
      
      do ilev = 2,nlevels

        ! Initialise the smoothers.
        select case (ismootherType)
        
        case (0:7,101:102,201:202)

          nullify(p_rsmoother)
        
          ! This is some kind of VANKA smoother. Initialise the correct one.
          select case (ismootherType)
          case (0)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
        
          case (1)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)

          case (2)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rnonlinearIteration%rprecSpecials%bneedVirtTransposedD = .true.

          case (3)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVSTDIRECT)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rnonlinearIteration%rprecSpecials%bneedVirtTransposedD = .true.

          case (4)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rnonlinearIteration%rprecSpecials%bneedVirtTransposedD = .true.

          case (5)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rnonlinearIteration%rprecSpecials%bneedVirtTransposedD = .true.

          case (6)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rnonlinearIteration%p_RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rnonlinearIteration%rprecSpecials%bneedVirtTransposedD = .true.

          case (7)
            call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DFNAVST)
            call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,&
                rnonlinearIteration%p_RfilterChain)

            ! We need virtually transposed B-matrices as D-matrices for this preconditioner.
            rnonlinearIteration%rprecSpecials%bneedVirtTransposedD = .true.

          ! --- NEW IMPLEMENTATION ---
          case (101)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_NAVST2D_DIAG)
          case (102)
            call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_NAVST2D_FULL)
          
          ! --- SP-SOR ---
          case (201)
            call linsol_initSPSOR (p_rsmoother,LINSOL_SPSOR_NAVST2D)
          case (202)
            call linsol_initSPSOR (p_rsmoother,LINSOL_SPSOR_NAVST2D_DIAG)

          case default
          
            call output_line ('Unknown smoother!', &
                OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
            call sys_halt()

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
          
          ! Set up the interlevel projection structure on all levels
          call linsol_initProjMultigrid2Level(p_rlevelInfo,&
              rnonlinearIteration%RcoreEquation(ilev+rnonlinearIteration%NLMIN-1)%&
              p_rprojection)
          
        case default
        
          call output_line ('Unknown smoother!', &
              OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
          call sys_halt()

        end select
      
      end do

      ! Get information about adaptive matrix generation from INI/DAT files
      call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
          'iAdaptiveMatrix', rnonlinearIteration%rprecSpecials%iadaptiveMatrices, 0)
                                
      call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
          'dAdMatThreshold', rnonlinearIteration%rprecSpecials%dAdMatThreshold, 20.0_DP)

    case default
    
      call output_line ('Unknown linear solver!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
      call sys_halt()
      
    end select

    ! Put the final solver node to the preconditioner structure.
    rnonlinearIteration%rpreconditioner%p_rsolverNode => p_rsolverNode

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneLinearSolver (rnonlinearIteration)
  
!<description>
  ! Releases information of the linear solver preconditioner from the
  ! structure rnonlinearIteration. Cleans up what was configured
  ! in cc_initLinearSolver.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure to be cleaned up.
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! Release the solver and all subsolvers.
    call linsol_releaseSolver(rnonlinearIteration%rpreconditioner%p_rsolverNode)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine cc_initPreconditioner (rproblem,rnonlinearIteration,rvector,rrhs)
  
!<description>
  ! This routine prepares the preconditioner we used during the
  ! nonlinear iteration. The structure rpreconditioner will be initialised
  ! based on the information in rproblem.
  ! Necessary variables will be added to the collection structure in
  ! rproblem\%rcollection to be available in the callback routines.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  type(t_vectorBlock), intent(in) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(in) :: rrhs
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! This is configured according to the preconditioner as specified in
  ! the DAT files.
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
    integer :: imaxmem
    character(LEN=PARLST_MLDATA) :: ssolverName,sstring

    ! At first, ask the parameters in the INI/DAT file which type of
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    call parlst_getvalue_int (rproblem%rparamList, 'CC2D-NONLINEAR', &
        'itypePreconditioning', &
        rnonlinearIteration%rpreconditioner%ctypePreconditioning, 1)
    
    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioner
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX
      
      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      call parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      if (sstring .ne. '') read (sstring,*) ssolverName
      if (ssolverName .eq. '') then
        call output_line ('No linear subsolver!', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initPreconditioner')
        call sys_halt()
      end if
      
      ! Initialise a standard interlevel projection structures.
      do i = NLMIN+1, NLMAX
      
        ! Allocate the projection structure
        allocate(rnonlinearIteration%RcoreEquation(i)%p_rprojection)
        
        ! Initialise the projection based on the discretisation.
        call mlprj_initProjectionDiscr (&
            rnonlinearIteration%RcoreEquation(i)%p_rprojection,&
            rnonlinearIteration%RcoreEquation(i)%p_rdiscretisation)
      
        ! Initialise the projection structure with data from the INI/DAT
        ! files. This allows to configure prolongation/restriction.
        call cc_getProlRest (rnonlinearIteration%RcoreEquation(i)%p_rprojection, &
            rproblem%RlevelInfo(i), rproblem%rparamList,  'CC-PROLREST')
      
      end do
      
      ! Initialise the linear solver as configured in the DAT file.
      call cc_initLinearSolver (rproblem,rnonlinearIteration,ssolverName)
      
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
            rnonlinearIteration%RcoreEquation(i)%p_rprojection,&
            rproblem%RlevelInfo(i-1)% &
              rdiscretisation%RspatialDiscr(1:rrhs%nblocks),&
            rproblem%RlevelInfo(i)% &
              rdiscretisation%RspatialDiscr(1:rrhs%nblocks)))
      end do
      
      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      allocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      call lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc,&
                                max(imaxmem,rrhs%NEQ),.false.)
      
      ! Set up a second temporary vector that we need for calculating
      ! the optimal defect correction.
      allocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      call lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2,&
                                rrhs%NEQ,.false.,ST_DOUBLE)
      
      ! Initialise the matrices.
      call cc_updatePreconditioner (rproblem,rnonlinearIteration,&
          rvector,rrhs,.true.,.true.)
      
    case default
      
      ! Unknown preconditioner
      call output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initPreconditioner')
      call sys_halt()
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_updatePreconditioner (rproblem,rnonlinearIteration,&
      rvector,rrhs,binit,bstructuralUpdate)
  
!<description>
  ! This routine has to be called whenever the system matrices change.
  ! It initialises (depending on the system matrices) the matrices of the
  ! preconditioner or performs an update of them.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  type(t_vectorBlock), intent(in) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(in) :: rrhs

  ! First initialisation.
  ! Has to be set to TRUE on the first call. Initialises the preconditioner
  ! with the structure of the matrices.
  logical, intent(in) :: binit

  ! Whether the structure of the system matrices is new.
  ! This variable has to be set to TRUE whenever there was a structure in
  ! the system matrices. This reinitialises the linear solver.
  logical, intent(in) :: bstructuralUpdate
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Preconditioner data is saved here.
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i,cmatrixType
    character(LEN=PARLST_MLDATA) :: sstring,snewton
    type(t_timer) :: rtimer

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A pointer to the matrix of the preconditioner
    type(t_matrixBlock), pointer :: p_rmatrixPreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    
    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioner
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
    
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Adjust the special preconditioner parameters to the current situation.
      call cc_adjustPrecSpecials (rproblem,rnonlinearIteration,&
          rnonlinearIteration%rprecSpecials)
      
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX
      
      ! Initialise the preconditioner matrices on all levels.
      do i=NLMIN,NLMAX
      
        ! Prepare the preconditioner matrices level i.
        if (binit .or. bstructuralUpdate) then

          ! What type of matrix do we have? Is it a 'primitive' matrix
          ! (with A11 and A22 being the same) or a 'full-tensor matrix'
          ! with Aij all independent from each other
          if ((rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
              CCPREC_NEWTON) .or. &
              (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
              CCPREC_NEWTONDYNAMIC)) then
            cmatrixType = CCMASM_MTP_FULLTENSOR
          else
            cmatrixType = CCMASM_MTP_AUTOMATIC
          end if

          ! Allocate memory for the matrix or release the existing matrix,
          ! if there is a structural update.
          if (binit) then
            allocate(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
          else
            call lsysbl_releaseMatrix(&
                rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
          end if
          
          ! Allocate memory for the basic submatrices.
          ! The routine also reserves memory for A(7,7) for the case it is needed,
          ! but switches that matrix off.
          call cc_allocPrecSystemMatrix (rproblem,rnonlinearIteration%rprecSpecials,&
              i,nlmin,nlmax,rproblem%RlevelInfo(i),cmatrixType,&
              rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
              
          ! Attach boundary conditions
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner%p_rdiscreteBC &
              => rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteBC
          
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner%p_rdiscreteBCfict &
              => rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteFBC
        end if
        
        ! On the current level, set up a global preconditioner matrix.
        p_rmatrixPreconditioner => &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner
      
        ! ----------------------------------------------------
        ! Should the linear solver use the adaptive Newton matrix?
            
        if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
          CCPREC_NEWTONDYNAMIC) then
          
          ! We have even the extended, dynamic Newton as preconditioner.
          ! Put the parameters for the extended Newton from the DAT file
          ! into the Adaptive-Newton configuration block.
          
          call parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                      'spreconditionerAdaptiveNewton', sstring, '')
          snewton = ''
          if (sstring .ne. '') read (sstring,*) snewton
          if (snewton .ne. '') then
            ! Initialise the parameters of the adaptive Newton
            call parlst_getvalue_int (rproblem%rparamList, snewton, &
                'nminFixPointIterations', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%nminFixPointIterations, 0)

            call parlst_getvalue_int (rproblem%rparamList, snewton, &
                'nmaxFixPointIterations', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%nmaxFixPointIterations, 999)

            call parlst_getvalue_double (rproblem%rparamList, snewton, &
                'depsAbsNewton', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%depsAbsNewton, 1E-5_DP)

            call parlst_getvalue_double (rproblem%rparamList, snewton, &
                'depsRelNewton', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%depsRelNewton, 1E99_DP)

            call parlst_getvalue_int (rproblem%rparamList, snewton, &
                'cinexactNewton', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%cinexactNewton, 1)

            call parlst_getvalue_double (rproblem%rparamList, snewton, &
                'dinexactNewtonEpsRel', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%dinexactNewtonEpsRel, 1.0E-2_DP)

            call parlst_getvalue_double (rproblem%rparamList, snewton, &
                'dinexactNewtonExponent', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%dinexactNewtonExponent, 2.0_DP)
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
        ! In case we have a pure-dirichlet problem, we activate the 7,7-submatrix
        ! It is probably needed for the preconditioner to make the pressure definite.
        if (bstructuralUpdate) then
          select case (rnonlinearIteration%rprecSpecials%isolverType)
          case (0)
            ! UMFPACK solver.
            if (i .eq. NLMAX) then

              if (rnonlinearIteration%rprecSpecials%bneedPressureDiagonalBlock) then
                ! Activate the 7,7-block, UMFPACK needs it in case the pressure
                ! is globally indefinite.
                p_rmatrixPreconditioner%RmatrixBlock(7,7)%dscaleFactor = 1.0_DP
              else
                p_rmatrixPreconditioner%RmatrixBlock(7,7)%dscaleFactor = 0.0_DP
              end if
            end if

          case (1)
            ! Multigrid solver. Treat the matrix at the coarse level if there is
            ! UMFPACK chosen as coarse grid solver.
            if (i .eq. NLMIN) then
            
              if ((rnonlinearIteration%rprecSpecials%icoarseGridSolverType .eq. 0) .or. &
                  (rnonlinearIteration%rprecSpecials%icoarseGridSolverType .eq. 4)) then
                
                if (rnonlinearIteration%rprecSpecials%bneedPressureDiagonalBlock) then
                  ! Activate the 7,7-block, UMFPACK needs it in case the pressure
                  ! is globally indefinite.
                  p_rmatrixPreconditioner%RmatrixBlock(7,7)%dscaleFactor = 1.0_DP
                else
                  p_rmatrixPreconditioner%RmatrixBlock(7,7)%dscaleFactor = 0.0_DP
                end if
                
              end if
              
            end if

          end select
          
        end if
        
      end do
      
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer)
      
      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
      allocate(Rmatrices(1:NLMAX))
      do i=NLMIN,NLMAX
        call lsysbl_duplicateMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
          Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end do
      
      call linsol_setMatrices(&
          rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(NLMIN:NLMAX))
          
      ! The solver got the matrices; clean up Rmatrices, it was only of temporary
      ! nature...
      do i=NLMIN,NLMAX
        call lsysbl_releaseMatrix (Rmatrices(i))
      end do
      deallocate(Rmatrices)
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      if (binit) then
      
        call linsol_initStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
            ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
        
      else if (bstructuralUpdate) then
      
        call linsol_updateStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
            ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
        
      end if
      
      ! Gather statistics
      call stat_stopTimer(rtimer)
      rproblem%rstatistics%dtimeLinearSolverFactorisation = &
        rproblem%rstatistics%dtimeLinearSolverFactorisation + rtimer%delapsedReal
      
    case DEFAULT
      
      ! Unknown preconditioner
      call output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_updatePreconditioner')
      call sys_halt()
      
    end select

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine cc_releasePreconditioner (rnonlinearIteration)
  
!<description>
  ! This routine releases the preconditioner for the nonlinear iteration
  ! which was prepared in cc_initPreconditioner. Memory is released
  ! from heap.
!</description>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! The preconditioner data is removed from that,
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Which preconditioner do we have?
    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CCPREC_NONE)
      ! No preconditioning
    case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the preconditioner matrix on every level
      do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        call lsysbl_releaseMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        deallocate(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        
      end do
      
      ! Release the temporary vector(s)
      call lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      deallocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      
      call lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      deallocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      
      ! Clean up data about the projection etc.
      do i = rnonlinearIteration%NLMIN+1, rnonlinearIteration%NLMAX
        call mlprj_doneProjection(rnonlinearIteration%RcoreEquation(i)%p_rprojection)
        deallocate(rnonlinearIteration%RcoreEquation(i)%p_rprojection)
      end do

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      call linsol_doneStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode)
      call cc_doneLinearSolver (rnonlinearIteration)
      
    case DEFAULT
      
      ! Unknown preconditioner
      call output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_releasePreconditioner')
      call sys_halt()
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getProlRest (rprojection, rlevelInfo, rparamList, sname)
  
!<description>
  ! Initialises an existing interlevel projection structure rprojection
  ! with parameters from the INI/DAT files. sname is the section in the
  ! parameter list containing parameters about prolongation restriction.
!</description>

!<input>
  ! The level info structure of the current level.
  type(t_problem_lvl), intent(in) :: rlevelInfo
  
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(in) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the prolongation/restriction.
  character(LEN=*), intent(in) :: sname
!</input>

!<output>
  ! An interlevel projection block structure containing an initial
  ! configuration of prolongation/restriction. The structure is modified
  ! according to the parameters in the INI/DAT file(s).
  type(t_interlevelProjectionBlock), intent(inout) :: rprojection
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
!   ###########################################################################
    ! Do we use matrix based projection for velocityies and displacement?
!   ###########################################################################
    call parlst_getvalue_int (p_rsection,'iprojTypeVelocity',i1,0)
    if(i1 .eq. 1) then
      
      ! Yes, so link the prolongation matrix to the projection structure,
      ! for both X- and Y-velocity.

! so what is the two numbers in RscalarProjection(1,3) ?
! 3: the block number
! 1: 1st FE space (for our special case we have just one FE space for each block)
!    however, in some cases, one may have more than one space for each block (say
! (mixture of 4-node quadraliteral and 3-node triangular elements) then
! one has to repeate the command 2 times one for RscalarProjection(1,3)
! and the other for RscalarProjection(2,3) where 1:quad. Elem. and 2: triangular elem.


! for uS1
      call mlprj_initMatrixProjection(&
          rprojection%RscalarProjection(1,1), &
          rlevelInfo%rasmTempl%rmatrixProlVelocity, &
          rlevelInfo%rasmTempl%rmatrixInterpVelocity)
! for uS2
      call mlprj_initMatrixProjection(&
          rprojection%RscalarProjection(1,2), &
          rlevelInfo%rasmTempl%rmatrixProlVelocity, &
          rlevelInfo%rasmTempl%rmatrixInterpVelocity)
! for vS1
      call mlprj_initMatrixProjection(&
          rprojection%RscalarProjection(1,3), &
          rlevelInfo%rasmTempl%rmatrixProlVelocity, &
          rlevelInfo%rasmTempl%rmatrixInterpVelocity)
! for vS2
      call mlprj_initMatrixProjection(&
          rprojection%RscalarProjection(1,4), &
          rlevelInfo%rasmTempl%rmatrixProlVelocity, &
          rlevelInfo%rasmTempl%rmatrixInterpVelocity)
! for vF1
      call mlprj_initMatrixProjection(&
          rprojection%RscalarProjection(1,5), &
          rlevelInfo%rasmTempl%rmatrixProlVelocity, &
          rlevelInfo%rasmTempl%rmatrixInterpVelocity)
! for vF2
      call mlprj_initMatrixProjection(&
          rprojection%RscalarProjection(1,6), &
          rlevelInfo%rasmTempl%rmatrixProlVelocity, &
          rlevelInfo%rasmTempl%rmatrixInterpVelocity)
!

          
    end if
    
! !  if the displacement and the 2 velocity have different FE spaces then
! ! we have to include three parameters in the discretisation.dat  for
! ! solid skeleton displacement and velocity (iprojTypeDisplacement, iprojTypeVel_S)
! ! and for pore fluid velocity (iprojTypeVel_F) and then do the follwoing:
! !   ###########################################################################
!     ! Do we use matrix based projection for solid skeleton displacement?
! !   ###########################################################################
!     call parlst_getvalue_int (p_rsection,'iprojTypeDisplacement',i1,0)
!     if(i1 .eq. 1) then
!       
!       ! Yes, so link the prolongation matrix to the projection structure,
!       ! for both X- and Y-velocity.
!       call mlprj_initMatrixProjection(&
!           rprojection%RscalarProjection(1,1), &
!           rlevelInfo%rasmTempl%rmatrixProlDisplacementSolid, &
!           rlevelInfo%rasmTempl%rmatrixInterpDisplacementSolid)
! !
!       call mlprj_initMatrixProjection(&
!           rprojection%RscalarProjection(1,2), &
!           rlevelInfo%rasmTempl%rmatrixProlDisplacementSolid, &
!           rlevelInfo%rasmTempl%rmatrixInterpDisplacementSolid)
!           
!     end if
! !   ###########################################################################
!     ! Do we use matrix based projection for solid skeleton velocity?
! !   ###########################################################################
!     call parlst_getvalue_int (p_rsection,'iprojTypeVelocity_S',i1,0)
!     if(i1 .eq. 1) then
!       
!       ! Yes, so link the prolongation matrix to the projection structure,
!       ! for both X- and Y-velocity.
!       call mlprj_initMatrixProjection(&
!           rprojection%RscalarProjection(1,3), &
!           rlevelInfo%rasmTempl%rmatrixProlVelocitySolid, &
!           rlevelInfo%rasmTempl%rmatrixInterpVelocitySolid)
!       call mlprj_initMatrixProjection(&
!           rprojection%RscalarProjection(1,4), &
!           rlevelInfo%rasmTempl%rmatrixProlVelocitySolid, &
!           rlevelInfo%rasmTempl%rmatrixInterpVelocitySolid)
!           
!     end if
! !   ###########################################################################
!     ! Do we use matrix based projection for pore fluid velocity?
! !   ###########################################################################
!     call parlst_getvalue_int (p_rsection,'iprojTypeVelocity_F',i1,0)
!     if(i1 .eq. 1) then
!       
!       ! Yes, so link the prolongation matrix to the projection structure,
!       ! for both X- and Y-velocity.
!       call mlprj_initMatrixProjection(&
!           rprojection%RscalarProjection(1,5), &
!           rlevelInfo%rasmTempl%rmatrixProlVelocityFluid, &
!           rlevelInfo%rasmTempl%rmatrixInterpVelocityFluid)
!       call mlprj_initMatrixProjection(&
!           rprojection%RscalarProjection(1,6), &
!           rlevelInfo%rasmTempl%rmatrixProlVelocityFluid, &
!           rlevelInfo%rasmTempl%rmatrixInterpVelocityFluid)
!           
!     end if
!   ###########################################################################
    ! Do we use matrix based projection for pressure?
!   ###########################################################################
    call parlst_getvalue_int (p_rsection,'iprojTypePressure',i1,0)
    if(i1 .eq. 1) then
      
      ! Yes, so link the prolongation matrix to the projection structure.
      call mlprj_initMatrixProjection(&
          rprojection%RscalarProjection(1,7), &
          rlevelInfo%rasmTempl%rmatrixProlPressure, &
          rlevelInfo%rasmTempl%rmatrixInterpPressure)
          
    end if
! what are the numbers in RscalarProjection(:,1:3*NDIM2D)
! : menas al the FE spaces,  1:3*NDIM2D: the blocks 1,2, ... 3*2

!   ###########################################################################
    ! Prolongation/restriction order for all the components except the pressure
!   ###########################################################################
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,1:3*NDIM2D)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,1:3*NDIM2D)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,1:3*NDIM2D)%iinterpolationOrder = i1
    end if


! ! in case of differen spaces for velocities and displacement:
! !    
! !   ###########################################################################
!     ! Prolongation/restriction order for solid skeleton displacement components
! !   ###########################################################################
!     call parlst_getvalue_int (p_rsection,'iinterpolationOrderDisp',i1,-1)
!     
!     if (i1 .ne. -1) then
!       ! Initialise order of prolongation/restriction for velocity components
!       rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
!       rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
!       rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
!     end if
! 
!     
! !   ###########################################################################
!     ! Prolongation/restriction order for solid skeleton velocity components
! !   ###########################################################################
!     call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel_S',i1,-1)
!     
!     if (i1 .ne. -1) then
!       ! Initialise order of prolongation/restriction for velocity components
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%iprolongationOrder  = i1
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%irestrictionOrder   = i1
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%iinterpolationOrder = i1
!     end if
! 
!     
! !   ###########################################################################
!     ! Prolongation/restriction order for pore fluid velocity NDIM components
! !   ###########################################################################
!     call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel_F',i1,-1)
!     
!     if (i1 .ne. -1) then
!       ! Initialise order of prolongation/restriction for velocity components
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%iprolongationOrder  = i1
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%irestrictionOrder   = i1
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%iinterpolationOrder = i1
!     end if

!   ###########################################################################
    ! Prolongation/restriction order for pressure
!   ###########################################################################
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,3*NDIM2D+1)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,3*NDIM2D+1)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,3*NDIM2D+1)%iinterpolationOrder = i1
    end if


!#############################################################
! Prolongation/restriction variant for displacement
! components in case of Q1~ discretisation
!#############################################################
    call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
    
    if (i1 .ne. -1) then
      rprojection%RscalarProjection(:,1:3*NDIM2D)%iprolVariant  = i1
      rprojection%RscalarProjection(:,1:3*NDIM2D)%irestVariant  = i1
    end if

! ! in case of differen spaces for velocities and displacement:
! 
! !#############################################################
! ! Prolongation/restriction variant for displacement
! ! components in case of Q1~ discretisation
! !#############################################################
!     call parlst_getvalue_int (p_rsection,'iinterpolationVariantDisp',i1,0)
!     
!     if (i1 .ne. -1) then
!       rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
!       rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1
!     end if
! !#############################################################
! ! Prolongation/restriction variant for solid skeleton velocity 
! !  components in case of Q1~ discretisation
! !#############################################################
!     call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel_S',i1,0)
!     
!     if (i1 .ne. -1) then
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%iprolVariant  = i1
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%irestVariant  = i1
!     end if
! !#############################################################
! ! Prolongation/restriction variant for pore fluid velocity 
! ! components in case of Q1~ discretisation
! !#############################################################
!     call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel_F',i1,0)
!     
!     if (i1 .ne. -1) then
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%iprolVariant  = i1
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%irestVariant  = i1
!     end if

! #############################################################
    ! Aspect-ratio indicator in case of Q1~ discretisation
    ! with extended prolongation/restriction, for u_S,v_S & v_F
! #############################################################
    call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
    
    if (i1 .ne. 1) then
      rprojection%RscalarProjection(:,1:3*NDIM2D)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,1:3*NDIM2D)%irestARIndicatorEX3Y  = i1
    end if

! ! in case of differen spaces for velocities and displacement:
! 
! ! ########################################################
!     ! Aspect-ratio indicator in case of Q1~ discretisation
!     ! with extended prolongation/restriction, for u_S
! ! ########################################################
!     call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YDisp',i1,1)
!     
!     if (i1 .ne. 1) then
!       rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
!       rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
!     end if
! ! ########################################################
!     ! Aspect-ratio indicator in case of Q1~ discretisation
!     ! with extended prolongation/restriction, for v_S
! ! ########################################################
!     call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel_S',i1,1)
!     
!     if (i1 .ne. 1) then
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%iprolARIndicatorEX3Y  = i1
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%irestARIndicatorEX3Y  = i1
!     end if
! ! ########################################################
!     ! Aspect-ratio indicator in case of Q1~ discretisation
!     ! with extended prolongation/restriction, for v_F
! ! ########################################################
!     call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel_F',i1,1)
!     
!     if (i1 .ne. 1) then
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%iprolARIndicatorEX3Y  = i1
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%irestARIndicatorEX3Y  = i1
!     end if


! ####################################################################### u_S
  ! Aspect-ratio bound for switching to constant prolongation/restriction
  ! in case of Q1~ discretisation with extended prolongation/restriction
! #######################################################################

    call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
    
    if (d1 .ne. 20.0_DP) then
      rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
    end if

! ! in case of differen spaces for velocities and displacement:
! 
! ! ####################################################################### u_S
!   ! Aspect-ratio bound for switching to constant prolongation/restriction
!   ! in case of Q1~ discretisation with extended prolongation/restriction
! ! #######################################################################
! 
!     call parlst_getvalue_double (p_rsection,'dintARboundEX3YDisp',d1,20.0_DP)
!     
!     if (d1 .ne. 20.0_DP) then
!       rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
!       rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
!     end if
! ! ####################################################################### v_S
!   ! Aspect-ratio bound for switching to constant prolongation/restriction
!   ! in case of Q1~ discretisation with extended prolongation/restriction
! ! #######################################################################
!     call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel_S',d1,20.0_DP)
!     
!     if (d1 .ne. 20.0_DP) then
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%dprolARboundEX3Y  = d1
!       rprojection%RscalarProjection(:,(NDIM2D+1):2*NDIM2D)%drestARboundEX3Y  = d1
!     end if
! ! ####################################################################### v_F
!   ! Aspect-ratio bound for switching to constant prolongation/restriction
!   ! in case of Q1~ discretisation with extended prolongation/restriction
! ! #######################################################################
!     call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel_F',d1,20.0_DP)
!     
!     if (d1 .ne. 20.0_DP) then
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%dprolARboundEX3Y  = d1
!       rprojection%RscalarProjection(:,(2*NDIM2D+1):3*NDIM2D)%drestARboundEX3Y  = d1
!     end if

  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine cc_adjustPrecSpecials (rproblem,rnonlinearIteration,rprecSpecials)
  
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
  type(t_problem), intent(inout), target :: rproblem

  ! Nonlinar iteration structure saving data about the actual configuration
  ! of the core equation.
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration
!</inputoutput>

!<inputoutput>
  ! Nonlinear iteration structure.
  ! The t_ccPreconditionerSpecials substructure that receives information how to
  ! finally assembly the matrices such that everything in the callback routines
  ! will work.
  type(t_ccPreconditionerSpecials), intent(inout) :: rprecSpecials
!</inputoutput>

!</subroutine>

    ! Check if we have Neumann boundary components. If not, the matrices
    ! may have to be changed, depending on the solver.
    rprecSpecials%bneedPressureDiagonalBlock = &
      .not. rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo%bhasNeumannBoundary

    ! If there are no Neumann BC`s, the pressure is indefinite.
    rprecSpecials%bpressureIndefinite = &
      .not. rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo%bhasNeumannBoundary

  end subroutine

end module
