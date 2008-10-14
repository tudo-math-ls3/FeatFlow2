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
!# 1.) cc_allocSystemMatrix
!#     -> Allocates memory for the basic system matrix representing the 
!#        core equation.
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

MODULE ccnonlinearcoreinit

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  USE linearsolverautoinitialise
  USE matrixrestriction
  USE statistics
  
  USE collection
  USE convection
    
  USE ccbasic
  USE ccnonlinearcore
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_allocSystemMatrix (rproblem,rlevelInfo,cmatrixType,rmatrix)
  
!<description>
  ! Allocates memory for the system matrix. rlevelInfo provides information
  ! about the level where the system matrix should be created.
  !
  ! Before this routine is called, the structure of all matrices in
  ! rlevelInfo must have been set up!
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(IN) :: rproblem

  ! Type of matrix.
  ! =CCMASM_MTP_AUTOMATIC: standard matrix, A11=A22
  ! =CCMASM_MTP_DECOUPLED: Decoupled velocity matrices A11 and A22
  ! =CCMASM_MTP_FULLTENSOR: Full-tensor matrix with A11,A12,A21,A22 independent
  INTEGER, INTENT(IN) :: cmatrixType

  ! A level-info structure specifying the matrices of the problem.
  TYPE(t_problem_lvl), INTENT(IN), TARGET :: rlevelInfo
!</input>

!<output>
  ! A block matrix that receives the basic system matrix.
  TYPE(t_matrixBlock), INTENT(OUT) :: rmatrix
!</output>

!</subroutine>

    ! local variables
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
  
    ! A t_nonlinearCCMatrix used for defining the matrix structure
    TYPE(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%rdiscretisation
    
    ! Get a pointer to the template FEM matrix.
    p_rmatrixTemplateFEM => rlevelInfo%rmatrixTemplateFEM

    ! In the global system, there are two gradient matrices B1 and B2.
    ! Get a pointer to the template structure for these.
    p_rmatrixTemplateGradient => rlevelInfo%rmatrixTemplateGradient

    ! Initialise the block matrix with default values based on
    ! the discretisation.
    CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)    
      
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
    ! For this purpose, we have to initialise a t_nonlinearCCMatrix structure
    ! which defines the shape of the matrix. We simply set the parameters
    ! of those terms which wshould appear in the matrix to a value <> 0,
    ! that's enough for the memory allocation.
    
    rnonlinearCCMatrix%dtheta = 1.0_DP   ! A velocity block
    rnonlinearCCMatrix%deta = 1.0_DP     ! A gradient block
    rnonlinearCCMatrix%dtau = 1.0_DP     ! A divergence block
    rnonlinearCCMatrix%p_rdiscretisation => rlevelInfo%rdiscretisation
    rnonlinearCCMatrix%p_rmatrixTemplateFEM => rlevelInfo%rmatrixTemplateFEM
    rnonlinearCCMatrix%p_rmatrixTemplateGradient => rlevelInfo%rmatrixTemplateGradient
    rnonlinearCCMatrix%p_rdiscretisation => rlevelInfo%rdiscretisation
    rnonlinearCCMatrix%p_rmatrixStokes => rlevelInfo%rmatrixStokes
    rnonlinearCCMatrix%p_rmatrixB1 => rlevelInfo%rmatrixB1
    rnonlinearCCMatrix%p_rmatrixB2 => rlevelInfo%rmatrixB2
    rnonlinearCCMatrix%p_rmatrixD1 => rlevelInfo%rmatrixD1
    rnonlinearCCMatrix%p_rmatrixD2 => rlevelInfo%rmatrixD2
    rnonlinearCCMatrix%p_rmatrixMass => rlevelInfo%rmatrixMass

    CALL cc_assembleMatrix (CCMASM_ALLOCMEM,cmatrixType,&
      rmatrix,rnonlinearCCMatrix)
                                  
    ! That's it, all submatrices are set up.
      
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initNonlinearLoop (rproblem,nlmin,nlmax,rvector,rrhs,&
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
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! Minimum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners.
  INTEGER, INTENT(IN) :: nlmin
  
  ! Maximum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners. This level must correspond to rvector and rrhs.
  INTEGER, INTENT(IN) :: nlmax

  ! The solution vector which is modified later during the nonlinear iteration.
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rvector

  ! The right-hand-side vector to use in the equation
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rrhs

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  CHARACTER(LEN=*), INTENT(IN) :: sname
!</input>

!<output>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Is filled with data.
  TYPE(t_ccnonlinearIteration), INTENT(OUT) :: rnonlinearIteration
!</output>

!</subroutine>

    ! local variables
    INTEGER :: ilevel
    LOGICAL :: bneumann
    TYPE(t_parlstSection), POINTER :: p_rsection

    rnonlinearIteration%NLMIN = NLMIN
    rnonlinearIteration%NLMAX = NLMAX

    ! Initialise the matrix pointers on all levels that we have to maintain.
    ALLOCATE(rnonlinearIteration%RcoreEquation(NLMIN:NLMAX))

    rnonlinearIteration%MT_OutputLevel = rproblem%MT_OutputLevel
    
    ! Get the minimum/maximum damping parameter from the parameter list, save
    ! them to the nonlinear iteration structure (which is now initialised).
    CALL parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMin', rnonlinearIteration%domegaMin, 0.0_DP)
                              
    CALL parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
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
    DO ilevel = nlmin,nlmax
      rnonlinearIteration%RcoreEquation(ilevel)%p_rdiscretisation => &
        rproblem%RlevelInfo(ilevel)%rdiscretisation

      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixStokes => &
        rproblem%RlevelInfo(ilevel)%rmatrixStokes
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB1 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB1
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB2 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB2

      ! The D1/D2 matrices are created by hand and not as pointers.
      ! This allows us to modify them if necessary.
      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilevel)%rmatrixD1,&
          rnonlinearIteration%RcoreEquation(ilevel)%rmatrixD1,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilevel)%rmatrixD2,&
          rnonlinearIteration%RcoreEquation(ilevel)%rmatrixD2,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixMass => &
        rproblem%RlevelInfo(ilevel)%rmatrixMass

      rnonlinearIteration%RcoreEquation(ilevel)%p_rtempVector => &
        rproblem%RlevelInfo(ilevel)%rtempVector

    END DO
      
    ! Clear auxiliary variables for the nonlinear iteration
    rnonlinearIteration%DresidualInit = 0.0_DP
    rnonlinearIteration%DresidualOld  = 0.0_DP
    
    CALL parlst_querysection(rproblem%rparamList, sname, p_rsection) 

    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      CALL output_line ('Cannot create nonlinear solver; no section '''//&
          TRIM(sname)//'''!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initNonlinearLoop')
      CALL sys_halt()
    END IF

    ! Get stopping criteria of the nonlinear iteration
    CALL parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnonlinearIteration%DepsNL(1), 0.1_DP)

    CALL parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnonlinearIteration%DepsNL(2), 0.1_DP)

    CALL parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnonlinearIteration%DepsNL(3), 0.1_DP)

    CALL parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnonlinearIteration%DepsNL(4), 0.1_DP)

    CALL parlst_getvalue_double (p_rsection, 'ddmpD', &
                                 rnonlinearIteration%DepsNL(5), 0.1_DP)
      
    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions. This filter chain is applied to each
    ! defect vector during the linear and nonlinear iteration.
    ALLOCATE(rnonlinearIteration%p_RfilterChain(3))
    
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter for defect vectors:
    rnonlinearIteration%p_RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! The second filter filters for boundary conditions of fictitious boundary
    ! components
    rnonlinearIteration%p_RfilterChain(2)%ifilterType = FILTER_DISCBCDEFFICT
    
    ! Do we have Neumann boundary?
    !
    ! The bhasNeumannBoundary flag of the higher level decides about that...
    bneumann = rproblem%RlevelInfo(rproblem%NLMAX)%bhasNeumannBoundary
    rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_DONOTHING
    IF (.NOT. bneumann) THEN
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_TOL20
      rnonlinearIteration%p_RfilterChain(3)%itoL20component = NDIM2D+1
    END IF
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneNonlinearLoop (rnonlinearIteration)
  
!<description>
  ! Releases memory allocated in cc_initNonlinearLoop.
  ! The routine automatically calls cc_releaseNonlinearLoop to release
  ! internal parameters.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that should be cleaned up.
  TYPE(t_ccNonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>
    
    ! Release the filter chain for the defect vectors.
    IF (ASSOCIATED(rnonlinearIteration%p_RfilterChain)) &
      DEALLOCATE(rnonlinearIteration%p_RfilterChain)
      
    IF (ASSOCIATED(rnonlinearIteration%RcoreEquation)) &
      DEALLOCATE(rnonlinearIteration%RcoreEquation)

    rnonlinearIteration%NLMIN = 0
    rnonlinearIteration%NLMAX = 0

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initLinearSolver (rproblem,rnonlinearIteration,ssection)
  
!<description>
  ! This routine initialises a linear solver structure to be used
  ! for preconditioning of the nonlinear defect.
!</description>

!<input>
  ! Name of the section in the parameter list that contains the parameters
  ! of the linear solver.
  CHARACTER(LEN=*), INTENT(IN) :: ssection
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! The nonlinear iteration structure to which a preconditioner should
  ! be initialised. The parameters for the linear solver are written to
  ! the rpreconditioner substructure.
  ! The rprecSpecials structure is configured according to the linear
  ! solver.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_parlstSection), POINTER :: p_rsection
    TYPE(t_parlist), POINTER :: p_rparamList
    INTEGER :: nlevels, ilev, nsm
    
    INTEGER :: isolverType,ismootherType,icoarseGridSolverType
    CHARACTER(LEN=SYS_STRLEN) :: sstring,ssolverSection,ssmootherSection
    CHARACTER(LEN=SYS_STRLEN) :: scoarseGridSolverSection,spreconditionerSection
    TYPE(t_linsolMGLevelInfo2), POINTER :: p_rlevelInfo
    TYPE(t_linsolNode), POINTER :: p_rpreconditioner, p_rsmoother
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!
    
    p_rparamList => rproblem%rparamList
        
    CALL parlst_querysection(p_rparamList, ssection, p_rsection) 
    
    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      CALL output_line ('Cannot create linear solver; no section '''//TRIM(ssection)//&
                        '''!', OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
      CALL sys_halt()
    END IF
    
    ! Get the parameters that configure the solver type
    
    CALL parlst_getvalue_int (p_rsection, 'isolverType', isolverType, 1)
    CALL parlst_getvalue_int (p_rsection, 'ismootherType', ismootherType, 3)
    CALL parlst_getvalue_int (p_rsection, 'icoarseGridSolverType', &
        icoarseGridSolverType, 1)
        
    rnonlinearIteration%rprecSpecials%isolverType = isolverType
    rnonlinearIteration%rprecSpecials%ismootherType = ismootherType
    rnonlinearIteration%rprecSpecials%icoarseGridSolverType = icoarseGridSolverType

    CALL parlst_getvalue_string (p_rsection, 'ssolverSection', sstring,'')
    READ (sstring,*) ssolverSection
    CALL parlst_getvalue_string (p_rsection, 'ssmootherSection', sstring,'')
    READ (sstring,*) ssmootherSection
    CALL parlst_getvalue_string (p_rsection, 'scoarseGridSolverSection', sstring,'')
    READ (sstring,*) scoarseGridSolverSection
    
    ! Which type of solver do we have?
    
    SELECT CASE (isolverType)
    
    CASE (0)
    
      ! This is the UMFPACK solver. Very easy to initialise. No parameters at all.
      CALL linsol_initUMFPACK4 (p_rsolverNode)
    
    CASE (1)
    
      ! Multigrid solver. This is a little bit harder.
      !
      ! In a first step, initialise the main solver node for all our levels.
      nlevels = rnonlinearIteration%NLMAX - rnonlinearIteration%NLMIN + 1
      
      CALL linsol_initMultigrid2 (p_rsolverNode,nlevels,&
          rnonlinearIteration%p_RfilterChain)
      
      ! Manually trim the coarse grid correction in Multigrid to multiply the 
      ! pressure equation with -1. This (un)symmetrises the operator and gives
      ! much better convergence rates.
      CALL cgcor_release(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
      CALL cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,NDIM2D+1)
      p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(3) &
          = -1.0_DP
      
      ! Init standard solver parameters and extended multigrid parameters
      ! from the DAT file.
      CALL linsolinit_initParams (p_rsolverNode,p_rparamList,ssolverSection,&
          LINSOL_ALG_UNDEFINED)
      CALL linsolinit_initParams (p_rsolverNode,p_rparamList,ssolverSection,&
          LINSOL_ALG_MULTIGRID2)
          
      ! Ok, now we have to initialise all levels. First, we create a coarse
      ! grid solver and configure it.
      CALL linsol_getMultigridLevel2 (p_rsolverNode,1,p_rlevelInfo)
      
      SELECT CASE (icoarseGridSolverType)
      CASE (0)
        ! UMFPACK coarse grid solver. Easy.
        CALL linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
        
      CASE (1)
        ! Defect correction with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        CALL linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
        
        CALL parlst_getvalue_string_direct (p_rparamList, scoarseGridSolverSection, &
            'spreconditionerSection', sstring, '')
        READ (sstring,*) spreconditionerSection
        CALL linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        CALL linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        CALL linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rnonlinearIteration%p_RfilterChain)
        CALL linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        CALL linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rpreconditioner%calgorithm)
        
      CASE (2)
        ! BiCGSTab with diagonal VANKA preconditioning.
        !
        ! Create VANKA and initialise it with the parameters from the DAT file.
        CALL linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DNAVST)
        
        CALL parlst_getvalue_string (p_rparamList, scoarseGridSolverSection, &
           'spreconditionerSection', sstring, '')
        READ (sstring,*) spreconditionerSection
        CALL linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,LINSOL_ALG_UNDEFINED)
        CALL linsolinit_initParams (p_rpreconditioner,p_rparamList,&
            spreconditionerSection,p_rpreconditioner%calgorithm)
        
        ! Create the defect correction solver, attach VANKA as preconditioner.
        CALL linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditioner,&
            rnonlinearIteration%p_RfilterChain)
        CALL linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,LINSOL_ALG_UNDEFINED)
        CALL linsolinit_initParams (p_rlevelInfo%p_rcoarseGridSolver,p_rparamList,&
            scoarseGridSolverSection,p_rlevelInfo%p_rcoarseGridSolver%calgorithm)
        
      END SELECT
      
      ! Now after the coarse grid solver is done, we turn to the smoothers
      ! on all levels. Their initialisation is similar to the coarse grid
      ! solver. Note that we use the same smoother on all levels, for 
      ! presmoothing as well as for postsmoothing.
      
      DO ilev = 2,nlevels

        ! Initialise the smoothers.
        SELECT CASE (ismootherType)
        
        CASE (0:5)

          NULLIFY(p_rsmoother)
        
          ! This is some kind of VANKA smoother. Initialise the correct one.
          SELECT CASE (ismootherType)
          CASE (0)
            CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERAL)
          CASE (1)
            CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_GENERALDIRECT)
          CASE (2)
            CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)
          CASE (3)
            CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVSTDIRECT)
          CASE (4)
            CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVST)
          CASE (5)
            CALL linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DFNAVSTDIRECT)
          END SELECT
          
          ! Initialise the parameters -- if there are any.
          CALL linsolinit_initParams (p_rsmoother,p_rparamList,&
              ssmootherSection,LINSOL_ALG_UNDEFINED)
          CALL linsolinit_initParams (p_rsmoother,p_rparamList,&
              ssmootherSection,p_rsmoother%calgorithm)
          
          ! Convert to a smoother with a defined number of smoothing steps.
          CALL parlst_getvalue_int (p_rparamList, ssmootherSection, &
                    'nsmoothingSteps', nsm, 4)
          CALL linsol_convertToSmoother (p_rsmoother,nsm)
          
          ! Put the smoother into the level info structure as presmoother
          ! and postsmoother
          CALL linsol_getMultigridLevel2 (p_rsolverNode,ilev,p_rlevelInfo)
          p_rlevelInfo%p_rpresmoother => p_rsmoother
          p_rlevelInfo%p_rpostsmoother => p_rsmoother
          
          ! Set up the interlevel projection structure on all levels
          p_rlevelInfo%rinterlevelProjection = &
            rnonlinearIteration%rpreconditioner%p_rprojection
          
        END SELECT
      
      END DO

      ! Get information about adaptive matrix generation from INI/DAT files
      CALL parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
          'iAdaptiveMatrix', rnonlinearIteration%rprecSpecials%iadaptiveMatrices, 0)
                                
      CALL parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
          'dAdMatThreshold', rnonlinearIteration%rprecSpecials%dAdMatThreshold, 20.0_DP)

    END SELECT    

    ! Put the final solver node to the preconditioner structure.
    rnonlinearIteration%rpreconditioner%p_rsolverNode => p_rsolverNode

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneLinearSolver (rnonlinearIteration)
  
!<description>
  ! Releases information of the linear solver preconditioner from the
  ! structure rnonlinearIteration. Cleans up what was configured
  ! in cc_initLinearSolver.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure to be cleaned up.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! Release the solver and all subsolvers.
    CALL linsol_releaseSolver(rnonlinearIteration%rpreconditioner%p_rsolverNode)

  END SUBROUTINE


  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initPreconditioner (rproblem,rnonlinearIteration,rvector,rrhs)
  
!<description>
  ! This routine prepares the preconditioner that us used during the
  ! nonlinear iteration. The structure rpreconditioner will be initialised
  ! based on the information in rproblem.
  ! Necessary variables will be added to the collection structure in
  ! rproblem\%rcollection to be available in the callback routines.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! This is configured according to the preconditioner as specified in
  ! the DAT files.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: NLMIN,NLMAX
    INTEGER :: i
    INTEGER(PREC_VECIDX) :: imaxmem
    CHARACTER(LEN=PARLST_MLDATA) :: ssolverName,sstring

    ! At first, ask the parameters in the INI/DAT file which type of 
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC2D-NONLINEAR', &
        'itypePreconditioning', &
        rnonlinearIteration%rpreconditioner%ctypePreconditioning, 1)
    
    SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_NONE)
      ! No preconditioner
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX
      
      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      CALL parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      IF (sstring .NE. '') READ (sstring,*) ssolverName
      IF (ssolverName .EQ. '') THEN
        CALL output_line ('No linear subsolver!', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initPreconditioner')
        CALL sys_halt()
      END IF
                                    
      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      ALLOCATE(rnonlinearIteration%rpreconditioner%p_rprojection)
      CALL mlprj_initProjectionVec (&
          rnonlinearIteration%rpreconditioner%p_rprojection,rrhs)
      
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      CALL cc_getProlRest (rnonlinearIteration%rpreconditioner%p_rprojection, &
          rproblem%rparamList,  'CC-PROLREST')
      
      ! Initialise the linear solver as configured in the DAT file.
      CALL cc_initLinearSolver (rproblem,rnonlinearIteration,ssolverName)
      
      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...

      imaxmem = 0
      DO i=NLMIN+1,NLMAX
        ! Pass the system metrices on the coarse/fine grid to 
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = MAX(imaxmem,mlprj_getTempMemoryDirect (&
            rnonlinearIteration%rpreconditioner%p_rprojection,&
            rproblem%RlevelInfo(i-1)% &
              rdiscretisation%RspatialDiscr(1:rrhs%nblocks),&
            rproblem%RlevelInfo(i)% &
              rdiscretisation%RspatialDiscr(1:rrhs%nblocks)))
      END DO
      
      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      ALLOCATE(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      CALL lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc,&
                                MAX(imaxmem,rrhs%NEQ),.FALSE.)
      
      ! Set up a second temporary vector that we need for calculating
      ! the optimal defect correction.
      ALLOCATE(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      CALL lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2,&
                                rrhs%NEQ,.FALSE.,ST_DOUBLE)
      
      ! Initialise the matrices.
      CALL cc_updatePreconditioner (rproblem,rnonlinearIteration,&
          rvector,rrhs,.TRUE.,.TRUE.)
      
    CASE DEFAULT
      
      ! Unknown preconditioner
      CALL output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initPreconditioner')
      CALL sys_halt()
      
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_updatePreconditioner (rproblem,rnonlinearIteration,&
      rvector,rrhs,binit,bstructuralUpdate)
  
!<description>
  ! This routine has to be called whenever the system matrices change.
  ! It initialises (depending on the system matrices) the matrices of the
  ! preconditioner or performs an update of them.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs

  ! First initialisation.
  ! Has to be set to TRUE on the first call. Initialises the preconditioner
  ! with the structure of the matrices.
  LOGICAL, INTENT(IN) :: binit

  ! Whether the structure of the system matrices is new.
  ! This variable has to be set to TRUE whenever there was a structure in 
  ! the system matrices. This reinitialises the linear solver.
  LOGICAL, INTENT(IN) :: bstructuralUpdate
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Preconditioner data is saved here.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: NLMIN,NLMAX
    INTEGER :: i,cmatrixType
    CHARACTER(LEN=PARLST_MLDATA) :: sstring,snewton
    LOGICAL :: btranspose
    TYPE(t_timer) :: rtimer

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
  
    ! A pointer to the matrix of the preconditioner
    TYPE(t_matrixBlock), POINTER :: p_rmatrixPreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(:), POINTER :: Rmatrices
    
    ! Pointer to the template FEM matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixTempateFEM
    
    SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_NONE)
      ! No preconditioner
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
    
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Adjust the special preconditioner parameters to the current situation.
      CALL cc_adjustPrecSpecials (rproblem,rnonlinearIteration,&
          rnonlinearIteration%rprecSpecials)
      
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX
      
      ! Initialise the preconditioner matrices on all levels.
      DO i=NLMIN,NLMAX
      
        ! Prepare the preconditioner matrices level i. 
        IF (binit .OR. bstructuralUpdate) THEN

          ! What type of matrix do we have? Is it a 'primitive' matrix
          ! (with A11 and A22 being the same) or a 'full-tensor matrix'
          ! with Aij all independent from each other
          IF ((rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
              CCPREC_NEWTON) .OR. &
              (rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
              CCPREC_NEWTONDYNAMIC)) THEN
            cmatrixType = CCMASM_MTP_FULLTENSOR
          ELSE
            cmatrixType = CCMASM_MTP_AUTOMATIC
          END IF

          ! Allocate memory for the matrix or release the existing matrix,
          ! if there's a structural update.
          IF (binit) THEN
            ALLOCATE(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
          ELSE
            CALL lsysbl_releaseMatrix(&
                rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
          END IF
          
          ! Allocate memory for the basic submatrices.
          ! The B1^T and B2^T matrices are saved as 'virtual transpose' of B1 and
          ! B2 by default, we must change them later if necessary.
          CALL cc_allocSystemMatrix (rproblem,rproblem%RlevelInfo(i),cmatrixType,&
              rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
              
          ! Attach boundary conditions
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner%p_rdiscreteBC &
              => rproblem%RlevelInfo(i)%p_rdiscreteBC
          
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner%p_rdiscreteBCfict &
              => rproblem%RlevelInfo(i)%p_rdiscreteFBC
        END IF
        
        ! On the current level, set up a global preconditioner matrix.
        p_rmatrixPreconditioner => &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner
      
        ! ----------------------------------------------------
        ! Should the linear solver use the adaptive Newton matrix?
            
        IF (rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
          CCPREC_NEWTONDYNAMIC) THEN
          
          ! We have even the extended, dynamic Newton as preconditioner.
          ! Put the parameters for the extended Newton from the DAT file
          ! into the Adaptive-Newton configuration block.
          
          CALL parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                      'spreconditionerAdaptiveNewton', sstring, '')
          snewton = ''
          IF (sstring .NE. '') READ (sstring,*) snewton
          IF (snewton .NE. '') THEN
            ! Initialise the parameters of the adaptive Newton
            CALL parlst_getvalue_int (rproblem%rparamList, snewton, &
                'nminFixPointIterations', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%nminFixPointIterations, 0)

            CALL parlst_getvalue_int (rproblem%rparamList, snewton, &
                'nmaxFixPointIterations', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%nmaxFixPointIterations, 999)

            CALL parlst_getvalue_double (rproblem%rparamList, snewton, &
                'depsAbsNewton', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%depsAbsNewton, 1E-5_DP)

            CALL parlst_getvalue_double (rproblem%rparamList, snewton, &
                'depsRelNewton', rnonlinearIteration%rpreconditioner% &
                radaptiveNewton%depsRelNewton, 1E99_DP)
          END IF
          
        END IF
          
        ! We add a zero diagonal matrix to the pressure block. This matrix
        ! is only used under rare circumstances, e.g. if we have a pure
        ! Dirichlet problem on that level and the solver does not support
        ! filtering. In this case, this matrix is used to ensure definiteness.
        !
        ! We ignore the scaling factor here as we only want to ensure that there's
        ! space available for the matrix.
        IF (.NOT. lsysbl_isSubmatrixPresent(p_rmatrixPreconditioner,3,3,.TRUE.)) THEN
          CALL lsyssc_createDiagMatrixStruc (p_rmatrixPreconditioner%RmatrixBlock(3,3),&
              p_rmatrixPreconditioner%RmatrixBlock(1,3)%NCOLS,LSYSSC_MATRIX9)
          CALL lsyssc_allocEmptyMatrix(&
              p_rmatrixPreconditioner%RmatrixBlock(3,3),LSYSSC_SETM_ZERO)
          p_rmatrixPreconditioner%RmatrixBlock(3,3)%dscaleFactor = 0.0_DP
        END IF
        
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
        IF (bstructuralUpdate) THEN
          btranspose = .FALSE.
          SELECT CASE (rnonlinearIteration%rprecSpecials%isolverType)
          CASE (0)
            ! UMFPACK solver.
            IF (i .EQ. NLMAX) THEN
              ! UMFPACK solver. Tweak the matrix on the max. level.
              ! Ignore the other levels.
              btranspose = .TRUE.
              
              IF (rnonlinearIteration%rprecSpecials%bpressureGloballyIndefinite) THEN
                ! Activate the 3,3-block, UMFPACK needs it in case the pressure
                ! is globally indefinite.
                p_rmatrixPreconditioner%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
              END IF
            END IF

          CASE (1)
            ! Multigrid solver. Treat the matrix at the coarse level if there's
            ! UMFPACK chosen as coarse grid solver.
            IF (i .EQ. NLMIN) THEN
            
              IF (rnonlinearIteration%rprecSpecials%icoarseGridSolverType .EQ. 0) THEN
                btranspose = .TRUE.
                
                IF (rnonlinearIteration%rprecSpecials%bpressureGloballyIndefinite) THEN
                  ! Activate the 3,3-block, UMFPACK needs it in case the pressure
                  ! is globally indefinite.
                  p_rmatrixPreconditioner%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
                END IF
                
              ELSE
                
                ! Tweak the matrix if the preconditioner needs transposed matrices.
                !
                ! Currently not implemented, as there is no configuration where a
                ! preconditioner needs transposed matrices...
                
              END IF
              
            ELSE
            
              ! On the other levels, tweak the matrix if the general VANKA is
              ! chosen as smoother; it needs transposed matrices.
              IF ((rnonlinearIteration%rprecSpecials%ismootherType .EQ. 0) .OR. &
                  (rnonlinearIteration%rprecSpecials%ismootherType .EQ. 1)) THEN
                btranspose = .TRUE.
              END IF              
              
            END IF

          END SELECT
          
          IF (btranspose) THEN
            ! Check if B1^T and B2^T are virtually transposed saved. If that's
            ! the case, create real-transposed matrices from them.
            if (iand(rnonlinearIteration%RcoreEquation(i)%rmatrixD1%imatrixSpec,&
                LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
            
              ! Transpose back and re-transpose; this allocates memory.
              call lsyssc_transposeMatrixInSitu(&
                  rnonlinearIteration%RcoreEquation(i)%rmatrixD1,LSYSSC_TR_VIRTUAL)
              call lsyssc_transposeMatrixInSitu(&
                  rnonlinearIteration%RcoreEquation(i)%rmatrixD1,LSYSSC_TR_ALL)
            
            end if
            
            ! The same for B2:
            if (iand(rnonlinearIteration%RcoreEquation(i)%rmatrixD2%imatrixSpec,&
                LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
            
              ! Transpose back and re-transpose; this allocates memory.
              call lsyssc_transposeMatrixInSitu(&
                  rnonlinearIteration%RcoreEquation(i)%rmatrixD2,LSYSSC_TR_VIRTUAL)
              call lsyssc_transposeMatrixInSitu(&
                  rnonlinearIteration%RcoreEquation(i)%rmatrixD2,LSYSSC_TR_ALL)
            
            end if
            
            ! B1^T and B2^T have the same structure, release unnecessary memory.
            CALL lsyssc_duplicateMatrix (&
                rnonlinearIteration%RcoreEquation(i)%rmatrixD1,&
                rnonlinearIteration%RcoreEquation(i)%rmatrixD2,&
                LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
          END IF

        END IF
        
      END DO
      
      CALL stat_clearTimer(rtimer)
      CALL stat_startTimer(rtimer)
      
      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
      ALLOCATE(Rmatrices(1:NLMAX))
      DO i=NLMIN,NLMAX
        CALL lsysbl_duplicateMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
          Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      END DO
      
      CALL linsol_setMatrices(&
          rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(NLMIN:NLMAX))
          
      ! The solver got the matrices; clean up Rmatrices, it was only of temporary
      ! nature...
      DO i=NLMIN,NLMAX
        CALL lsysbl_releaseMatrix (Rmatrices(i))
      END DO
      DEALLOCATE(Rmatrices)
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      IF (binit) THEN
      
        CALL linsol_initStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
            ierror)
        IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
        
      ELSE IF (bstructuralUpdate) THEN
      
        CALL linsol_updateStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
            ierror)
        IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
        
      END IF
      
      ! Gather statistics
      CALL stat_stopTimer(rtimer)
      rproblem%rstatistics%dtimeLinearSolverFactorisation = &
        rproblem%rstatistics%dtimeLinearSolverFactorisation + rtimer%delapsedReal
      
    CASE DEFAULT
      
      ! Unknown preconditioner
      CALL output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_updatePreconditioner')
      CALL sys_halt()
      
    END SELECT

  END SUBROUTINE

! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_releasePreconditioner (rnonlinearIteration)
  
!<description>
  ! This routine releases the preconditioner for the nonlinear iteration
  ! which was prepared in cc_initPreconditioner. Memory is released
  ! from heap.
!</description>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! The preconditioner data is removed from that,
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i

    ! Which preconditioner do we have?    
    SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_NONE)
      ! No preconditioning
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the preconditioner matrix on every level
      DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        CALL lsysbl_releaseMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        DEALLOCATE(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        
        ! Also release B1^T and B2^T everywhere where they exist.
        CALL lsyssc_releaseMatrix (rnonlinearIteration%RcoreEquation(i)%rmatrixD1)
        CALL lsyssc_releaseMatrix (rnonlinearIteration%RcoreEquation(i)%rmatrixD2)
        
      END DO
      
      ! Release the temporary vector(s)
      CALL lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      DEALLOCATE(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      
      CALL lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      DEALLOCATE(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      
      ! Clean up data about the projection etc.
      CALL mlprj_doneProjection(rnonlinearIteration%rpreconditioner%p_rprojection)
      DEALLOCATE(rnonlinearIteration%rpreconditioner%p_rprojection)

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      CALL linsol_doneStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode)
      CALL cc_doneLinearSolver (rnonlinearIteration)
      
    CASE DEFAULT
      
      ! Unknown preconditioner
      CALL output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_releasePreconditioner')
      CALL sys_halt()
      
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_getProlRest (rprojection, rparamList, sname)
  
!<description>
  ! Initialises an existing interlevel projection structure rprojection
  ! with parameters from the INI/DAT files. sname is the section in the
  ! parameter list containing parameters about prolongation restriction.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  TYPE(t_parlist), INTENT(IN) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the prolongation/restriction.
  CHARACTER(LEN=*), INTENT(IN) :: sname
!</input>

!<output>
  ! An interlevel projection block structure containing an initial
  ! configuration of prolongation/restriction. The structure is modified
  ! according to the parameters in the INI/DAT file(s).
  TYPE(t_interlevelProjectionBlock), INTENT(INOUT) :: rprojection
!</output>

!</subroutine>

    ! local variables
    TYPE(t_parlstSection), POINTER :: p_rsection
    INTEGER :: i1
    REAL(DP) :: d1

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    CALL parlst_querysection(rparamList, sname, p_rsection) 

    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      ! We use the default configuration; stop here.
      RETURN
    END IF
    
    ! Now take a look which parameters appear in that section.

    ! Prolongation/restriction order for velocity components
    CALL parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
    
    IF (i1 .NE. -1) THEN
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
    END IF

    ! Prolongation/restriction order for pressure
    CALL parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
    
    IF (i1 .NE. -1) THEN
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,NDIM2D+1)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%iinterpolationOrder = i1
    END IF
    
    ! Prolongation/restriction variant for velocity components
    ! in case of Q1~ discretisation
    CALL parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
    
    IF (i1 .NE. -1) THEN
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1
    END IF
    
    ! Aspect-ratio indicator in case of Q1~ discretisation
    ! with extended prolongation/restriction
    CALL parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
    
    IF (i1 .NE. 1) THEN
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
    END IF

    ! Aspect-ratio bound for switching to constant prolongation/restriction
    ! in case of Q1~ discretisation with extended prolongation/restriction
    CALL parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
    
    IF (d1 .NE. 20.0_DP) THEN
      rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
    END IF

  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_adjustPrecSpecials (rproblem,rnonlinearIteration,rprecSpecials)
  
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
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! Nonlinar iteration structure saving data about the actual configuration
  ! of the core equation.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
!</inputoutput>

!<inputoutput>
  ! Nonlinear iteration structure. 
  ! The t_ccPreconditionerSpecials substructure that receives information how to 
  ! finally assembly the matrices such that everything in the callback routines 
  ! will work.
  TYPE(t_ccPreconditionerSpecials), INTENT(INOUT) :: rprecSpecials
!</inputoutput>

!</subroutine>

    ! Check if we have Neumann boundary components. If not, the matrices
    ! may have to be changed, depending on the solver.
    rprecSpecials%bpressureGloballyIndefinite = &
      .NOT. rproblem%RlevelInfo(rproblem%NLMAX)%bhasNeumannBoundary

  END SUBROUTINE

END MODULE
