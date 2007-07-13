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
!# 1.) c2d2_allocSystemMatrix
!#     -> Allocates memory for the system matrix representing the 
!#        core equation.
!#
!# 2.) c2d2_initPreconditioner
!#     -> Initialises a spatial preconditioner structure with parameters from
!#        the DAT file. This is needed for preconditioning in space.
!#     -> Extension to c2d2_createPreconditioner.
!#
!# 3.) c2d2_donePreconditioner
!#     -> Cleans up a spatial preconditioner structure initialised by
!#        c2d2_initPreconditioner. 
!#     -> Extension to c2d2_releasePreconditioner.
!#
!# 4.) c2d2_configPreconditioner
!#     -> Configures a preconditioner according to the actual situation
!#
!# 5.) c2d2_updatePreconditioner
!#     -> Updates the preconditioner if there was a change in the system matrices
!#
!# 7.) c2d2_getProlRest
!#     -> Auxiliary routine: Set up interlevel projection structure
!#        with information from INI/DAT files
!#
!# Auxiliary routines, not to be called from outside:
!#
!# 1.) c2d2_checkAssembly
!#     -> Checks if the system matrices are compatible to the preconditioner.
!#        Set up some situation dependent 'tweak' flags for the assembly
!#        of the matrices.
!#
!# 2.) c2d2_finaliseMatrices
!#     -> Makes the system matrices compatible to the preconditioner if
!#        necessary.
!#
!# 3.) c2d2_unfinaliseMatrices 
!#     -> Reverts the changes of c2d2_finaliseMatrices and brings matrices
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

MODULE cc2dmediumm2nonlinearcoreinit

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
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmediumm2nonlinearcore
  
  IMPLICIT NONE
  

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_allocSystemMatrix (rproblem,rlevelInfo,rmatrix)
  
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
  
    ! A t_ccmatrixComponents used for defining the matrix structure
    TYPE(t_ccmatrixComponents) :: rmatrixAssembly
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%p_rdiscretisation
    
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
    ! We allocate the system matrix using the c2d2_assembleMatrix routine.
    ! For this purpose, we have to initialise a t_ccmatrixComponents structure
    ! which defines the shape of the matrix. We simply set the parameters
    ! of those terms which wshould appear in the matrix to a value <> 0,
    ! that's enough for the memory allocation.
    
    rmatrixAssembly%dtheta1 = 1.0_DP   ! A velocity block
    rmatrixAssembly%deta1 = 1.0_DP     ! A gradient block
    rmatrixAssembly%dtau1 = 1.0_DP     ! A divergence block
    rmatrixAssembly%dkappa1 = 1.0_DP   ! Pressure block
    
    rmatrixAssembly%dtheta2 = 1.0_DP   ! A velocity block
    rmatrixAssembly%deta2 = 1.0_DP     ! A gradient block
    rmatrixAssembly%dtau2 = 1.0_DP     ! A divergence block
    rmatrixAssembly%dkappa2 = 1.0_DP   ! Pressure block
    
    rmatrixAssembly%p_rdiscretisation => rlevelInfo%p_rdiscretisation
    rmatrixAssembly%p_rmatrixTemplateFEM => rlevelInfo%rmatrixTemplateFEM
    rmatrixAssembly%p_rmatrixTemplateGradient => rlevelInfo%rmatrixTemplateGradient
    rmatrixAssembly%p_rmatrixStokes => rlevelInfo%rmatrixStokes
    rmatrixAssembly%p_rmatrixB1 => rlevelInfo%rmatrixB1
    rmatrixAssembly%p_rmatrixB2 => rlevelInfo%rmatrixB2
    rmatrixAssembly%p_rmatrixMass => rlevelInfo%rmatrixMass
    rmatrixAssembly%p_rmatrixIdentityPressure => rlevelInfo%rmatrixIdentityPressure

    IF (.NOT. rproblem%bdecoupledXY) THEN    
      CALL c2d2_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_AUTOMATIC,&
        rmatrix,rmatrixAssembly)
    ELSE
      CALL c2d2_assembleMatrix (CCMASM_ALLOCMEM,CCMASM_MTP_DECOUPLED,&
        rmatrix,rmatrixAssembly)
    END IF
                                  
    ! That's it, all submatrices are set up.
      
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initPreconditioner (rproblem,nlmin,nlmax,rpreconditioner)
  
!<description>
  ! Initialises the given spatial preconditioner structure rpreconditioner.
  ! Creates the structure with c2d2_createPreconditioner and saves all
  ! problem dependent parameters and matrices in it.
  !
  ! This routine initialises only the basic structure. However, it does
  ! not set/initialise the type of preconditioner (Defect corection,
  ! Newton,´...).
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! Minimum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioners.
  INTEGER, INTENT(IN) :: nlmin
  
  ! Maximum refinement level in the rproblem structure that is allowed to be used
  ! by the preconditioner. This is the level where the preconditioner is to be
  ! applied!
  INTEGER, INTENT(IN) :: nlmax
!</input>

!<output>
  ! A spatial preconditioner structure to be initialised.
  TYPE(t_ccspatialPreconditioner), INTENT(OUT) :: rpreconditioner
!</output>

!</subroutine>

    ! local variables
    INTEGER :: ilevel
    LOGICAL :: bneumann
    TYPE(t_parlstSection), POINTER :: p_rsection

    ! Basic initialisation of the nonlinenar iteration structure.
    CALL c2d2_createPreconditioner (rpreconditioner,nlmin,nlmax)    
    
    ! Assign the matrix pointers in the nonlinear iteration structure to
    ! all our matrices that we want to use.
    DO ilevel = nlmin,nlmax
      rpreconditioner%RcoreEquation(ilevel)%p_rmatrix => &
        rproblem%RlevelInfo(ilevel)%rpreallocatedSystemMatrix
        
      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixStokes => &
        rproblem%RlevelInfo(ilevel)%rmatrixStokes
        
      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixB1 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB1
        
      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixB2 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB2

      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixMass => &
        rproblem%RlevelInfo(ilevel)%rmatrixMass

      rpreconditioner%RcoreEquation(ilevel)%p_rtempVector => &
        rproblem%RlevelInfo(ilevel)%rtempVector

      rpreconditioner%RcoreEquation(ilevel)%p_rmatrixIdentityPressure => &
        rproblem%RlevelInfo(ilevel)%rmatrixIdentityPressure

    END DO
      
    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions.
    ALLOCATE(rpreconditioner%p_RfilterChain(4))
    
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter for defect vectors:
    rpreconditioner%p_RfilterChain(1)%ifilterType = &
        FILTER_DISCBCDEFREAL

    ! The second filter filters for boundary conditions of fictitious boundary
    ! components
    rpreconditioner%p_RfilterChain(2)%ifilterType = &
        FILTER_DISCBCDEFFICT
    
    ! Do we have Neumann boundary?
    bneumann = collct_getvalue_int (rproblem%rcollection, 'INEUMANN') .EQ. YES
    rpreconditioner%p_RfilterChain(3)%ifilterType = FILTER_DONOTHING
    IF (.NOT. bneumann) THEN
      ! Pure Dirichlet problem -- Neumann boundary for the pressure.
      ! Filter the pressure to avoid indefiniteness.
      rpreconditioner%p_RfilterChain(3)%ifilterType = FILTER_TOL20
      rpreconditioner%p_RfilterChain(3)%itoL20component = NDIM2D+1

      rpreconditioner%p_RfilterChain(4)%ifilterType = FILTER_TOL20
      rpreconditioner%p_RfilterChain(4)%itoL20component = 2*(NDIM2D+1)
    END IF
      
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_donePreconditioner (rpreconditioner)
  
!<description>
  ! Releases memory allocated in c2d2_initPrecoditioner..
  ! The routine automatically calls c2d2_releasePreconditioner to release
  ! internal parameters.
!</description>

!<inputoutput>
  ! A spatial preconditioner structure to be cleaned up.
  TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>
    
    ! local variables
    INTEGER :: i

    ! Which preconditioner do we have?    
    SELECT CASE (rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_NONE)
      ! No preconditioning
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the preconditioner matrix on every level
      DO i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
        CALL lsysbl_releaseMatrix ( &
          rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner)
        DEALLOCATE(rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner)
      END DO
      
      ! Release the temporary vector(s)
      CALL lsyssc_releaseVector (rpreconditioner%p_rtempVectorSc)
      DEALLOCATE(rpreconditioner%p_rtempVectorSc)
      
      ! Clean up data about the projection etc.
      CALL mlprj_doneProjection(rpreconditioner%p_rprojection)
      DEALLOCATE(rpreconditioner%p_rprojection)

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      CALL linsol_doneStructure (rpreconditioner%p_rsolverNode)
      CALL linsol_releaseSolver (rpreconditioner%p_rsolverNode)
      
    CASE DEFAULT
      
      ! Unknown preconditioner
      PRINT *,'Unknown preconditioner for nonlinear iteration!'
      STOP
      
    END SELECT

    ! Release the filter chain for the defect vectors.
    IF (ASSOCIATED(rpreconditioner%p_RfilterChain)) &
      DEALLOCATE(rpreconditioner%p_RfilterChain)
      
    CALL c2d2_releasePreconditioner (rpreconditioner)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_configPreconditioner (rproblem,rpreconditioner)
  
!<description>
  ! This routine prepares the preconditioner by means of the parameters
  ! in the DAT files and on the information in rproblem.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</input>

!<inputoutput>
  ! A spatial preconditioner structure to be initialised. Must have been
  ! created previously with initPreconditioner.
  ! This is configured according to the preconditioner as specified in
  ! the DAT files.
  TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: NLMIN,NLMAX
    INTEGER :: i,icomponents
    INTEGER(PREC_VECIDX) :: imaxmem
    CHARACTER(LEN=PARLST_MLDATA) :: ssolverName,sstring
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! At first, ask the parameters in the INI/DAT file which type of 
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC2D-NONLINEAR', &
        'itypePreconditioning', &
        rpreconditioner%ctypePreconditioning, 1)
    
    SELECT CASE (rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_NONE)
      ! No preconditioner
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rpreconditioner%NLMIN
      NLMAX = rpreconditioner%NLMAX
      
      ! Get a pointer to the discretsation structure on the level
      ! where the preconditioner should act
      p_rdiscretisation => rproblem%RlevelInfo(NLMAX)%p_rdiscretisation
      
      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      CALL parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      IF (sstring .NE. '') READ (sstring,*) ssolverName
      IF (ssolverName .EQ. '') THEN
        PRINT *,'No linear subsolver!'
        STOP
      END IF
                                    
      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      ALLOCATE(rpreconditioner%p_rprojection)
      CALL mlprj_initProjectionDiscr (rpreconditioner%p_rprojection,p_rdiscretisation)
      
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      CALL c2d2_getProlRest (rpreconditioner%p_rprojection, &
          rproblem%rparamList,  'CC-PROLREST')
      
      ! Initialise the linear subsolver using the parameters from the INI/DAT
      ! files, the prepared filter chain and the interlevel projection structure.
      ! This gives us the linear solver node rpreconditioner%p_rsolverNode
      ! which identifies the linear solver.
      CALL linsolinit_initFromFile (rpreconditioner%p_rsolverNode,&
                                    rproblem%rparamList,ssolverName,&
                                    NLMAX-NLMIN+1,&
                                    rpreconditioner%p_RfilterChain,&
                                    rpreconditioner%p_rprojection)
      
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
            rpreconditioner%p_rprojection,&
            rproblem%RlevelInfo(i-1)% &
              p_rdiscretisation%RspatialDiscretisation(1:p_rdiscretisation%ncomponents),&
            rproblem%RlevelInfo(i)% &
              p_rdiscretisation%RspatialDiscretisation(1:p_rdiscretisation%ncomponents)))
      END DO
      
      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      ALLOCATE(rpreconditioner%p_rtempVectorSc)
      CALL lsyssc_createVector (rpreconditioner%p_rtempVectorSc,&
        MAX(imaxmem,dof_igetNDofGlobBlock(p_rdiscretisation)),.FALSE.)
      
      ! Initialise the matrices.
      CALL c2d2_updatePreconditioner (rproblem,rpreconditioner,.TRUE.,.TRUE.)

    CASE DEFAULT
      
      ! Unknown preconditioner
      PRINT *,'Unknown preconditioner for nonlinear iteration!'
      STOP
      
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_updatePreconditioner (rproblem,rpreconditioner,&
      binit,bstructuralUpdate)
  
!<description>
  ! This routine has to be called whenever the system matrices change.
  ! It initialises (depending on the system matrices) the matrices of the
  ! preconditioner or performs an update of them.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

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
  ! A spatial preconditioner structure to be úpdated. Must have been
  ! created previously with initPreconditioner and configPreconditioner.
  TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: NLMIN,NLMAX
    INTEGER :: i
    CHARACTER(LEN=PARLST_MLDATA) :: sstring,snewton

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
  
    ! A pointer to the matrix of the preconditioner
    TYPE(t_matrixBlock), POINTER :: p_rmatrixPreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(NNLEV) :: Rmatrices
    
    ! Pointer to the template FEM matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixTempateFEM
    
    ! At first, ask the parameters in the INI/DAT file which type of 
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC2D-NONLINEAR', &
        'itypePreconditioning', &
        rpreconditioner%ctypePreconditioning, 1)
    
    SELECT CASE (rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_NONE)
      ! No preconditioner
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)

      IF (.NOT. binit) THEN
        ! Restore the standard matrix structure in case the matrices had been
        ! modified by c2d2_finaliseMatrices for the preconditioner -- i.e. 
        ! temporarily switch to the matrix structure that is compatible to 
        ! the discretisation.
        CALL c2d2_unfinaliseMatrices (rpreconditioner,.FALSE.)
      END IF
    
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rpreconditioner%NLMIN
      NLMAX = rpreconditioner%NLMAX
      
      IF (binit .OR. bstructuralUpdate) THEN
      
        ! Initialise the preconditioner matrices on all levels.
        DO i=NLMIN,NLMAX
        
          ! Prepare the preconditioner matrices level i. This is
          ! basically the system matrix.
          ! Create an empty matrix structure or clean up the old one.
          IF (.NOT. ASSOCIATED(rpreconditioner%RcoreEquation(i)%&
              p_rmatrixPreconditioner)) THEN
            ALLOCATE(rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner)
          ELSE
            CALL lsysbl_releaseMatrix (rpreconditioner%RcoreEquation(i)%&
                p_rmatrixPreconditioner)
          END IF
          p_rmatrixPreconditioner => &
              rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner
        
          CALL lsysbl_duplicateMatrix (&
              rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix,&
              p_rmatrixPreconditioner,&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
              
          ! ----------------------------------------------------
          ! Should the linear solver use the Newton matrix?
          IF (rpreconditioner%ctypePreconditioning .EQ. CCPREC_NEWTON) THEN
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
            
            ! If we have a Stokes problem, A12 and A21 don't exist.
            IF (rproblem%iequation .EQ. 0) THEN
            
              IF (p_rmatrixPreconditioner%RmatrixBlock(1,2)%cmatrixFormat &
                  .EQ. LSYSSC_MATRIXUNDEFINED) THEN
                  
                CALL lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                  p_rmatrixPreconditioner%RmatrixBlock(1,2), &
                  LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
                ! Allocate memory for the entries; don't initialise the memory.
                CALL lsyssc_allocEmptyMatrix (&
                    p_rmatrixPreconditioner%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
                  
              END IF

              IF (p_rmatrixPreconditioner%RmatrixBlock(2,1)%cmatrixFormat &
                  .EQ. LSYSSC_MATRIXUNDEFINED) THEN
                  
                ! Create a new matrix A21 in memory. create a new matrix
                ! using the template FEM matrix...
                CALL lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                  p_rmatrixPreconditioner%RmatrixBlock(2,1), &
                  LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
                ! Allocate memory for the entries; don't initialise the memory.
                CALL lsyssc_allocEmptyMatrix (&
                    p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)

              END IF               

            END IF

            ! A22 may share its entries with A11. If that's the case,
            ! allocate additional memory for A22, as the Newton matrix
            ! requires a separate A22!
            IF (lsyssc_isMatrixContentShared( &
                p_rmatrixPreconditioner%RmatrixBlock(1,1), &
                p_rmatrixPreconditioner%RmatrixBlock(2,2)) ) THEN
              ! Release the matrix structure. As the matrix is a copy
              ! of another one, this will clean up the structure but
              ! not release any memory.
              CALL lsyssc_releaseMatrix ( &
                  p_rmatrixPreconditioner%RmatrixBlock(2,2))

              ! Create a new matrix A21 in memory. create a new matrix
              ! using the template FEM matrix...
              CALL lsyssc_duplicateMatrix ( &
                p_rmatrixPreconditioner%RmatrixBlock(1,1), &
                p_rmatrixPreconditioner%RmatrixBlock(2,2),&
                LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
              ! ... then allocate memory for the entries; 
              ! don't initialise the memory.
              CALL lsyssc_allocEmptyMatrix (&
                  p_rmatrixPreconditioner%RmatrixBlock(2,2),&
                  LSYSSC_SETM_UNDEFINED)
            END IF
            
          END IF
        END DO
        
      END IF
      
      IF (binit) THEN
        ! Check the matrices if they are compatible to our
        ! preconditioner. If not, we later have to modify the matrices a little
        ! bit to make it compatible. 
        ! The result of this matrix analysis is saved to the rfinalAssembly structure 
        ! in rpreconditioner and allows us later to switch between these two
        ! matrix representations: Compatibility to the discretisation routines
        ! and compatibity to the preconditioner.
        ! The c2d2_checkAssembly routine below uses this information to perform
        ! the actual modification in the matrices.
        CALL c2d2_checkAssembly (rproblem,rpreconditioner)
      END IF
      ! Otherwise, checkAssembly was already called and does not have to be 
      ! called again.

      ! Using rfinalAssembly as computed above, make the matrices compatible 
      ! to our preconditioner if they are not.
      CALL c2d2_finaliseMatrices (rpreconditioner)

      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
      DO i=NLMIN,NLMAX
        CALL lsysbl_duplicateMatrix ( &
          rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner, &
          Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      END DO
      
      CALL linsol_setMatrices(&
          rpreconditioner%p_rsolverNode,Rmatrices(NLMIN:NLMAX))
          
      ! The solver got the matrices; clean up Rmatrices, it was only of temporary
      ! nature...
      DO i=NLMIN,NLMAX
        CALL lsysbl_releaseMatrix (Rmatrices(i))
      END DO
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      IF (binit) THEN
        CALL linsol_initStructure (rpreconditioner%p_rsolverNode,ierror)
        IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
      ELSE IF (bstructuralUpdate) THEN
        CALL linsol_updateStructure (rpreconditioner%p_rsolverNode,ierror)
        IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
      END IF
      
    CASE DEFAULT
      
      ! Unknown preconditioner
      PRINT *,'Unknown preconditioner for nonlinear iteration!'
      STOP
      
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_getProlRest (rprojection, rparamList, sname)
  
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

  SUBROUTINE c2d2_checkAssembly (rproblem,rpreconditioner)
  
!<description>
  ! This routine checks the matrices against an existing preconditioner.
  ! It may happen that e.g. VANCA does not like our matrices (probably
  ! they have to be saved transposed or whatever). Information about
  ! which things must be changed in the assembly to make 'everything proper'
  ! is saved to rfinalAssembly in rpreconditioner.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A spatial preconditioner structure to be úpdated.
  TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: NLMIN,NLMAX,ccompatible,iprecType
  CHARACTER(LEN=PARLST_MLDATA) :: ssolverName,sstring
  TYPE(t_interlevelProjectionBlock) :: rprojection

  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  TYPE(t_matrixBlock), DIMENSION(NNLEV) :: Rmatrices
  TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    ! At first, ask the parameters in the INI/DAT file which type of 
    ! preconditioner is to be used. 
    CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                     'itypePreconditioning', &
                                     iprecType, 1)

    SELECT CASE (iprecType)
    CASE (CCPREC_NONE)
      ! No preconditioning
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! That preconditioner is a solver for a linear system.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rpreconditioner%NLMIN
      NLMAX = rpreconditioner%NLMAX

      ! Temporarily set up the solver node for the linear subsolver.
      !
      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      CALL parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      IF (sstring .NE. '') READ (sstring,*) ssolverName
      IF (ssolverName .EQ. '') THEN
        PRINT *,'No linear subsolver!'
        STOP
      END IF
                                    
      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      CALL mlprj_initProjectionDiscr (rprojection,&
          rproblem%RlevelInfo(NLMAX)%p_rdiscretisation)
      
      ! Initialise the linear subsolver using the parameters from the INI/DAT
      ! files, the prepared filter chain and the interlevel projection structure.
      ! This gives us the linear solver node p_rsolverNode
      ! which identifies the linear solver.
      CALL linsolinit_initFromFile (p_rsolverNode,&
                                    rproblem%rparamList,ssolverName,&
                                    NLMAX-NLMIN+1,&
                                    rinterlevelProjection=rprojection)
      
      ! Check the matrices.
      !
      ! We copy our matrices to a big matrix array and transfer that
      ! to the setMatrices routines. This intitialises then the matrices
      ! on all levels according to that array.
      Rmatrices(NLMIN:NLMAX) = rproblem%RlevelInfo(NLMIN:NLMAX)%rpreallocatedSystemMatrix
      CALL linsol_matricesCompatible(p_rsolverNode, &
          Rmatrices(NLMIN:NLMAX),ccompatible)
      
      SELECT CASE (ccompatible)
      CASE (LINSOL_COMP_OK) ! nothing to do
        rpreconditioner%rfinalAssembly%iBmatricesTransposed = NO
      CASE (LINSOL_COMP_ERRTRANSPOSED)
        ! The B-matrices must be assembled in a transposed way. Remember that.
        rpreconditioner%rfinalAssembly%iBmatricesTransposed = YES
      CASE DEFAULT
        PRINT *,'Preconditioner incompatible to the matrices. Don''t know why!?!'
        STOP
      END SELECT
      
      ! Release the solver node again.
      CALL linsol_releaseSolver (p_rsolverNode)
      
      ! We also don't need the temporary projection structure anymore.
      CALL mlprj_doneProjection (rprojection)
      
    END SELECT

    ! Add information about adaptive matrix generation from INI/DAT files
    ! to the collection.
    CALL parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'iAdaptiveMatrix', rpreconditioner%rfinalAssembly%iadaptiveMatrices, 0)
                              
    CALL parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'dAdMatThreshold', rpreconditioner%rfinalAssembly%dAdMatThreshold, 20.0_DP)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_finaliseMatrices (rpreconditioner)
  
!<description>
  ! This routine performs final assembly tasks to the matrices such that they
  ! are compatible to the preconditioner.
  ! It may happen that e.g. VANCA does not like our matrices (probably
  ! they have to be saved transposed or whatever). In that case, we
  ! have to make slight modifications to our matrices in order to 
  ! make them compatible.
!</description>

!<inputoutput>
  ! A spatial preconditioner structure that contains pointers to the
  ! matrices which are to be used.
  TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ilev,NLMIN,NLMAX
    TYPE(t_matrixBlock), POINTER :: p_rmatrix

    IF (rpreconditioner%rfinalAssembly%iBmatricesTransposed .EQ. YES) THEN
      ! There is usually a VANCA subsolver in the main linear solver which
      ! cannot deal with our virtually transposed matrices. So we make
      ! a copy of B1/B2 on every level and really transpose them.

      ! Which levels have we to take care of during the solution process?
      NLMIN = rpreconditioner%NLMIN
      NLMAX = rpreconditioner%NLMAX

      ! Loop through the levels, transpose the B-matrices  
      DO ilev=NLMIN,NLMAX
        ! Get the matrix of the preconditioner
        p_rmatrix => rpreconditioner%RcoreEquation(ilev)%p_rmatrix
        
        ! Release the old B1/B2 matrices from the system matrix. This
        ! does not release any memory, as the content of the matrix
        ! is saved elsewhere.
        CALL lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,1))
        CALL lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,2))
        
        ! Transpose B1/B2, write result to the system matrix.
        CALL lsyssc_transposeMatrix (&
            rpreconditioner%RcoreEquation(ilev)%p_rmatrixB1, &
            p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_ALL)

        CALL lsyssc_transposeMatrix (&
            rpreconditioner%RcoreEquation(ilev)%p_rmatrixB2, &
            p_rmatrix%RmatrixBlock(3,2),LSYSSC_TR_ALL)
                                    
        ! Release the memory that was allocated for the B2 structure by
        ! the matrix-transpose routine. 
        ! Replace the structure of B2 by that of B1; more precisely,
        ! share the structure. We can do this as we know that B1 and B2
        ! have exactly the same structure!
        CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_REMOVE,LSYSSC_DUP_IGNORE)
        CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
      
        ! -----------------------------------------------------------
        ! Optimal control problem extension
        ! -----------------------------------------------------------
        ! Put a reference to the transposed B-matrices also to the
        ! dual equation part.
        CALL lsyssc_duplicateMatrix ( &
            p_rmatrix%RmatrixBlock(3,1),p_rmatrix%RmatrixBlock(6,4), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        CALL lsyssc_duplicateMatrix ( &
            p_rmatrix%RmatrixBlock(3,2),p_rmatrix%RmatrixBlock(6,5), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      END DO
      
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_unfinaliseMatrices (rpreconditioner,bdata)
  
!<description>
  ! Reverts the changes that were done by c2d2_finaliseMatrices, brings
  ! matrices back to their 'standard' form. Necessary if multiple simulations
  ! are performed, e.g. in a time-dependent simulation: The assembly
  ! routines expect matrices in the standard form and cannot work with
  ! 'specially tuned' matrices which are used in the actual solving/
  ! preconditioning process.
!</description>

!<inputoutput>
  ! A spatial preconditioner structure that contains pointers to the
  ! matrices which are to be used.
  TYPE(t_ccspatialPreconditioner), INTENT(INOUT) :: rpreconditioner
!</inputoutput>

!<input>
  ! TRUE  = restore the original matrices in their whole -- structure and data
  ! FALSE = Ignore the data, only restore the original matrix structure.
  !         Used to save time if the matrix content is thrown away in the
  !         reassembling process anyway.
  LOGICAL, INTENT(IN) :: bdata
!</input>

!</subroutine>

    ! local variables
    INTEGER :: ilev,NLMIN,NLMAX
    TYPE(t_matrixBlock), POINTER :: p_rmatrix

    ! Which levels have we to take care of during the solution process?
    NLMIN = rpreconditioner%NLMIN
    NLMAX = rpreconditioner%NLMAX

    ! Loop through the levels, transpose the B-matrices  
    DO ilev=NLMIN,NLMAX
      ! Get the matrix of the preconditioner
      p_rmatrix => rpreconditioner%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
      IF ((rpreconditioner%rfinalAssembly%iBmatricesTransposed .EQ. YES) .AND. &
          (IAND(p_rmatrix%RmatrixBlock(3,1)%imatrixSpec,&
                LSYSSC_MSPEC_TRANSPOSED) .NE. 0)) THEN
                
        ! There is usually a VANCA subsolver in the main linear solver which
        ! cannot deal with our virtually transposed matrices. 
        ! The B1/B2 matrices are transposed -- so we re-transpose them
        ! to get the original matrices.

        ! Release the old B1/B2 matrices from the system matrix. This
        ! does not release any memory, as the content of the matrix
        ! is saved elsewhere.
        CALL lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,1))
        CALL lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,2))
        
        ! Transpose B1/B2, write result to the system matrix.
        IF (bdata) THEN
          ! Structure and data
          CALL lsyssc_transposeMatrix (&
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixB1, &
              p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_ALL)

          CALL lsyssc_transposeMatrix (&
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixB2, &
              p_rmatrix%RmatrixBlock(3,2),LSYSSC_TR_ALL)
        ELSE
          ! Only the structure; content gets invalid.
          CALL lsyssc_transposeMatrix (&
              rpreconditioner%RcoreEquation(ilev)%p_rmatrixB1, &
              p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_STRUCTURE)

          ! No change has to be done to the B2^T block in the global system 
          ! matrix since that one will later share the same structure as B1^T!
          ! So the following command is commented out and should 
          ! not be commented in!
          ! CALL lsyssc_transposeMatrix (&
          !     rpreconditioner%RcoreEquation(ilev)%p_rmatrixB2, &
          !     p_rmatrix%RmatrixBlock(3,2),LSYSSC_TR_STRUCTURE)
        END IF
                                    
        ! Release the memory that was allocated for the B2 structure by
        ! the matrix-transpose routine. 
        ! Replace the structure of B2 by that of B1; more precisely,
        ! share the structure. We can do this as we know that B1 and B2
        ! have exactly the same structure!
        CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_REMOVE,LSYSSC_DUP_IGNORE)
        CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(3,1),&
                p_rmatrix%RmatrixBlock(3,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
      
        ! -----------------------------------------------------------
        ! Optimal control problem extension
        ! -----------------------------------------------------------
        ! Put a reference to the non-transposed B-matrices also to the
        ! dual equation part.
        CALL lsyssc_duplicateMatrix ( &
            p_rmatrix%RmatrixBlock(3,1),p_rmatrix%RmatrixBlock(6,4), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        CALL lsyssc_duplicateMatrix ( &
            p_rmatrix%RmatrixBlock(3,2),p_rmatrix%RmatrixBlock(6,5), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      END IF

    END DO
      
  END SUBROUTINE

END MODULE
