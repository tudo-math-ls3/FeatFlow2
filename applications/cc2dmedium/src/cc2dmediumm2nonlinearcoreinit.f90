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
!# 1.) c2d2_initNonlinearLoop
!#     -> Initialises a 'nonlinear iteration structure' with parameters from
!#        the DAT file. This is needed for solving the core equation.
!#
!# 2.) c2d2_preparePreconditioner
!#     -> Prepare preconditioner of nonlinear iteration
!#
!# 3.) c2d2_releasePreconditioner
!#     -> Clean up preconditioner of nonlinear iteration
!#
!# 4.) c2d2_getProlRest
!#     -> Auxiliary routine: Set up interlevel projection structure
!#        with information from INI/DAT files
!#
!# 5.) c2d2_checkAssembly
!#     -> Checks if the system matrices are compatible to the preconditioner.
!#        Set up some situation dependent 'tweak' flags for the assembly
!#        of the matrices.
!#
!# 6.) c2d2_finaliseMatrices
!#     -> Makes the system matrices compatible to the preconditioner if
!#        necessary.
!#
!# 7.) c2d2_unfinaliseMatrices 
!#     -> Reverts the changes of c2d2_finaliseMatrices and brings matrices
!#        into their original form.
!#
!# The module works in tight relationship to cc2dmediumm2nonlinearcode.
!# cc2dmediumm2nonlinearcodeinit provides the routines to initialise
!# preconditioner and important structures using the problem related
!# structure. This module cc2dmediumm2nonlinearcore on the other hand
!# contains only the 'main' worker routines that do the work of the
!# nonlinear iteration.
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
  USE cc2dmediumm2discretisation
  
  IMPLICIT NONE
  

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initNonlinearLoop (rproblem,rvector,rrhs,rnonlinearIteration,&
      sname)
  
!<description>
  ! Initialises the given nonlinear iteration structure rnonlinearIteration.
  ! Creates the structure with c2d2_createNonlinearLoop and saves all
  ! problem dependent parameters and matrices in it.
  ! Note: This does not initialise the preconditioner in that structure!
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(IN), TARGET :: rproblem

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
    TYPE(t_parlstSection), POINTER :: p_rsection

    ! Basic initialisation of the nonlinenar iteration structure.
    CALL c2d2_createNonlinearLoop (rnonlinearIteration,rproblem%NLMIN,rproblem%NLMAX)    
    
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
    rnonlinearIteration%rfinalAssembly%iBmatricesTransposed = NO
    rnonlinearIteration%rfinalAssembly%iadaptiveMatrices = 0
    
    ! Assign the matrix pointers in the nonlinear iteration structure to
    ! all our matrices that we want to use.
    DO ilevel = rproblem%NLMIN,rproblem%NLMAX
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrix => &
        rproblem%RlevelInfo(ilevel)%rmatrix
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixStokes => &
        rproblem%RlevelInfo(ilevel)%rmatrixStokes
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB1 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB1
        
      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB2 => &
        rproblem%RlevelInfo(ilevel)%rmatrixB2

      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixMass => &
        rproblem%RlevelInfo(ilevel)%rmatrixMass

      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixVelocityCoupling12 => &
        rproblem%RlevelInfo(ilevel)%rmatrixVelocityCoupling12

      rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixVelocityCoupling21 => &
        rproblem%RlevelInfo(ilevel)%rmatrixVelocityCoupling21
        
    END DO
      
    ! Clear auxiliary variables for the nonlinear iteration
    rnonlinearIteration%DresidualInit = 0.0_DP
    rnonlinearIteration%DresidualOld  = 0.0_DP
    
    CALL parlst_querysection(rproblem%rparamList, sname, p_rsection) 

    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      PRINT *,'Nonlinear solver not available; no section '''&
              //TRIM(sname)//'''!'
      STOP
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
      
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_preparePreconditioner (rproblem,rnonlinearIteration,&
      rvector,rrhs)
  
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
    CHARACTER(LEN=PARLST_MLDATA) :: ssolverName,sstring,snewton
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscr
    LOGICAL :: bneumann

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
  
    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix

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
        rnonlinearIteration%rpreconditioner%ctypePreconditioning, 1)
    
    SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rproblem%NLMIN
      NLMAX = rproblem%NLMAX
      
      ! Get our right hand side / solution / matrix on the finest
      ! level from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(NLMAX)%rmatrix
      
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
      ALLOCATE(rnonlinearIteration%rpreconditioner%p_rprojection)
      CALL mlprj_initProjectionVec (&
          rnonlinearIteration%rpreconditioner%p_rprojection,rrhs)
      
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      CALL c2d2_getProlRest (rnonlinearIteration%rpreconditioner%p_rprojection, &
          rproblem%rparamList,  'CC-PROLREST')
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ALLOCATE(rnonlinearIteration%rpreconditioner%p_RfilterChain(3))
      
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      rnonlinearIteration%rpreconditioner%p_RfilterChain(1)%ifilterType = &
          FILTER_DISCBCDEFREAL

      ! The second filter filters for boundary conditions of fictitious boundary
      ! components
      rnonlinearIteration%rpreconditioner%p_RfilterChain(2)%ifilterType = &
          FILTER_DISCBCDEFFICT
      
      ! Do we have Neumann boundary?
      bneumann = collct_getvalue_int (rproblem%rcollection, 'INEUMANN') .EQ. YES
      rnonlinearIteration%rpreconditioner%p_RfilterChain(3)%ifilterType = &
          FILTER_DONOTHING
      IF (.NOT. bneumann) THEN
        ! Pure Dirichlet problem -- Neumann boundary for the pressure.
        ! Filter the pressure to avoid indefiniteness.
        rnonlinearIteration%rpreconditioner%p_RfilterChain(3)%ifilterType = &
            FILTER_TOL20
        rnonlinearIteration%rpreconditioner%p_RfilterChain(3)%itoL20component = &
            NDIM2D+1
      END IF
      
      ! Initialise the linear subsolver using the parameters from the INI/DAT
      ! files, the prepared filter chain and the interlevel projection structure.
      ! This gives us the linear solver node rpreconditioner%p_rsolverNode
      ! which identifies the linear solver.
      CALL linsolinit_initFromFile (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
                                    rproblem%rparamList,ssolverName,&
                                    NLMAX-NLMIN+1,&
                                    rnonlinearIteration%rpreconditioner%p_RfilterChain,&
                                    rnonlinearIteration%rpreconditioner%p_rprojection)
      
      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...

      imaxmem = 0
      DO i=NLMIN+1,NLMAX
        ! Pass the system metrices on the coarse/fine grid to 
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = MAX(imaxmem,mlprj_getTempMemoryMat (&
            rnonlinearIteration%rpreconditioner%p_rprojection,&
            rproblem%RlevelInfo(i-1)%rmatrix,&
            rproblem%RlevelInfo(i)%rmatrix))
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
      
      ! Initialise the preconditioner matrices on all levels.
      DO i=NLMIN,NLMAX
      
        ! Prepare the preconditioner matrices level i. This is
        ! basically the system matrix...
        ALLOCATE(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        p_rmatrixPreconditioner => &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner
      
        CALL lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
            p_rmatrixPreconditioner,&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
            
        ! ----------------------------------------------------
        ! Should the linear solver use the Newton matrix?
        IF ((rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
            CCPREC_NEWTON) .OR. &
            (rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
            CCPREC_NEWTONDYNAMIC)) THEN
          ! That mens, our preconditioner matrix must look like
          !
          !  A11  A12  B1
          !  A21  A22  B2
          !  B1^T B2^T 0
          !
          ! With A12, A21, A11, A22 independent of each other!
          ! Do we have that case? If not, we have to allocate memory 
          ! for these matrices.
          p_rmatrixTempateFEM => rproblem%RlevelInfo(i)%rmatrixTemplateFEM
          
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
             
          ELSE
          
            ! A21 may share its entries with A12. If that's the case,
            ! allocate additional memory for A21!
            IF (lsyssc_isMatrixContentShared( &
                p_rmatrixPreconditioner%RmatrixBlock(1,2), &
                p_rmatrixPreconditioner%RmatrixBlock(2,1)) ) THEN

              ! Release the matrix structure. As the matrix is a copy
              ! of another one, this will clean up the structure but
              ! not release any memory.
              CALL lsyssc_releaseMatrix ( &
                  p_rmatrixPreconditioner%RmatrixBlock(2,1))

              ! Create a new matrix A21 in memory. create a new matrix
              ! using the template FEM matrix...
              CALL lsyssc_duplicateMatrix (p_rmatrixTempateFEM, &
                p_rmatrixPreconditioner%RmatrixBlock(2,1),&
                LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
              ! ... then allocate memory for the entries; 
              ! don't initialise the memory.
              CALL lsyssc_allocEmptyMatrix (&
                  p_rmatrixPreconditioner%RmatrixBlock(2,1),&
                  LSYSSC_SETM_UNDEFINED)
                
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
          
          ! ----------------------------------------------------
          ! Should we even use the adaptive Newton?
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
          
        END IF
      END DO
      
      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
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
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      CALL linsol_initStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
          ierror)
      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
      
      ! Switch off adaptive matrix generation if our discretisation is not a 
      ! uniform Q1~ discretisation - the matrix restriction does not support
      ! other cases.
      p_rdiscr => rproblem%RlevelInfo(NLMAX)%p_rdiscretisation%RspatialDiscretisation(1)
      IF ((p_rdiscr%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
          (elem_getPrimaryElement(p_rdiscr%RelementDistribution(1)%itrialElement) &
              .NE. EL_Q1T)) THEN
        i = 0
      END IF
      
    CASE DEFAULT
      
      ! Unknown preconditioner
      PRINT *,'Unknown preconditioner for nonlinear iteration!'
      STOP
      
    END SELECT

  END SUBROUTINE

! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_releasePreconditioner (rnonlinearIteration)
  
!<description>
  ! This routine releases the preconditioner for the nonlinear iteration
  ! which was prepared in c2d2_preparePreconditioner. Memory is released
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
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the preconditioner matrix on every level
      DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        CALL lsysbl_releaseMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        DEALLOCATE(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
      END DO
      
      ! Release the temporary vector(s)
      CALL lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      DEALLOCATE(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      
      CALL lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      DEALLOCATE(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      
      ! Clean up data about the projection etc.
      DEALLOCATE(rnonlinearIteration%rpreconditioner%p_RfilterChain)
      CALL mlprj_doneProjection(rnonlinearIteration%rpreconditioner%p_rprojection)
      DEALLOCATE(rnonlinearIteration%rpreconditioner%p_rprojection)

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      CALL linsol_doneStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode)
      CALL linsol_releaseSolver (rnonlinearIteration%rpreconditioner%p_rsolverNode)
      
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
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolEX3Yvariant  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestEX3Yvariant  = i1
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

  SUBROUTINE c2d2_checkAssembly (rproblem,rrhs,rfinalAssembly)
  
!<description>
  ! This routine checks the matrices against an existing preconditioner.
  ! It may happen that e.g. VANCA does not like our matrices (probably
  ! they have to be saved transposed or whatever). Information about
  ! which things must be changed in the assembly to make 'everything proper'
  ! is saved to rfinalAssembly.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!<input>
  ! The right-hand-side vector to use in the equation
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
!</input>

!<inputoutput>
  ! Nonlinear iteration structure. 
  ! The t_ccFinalAssemblyInfo substructure that receives information how to 
  ! finally assembly the matrices such that everything in the callback routines 
  ! will work.
  TYPE(t_ccFinalAssemblyInfo), INTENT(INOUT) :: rfinalAssembly
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
    CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
      ! That preconditioner is a solver for a linear system.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rproblem%NLMIN
      NLMAX = rproblem%NLMAX

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
      CALL mlprj_initProjectionVec (rprojection,rrhs)
      
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
      Rmatrices(NLMIN:NLMAX) = rproblem%RlevelInfo(NLMIN:NLMAX)%rmatrix
      CALL linsol_matricesCompatible(p_rsolverNode, &
          Rmatrices(NLMIN:NLMAX),ccompatible)
      
      SELECT CASE (ccompatible)
      CASE (LINSOL_COMP_OK) ! nothing to do
        rfinalAssembly%iBmatricesTransposed = NO
      CASE (LINSOL_COMP_ERRTRANSPOSED)
        ! The B-matrices must be assembled in a transposed way. Remember that.
        rfinalAssembly%iBmatricesTransposed = YES
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
        'iAdaptiveMatrix', rfinalAssembly%iadaptiveMatrices, 0)
                              
    CALL parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'dAdMatThreshold', rfinalAssembly%dAdMatThreshold, 20.0_DP)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_finaliseMatrices (rnonlinearIteration)
  
!<description>
  ! This routine performs final assembly tasks to the matrices such that they
  ! are compatible to the preconditioner.
  ! It may happen that e.g. VANCA does not like our matrices (probably
  ! they have to be saved transposed or whatever). In that case, we
  ! have to make slight modifications to our matrices in order to 
  ! make them compatible.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that contains pointers to the
  ! matrices which are to be used.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT), TARGET :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ilev,NLMIN,NLMAX
    TYPE(t_matrixBlock), POINTER :: p_rmatrix

    IF (rnonlinearIteration%rfinalAssembly%iBmatricesTransposed .EQ. YES) THEN
      ! There is usually a VANCA subsolver in the main linear solver which
      ! cannot deal with our virtually transposed matrices. So we make
      ! a copy of B1/B2 on every level and really transpose them.

      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX

      ! Loop through the levels, transpose the B-matrices  
      DO ilev=NLMIN,NLMAX
        ! Get the matrix of the preconditioner
        p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrix
        
        ! Release the old B1/B2 matrices from the system matrix. This
        ! does not release any memory, as the content of the matrix
        ! is saved elsewhere.
        CALL lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,1))
        CALL lsyssc_releaseMatrix (p_rmatrix%RmatrixBlock(3,2))
        
        ! Transpose B1/B2, write result to the system matrix.
        CALL lsyssc_transposeMatrix (&
            rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1, &
            p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_ALL)

        CALL lsyssc_transposeMatrix (&
            rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2, &
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
      
      END DO
      
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_unfinaliseMatrices (rnonlinearIteration,rfinalAssembly,bdata)
  
!<description>
  ! Reverts the changes that were done by c2d2_finaliseMatrices, brings
  ! matrices back to their 'standard' form. Necessary if multiple simulations
  ! are performed, e.g. in a time-dependent simulation: The assembly
  ! routines expect matrices in the standard form and cannot work with
  ! 'specially tuned' matrices which are used in the actual solving/
  ! preconditioning process.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that contains pointers to the
  ! matrices which are to be used.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT), TARGET :: rnonlinearIteration
!</inputoutput>

!<input>
  ! The t_ccFinalAssemblyInfo structure that receives information how to finally 
  ! assembly the matrices such that everything in the callback routines will 
  ! work. Must be set up with c2d2_checkAssembly above.
  TYPE(t_ccFinalAssemblyInfo), INTENT(IN) :: rfinalAssembly
  
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
    NLMIN = rnonlinearIteration%NLMIN
    NLMAX = rnonlinearIteration%NLMAX

    ! Loop through the levels, transpose the B-matrices  
    DO ilev=NLMIN,NLMAX
      ! Get the matrix of the preconditioner
      p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
      IF ((rfinalAssembly%iBmatricesTransposed .EQ. YES) .AND. &
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
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1, &
              p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_ALL)

          CALL lsyssc_transposeMatrix (&
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2, &
              p_rmatrix%RmatrixBlock(3,2),LSYSSC_TR_ALL)
        ELSE
          ! Only the structure; content gets invalid.
          CALL lsyssc_transposeMatrix (&
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1, &
              p_rmatrix%RmatrixBlock(3,1),LSYSSC_TR_STRUCTURE)

          ! No change has to be done to the B2^T block in the global system 
          ! matrix since that one will later share the same structure as B1^T!
          ! So the following command is commented out and should 
          ! not be commented in!
          ! CALL lsyssc_transposeMatrix (&
          !     rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2, &
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
      
      END IF

    END DO
      
  END SUBROUTINE

END MODULE
