!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2boundary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the definition of the analytic boundary conditions
!# as well as discretisation routines for the boundary.
!#
!# The following files can be found here:
!#
!# 1.) c2d2_initAnalyticBC
!#     -> Initialise analytic boundary conditions
!#
!# 2.) c2d2_initDiscreteBC
!#     -> Discretise analytic boundary conditions, create discrete
!#        boundary conditions
!#
!# 3.) c2d2_updateDiscreteBC
!#     -> Reassemble discrete boundary conditions
!#
!# 4.) c2d2_doneBC
!#     -> Release discrete and analytic boundary conditions
!#
!# 5.) c2d2_implementBC
!#     -> Implement discrete boundary conditions into solution/RHS vectors
!#        and matrix on finest level
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2boundary

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
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmediumm2boundarydef
  USE cc2dmedium_callback
  
  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initAnalyticBC (rproblem)
  
!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    INTEGER :: i
  
    ! Initialise the boundary condition by parsing the parameter files.
    ! This initialises rproblem%p_rboundaryConditions by parsing DAT-file 
    ! parameters in the parameter list rproblem%rparamList
    CALL c2d2_parseBDconditions (rproblem)

    ! Initialise the boundary conditions of fictitious boundary components
    CALL c2d2_parseFBDconditions (rproblem)    

    ! Install the analytic boundary conditions into all discretisation
    ! structures on all levels.
    DO i=rproblem%NLMIN,rproblem%NLMAX
      
      ! Ask the problem structure to give us the discretisation structure...
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! and inform the discretisation which analytic boundary conditions to use:
      p_rdiscretisation%p_rboundaryConditions => rproblem%p_rboundaryConditions

    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initDiscreteBC (rproblem,rvector,rrhs)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! A vector structure for the solution vector. The discrete BC structures are 
  ! attached to that.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! A vector structure for the RHS vector. The discrete BC structures are 
  ! attached to that.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

  ! A pointer to the system matrix and the RHS vector as well as 
  ! the discretisation
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

  ! Pointer to structure for saving discrete BC's:
  TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  TYPE(t_discreteFBC), POINTER :: p_rdiscreteFBC
    
    DO i=rproblem%NLMIN,rproblem%NLMAX
    
      ! Get our velocity matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscretisation
      
      ! For the discrete problem, we need a discrete version of the above
      ! boundary conditions. So we have to discretise them.
      ! The following routine gives back p_rdiscreteBC, a pointer to a
      ! discrete version of the boundary conditions. Remark that
      ! the pointer has to be nullified before calling the routine,
      ! otherwise, the routine tries to update the boundary conditions
      ! in p_rdiscreteBC!
      ! getBoundaryValues is a callback routine that specifies the
      ! values on the boundary. We pass our collection structure as well
      ! to this routine, so the callback routine has access to everything what is
      ! in the collection.
      !
      ! On maximum level, discretrise everything. On lower level, discretise
      ! only for the implementation into the matrices and defect vector. 
      ! That's enough, as the lower levels are only used for preconditioning 
      ! of defect vectors.
      
      NULLIFY(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      IF (i .EQ. rproblem%NLMAX) THEN
        CALL bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                .FALSE.,getBoundaryValues, &
                                rproblem%rcollection)
      ELSE
        CALL bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                .FALSE.,getBoundaryValues, &
                                rproblem%rcollection,BCASM_DISCFORDEFMAT)
      END IF
                
      ! The same way, discretise the fictitious boundary conditions and hang
      ! them in into the matrix/vectors.
      NULLIFY(rproblem%RlevelInfo(i)%p_rdiscreteFBC)
      IF (i .EQ. rproblem%NLMAX) THEN
        CALL bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,.FALSE., &
                                  getBoundaryValuesFBC,rproblem%rcollection)
      ELSE
        CALL bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,.FALSE., &
                                  getBoundaryValuesFBC,rproblem%rcollection,&
                                  BCASM_DISCFORDEFMAT)
      END IF

      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
      ! Also hang in the boundary conditions into the temporary vector that is
      ! used for the creation of solutions on lower levels.
      ! This allows us to filter this vector when we create it.
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBC => p_rdiscreteBC
      
      ! The same for the fictitious boudary boundary conditions.
      p_rdiscreteFBC => rproblem%RlevelInfo(i)%p_rdiscreteFBC
      p_rmatrix%p_rdiscreteBCfict => p_rdiscreteFBC
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBCfict => p_rdiscreteFBC
      
    END DO

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteBC
    
    rrhs%p_rdiscreteBC => p_rdiscreteBC
    rvector%p_rdiscreteBC => p_rdiscreteBC

    ! The same with the fictitious boundary BC's
    p_rdiscreteFBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteFBC
    rrhs%p_rdiscreteBCfict => p_rdiscreteFBC
    rvector%p_rdiscreteBCfict => p_rdiscreteFBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_updateDiscreteBC (rproblem, bforce)
  
!<description>
  ! This updates the discrete version of the boundary conditions. The BC's
  ! are reassembled according to the current situation of the simulation.
!</description>

!<input>
  ! Normally, only those BC's are update that are marked as 'nonstatic'
  ! (e.g. for time-dependent simulation). When bforce=TRUE, all discrete
  ! BC's are rebuild, which might be necessary in optimisation environments
  ! when testing multiple stationary configurations.
  LOGICAL, INTENT(IN) :: bforce
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    DO i=rproblem%NLMIN,rproblem%NLMAX
    
      ! Get our velocity matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscretisation
      
      ! For the discrete problem, we need a discrete version of the above
      ! boundary conditions. So we have to discretise them.
      ! The following routine gives back p_rdiscreteBC, a pointer to a
      ! discrete version of the boundary conditions. Remark that
      ! the pointer has to be nullified before calling the routine,
      ! otherwise, the routine tries to update the boundary conditions
      ! in p_rdiscreteBC!
      ! getBoundaryValues is a callback routine that specifies the
      ! values on the boundary. We pass our collection structure as well
      ! to this routine, so the callback routine has access to everything what is
      ! in the collection.
      !
      ! On maximum level, discretrise everything. On lower level, discretise
      ! only for the implementation into the matrices and defect vector. 
      ! That's enough, as the lower levels are only used for preconditioning 
      ! of defect vectors.
      !
      ! Note that we do not NULLIFY the pointer to the discrete BC's here.
      ! That's the only difference to the initDiscreteBC routine above,
      ! as this instructs bcasm_discretiseBC to update nonstatic BC's.
      
      IF (i .EQ. rproblem%NLMAX) THEN
        CALL bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                bforce,getBoundaryValues, &
                                rproblem%rcollection)
      ELSE
        CALL bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                bforce,getBoundaryValues, &
                                rproblem%rcollection,BCASM_DISCFORDEFMAT)
      END IF
                
      ! The same way, discretise the fictitious boundary conditions and hang
      ! them in into the matrix/vectors.
      IF (i .EQ. rproblem%NLMAX) THEN
        CALL bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,bforce, &
                                  getBoundaryValuesFBC,rproblem%rcollection)
      ELSE
        CALL bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,bforce, &
                                  getBoundaryValuesFBC,rproblem%rcollection,&
                                  BCASM_DISCFORDEFMAT)
      END IF

    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_implementBC (rproblem,rvector,rrhs,bmatrices,&
      bsolvector,brhsvector)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
  ! Implements boundary conditions into the system matrix on every level
  ! defined by rproblem.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A vector structure for the solution vector. The discrete BC's are implemented
  ! into that.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! A vector structure for the RHS vector. The discrete BC's are implamented into that.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
  
  ! Whether to implement the BC's into the system matrix/matrices or not.
  LOGICAL, INTENT(IN) :: bmatrices

  ! Whether to implement the BC's into the solution vector or not
  LOGICAL, INTENT(IN) :: bsolvector

  ! Whether to implement the BC's into the solution vector or not
  LOGICAL, INTENT(IN) :: brhsvector
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,ilvmax
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  
  ! A pointer to the system matrix and the RHS vector as well as 
  ! the discretisation
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%NLMAX
    
    IF (bsolvector) THEN
    
      ! Implement discrete boundary conditions into solution vector by
      ! filtering the vector.
      CALL vecfil_discreteBCsol (rvector)
      
      ! Implement discrete boundary conditions of fictitioous boundary comnponents
      ! into solution vector by filtering the vector.
      CALL vecfil_discreteFBCsol (rvector)
      
    END IF
    
    IF (brhsvector) THEN
    
      ! Implement pressure drop boundary conditions into RHS vector
      ! if there are any.
      CALL vecfil_discreteNLPDropBCrhs (rrhs)
      
      ! Implement discrete boundary conditions into RHS vector by 
      ! filtering the vector.
      CALL vecfil_discreteBCrhs (rrhs)

      ! Implement discrete boundary conditions of fictitious boundary components
      ! into RHS vector by filtering the vector.
      CALL vecfil_discreteFBCrhs (rrhs)

    END IF
    
    IF (bmatrices) THEN
    
      ! Implement discrete boundary conditions into the matrices on all 
      ! levels, too.
      ! In fact, this modifies the B-matrices. The A-matrices are overwritten
      ! later and must then be modified again!
      DO i=rproblem%NLMIN ,rproblem%NLMAX
        p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
        CALL matfil_discreteBC (p_rmatrix)  ! standard boundary conditions
        CALL matfil_discreteFBC (p_rmatrix)  ! fictitious boundary boundary conditions
      END DO
      
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  TYPE(t_fparser), POINTER :: p_rparser

    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Release our discrete version of the boundary conditions
      CALL bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%p_rdiscreteBC)
      
      ! as well as the discrete version of the BC's for fictitious boundaries
      CALL bcasm_releaseDiscreteFBC (rproblem%RlevelInfo(i)%p_rdiscreteFBC)

      ! ...and also the corresponding analytic description.
      CALL bcond_doneBC (rproblem%p_rboundaryConditions)
    END DO
    
    ! Remove the Neumann flag from the collection
    CALL collct_deletevalue (rproblem%rcollection, 'INEUMANN')
    
    ! Remove the identifier that we should use BDC from the DAT file
    CALL collct_deletevalue (rproblem%rcollection, 'BDCFROMDAT')

    ! Release the parser object with all the expressions to be evaluated
    ! on the boundary.
    p_rparser => collct_getvalue_pars (rproblem%rcollection, BDC_BDPARSER, &
                                       0, SEC_SBDEXPRESSIONS)
    CALL fparser_release (p_rparser)
    DEALLOCATE (p_rparser)

    ! Remove the boundary value parser from the collection
    CALL collct_deletevalue (rproblem%rcollection, BDC_BDPARSER)
    
    ! Remove the boundary-expression section we added earlier,
    ! with all their content.
    CALL collct_deletesection(rproblem%rcollection,SEC_SBDEXPRESSIONS)
    
  END SUBROUTINE

END MODULE
