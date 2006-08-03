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
!# 3.) c2d2_doneBC
!#     -> Release discrete and analytic boundary conditions
!#
!# 4.) c2d2_implementBC
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
  USE cc2dmedium_callback
  
  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

  SUBROUTINE c2d2_parseBDconditions (rproblem)

!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
  !
  ! For this purpose, the parameters in the [BDEXPRESSIONS] and [BDCONDITIONS]
  ! sections of the DAT files (saved in rproblem\%rparamList) are evaluated.
!</description>
  
!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    LOGICAL :: bNeumann
    INTEGER :: i,ityp,ivalue,ibdComponent,isegment,iintervalEnds,ibctyp,icount
    INTEGER :: iexptyp
    INTEGER, DIMENSION(2) :: IminIndex,imaxIndex
    REAL(DP) :: dvalue,dpar1,dpar2
    !INTEGER, DIMENSION(2) :: IvelComp
    CHARACTER(LEN=PARLST_MLDATA) :: cstr,cexpr,sbdex1,sbdex2
    CHARACTER(LEN=PARLST_MLNAME) :: cname
    
    ! A local collection we use for storing named parameters during the
    ! parsing process.
    TYPE(t_collection) :: rcoll
    
    ! Triangulation on currently highest level.
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! A pointer to the section with the expressions and the boundary conditions
    TYPE(t_parlstSection), POINTER :: p_rsection,p_rbdcond
    
    ! Get the domain from the problem structure
    p_rboundary => rproblem%p_rboundary
    
    ! Get the triangulation on the highest level
    p_rtriangulation => rproblem%RlevelInfo(rproblem%NLMAX)%p_rtriangulation
    
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, we need the analytic description of the boundary conditions.
    ! Initialise a structure for boundary conditions, which accepts this,
    ! on the heap.
    !
    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
    !
    ! Set p_rboundaryConditions to NULL() to create a new structure on the heap.
    NULLIFY (rproblem%p_rboundaryConditions)
    CALL bcond_initBC (rproblem%p_rboundaryConditions,p_rboundary)
    
    ! Get the expression/bc sections from the bondary condition block
    CALL parlst_querysection(rproblem%rparamList, 'BDEXPRESSIONS', p_rsection) 
    CALL parlst_querysection(rproblem%rparamList, 'BDCONDITIONS', p_rbdcond) 
    
    ! Add a section to the collection that accepts the boundary expressions
    CALL collct_addsection (rproblem%rcollection, 'BDEXPRESSIONS')
    
    ! For intermediate storing of expression types, we use a local collection
    CALL collct_init (rcoll)
    
    ! Add the boundary expressions to the collection into the
    ! specified section.
    DO i=1,parlst_querysubstrings_indir (p_rsection, 'bdExpressions')
    
      CALL parlst_getvalue_string_indir (p_rsection, 'bdExpressions', cstr, '', i)
      
      ! Get the type and decide on the identifier how to save the expression.
      READ(cstr,*) cname,ityp
      
      SELECT CASE (ityp)
      CASE (0)
        ! Real-value
        READ(cstr,*) cname,ityp,dvalue
        CALL collct_setvalue_real (rproblem%rcollection, cname, dvalue, .TRUE., &
                                   0, 'BDEXPRESSIONS') 
      CASE (1)
        ! Integer-value
        READ(cstr,*) cname,ityp,ivalue
        CALL collct_setvalue_int (rproblem%rcollection, cname, ivalue, .TRUE., &
                                   0, 'BDEXPRESSIONS')
                
      CASE (2)
        ! Parabolic profile with specified maximum velocity
        READ(cstr,*) cname,ityp,dvalue
        CALL collct_setvalue_real (rproblem%rcollection, cname, dvalue, .TRUE., &
                                   0, 'BDEXPRESSIONS') 
                                   
      CASE DEFAULT
        PRINT *,'Expressions not implemented!'
        STOP

      END SELECT
      
      ! Put the type of the expression to the temporary collection section
      CALL collct_setvalue_int (rcoll, cname, ityp, .TRUE.) 
      
    END DO
    
    bNeumann = .FALSE.
    
    ! Now to the actual boundary conditions.
    ! Loop through all boundary components we have.
    DO ibdComponent = 1,boundary_igetNBoundComp(p_rboundary)
      
      ! Parse the parameter 'bdComponentX'
      WRITE (cexpr,'(I10)') ibdComponent
      cstr = 'bdComponent' // ADJUSTL(cexpr)
      
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      i = parlst_queryvalue_indir (p_rbdcond, cstr)
      IF (i .NE. 0) THEN
         ! Parameter exists. Get the values in there.
         DO isegment = 1,parlst_querysubstrings_indir (p_rbdcond, cstr)
           
           CALL parlst_getvalue_string_fetch (p_rbdcond, &
                                              i, cstr, isubstring=isegment)
           ! Which type is that?
           READ(cstr,*) ityp
           
           SELECT CASE (ityp)
           CASE (0)
             READ(cstr,*) ityp,dpar2,iintervalEnds,ibctyp
           CASE (1)
             READ(cstr,*) ityp,ivalue,iintervalEnds,ibctyp
             dpar2 = REAL(ivalue,DP)
           END SELECT                                   
                         
           ! Form a boundary condition segment that covers that boundary part
           IF (dpar2 .GE. dpar1) THEN
             
             rboundaryRegion%dminParam = dpar1
             rboundaryRegion%dmaxParam = dpar2
             rboundaryRegion%iboundCompIdx = ibdComponent
             rboundaryRegion%dmaxParamBC = &
               boundary_dgetMaxParVal(p_rboundary, ibdComponent)
             rboundaryRegion%iproperties = iintervalEnds
             
             ! Now, which type of BC is to be created?
             SELECT CASE (ibctyp)
             
             CASE (0)
               ! Usually there's Neumann boundary in this region, but we can't be 
               ! sure. Check if, on the highest level, there's at least one edge
               ! of the triangulation belonging to the boundary. If yes, we
               ! have found Neumann boundary. If no, the segment is just too
               ! small to be considered as Neumann boundary.
               
               CALL bcasm_getEdgesInBCregion (p_rtriangulation,p_rboundary,rboundaryRegion, &
                                              IminIndex,ImaxIndex,icount)
               IF (icount .GT. 0) bNeumann = .TRUE.
             
             CASE (1)
               ! Simple Dirichlet boundary
               ! Read the line again, get the expressions for X- and Y-velocity
               READ(cstr,*) ityp,dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
               
               ! For any string <> '', create the appropriate Dirichlet boundary
               ! condition and add it to the list of boundary conditions.
               ! Set the string tag of the boundary condition to the name of
               ! the expression to evaluate.
               !
               ! Set the integer tag of the structure to the type of the expression.
               ! That way the callback routine for discretising the BC's can quickly
               ! check wht is to evaluate.
               ! If the type is a double precision value, set the double precision
               ! tag to that value so it can be evaluated quicker that taking the
               ! value from the collection.
               IF (sbdex1 .NE. '') THEN
                 CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,&
                                                    1,rboundaryRegion,p_rbcRegion)
                 p_rbcRegion%stag = sbdex1

                 iexptyp = collct_getvalue_int (rcoll, sbdex1)
                 p_rbcRegion%itag = iexptyp

                 IF ((iexptyp .EQ. 0) .OR. (iexptyp .EQ. 2)) THEN
                   p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                      sbdex1, 0, 'BDEXPRESSIONS')
                 END IF

               END IF
               
               IF (sbdex2 .NE. '') THEN
                 CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,&
                                                    2,rboundaryRegion,p_rbcRegion)
                 p_rbcRegion%stag = sbdex2

                 iexptyp = collct_getvalue_int (rcoll, sbdex2)
                 p_rbcRegion%itag = iexptyp

                 IF ((iexptyp .EQ. 0) .OR. (iexptyp .EQ. 2)) THEN
                   p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                      sbdex2, 0, 'BDEXPRESSIONS')
                 END IF

               END IF
               
               ! If we have no-slip boundary conditions, the X- and Y-velocity
               ! matrices are 'decoupled' as they are modified differently
               ! by the boundary-conditions implementation filter!
               ! In this case, change rproblem%bdecoupledXY from FALSE to TRUE
               ! to indicate that.
               IF ( ((sbdex1 .EQ. '') .AND. (sbdex2 .NE. '')) .OR.&
                    ((sbdex1 .NE. '') .AND. (sbdex2 .EQ. '')) ) THEN
                 rproblem%bdecoupledXY = .TRUE.
               END IF
             
             CASE DEFAULT
               PRINT *,'Unknown boundary condition'
               STOP
             END SELECT
             
             ! Move on to the next parameter value
             dpar1 = dpar2
             
           END IF
                                              
         END DO
      
      ELSE
        ! There is no parameter configuring the boundary condition on that 
        ! component - so we have Neumann boundary there.
        bNeumann = .TRUE.
      END IF
      
    END DO

    ! Add to the collection whether there is Neumann boundary or not.    
    IF (bNeumann) THEN
      CALL collct_setvalue_int (rproblem%rcollection, 'INEUMANN', YES, .TRUE.)
    ELSE
      CALL collct_setvalue_int (rproblem%rcollection, 'INEUMANN', NO, .TRUE.)
    END IF
    
    ! Remove the temporary collection from memory.
    CALL collct_done (rcoll)

  END SUBROUTINE  
  
  ! ***************************************************************************

  SUBROUTINE c2d2_parseFBDconditions (rproblem)

!<description>
  ! This initialises the analytic boundary conditions for fictitious boundary
  ! components in the problem and saves them to the problem structure.
!</description>
  
!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! An identifier array for the equations to be tackled when discretising
    ! boundary conditions on fictitious boundary components
    INTEGER, DIMENSION(2) :: Iequations
    
    ! A structure identifying the fictitious boundary component
    TYPE(t_fictBoundaryRegion) :: rfictBoundaryRegion

    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    
    ! Add a new fictitious boundary object It should impose Dirichlet
    ! boundary conditions in the domain in the one and only solution component.
    ! We use the default initialisation of rfictBoundaryRegion and only
    ! change the name of the component.
    rfictBoundaryRegion%sname = 'CIRCLE'
    Iequations = (/1,2/)    ! 1=x, 2=y-velocity
    CALL bcond_newDirichletBConFictBD (rproblem%p_rboundaryConditions,Iequations,&
                                       rfictBoundaryRegion,p_rbcRegion)    

  END SUBROUTINE  
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initAnalyticBC (rproblem)
  
!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
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

  SUBROUTINE c2d2_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  INTEGER :: i

  ! A pointer to the system matrix and the RHS vector as well as 
  ! the discretisation
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
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
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC

    ! The same with the fictitious boundary BC's
    p_rdiscreteFBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteFBC
    p_rrhs%p_rdiscreteBCfict => p_rdiscreteFBC
    p_rvector%p_rdiscreteBCfict => p_rdiscreteFBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  INTEGER :: i,ilvmax
  
    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%NLMAX
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    ! Implement discrete boundary conditions into RHS vector by 
    ! filtering the vector.
    CALL vecfil_discreteBCrhs (p_rrhs)

    ! Implement discrete boundary conditions into solution vector by
    ! filtering the vector.
    CALL vecfil_discreteBCsol (p_rvector)
    
    ! Implement discrete boundary conditions of fictitious boundary components
    ! into RHS vector by filtering the vector.
    CALL vecfil_discreteFBCrhs (p_rrhs)

    ! Implement discrete boundary conditions of fictitioous boundary comnponents
    ! into solution vector by filtering the vector.
    CALL vecfil_discreteFBCsol (p_rvector)
    
    ! Implement discrete boundary conditions into the matrices on all 
    ! levels, too.
    ! In fact, this modifies the B-matrices. The A-matrices are overwritten
    ! later and must then be modified again!
    DO i=rproblem%NLMIN ,rproblem%NLMAX
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      CALL matfil_discreteBC (p_rmatrix)  ! standard boundary conditions
      CALL matfil_discreteFBC (p_rmatrix)  ! fictitious boundary boundary conditions
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

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
    
    ! Remove the 'BDEXPRESSIONS' we added earlier,
    ! with all their content.
    CALL collct_deletesection(rproblem%rcollection,'BDEXPRESSIONS')
    
  END SUBROUTINE

END MODULE
