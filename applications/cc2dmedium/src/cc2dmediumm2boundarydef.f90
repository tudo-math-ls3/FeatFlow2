!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2boundarydef </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for parsing boundary conditions from a DAT
!# file
!#
!# The following files can be found here:
!#
!# 1.) c2d2_parseBDconditions
!#     -> Parses the definition of the BC's from sections given by DAT files
!#        and sets up an analytical boundary condition description
!#
!# 2.) c2d2_parseFBDconditions
!#     -> Parses the definition of the BC's given by fictitious boundaries
!#        by evaluating sections given by DAT files
!#        and sets up an analytical boundary condition description
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2boundarydef

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
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Names of the sections used in the main collection">

  ! Section name of the section saving the boundary expressions, i.e. the
  ! expressions that are to be evaluated in each point on the boundary.
  CHARACTER(LEN=COLLCT_MLSECTION), PARAMETER :: SEC_SBDEXPRESSIONS = 'BDEXPRESSIONS'

  ! Name of the parser object for boundary value expressions
  CHARACTER(LEN=COLLCT_MLSECTION), PARAMETER :: BDC_BDPARSER = 'BDEXPRPARSER'
!</constantblock>

!<constantblock description="The different types of boundary conditions supported by this module">

  ! User defined expression, evaluated by calling the callback routine
  INTEGER, PARAMETER :: BDC_USERDEF = -2
  
  ! Text expression, evaluated by parsing
  INTEGER, PARAMETER :: BDC_EXPRESSION = -1
  
  ! Fixed double precision value 
  INTEGER, PARAMETER :: BDC_VALDOUBLE = 0

  ! Fixed integer value 
  INTEGER, PARAMETER :: BDC_VALINT    = 1

  ! Parabolic profile with prescribed maximum value
  INTEGER, PARAMETER :: BDC_VALPARPROFILE = 2

!<constantblock description="Variables in expressions">

  ! Basic variables that are allowed in expressions.
  ! Variables that are not defined in the actual situation are set to 0.
  !
  ! X,Y,Z - coordinate of a point (z=0 in 2D case), \\
  ! L     - local parameter value in the range [0..1],\\
  ! R     - parameter value of a boundary point, 0-1 parametrisation, \\
  ! S     - parameter value of a boundary point, arc length parametrisation, \\
  ! TIME  - current simulation time (=0 in stationary simulation)
  !
  ! Depending on the situation, this list may be extended by situation
  ! specific variables or variables that are only available at runtime.
  CHARACTER(LEN=10), DIMENSION(7), PARAMETER :: SEC_EXPRVARIABLES = &
    (/'X    ','Y    ','Z    ','L    ','R    ','S    ','TIME '/)

!</constantblock>

!</constants>

CONTAINS
  
  ! ***************************************************************************

  SUBROUTINE c2d2_parseBDconditions (rproblem)

!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
  !
  ! For this purpose, the parameters in the [BDEXPRESSIONS] and [BDCONDITIONS]
  ! sections of the DAT files (saved in rproblem\%rparamList) are evaluated.
  !
  ! A flag 'INEUMANN'==YES/NO is added to the collection that tells about
  ! whether there is Neumann boundary somewhere.
!</description>
  
!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    LOGICAL :: bNeumann
    INTEGER :: i,ityp,ivalue,ibdComponent,isegment,iintervalEnds
    INTEGER :: ibctyp,icount,iexptyp
    INTEGER, DIMENSION(2) :: IminIndex,imaxIndex
    INTEGER, DIMENSION(NDIM2D) :: IvelComp
    REAL(DP) :: dvalue,dpar1,dpar2
    CHARACTER(LEN=PARLST_MLDATA) :: cstr,cexpr,sbdex1,sbdex2
    CHARACTER(LEN=PARLST_MLNAME) :: cname
    INTEGER, DIMENSION(NDIM2D) :: IvelEqns
    
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
    
    ! A compiled expression for evaluation at runtime
    TYPE(t_fparser), POINTER :: p_rparser
    
    ! Get the domain from the problem structure
    p_rboundary => rproblem%p_rboundary
    
    ! Get the triangulation on the highest level
    p_rtriangulation => rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation
    
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
    CALL collct_addsection (rproblem%rcollection, SEC_SBDEXPRESSIONS)
    
    ! For intermediate storing of expression types, we use a local collection
    CALL collct_init (rcoll)
    
    ! Create a parser structure for as many expressions as configured
    ALLOCATE (p_rparser)
    CALL fparser_create (p_rparser,&
         parlst_querysubstrings_indir (p_rsection, 'bdExpressions'))
    
    ! Add the parser to the collection
    CALL collct_setvalue_pars (rproblem%rcollection, BDC_BDPARSER, p_rparser, &
                                .TRUE., 0, SEC_SBDEXPRESSIONS)
    
    ! Add the boundary expressions to the collection into the
    ! specified section.
    DO i=1,parlst_querysubstrings_indir (p_rsection, 'bdExpressions')
    
      CALL parlst_getvalue_string_indir (p_rsection, 'bdExpressions', cstr, '', i)
      
      ! Get the type and decide on the identifier how to save the expression.
      READ(cstr,*) cname,ityp
      
      SELECT CASE (ityp)
      CASE (BDC_USERDEF)
        ! Name of a hardcoded expression realised in the callback routines
        READ(cstr,*) cname,ityp,cexpr
        CALL collct_setvalue_string (rproblem%rcollection, cname, cexpr, .TRUE., &
                                     0, SEC_SBDEXPRESSIONS)

      CASE (BDC_EXPRESSION)
        ! General expression; not implemented yet
        READ(cstr,*) cname,ityp,cexpr
        
        ! Compile the expression; the expression gets number i
        CALL fparser_parseFunction (p_rparser, i, cexpr, SEC_EXPRVARIABLES)
        
        ! Add the number of the function in the parser object to
        ! the collection with the name of the expression
        CALL collct_setvalue_int (rproblem%rcollection, cname, i, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS)
        
      CASE (BDC_VALDOUBLE)
        ! Real-value
        READ(cstr,*) cname,ityp,dvalue
        CALL collct_setvalue_real (rproblem%rcollection, cname, dvalue, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS) 
      CASE (BDC_VALINT)
        ! Integer-value
        READ(cstr,*) cname,ityp,ivalue
        CALL collct_setvalue_int (rproblem%rcollection, cname, ivalue, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS)
                
      CASE (BDC_VALPARPROFILE)
        ! Parabolic profile with specified maximum velocity
        READ(cstr,*) cname,ityp,dvalue
        CALL collct_setvalue_real (rproblem%rcollection, cname, dvalue, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS) 
                                   
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
              
              CALL bcasm_getEdgesInBCregion (p_rtriangulation,p_rboundary,&
                                            rboundaryRegion, &
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
                p_rbcRegion%ibdrexprtype = iexptyp

                SELECT CASE (iexptyp)
                CASE (0,2)
                  ! Constant or parabolic profile
                  p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                    sbdex1, 0, SEC_SBDEXPRESSIONS)
                CASE (-1)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  p_rbcRegion%itag = collct_getvalue_int (rproblem%rcollection, &
                                     sbdex1, 0, SEC_SBDEXPRESSIONS)
                END SELECT

              END IF
              
              IF (sbdex2 .NE. '') THEN
                CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,&
                                                  2,rboundaryRegion,p_rbcRegion)
                p_rbcRegion%stag = sbdex2

                iexptyp = collct_getvalue_int (rcoll, sbdex2)
                p_rbcRegion%ibdrexprtype = iexptyp

                SELECT CASE (iexptyp)
                CASE (0,2)
                  p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                     sbdex2, 0, SEC_SBDEXPRESSIONS)
                CASE (-1)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  p_rbcRegion%itag = collct_getvalue_int (rproblem%rcollection, &
                                     sbdex2, 0, SEC_SBDEXPRESSIONS)
                END SELECT

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
              
            CASE (2)
            
              ! Pressure drop boundary conditions.
              ! Read the line again to get the actual parameters
              READ(cstr,*) ityp,dvalue,iintervalEnds,ibctyp,sbdex1
              
              ! For any string <> '', create the appropriate pressure drop boundary
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
                IvelEqns = (/1,2/)
                CALL bcond_newPressureDropBConRealBD (rproblem%p_rboundaryConditions,&
                                              IvelEqns,&
                                              rboundaryRegion,p_rbcRegion)          
                p_rbcRegion%stag = sbdex1

                iexptyp = collct_getvalue_int (rcoll, sbdex1)
                p_rbcRegion%ibdrexprtype = iexptyp

                SELECT CASE (iexptyp)
                CASE (0,2)
                  ! Constant or parabolic profile
                  p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                    sbdex1, 0, SEC_SBDEXPRESSIONS)
                CASE (-1)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  p_rbcRegion%itag = collct_getvalue_int (rproblem%rcollection, &
                                     sbdex1, 0, SEC_SBDEXPRESSIONS)
                END SELECT

              END IF
              
              
            CASE (3)
              
              ! Nonlinear slip boundary conditions.
              IvelComp = (/1,2/)
              CALL bcond_newSlipBConRealBD (rproblem%p_rboundaryConditions,&
                                            IvelComp(1:NDIM2D),&
                                            rboundaryRegion,p_rbcRegion)
            
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
    !CALL bcond_newDirichletBConFictBD (rproblem%p_rboundaryConditions,Iequations,&
    !                                   rfictBoundaryRegion,p_rbcRegion)    

  END SUBROUTINE  
  
END MODULE
