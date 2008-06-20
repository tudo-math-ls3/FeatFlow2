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
!# 1.) cc_assembleBDconditions
!#     -> Assembles the definition of the BC's from sections given by DAT files
!#        and sets up an analytical boundary condition description
!#
!# 2.) cc_assembleFBDconditions
!#     -> Assembles the definition of the BC's given by fictitious boundaries
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
  USE cc2dmedium_callback
  
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
!</constantblock>

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

!<subroutine>

  SUBROUTINE cc_assembleBDconditions (rproblem,rdiscretisation,rdiscreteBC,rcollection)

!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
  !
  ! For this purpose, the parameters in the [BDEXPRESSIONS] and [BDCONDITIONS]
  ! sections of the DAT files (saved in rproblem\%rparamList) are evaluated.
  !
  ! The bhasNeumannBoudary flag on every level is initialised according to
  ! whether there exist Neumann boundary components on the boundary or not.
!</description>
  
!<input>
  ! A discretisation structure defining the discretisation of the current
  ! level.
  TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscretisation
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! Collection structure to be passed to callback routines
  TYPE(t_collection), INTENT(INOUT), TARGET :: rcollection
  
  ! A t_discreteBC structure that receives a discretised version
  ! of the boundary boundary conditions. The structure should
  ! be empty; new BC's are simply added to the structure.
  TYPE(t_discreteBC), INTENT(INOUT) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    LOGICAL :: bNeumann
    INTEGER :: i,ityp,ivalue,ibdComponent,isegment,iintervalEnds
    INTEGER :: ibctyp,icount,iexptyp
    INTEGER, DIMENSION(2) :: IminIndex,imaxIndex
    INTEGER, DIMENSION(NDIM2D) :: IvelComp
    REAL(DP) :: dvalue,dpar1,dpar2
    CHARACTER(LEN=PARLST_MLDATA) :: cstr,cexpr,sbdex1,sbdex2,sbdex3,sbdex4
    CHARACTER(LEN=PARLST_MLNAME) :: cname
    INTEGER, DIMENSION(NDIM2D) :: IvelEqns
    TYPE(t_collection) :: rlocalCollection
    
    ! A local collection we use for storing named parameters during the
    ! parsing process.
    TYPE(t_collection) :: rcoll
    
    ! Triangulation on currently highest level.
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! A pointer to the section with the expressions and the boundary conditions
    TYPE(t_parlstSection), POINTER :: p_rsection,p_rbdcond
    
    ! A compiled expression for evaluation at runtime
    TYPE(t_fparser), TARGET :: rparser
    
    ! Get the domain from the problem structure
    p_rboundary => rproblem%rboundary
    
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
    ! Get the expression/bc sections from the bondary condition block
    CALL parlst_querysection(rproblem%rparamList, 'BDEXPRESSIONS', p_rsection) 
    CALL parlst_querysection(rproblem%rparamList, 'BDCONDITIONS', p_rbdcond) 
    
    ! For intermediate storing of expression types, we use a local collection
    CALL collct_init (rcoll)
    
    ! Add a section to the collection that accepts the boundary expressions
    CALL collct_addsection (rcoll, SEC_SBDEXPRESSIONS)
    
    ! Create a parser structure for as many expressions as configured
    CALL fparser_create (rparser,&
         parlst_querysubstrings_indir (p_rsection, 'bdExpressions'))
    
    ! Add the parser to the collection
    CALL collct_setvalue_pars (rcoll, BDC_BDPARSER, rparser, &
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
        CALL collct_setvalue_string (rcoll, cname, cexpr, .TRUE., &
                                     0, SEC_SBDEXPRESSIONS)

      CASE (BDC_EXPRESSION)
        ! General expression; not implemented yet
        READ(cstr,*) cname,ityp,cexpr
        
        ! Compile the expression; the expression gets number i
        CALL fparser_parseFunction (rparser, i, cexpr, SEC_EXPRVARIABLES)
        
        ! Add the number of the function in the parser object to
        ! the collection with the name of the expression
        CALL collct_setvalue_int (rcoll, cname, i, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS)
        
      CASE (BDC_VALDOUBLE)
        ! Real-value
        READ(cstr,*) cname,ityp,dvalue
        CALL collct_setvalue_real (rcoll, cname, dvalue, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS) 
      CASE (BDC_VALINT)
        ! Integer-value
        READ(cstr,*) cname,ityp,ivalue
        CALL collct_setvalue_int (rcoll, cname, ivalue, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS)
                
      CASE (BDC_VALPARPROFILE)
        ! Parabolic profile with specified maximum velocity
        READ(cstr,*) cname,ityp,dvalue
        CALL collct_setvalue_real (rcoll, cname, dvalue, .TRUE., &
                                   0, SEC_SBDEXPRESSIONS) 
                                   
      CASE DEFAULT
        CALL output_line ('Expressions not implemented!', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
        CALL sys_halt()

      END SELECT
      
      ! Put the type of the expression to the temporary collection section
      CALL collct_setvalue_int (rcoll, cname, ityp, .TRUE.) 
      
    END DO
    
    bNeumann = .FALSE.
    
    ! Now to the actual boundary conditions.
    !    
    ! Put some information to the quick access arrays for access
    ! in the callback routine.
    rcoll%Iquickaccess(1) = rproblem%itimedependence
    SELECT CASE (rproblem%itimedependence)
    CASE (0)
      ! Stationary simulation
      rcoll%Dquickaccess(1) = 0.0_DP
      rcoll%Dquickaccess(2) = 0.0_DP
      rcoll%Dquickaccess(3) = 0.0_DP
    CASE (1)
      rcoll%Dquickaccess(1) = rproblem%rtimedependence%dtime
      rcoll%Dquickaccess(2) = rproblem%rtimedependence%dtimeInit
      rcoll%Dquickaccess(3) = rproblem%rtimedependence%dtimeMax
    END SELECT
    
    ! DquickAccess(4:) is reserved for BC specific information.
    !    
    ! Put a link to the previous collection into that local collection.
    ! That allows us to access it or to pass it to user defined callback
    ! functions.
    rcoll%p_rnextCollection => rcollection

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
          ! Read the segment parameters
          READ(cstr,*) dpar2,iintervalEnds,ibctyp
          
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
              ! Simple Dirichlet boundary.
              ! Prescribed primal velocity; dual velocity os always zero.
              ! Read the line again, get the expressions for X- and Y-velocity
              READ(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
              
              ! For any string <> '', create the appropriate Dirichlet boundary
              ! condition and add it to the list of boundary conditions.
              !
              ! The IquickAccess array is set up as follows:
              !  IquickAccess(1) = Type of boundary condition
              !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
              !  IquickAccess(3) = expression type
              !  IquickAccess(4) = expression identifier
              !
              ! The SquickAccess array is set up as follows:
              !  SquickAccess(1) = Name of the expression
              !
              rcoll%IquickAccess(1) = ibctyp
              
              ! If the type is a double precision value, set the DquickAccess(4)
              ! to that value so it can quickly be accessed.
              IF (sbdex1 .NE. '') THEN
                ! X-velocity
                !
                ! The 2nd element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 1
                
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex1)
                rcoll%IquickAccess(3) = iexptyp
                
                ! The 1st element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex1
                
                ! Dquickaccess(3) / IquickAccess(3) saves information
                ! about the expression.
                SELECT CASE (iexptyp)
                CASE (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll, sbdex1, 0, SEC_SBDEXPRESSIONS)
                CASE (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll, sbdex1, 0, SEC_SBDEXPRESSIONS)
                END SELECT
              
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,1,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)
                    
              END IF
              
              IF (sbdex2 .NE. '') THEN
              
                ! Y-velocity
                !
                ! The 1st element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 2
                
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex2)
                rcoll%IquickAccess(3) = iexptyp
                
                ! The 1st element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex2
                
                ! Dquickaccess(4) / IquickAccess(4) saves information
                ! about the expression.
                SELECT CASE (iexptyp)
                CASE (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll,sbdex2, 0, SEC_SBDEXPRESSIONS)
                CASE (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll,sbdex2, 0, SEC_SBDEXPRESSIONS)
                END SELECT
              
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,2,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)

              END IF
              
              ! Now the same thing again, this time separately for primal and dual
              ! variables.
              ! If a primal velocity is specified, Dirichlet-0-boundary conditions are
              ! assumed for the corresponding dual.
              IF (sbdex1 .NE. '') THEN
              
                ! dual X-velocity if primal X-velocity exists
                !
                ! The 2nd element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 4
                
                ! Ditichlet-0 boundary
                iexptyp = BDC_VALDOUBLE
                rcoll%Dquickaccess(4) = 0.0_DP
                rcoll%IquickAccess(3) = iexptyp
                rcoll%SquickAccess(1) = ''
                
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,4,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)
                    
              END IF
              
              IF (sbdex2 .NE. '') THEN
              
                ! dual Y-velocity if primal Y-velocity exists
                !
                ! The 2nd element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 5
                
                ! Ditichlet-0 boundary
                iexptyp = BDC_VALDOUBLE
                rcoll%Dquickaccess(4) = 0.0_DP
                rcoll%IquickAccess(3) = iexptyp
                rcoll%SquickAccess(1) = ''
                
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,5,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)
                    
              END IF
              
            CASE (2)
            
              ! Pressure drop boundary conditions.
              ! Read the line again to get the actual parameters
              READ(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1
              
              ! For any string <> '', create the appropriate pressure drop boundary
              ! condition and add it to the list of boundary conditions.
              !
              ! The IquickAccess array is set up as follows:
              !  IquickAccess(1) = Type of boundary condition
              !  IquickAccess(2) = 0 (undefined)
              !  IquickAccess(3) = expression type
              !  IquickAccess(4) = expression identifier
              !
              ! The SquickAccess array is set up as follows:
              !  SquickAccess(1) = Name of the expression
              !
              IF (sbdex1 .NE. '') THEN
              
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex1)
                rcoll%IquickAccess(3) = iexptyp

                ! The first element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex1
                
                ! Dquickaccess(3) / IquickAccess(2) saves information
                ! about the expression.
                SELECT CASE (iexptyp)
                CASE (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll,sbdex1, 0, SEC_SBDEXPRESSIONS)
                CASE (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll,sbdex1, 0, SEC_SBDEXPRESSIONS)
                END SELECT
              
                IvelEqns = (/1,2/)
                CALL bcasm_newPdropBConRealBd (&
                    rdiscretisation,IvelEqns,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)     
              END IF
              
              
            CASE (3)
              
              ! Nonlinear slip boundary conditions.
              IvelComp = (/1,2/)
              CALL bcasm_newSlipBConRealBd (&
                  rdiscretisation,IvelEqns(1:NDIM2D),rboundaryRegion,rdiscreteBC)     
            
            CASE (4)
              ! Simple Dirichlet boundary, separate definition for all
              ! solution components.
              ! Read the line again, get the expressions for X- and Y-velocity
              READ(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2,sbdex3,sbdex4
              
              ! For any string <> '', create the appropriate Dirichlet boundary
              ! condition and add it to the list of boundary conditions.
              !
              ! The IquickAccess array is set up as follows:
              !  IquickAccess(1) = Type of boundary condition
              !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
              !  IquickAccess(3) = expression type
              !  IquickAccess(4) = expression identifier
              !
              ! The SquickAccess array is set up as follows:
              !  SquickAccess(1) = Name of the expression
              !
              rcoll%IquickAccess(1) = ibctyp
              
              ! If the type is a double precision value, set the DquickAccess(4)
              ! to that value so it can quickly be accessed.
              IF (sbdex1 .NE. '') THEN
                ! X-velocity
                !
                ! The 2nd element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 1
                
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex1)
                rcoll%IquickAccess(3) = iexptyp
                
                ! The 1st element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex1
                
                ! Dquickaccess(3) / IquickAccess(3) saves information
                ! about the expression.
                SELECT CASE (iexptyp)
                CASE (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll, sbdex1, 0, SEC_SBDEXPRESSIONS)
                CASE (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll, sbdex1, 0, SEC_SBDEXPRESSIONS)
                END SELECT
              
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,1,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)
                    
              END IF
              
              IF (sbdex2 .NE. '') THEN
              
                ! Y-velocity
                !
                ! The 1st element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 2
                
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex2)
                rcoll%IquickAccess(3) = iexptyp
                
                ! The 1st element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex2
                
                ! Dquickaccess(4) / IquickAccess(4) saves information
                ! about the expression.
                SELECT CASE (iexptyp)
                CASE (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll,sbdex2, 0, SEC_SBDEXPRESSIONS)
                CASE (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll,sbdex2, 0, SEC_SBDEXPRESSIONS)
                END SELECT
              
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,2,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)

              END IF
              
              ! Now the same thing again, this time separately for primal and dual
              ! variables.
              ! If a velocity is not specified, Dirichlet-0-boundary conditions are
              ! assumed.
              
              IF (sbdex3 .NE. '') THEN
                ! X-velocity
                !
                ! The 2nd element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 4
                
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex3)
                rcoll%IquickAccess(3) = iexptyp
                
                ! The 1st element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex3
                
                ! Dquickaccess(3) / IquickAccess(3) saves information
                ! about the expression.
                SELECT CASE (iexptyp)
                CASE (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll, sbdex3, 0, SEC_SBDEXPRESSIONS)
                CASE (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll, sbdex3, 0, SEC_SBDEXPRESSIONS)
                END SELECT
              
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,4,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)
                    
              END IF
              
              IF (sbdex4 .NE. '') THEN
              
                ! Y-velocity
                !
                ! The 1st element in IquickAccess saves the component number.
                rcoll%IquickAccess(2) = 5
                
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex4)
                rcoll%IquickAccess(3) = iexptyp
                
                ! The 1st element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex4
                
                ! Dquickaccess(4) / IquickAccess(4) saves information
                ! about the expression.
                SELECT CASE (iexptyp)
                CASE (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll,sbdex4, 0, SEC_SBDEXPRESSIONS)
                CASE (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll,sbdex4, 0, SEC_SBDEXPRESSIONS)
                END SELECT
              
                ! Assemble the BC's.
                CALL bcasm_newDirichletBConRealBD (&
                    rdiscretisation,5,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll)
                    
              END IF

            CASE DEFAULT
              CALL output_line ('Unknown boundary condition!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
              CALL sys_halt()
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

    ! Release the parser object with all the expressions to be evaluated
    ! on the boundary.
    CALL fparser_release (rparser)

    ! Remove the boundary value parser from the collection
    CALL collct_deletevalue (rcoll, BDC_BDPARSER)
    
    ! Remove the boundary-expression section we added earlier,
    ! with all their content.
    CALL collct_deletesection(rcoll,SEC_SBDEXPRESSIONS)

    ! Remove the temporary collection from memory.
    CALL collct_done (rcoll)

    ! The setting in the DAT file may overwrite our guess about Neumann boundaries.
    CALL parlst_getvalue_int_indir (p_rbdcond, 'ineumannBoundary', i, -1)
    SELECT CASE (i)
    CASE (0)
      bNeumann = .FALSE.
    CASE (1)
      bNeumann = .TRUE.
    END SELECT
    
    ! Remember whether there is Neumann boundary or not.
    !
    ! Note: actually the bhasNeumannBoundary flag signales whether there
    ! are Neumann boundary components discretely *visible* on that level or not!
    ! We just initialise it here according to whether there are analytically
    ! visible or not. This is still a lack in the design and has to be
    ! somehow fixed later!
    rproblem%RlevelInfo(rproblem%NLMIN:rproblem%NLMAX)%bhasNeumannBoundary = &
        bNeumann
        
  END SUBROUTINE  

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_getBDconditions (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                 cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  USE collection
  USE spatialdiscretisation
  USE discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the solution component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...,
  !   3=3rd solution component, e.g. pressure)
  ! For pressure drop boundary / normal stress:
  !   Velocity components that are affected by the normal stress
  !   (usually "1 2" for x- and y-velocity while returned value must specify
  !   the pressure at the boundary)
  INTEGER, DIMENSION(:), INTENT(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  TYPE(t_boundaryRegion), INTENT(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  INTEGER(I32), INTENT(IN)                                    :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  INTEGER, INTENT(IN)                                         :: cinfoNeeded
  
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDFUNCMID : 
  !   iwhere = number of the edge in which midpoint the value
  !            should be computed
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
  !   iwhere = Number of the edge where the normal stress should be computed.
  INTEGER(I32), INTENT(IN)                                    :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
  !   dwhere = parameter value of the point on edge iwhere where the normal
  !            stress should be computed.
  REAL(DP), INTENT(IN)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(IN), OPTIONAL      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  !
  ! The function may return SYS_INFINITY as a value. This indicates the
  ! framework to ignore the node and treat it as 'natural boundary condition'
  ! node.
  REAL(DP), DIMENSION(:), INTENT(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

    INTEGER :: icomponent,iexprtyp
    
    REAL(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    dtime = rcollection%Dquickaccess(1)

    ! Use boundary conditions from DAT files.
    SELECT CASE (cinfoNeeded)
    
    CASE (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
          DISCBC_NEEDINTMEAN,DISCBC_NEEDNORMALSTRESS)
      
      ! Dirichlet boundary conditions
    
      ! The IquickAccess array is set up as follows:
      !  IquickAccess(1) = Type of boundary condition
      !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
      !  IquickAccess(3) = expression type
      !  IquickAccess(4) = expression identifier
      !
      ! Get from the current component of the PDE we are discretising:
      icomponent = rcollection%IquickAccess(2)
      
      ! -> 1=X-velocity, 2=Y-velocity.
      
      ! Return zero Dirichlet boundary values for all situations by default.
      Dvalues(1) = 0.0_DP

      ! Now, depending on the problem, calculate the return value.
      
      ! The IquickAccess array is set up as follows:
      !  IquickAccess(1) = Type of boundary condition
      !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
      !  IquickAccess(3) = expression type
      !  IquickAccess(4) = expression identifier
      !
      ! Get the type of the expression to evaluate from the 
      ! integer tag of the BC-region - if there is an expression to evaluate
      ! at all.
      iexprtyp = rcollection%IquickAccess(3)
            
      ! Now, which boundary condition do we have here?
      !
      ! Get the information from evalBoundary.
      ! Note: The information about the BC's can be retrieved from the
      ! quick-access arrays in the collection as initialised above.
      SELECT CASE (rcollection%IquickAccess(1))
      CASE (1)
        ! Simple Dirichlet BC's. Evaluate the expression iexprtyp.
        Dvalues(1) = evalBoundary (icomponent,rdiscretisation, rboundaryRegion, &
            iexprtyp, rcollection%IquickAccess(4), rcollection%DquickAccess(4), &
            dwhere, rcollection%SquickAccess(1),dtime,&
            rcollection)
    
      CASE (2)
        ! Normal stress / pressure drop. Evaluate the  expression iexprtyp.
        Dvalues(1) = evalBoundary (icomponent,rdiscretisation, rboundaryRegion, &
            iexprtyp, rcollection%IquickAccess(4), rcollection%DquickAccess(4), &
            dwhere, rcollection%SquickAccess(1),dtime,&
            rcollection)
            
      END SELECT
      
    END SELECT
  
  CONTAINS
  
    ! Auxiliary function: Evaluate a scalar expression on the boundary.
    
    REAL(DP) FUNCTION evalBoundary (icomponent,rdiscretisation, rboundaryRegion, &
                                    ityp, ivalue, dvalue, dpar, stag, dtime, rcollection)
    
    ! Solution component for which the expression is evaluated.
    ! 1 = X-velocity, 2 = y-velocity,...
    INTEGER, INTENT(IN) :: icomponent
    
    ! Discretisation structure of the underlying discretisation
    TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation
    
    ! Current boundary region
    TYPE(t_boundaryRegion), INTENT(IN) :: rboundaryRegion
    
    ! Type of expression to evaluate.
    ! One of the BDC_xxxx constants from ccboundaryconditionparser.f90.
    INTEGER, INTENT(IN) :: ityp
    
    ! Integer tag. If ityp=BDC_EXPRESSION, this must specify the number of
    ! the expression in the expression object to evaluate.
    ! Otherwise unused.
    INTEGER, INTENT(IN) :: ivalue
    
    ! Double precision parameter for simple expressions
    REAL(DP), INTENT(IN) :: dvalue
    
    ! Current parameter value of the point on the boundary.
    ! 0-1-parametrisation.
    REAL(DP), INTENT(IN) :: dpar

    ! String tag that defines more complicated BC's.
    CHARACTER(LEN=*), INTENT(IN) :: stag
    
    ! For nonstationary simulation: Simulation time.
    ! =0 for stationary simulations.
    REAL(DP), INTENT(IN) :: dtime
    
    ! A compiled expression for evaluation at runtime
    TYPE(t_fparser), POINTER :: p_rparser
    
    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), OPTIONAL                  :: rcollection

      ! local variables
      REAL(DP) :: d,dx,dy
      CHARACTER(LEN=PARLST_MLDATA) :: sexpr
      REAL(DP), DIMENSION(SIZE(SEC_EXPRVARIABLES)) :: Rval
      
      SELECT CASE (ityp)
      CASE (BDC_USERDEF)
        ! This is a hardcoded, user-defined identifier.
        ! In stag, the name of the identifier is noted.
        ! Get the identifier itself from the collection.
        CALL collct_getvalue_string (rcollection, stag, sexpr, &
                                     0, SEC_SBDEXPRESSIONS)
                                
        ! Call the user defined callback routine to evaluate the expression.
        !
        ! As collection, we pass rcollection%p_rcollection here; this is a pointer
        ! to the application specific, global collection that may be of interest for
        ! callback routines. rcollection itself is actually a 'local' collection,
        ! a 'wrapper' for the rcollection of the application!
        CALL getBoundaryValues (&
            stag,icomponent,rdiscretisation,rboundaryRegion,&
            dpar, d, rcollection%p_rnextCollection)     
        
        evalBoundary = d
      
      CASE (BDC_VALDOUBLE)
        ! A simple constant, given by dvalue
        evalBoundary = dvalue

      CASE (BDC_EXPRESSION)
        ! A complex expression.
        ! Get the expression object from the collection.
        
        p_rparser => collct_getvalue_pars (rcollection, BDC_BDPARSER, &
                                   0, SEC_SBDEXPRESSIONS)
                                   
        ! Set up an array with variables for evaluating the expression.
        ! Give the values in exactly the same order as specified
        ! by SEC_EXPRVARIABLES!
        Rval = 0.0_DP
        
        CALL boundary_getCoords(rdiscretisation%p_rboundary, &
                                rboundaryRegion%iboundCompIdx, &
                                dpar, dx, dy)
        
        ! Get the local parameter value 0 <= d <= 1.
        ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
        ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
        ! although 0 <= dminpar <= max.par 
        !      and 0 <= dmaxpar <= max.par!
        d = dpar 
        IF (d .LT. rboundaryRegion%dminParam) &
          d = d + boundary_dgetMaxParVal(rdiscretisation%p_rboundary,&
                                         rboundaryRegion%iboundCompIdx)
        d = d - rboundaryRegion%dminParam
        
        ! Normalise to 0..1 using the length of the parameter region.
        ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
        d = d / (rboundaryRegion%dmaxParam - rboundaryRegion%dminParam)
        
        Rval(1) = dx
        Rval(2) = dy
        ! Rval(3) = .
        Rval(4) = d
        Rval(5) = dpar
        Rval(6) = boundary_convertParameter(rdiscretisation%p_rboundary, &
                                            rboundaryRegion%iboundCompIdx, dpar, &
                                            BDR_PAR_01, BDR_PAR_LENGTH) 
        Rval(7) = dtime
        
        ! Evaluate the expression. ivalue is the number of
        ! the expression to evaluate.
        CALL fparser_evalFunction (p_rparser, ivalue, Rval, evalBoundary)
        
      CASE (BDC_VALPARPROFILE)
        ! A parabolic profile. dvalue expresses the
        ! maximum value of the profile. 
        !
        ! Get the local parameter value 0 <= d <= 1.
        ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
        ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
        ! although 0 <= dminpar <= max.par 
        !      and 0 <= dmaxpar <= max.par!
        d = dpar 
        IF (d .LT. rboundaryRegion%dminParam) &
          d = d + boundary_dgetMaxParVal(rdiscretisation%p_rboundary,&
                                         rboundaryRegion%iboundCompIdx)
        d = d - rboundaryRegion%dminParam
    
        ! Normalise to 0..1 using the length of the parameter region.
        ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
        d = d / (rboundaryRegion%dmaxParam - rboundaryRegion%dminParam)
        
        evalBoundary = mprim_getParabolicProfile (d,1.0_DP,dvalue) 
      END SELECT
    
    END FUNCTION

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_assembleFBDconditions (rproblem,rdiscretisation,rdiscreteFBC,rcollection)

!<description>
  ! This parses the boundary conditions for fictitious boundary
  ! components in the problem and assembles them into rdiscreteFBC.
!</description>
  
!<input>
  ! A discretisation structure defining the discretisation of the current
  ! level.
  TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscretisation
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! Collection structure to be passed to callback routines
  TYPE(t_collection), INTENT(INOUT), TARGET :: rcollection
  
  ! A t_discreteFBC structure that receives a discretised version
  ! of the fictitious boundary boundary conditions. The structure should
  ! be empty; new BC's are simply added to the structure.
  TYPE(t_discreteFBC), INTENT(INOUT) :: rdiscreteFBC
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
    ! CALL bcasm_newDirichletBConFBD (rdiscretisation,Iequations,rdiscreteFBC,&
    !     getBoundaryValuesFBC,rcollection)

  END SUBROUTINE  
  
END MODULE
