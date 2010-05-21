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

module cc2dmediumm2boundarydef

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
  use boundarycondition

  use collection
  use convection
  use basicgeometry
    
  use cc2dmediumm2basic
  
  implicit none

!<constants>

!<constantblock description="Names of the sections used in the main collection">

  ! Section name of the section saving the boundary expressions, i.e. the
  ! expressions that are to be evaluated in each point on the boundary.
  character(LEN=COLLCT_MLSECTION), parameter :: SEC_SBDEXPRESSIONS = 'BDEXPRESSIONS'

  ! Name of the parser object for boundary value expressions
  character(LEN=COLLCT_MLSECTION), parameter :: BDC_BDPARSER = 'BDEXPRPARSER'
!</constantblock>

!<constantblock description="The different types of boundary conditions supported by this module">

  ! User defined expression, evaluated by calling the callback routine
  integer, parameter :: BDC_USERDEF = -2
  
  ! Text expression, evaluated by parsing
  integer, parameter :: BDC_EXPRESSION = -1
  
  ! Fixed double precision value 
  integer, parameter :: BDC_VALDOUBLE = 0

  ! Fixed integer value 
  integer, parameter :: BDC_VALINT    = 1

  ! Parabolic profile with prescribed maximum value
  integer, parameter :: BDC_VALPARPROFILE = 2
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
  character(LEN=10), dimension(7), parameter :: SEC_EXPRVARIABLES = &
    (/'X    ','Y    ','Z    ','L    ','R    ','S    ','TIME '/)

!</constantblock>

!</constants>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_parseBDconditions (rproblem)

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
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    logical :: bNeumann
    integer :: i,ityp,ivalue,ibdComponent,isegment,iintervalEnds
    integer :: ibctyp,icount,iexptyp
    integer, dimension(2) :: IminIndex,imaxIndex
    integer, dimension(NDIM2D) :: IvelComp
    real(DP) :: dvalue,dpar1,dpar2
    character(LEN=PARLST_MLDATA) :: cstr,cexpr,sbdex1,sbdex2
    character(LEN=PARLST_MLNAME) :: cname
    integer, dimension(NDIM2D) :: IvelEqns
    
    ! A local collection we use for storing named parameters during the
    ! parsing process.
    type(t_collection) :: rcoll
    
    ! Triangulation on currently highest level.
    type(t_triangulation), pointer :: p_rtriangulation

    ! A set of variables describing the analytic boundary conditions.    
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_bcRegion), pointer :: p_rbcRegion
    
    ! A pointer to the domain
    type(t_boundary), pointer :: p_rboundary
    
    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection,p_rbdcond
    
    ! A compiled expression for evaluation at runtime
    type(t_fparser), pointer :: p_rparser
    
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
    nullify (rproblem%p_rboundaryConditions)
    call bcond_initBC (rproblem%p_rboundaryConditions,p_rboundary)
    
    ! Get the expression/bc sections from the bondary condition block
    call parlst_querysection(rproblem%rparamList, 'BDEXPRESSIONS', p_rsection) 
    call parlst_querysection(rproblem%rparamList, 'BDCONDITIONS', p_rbdcond) 
    
    ! Add a section to the collection that accepts the boundary expressions
    call collct_addsection (rproblem%rcollection, SEC_SBDEXPRESSIONS)
    
    ! For intermediate storing of expression types, we use a local collection
    call collct_init (rcoll)
    
    ! Create a parser structure for as many expressions as configured
    allocate (p_rparser)
    call fparser_create (p_rparser,&
         parlst_querysubstrings (p_rsection, 'bdExpressions'))
    
    ! Add the parser to the collection
    call collct_setvalue_pars (rproblem%rcollection, BDC_BDPARSER, p_rparser, &
                                .true., 0, SEC_SBDEXPRESSIONS)
    
    ! Add the boundary expressions to the collection into the
    ! specified section.
    do i=1,parlst_querysubstrings (p_rsection, 'bdExpressions')
    
      call parlst_getvalue_string (p_rsection, 'bdExpressions', cstr, '', i)
      
      ! Get the type and decide on the identifier how to save the expression.
      read(cstr,*) cname,ityp
      
      select case (ityp)
      case (BDC_USERDEF)
        ! Name of a hardcoded expression realised in the callback routines
        read(cstr,*) cname,ityp,cexpr
        call collct_setvalue_string (rproblem%rcollection, cname, cexpr, .true., &
                                     0, SEC_SBDEXPRESSIONS)

      case (BDC_EXPRESSION)
        ! General expression; not implemented yet
        read(cstr,*) cname,ityp,cexpr
        
        ! Compile the expression; the expression gets number i
        call fparser_parseFunction (p_rparser, i, cexpr, SEC_EXPRVARIABLES)
        
        ! Add the number of the function in the parser object to
        ! the collection with the name of the expression
        call collct_setvalue_int (rproblem%rcollection, cname, i, .true., &
                                   0, SEC_SBDEXPRESSIONS)
        
      case (BDC_VALDOUBLE)
        ! Real-value
        read(cstr,*) cname,ityp,dvalue
        call collct_setvalue_real (rproblem%rcollection, cname, dvalue, .true., &
                                   0, SEC_SBDEXPRESSIONS) 
      case (BDC_VALINT)
        ! Integer-value
        read(cstr,*) cname,ityp,ivalue
        call collct_setvalue_int (rproblem%rcollection, cname, ivalue, .true., &
                                   0, SEC_SBDEXPRESSIONS)
                
      case (BDC_VALPARPROFILE)
        ! Parabolic profile with specified maximum velocity
        read(cstr,*) cname,ityp,dvalue
        call collct_setvalue_real (rproblem%rcollection, cname, dvalue, .true., &
                                   0, SEC_SBDEXPRESSIONS) 
                                   
      case DEFAULT
        print *,'Expressions not implemented!'
        stop

      end select
      
      ! Put the type of the expression to the temporary collection section
      call collct_setvalue_int (rcoll, cname, ityp, .true.) 
      
    end do
    
    bNeumann = .false.
    
    ! Now to the actual boundary conditions.
    ! Loop through all boundary components we have.
    do ibdComponent = 1,boundary_igetNBoundComp(p_rboundary)
      
      ! Parse the parameter 'bdComponentX'
      write (cexpr,'(I10)') ibdComponent
      cstr = 'bdComponent' // adjustl(cexpr)
      
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      i = parlst_queryvalue (p_rbdcond, cstr)
      if (i .ne. 0) then
        ! Parameter exists. Get the values in there.
        do isegment = 1,parlst_querysubstrings (p_rbdcond, cstr)
          
          call parlst_getvalue_string_fetch (p_rbdcond, &
                                            i, cstr, isubstring=isegment)
          ! Which type is that?
          read(cstr,*) ityp
          
          select case (ityp)
          case (0)
            read(cstr,*) ityp,dpar2,iintervalEnds,ibctyp
          case (1)
            read(cstr,*) ityp,ivalue,iintervalEnds,ibctyp
            dpar2 = real(ivalue,DP)
          end select                                   
                        
          ! Form a boundary condition segment that covers that boundary part
          if (dpar2 .ge. dpar1) then
            
            rboundaryRegion%dminParam = dpar1
            rboundaryRegion%dmaxParam = dpar2
            rboundaryRegion%iboundCompIdx = ibdComponent
            rboundaryRegion%dmaxParamBC = &
              boundary_dgetMaxParVal(p_rboundary, ibdComponent)
            rboundaryRegion%iproperties = iintervalEnds
            
            ! Now, which type of BC is to be created?
            select case (ibctyp)
            
            case (0)
              ! Usually there's Neumann boundary in this region, but we can't be 
              ! sure. Check if, on the highest level, there's at least one edge
              ! of the triangulation belonging to the boundary. If yes, we
              ! have found Neumann boundary. If no, the segment is just too
              ! small to be considered as Neumann boundary.
              
              call bcasm_getEdgesInBCregion (p_rtriangulation,p_rboundary,&
                                            rboundaryRegion, &
                                            IminIndex,ImaxIndex,icount)
              if (icount .gt. 0) bNeumann = .true.
            
            case (1)
              ! Simple Dirichlet boundary
              ! Read the line again, get the expressions for X- and Y-velocity
              read(cstr,*) ityp,dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
              
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
              if (sbdex1 .ne. '') then
                call bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,&
                                                  1,rboundaryRegion,p_rbcRegion)
                p_rbcRegion%stag = sbdex1

                iexptyp = collct_getvalue_int (rcoll, sbdex1)
                p_rbcRegion%ibdrexprtype = iexptyp

                select case (iexptyp)
                case (0,2)
                  ! Constant or parabolic profile
                  p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                    sbdex1, 0, SEC_SBDEXPRESSIONS)
                case (-1)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  p_rbcRegion%itag = collct_getvalue_int (rproblem%rcollection, &
                                     sbdex1, 0, SEC_SBDEXPRESSIONS)
                end select

              end if
              
              if (sbdex2 .ne. '') then
                call bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,&
                                                  2,rboundaryRegion,p_rbcRegion)
                p_rbcRegion%stag = sbdex2

                iexptyp = collct_getvalue_int (rcoll, sbdex2)
                p_rbcRegion%ibdrexprtype = iexptyp

                select case (iexptyp)
                case (0,2)
                  p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                     sbdex2, 0, SEC_SBDEXPRESSIONS)
                case (-1)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  p_rbcRegion%itag = collct_getvalue_int (rproblem%rcollection, &
                                     sbdex2, 0, SEC_SBDEXPRESSIONS)
                end select

              end if
              
              ! If we have no-slip boundary conditions, the X- and Y-velocity
              ! matrices are 'decoupled' as they are modified differently
              ! by the boundary-conditions implementation filter!
              ! In this case, change rproblem%bdecoupledXY from FALSE to TRUE
              ! to indicate that.
              if ( ((sbdex1 .eq. '') .and. (sbdex2 .ne. '')) .or.&
                   ((sbdex1 .ne. '') .and. (sbdex2 .eq. '')) ) then
                rproblem%bdecoupledXY = .true.
              end if
              
            case (2)
            
              ! Pressure drop boundary conditions.
              ! Read the line again to get the actual parameters
              read(cstr,*) ityp,dvalue,iintervalEnds,ibctyp,sbdex1
              
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
              if (sbdex1 .ne. '') then
                IvelEqns = (/1,2/)
                call bcond_newPressureDropBConRealBD (rproblem%p_rboundaryConditions,&
                                              IvelEqns,&
                                              rboundaryRegion,p_rbcRegion)          
                p_rbcRegion%stag = sbdex1

                iexptyp = collct_getvalue_int (rcoll, sbdex1)
                p_rbcRegion%ibdrexprtype = iexptyp

                select case (iexptyp)
                case (0,2)
                  ! Constant or parabolic profile
                  p_rbcRegion%dtag = collct_getvalue_real (rproblem%rcollection, &
                                    sbdex1, 0, SEC_SBDEXPRESSIONS)
                case (-1)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  p_rbcRegion%itag = collct_getvalue_int (rproblem%rcollection, &
                                     sbdex1, 0, SEC_SBDEXPRESSIONS)
                end select

              end if
              
              
            case (3)
              
              ! Nonlinear slip boundary conditions.
              IvelComp = (/1,2/)
              call bcond_newSlipBConRealBD (rproblem%p_rboundaryConditions,&
                                            IvelComp(1:NDIM2D),&
                                            rboundaryRegion,p_rbcRegion)
            
            case DEFAULT
              print *,'Unknown boundary condition'
              stop
            end select
            
            ! Move on to the next parameter value
            dpar1 = dpar2
            
          end if
                                            
        end do
      
      else
        ! There is no parameter configuring the boundary condition on that 
        ! component - so we have Neumann boundary there.
        bNeumann = .true.
      end if
      
    end do

    ! Remember whether there is Neumann boundary or not.
    !
    ! Note: actually the bhasNeumannBoundary flag signales whether there
    ! are Neumann boundary components discretely *visible* on that level or not!
    ! We just initialise it here according to whether there are analytically
    ! visible or not. This is still a lack in the design and has to be
    ! somehow fixed later!
    rproblem%RlevelInfo(rproblem%NLMIN:rproblem%NLMAX)%bhasNeumannBoundary = &
        bNeumann
    
    ! Remove the temporary collection from memory.
    call collct_done (rcoll)

  end subroutine  
  
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_parseFBDconditions (rproblem)

!<description>
  ! This initialises the analytic boundary conditions for fictitious boundary
  ! components in the problem and saves them to the problem structure.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    ! An identifier array for the equations to be tackled when discretising
    ! boundary conditions on fictitious boundary components
    integer, dimension(2) :: Iequations
    
    ! A structure identifying the fictitious boundary component
    type(t_fictBoundaryRegion) :: rfictBoundaryRegion

    ! A set of variables describing the analytic boundary conditions.    
    type(t_bcRegion), pointer :: p_rbcRegion
    
    ! Add a new fictitious boundary object It should impose Dirichlet
    ! boundary conditions in the domain in the one and only solution component.
    ! We use the default initialisation of rfictBoundaryRegion and only
    ! change the name of the component.
    rfictBoundaryRegion%sname = 'CIRCLE'
    Iequations = (/1,2/)    ! 1=x, 2=y-velocity
    !CALL bcond_newDirichletBConFictBD (rproblem%p_rboundaryConditions,Iequations,&
    !                                   rfictBoundaryRegion,p_rbcRegion)    

  end subroutine  
  
end module
