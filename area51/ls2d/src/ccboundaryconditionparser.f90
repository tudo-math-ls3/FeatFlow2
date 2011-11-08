!##############################################################################
!# ****************************************************************************
!# <name> ccboundaryconditionparser </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for parsing boundary conditions from a DAT
!# file.
!#
!# The following files can be found here:
!#
!# 1.) cc_assembleBDconditions
!#     -> Parses the definition of the BC's from sections given by DAT files
!#        and sets up an discrete boundary condition description
!#
!# 2.) cc_assembleFBDconditions
!#     -> Assembles fictitious boundary boundary conditions.
!#
!# </purpose>
!##############################################################################

module ccboundaryconditionparser

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use basicgeometry
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use bcassemblybase
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use bcassembly
  use fparser
  use paramlist
  use bcassembly
  
  use collection
  use convection
    
  use ccbasic
  use cccallback
  
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

  subroutine cc_assembleBDconditions (rproblem,rdiscretisation,rdiscreteBC,&
    rcollection,bforPostprocessing)

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
  type(t_blockDiscretisation), intent(IN) :: rdiscretisation
  
  ! OPTIONAL: If this flag ist set to TRUE, the boundary conditions are
  ! assembled for postprocessing of a solution vector. When being set to FALSE
  ! or not present, the boundary condition are assembled for computation.
  logical, intent(IN), optional :: bforPostprocessing
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! Collection structure to be passed to callback routines
  type(t_collection), intent(INOUT), target :: rcollection
  
  ! A t_discreteBC structure that receives a discretised version
  ! of the boundary boundary conditions. The structure should
  ! be empty; new BC's are simply added to the structure.
  type(t_discreteBC), intent(INOUT) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    logical :: bNeumann
    integer :: i,ityp,ivalue,ibdComponent,isegment,iintervalEnds
    integer :: ibctyp,icount,iexptyp
    integer, dimension(2) :: IminIndex,imaxIndex
    real(DP) :: dvalue,dpar1,dpar2
    character(LEN=PARLST_MLDATA) :: cstr,cexpr,sbdex1,sbdex2
    character(LEN=PARLST_MLNAME) :: cname
    integer, dimension(NDIM2D) :: IvelEqns
    type(t_collection) :: rlocalCollection
    integer(I32) :: casmComplexity
    
    ! A local collection we use for storing named parameters during the
    ! parsing process.
    type(t_collection) :: rcoll
    
    ! Triangulation on currently highest level.
    type(t_triangulation), pointer :: p_rtriangulation

    ! A set of variables describing the analytic boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    
    ! A pointer to the domain
    type(t_boundary), pointer :: p_rboundary
    
    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection,p_rbdcond
    
    ! A compiled expression for evaluation at runtime
    type(t_fparser), target :: rparser
    
    ! Determine what to assemble
    casmComplexity = BCASM_DISCFORALL
    if (present(bforPostprocessing)) then
      ! Assemble only for the solution vector
      if (bforPostprocessing) casmComplexity = BCASM_DISCFORSOL
    endif
    
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
    call parlst_querysection(rproblem%rparamList, 'BDEXPRESSIONS', p_rsection)
    call parlst_querysection(rproblem%rparamList, 'BDCONDITIONS', p_rbdcond)
    
    ! For intermediate storing of expression types, we use a local collection
    call collct_init (rcoll)
    
    ! Add a section to the collection that accepts the boundary expressions
    call collct_addsection (rcoll, SEC_SBDEXPRESSIONS)
    
    ! Create a parser structure for as many expressions as configured
    call fparser_create (rparser,&
         parlst_querysubstrings (p_rsection, 'bdExpressions'))
    
    ! Add the parser to the collection
    call collct_setvalue_pars (rcoll, BDC_BDPARSER, rparser, &
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
        call collct_setvalue_string (rcoll, cname, cexpr, .true., &
                                     0, SEC_SBDEXPRESSIONS)

      case (BDC_EXPRESSION)
        ! General expression; not implemented yet
        read(cstr,*) cname,ityp,cexpr
        
        ! Compile the expression; the expression gets number i
        call fparser_parseFunction (rparser, i, cexpr, SEC_EXPRVARIABLES)
        
        ! Add the number of the function in the parser object to
        ! the collection with the name of the expression
        call collct_setvalue_int (rcoll, cname, i, .true., &
                                   0, SEC_SBDEXPRESSIONS)
        
      case (BDC_VALDOUBLE)
        ! Real-value
        read(cstr,*) cname,ityp,dvalue
        call collct_setvalue_real (rcoll, cname, dvalue, .true., &
                                   0, SEC_SBDEXPRESSIONS)
      case (BDC_VALINT)
        ! Integer-value
        read(cstr,*) cname,ityp,ivalue
        call collct_setvalue_int (rcoll, cname, ivalue, .true., &
                                   0, SEC_SBDEXPRESSIONS)
                
      case (BDC_VALPARPROFILE)
        ! Parabolic profile with specified maximum velocity
        read(cstr,*) cname,ityp,dvalue
        call collct_setvalue_real (rcoll, cname, dvalue, .true., &
                                   0, SEC_SBDEXPRESSIONS)
                                   
      case DEFAULT
        call output_line ('Expressions not implemented!', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
        call sys_halt()

      end select
      
      ! Put the type of the expression to the temporary collection section
      call collct_setvalue_int (rcoll, cname, ityp, .true.)
      
    end do
    
    bNeumann = .false.
    
    ! Now to the actual boundary conditions.
    !
    ! Put some information to the quick access arrays for access
    ! in the callback routine.
    rcoll%Iquickaccess(1) = rproblem%itimedependence
    select case (rproblem%itimedependence)
    case (0)
      ! Stationary simulation
      rcoll%Dquickaccess(1) = 0.0_DP
      rcoll%Dquickaccess(2) = 0.0_DP
      rcoll%Dquickaccess(3) = 0.0_DP
    case (1)
      rcoll%Dquickaccess(1) = rproblem%rtimedependence%dtime
      rcoll%Dquickaccess(2) = rproblem%rtimedependence%dtimeInit
      rcoll%Dquickaccess(3) = rproblem%rtimedependence%dtimeMax
    end select
    
    ! DquickAccess(4:) is reserved for BC specific information.
    !
    ! Put a link to the previous collection into that local collection.
    ! That allows us to access it or to pass it to user defined callback
    ! functions.
    rcoll%p_rnextCollection => rcollection

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
          
          call parlst_getvalue_string (p_rbdcond, i, cstr, isubstring=isegment)
          ! Read the segment parameters
          read(cstr,*) dpar2,iintervalEnds,ibctyp
          
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
              read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
              
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
              if (sbdex1 .ne. '') then
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
                select case (iexptyp)
                case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll, sbdex1, 0, SEC_SBDEXPRESSIONS)
                case (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll, sbdex1, 0, SEC_SBDEXPRESSIONS)
                end select
              
                ! Assemble the BC's.
                call bcasm_newDirichletBConRealBD (&
                    rdiscretisation,1,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll,casmComplexity)
                    
              end if
              
              if (sbdex2 .ne. '') then
              
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
                select case (iexptyp)
                case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll,sbdex2, 0, SEC_SBDEXPRESSIONS)
                case (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll,sbdex2, 0, SEC_SBDEXPRESSIONS)
                end select
              
                ! Assemble the BC's.
                call bcasm_newDirichletBConRealBD (&
                    rdiscretisation,2,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll,casmComplexity)

              end if
              
            case (2)
            
              ! Pressure drop boundary conditions.
              ! Read the line again to get the actual parameters
              read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1
              
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
              if (sbdex1 .ne. '') then
              
                ! IquickAccess(3) saves the type of the expression
                iexptyp = collct_getvalue_int (rcoll, sbdex1)
                rcoll%IquickAccess(3) = iexptyp

                ! The first element in the sting quick access array is
                ! the name of the expression to evaluate.
                rcoll%SquickAccess(1) = sbdex1
                
                ! Dquickaccess(3) / IquickAccess(2) saves information
                ! about the expression.
                select case (iexptyp)
                case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                  ! Constant or parabolic profile
                  rcoll%Dquickaccess(4) = &
                      collct_getvalue_real (rcoll,sbdex1, 0, SEC_SBDEXPRESSIONS)
                case (BDC_EXPRESSION)
                  ! Expression. Write the identifier for the expression
                  ! as itag into the boundary condition structure.
                  rcoll%IquickAccess(4) = &
                      collct_getvalue_int (rcoll,sbdex1, 0, SEC_SBDEXPRESSIONS)
                end select
              
                IvelEqns = (/1,2/)
                call bcasm_newPdropBConRealBd (&
                    rdiscretisation,IvelEqns,rboundaryRegion,rdiscreteBC,&
                    cc_getBDconditions,rcoll,casmComplexity)
              end if
              
              
            case (3)
              
              ! Nonlinear slip boundary conditions.
              IvelEqns = (/1,2/)
              call bcasm_newSlipBConRealBd (&
                  rdiscretisation,IvelEqns(1:NDIM2D),rboundaryRegion,&
                  rdiscreteBC,casmComplexity)
            
            case DEFAULT
              call output_line ('Unknown boundary condition!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
              call sys_halt()
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

    ! Release the parser object with all the expressions to be evaluated
    ! on the boundary.
    call fparser_release (rparser)

    ! Remove the boundary value parser from the collection
    call collct_deletevalue (rcoll, BDC_BDPARSER)
    
    ! Remove the boundary-expression section we added earlier,
    ! with all their content.
    call collct_deletesection(rcoll,SEC_SBDEXPRESSIONS)

    ! Remove the temporary collection from memory.
    call collct_done (rcoll)

    ! The setting in the DAT file may overwrite our guess about Neumann boundaries.
    call parlst_getvalue_int (p_rbdcond, 'ineumannBoundary', i, -1)
    select case (i)
    case (0)
      bNeumann = .false.
    case (1)
      bNeumann = .true.
    end select
    
    ! Remember whether there is Neumann boundary or not.
    !
    ! Note: actually the bhasNeumannBoundary flag signales whether there
    ! are Neumann boundary components discretely *visible* on that level or not!
    ! We just initialise it here according to whether there are analytically
    ! visible or not. This is still a lack in the design and has to be
    ! somehow fixed later!
    rproblem%RlevelInfo(rproblem%NLMIN:rproblem%NLMAX)%bhasNeumannBoundary = &
        bNeumann
        
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_getBDconditions (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                 cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(IN)                                         :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(IN)                                         :: cinfoNeeded
  
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
  integer, intent(IN)                                         :: iwhere

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
  real(DP), intent(IN)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  !
  ! The function may return SYS_INFINITY_DP as a value. This indicates the
  ! framework to ignore the node and treat it as 'natural boundary condition'
  ! node.
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

    integer :: icomponent,iexprtyp
    
    real(DP) :: dtime
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    dtime = rcollection%Dquickaccess(1)

    ! Use boundary conditions from DAT files.
    select case (cinfoNeeded)
    
    case (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
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
      select case (rcollection%IquickAccess(1))
      case (1)
        ! Simple Dirichlet BC's. Evaluate the expression iexprtyp.
        Dvalues(1) = evalBoundary (icomponent,rdiscretisation, rboundaryRegion, &
            iexprtyp, rcollection%IquickAccess(4), rcollection%DquickAccess(4), &
            dwhere, rcollection%SquickAccess(1),dtime,&
            rcollection)
    
      case (2)
        ! Normal stress / pressure drop. Evaluate the  expression iexprtyp.
        Dvalues(1) = evalBoundary (icomponent,rdiscretisation, rboundaryRegion, &
            iexprtyp, rcollection%IquickAccess(4), rcollection%DquickAccess(4), &
            dwhere, rcollection%SquickAccess(1),dtime,&
            rcollection)
            
      end select
      
    end select
  
  contains
  
    ! Auxiliary function: Evaluate a scalar expression on the boundary.
    
    real(DP) function evalBoundary (icomponent,rdiscretisation, rboundaryRegion, &
                                    ityp, ivalue, dvalue, dpar, stag, dtime, rcollection)
    
    ! Solution component for which the expression is evaluated.
    ! 1 = X-velocity, 2 = y-velocity,...
    integer, intent(IN) :: icomponent
    
    ! Discretisation structure of the underlying discretisation
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! Current boundary region
    type(t_boundaryRegion), intent(IN) :: rboundaryRegion
    
    ! Type of expression to evaluate.
    ! One of the BDC_xxxx constants from ccboundaryconditionparser.f90.
    integer, intent(IN) :: ityp
    
    ! Integer tag. If ityp=BDC_EXPRESSION, this must specify the number of
    ! the expression in the expression object to evaluate.
    ! Otherwise unused.
    integer, intent(IN) :: ivalue
    
    ! Double precision parameter for simple expressions
    real(DP), intent(IN) :: dvalue
    
    ! Current parameter value of the point on the boundary.
    ! 0-1-parametrisation.
    real(DP), intent(IN) :: dpar

    ! String tag that defines more complicated BC's.
    character(LEN=*), intent(IN) :: stag
    
    ! For nonstationary simulation: Simulation time.
    ! =0 for stationary simulations.
    real(DP), intent(IN) :: dtime
    
    ! A compiled expression for evaluation at runtime
    type(t_fparser), pointer :: p_rparser
    
    ! A pointer to a collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), optional                  :: rcollection

      ! local variables
      real(DP) :: d,dx,dy
      character(LEN=PARLST_MLDATA) :: sexpr
      real(DP), dimension(size(SEC_EXPRVARIABLES)) :: Rval
      
      select case (ityp)
      case (BDC_USERDEF)
        ! This is a hardcoded, user-defined identifier.
        ! In stag, the name of the identifier is noted.
        ! Get the identifier itself from the collection.
        call collct_getvalue_string (rcollection, stag, sexpr, &
                                     0, SEC_SBDEXPRESSIONS)
                                
        ! Call the user defined callback routine to evaluate the expression.
        !
        ! As collection, we pass rcollection%p_rcollection here; this is a pointer
        ! to the application specific, global collection that may be of interest for
        ! callback routines. rcollection itself is actually a 'local' collection,
        ! a 'wrapper' for the rcollection of the application!
        call getBoundaryValues (&
            stag,icomponent,rdiscretisation,rboundaryRegion,&
            dpar, d, rcollection%p_rnextCollection)
        
        evalBoundary = d
      
      case (BDC_VALDOUBLE)
        ! A simple constant, given by dvalue
        evalBoundary = dvalue

      case (BDC_EXPRESSION)
        ! A complex expression.
        ! Get the expression object from the collection.
        
        p_rparser => collct_getvalue_pars (rcollection, BDC_BDPARSER, &
                                   0, SEC_SBDEXPRESSIONS)
                                   
        ! Set up an array with variables for evaluating the expression.
        ! Give the values in exactly the same order as specified
        ! by SEC_EXPRVARIABLES!
        Rval = 0.0_DP
        
        call boundary_getCoords(rdiscretisation%p_rboundary, &
                                rboundaryRegion%iboundCompIdx, &
                                dpar, dx, dy)
        
        ! Get the local parameter value 0 <= d <= 1 in the boundary region.
        ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
        ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
        ! although 0 <= dminpar <= max.par
        !      and 0 <= dmaxpar <= max.par!
        d = dpar
        if (d .lt. rboundaryRegion%dminParam) &
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
        call fparser_evalFunction (p_rparser, ivalue, Rval, evalBoundary)
        
      case (BDC_VALPARPROFILE)
        ! A parabolic profile. dvalue expresses the
        ! maximum value of the profile.
        !
        ! Get the local parameter value 0 <= d <= 1.
        ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
        ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
        ! although 0 <= dminpar <= max.par
        !      and 0 <= dmaxpar <= max.par!
        d = dpar
        if (d .lt. rboundaryRegion%dminParam) &
          d = d + boundary_dgetMaxParVal(rdiscretisation%p_rboundary,&
                                         rboundaryRegion%iboundCompIdx)
        d = d - rboundaryRegion%dminParam
        
        ! Normalise to 0..1 using the length of the parameter region.
        ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
        d = d / (rboundaryRegion%dmaxParam - rboundaryRegion%dminParam)
    
        evalBoundary = mprim_getParabolicProfile (d,1.0_DP,dvalue)
      end select
    
    end function

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleFBDconditions (rproblem,rdiscretisation,rdiscreteFBC,rcollection)

!<description>
  ! This parses the boundary conditions for fictitious boundary
  ! components in the problem and assembles them into rdiscreteFBC.
!</description>
  
!<input>
  ! A discretisation structure defining the discretisation of the current
  ! level.
  type(t_blockDiscretisation), intent(IN) :: rdiscretisation
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! Collection structure to be passed to callback routines
  type(t_collection), intent(INOUT), target :: rcollection
  
  ! A t_discreteFBC structure that receives a discretised version
  ! of the fictitious boundary boundary conditions. The structure should
  ! be empty; new BC's are simply added to the structure.
  type(t_discreteFBC), intent(INOUT) :: rdiscreteFBC
!</inputoutput>

!</subroutine>

    ! An identifier array for the equations to be tackled when discretising
    ! boundary conditions on fictitious boundary components
    integer, dimension(2) :: Iequations
    
    ! Add a new fictitious boundary object It should impose Dirichlet
    ! boundary conditions in the domain in the one and only solution component.
    ! We use the default initialisation of rfictBoundaryRegion and only
    ! change the name of the component.
    Iequations = (/1,2/)    ! 1=x, 2=y-velocity
    ! CALL bcasm_newDirichletBConFBD (rdiscretisation,Iequations,rdiscreteFBC,&
    !     getBoundaryValuesFBC,rcollection)

  end subroutine

end module
