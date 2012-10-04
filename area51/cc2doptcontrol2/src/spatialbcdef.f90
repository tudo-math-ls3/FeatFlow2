!##############################################################################
!# ****************************************************************************
!# <name> spatialbcdef </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for parsing boundary conditions from a DAT
!# file
!#
!# The following files can be found here:
!#
!# 1.) sbc_assembleBDconditions
!#     -> Assembles the definition of the BC's from sections given by DAT files
!#        and sets up an analytical boundary condition description
!#
!# 2.) sbc_assembleFBDconditions
!#     -> Assembles the definition of the BC's given by fictitious boundaries
!#        by evaluating sections given by DAT files
!#        and sets up an analytical boundary condition description
!#
!# </purpose>
!##############################################################################

module spatialbcdef

  use fsystem
  use genoutput
  use storage
  use boundary
  use basicgeometry
  use discretebc
  use discretefbc
  use bcassembly
  use bcassemblybase
  use triangulation
  use spatialdiscretisation
  use paramlist
  use mprimitives
  use derivatives
  use cubature
  use fparser
  use linearsystemblock
  use scalarpde
  use linearformevaluation
  use linearsystemscalar
  
  use collection
  use convection
    
  use constantsoptc
  use structuresoptc
  
  use timediscretisation

  use user_callback
  
  implicit none
  
  private
  
!<constants>

!<constantblock description="Names of the sections used in the main collection">

  ! Section name of the section saving the boundary expressions, i.e. the
  ! expressions that are to be evaluated in each point on the boundary.
  character(LEN=COLLCT_MLSECTION), parameter, public :: SEC_SBDEXPRESSIONS = "BDEXPRESSIONS"

  ! Section name of the section saving the boundary conditiuon definitions
  ! for all boundary segments.
  character(LEN=COLLCT_MLSECTION), parameter, public :: SEC_SBDCONDITIONS = "BDCONDITIONS"

  ! Name of the parser object for boundary value expressions
  character(LEN=COLLCT_MLSECTION), parameter, public :: BDC_BDPARSER = "BDEXPRPARSER"
!</constantblock>

!<constantblock description="The different types of boundary conditions supported by this module">

  ! Automatic, analytical. Defined via analytical reference function.
  integer, parameter, public :: BDC_USERDEFID = -3

  ! User defined expression, evaluated by calling the callback routine
  integer, parameter, public :: BDC_USERDEF = -2
  
  ! Text expression, evaluated by parsing
  integer, parameter, public :: BDC_EXPRESSION = -1
  
  ! Fixed double precision value
  integer, parameter, public :: BDC_VALDOUBLE = 0

  ! Fixed integer value
  integer, parameter, public :: BDC_VALINT    = 1

  ! Parabolic profile with prescribed maximum value
  integer, parameter, public :: BDC_VALPARPROFILE = 2
!</constantblock>


!<constantblock description="Assembly mode for boundary conditions">
  ! Assembles only standard boundary conditions
  integer, parameter, public :: SBC_BDC = 1

  ! Assembles Neumann boundary conditions.
  integer, parameter, public :: SBC_NEUMANN = 2

  ! Assembles Dirichlet boudnary control boundary conditions.
  integer, parameter, public :: SBC_DIRICHLETBCC = 4
  
  ! Assembles all boundary conditions.
  integer, parameter, public :: SBC_ALL = SBC_BDC + SBC_NEUMANN + SBC_DIRICHLETBCC
!</constantblock>


!<constantblock description="Variables in expressions">

  ! Basic variables that are allowed in expressions.
  ! Variables that are not defined in the actual situation are set to 0.
  !
  ! X,Y,Z - coordinate of a point (z=0 in 2D case),
  ! L     - local parameter value in the range [0..1],
  ! R     - parameter value of a boundary point, 0-1 parametrisation,
  ! S     - parameter value of a boundary point, arc length parametrisation,
  ! TIME  - current simulation time (=0 in stationary simulation)
  !
  ! Depending on the situation, this list may be extended by situation
  ! specific variables or variables that are only available at runtime.
  character(LEN=10), dimension(7), parameter, public :: SEC_EXPRVARIABLES = &
    (/'X    ','Y    ','Z    ','L    ','R    ','S    ','TIME '/)

!</constantblock>

!</constants>

!<types>

!<typeblock>
  
  ! Encapsules a region on the boundary where Neumann boundary conditions
  ! are present.
  type t_bdRegionEntry
  
    ! The boundary region where there are Neumann boudary conditions
    type(t_boundaryRegion) :: rboundaryRegion
    
    ! Pointer to the next Neumann boundaty region
    type(t_bdRegionEntry), pointer :: p_nextBdRegion
  
  end type
  
!</typeblock>

!<typeblock>

  ! Encapsules all Neumann boundary parts in the domain
  type t_boundaryRegionList
  
    ! Number of Neumann boundary regions in the primal equation
    integer :: nregionsPrimal = 0

    ! Number of Neumann boundary regions in the dual equation
    integer :: nregionsDual = 0
  
    ! Pointer to the head of a list of Neumann boundary regions, primal equation
    type(t_bdRegionEntry), pointer :: p_rprimalBdHead => null()

    ! Pointer to the head of a list of Neumann boundary regions, dual equation
    type(t_bdRegionEntry), pointer :: p_rdualBdHead => null()

    ! Pointer to the tail of a list of Neumann boundary regions, primal equation
    type(t_bdRegionEntry), pointer :: p_rprimalBdTail => null()

    ! Pointer to the tail of a list of Neumann boundary regions, dual equation
    type(t_bdRegionEntry), pointer :: p_rdualBdTail => null()
  
  end type

!</typeblock>

!<types>

  public :: t_bdRegionEntry
  public :: t_boundaryRegionList
  public :: sbc_assembleBDconditions
  public :: sbc_assembleFBDconditions
  public :: sbc_releaseBoundaryList

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_addBoundaryRegion (rboundaryRegion,rneumannBoundary,cprimaldual)
  
!<description>
  ! Adds a Neumann boundary region to the list of Neumann boundary regions.
!</desctiption>

!<input>
  ! Region to attach.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
  
  ! Whether the region is to be attached to the primal or dual Neumann boudary.
  ! =1: primal space, =2: dual space.
  integer, intent(in) :: cprimaldual
!</input>

!<inputoutput>
  ! List of Neumann boundary regions
  type(t_boundaryRegionList), intent(inout), target :: rneumannBoundary
!</inputoutput>

!</subroutine>

    select case (cprimaldual)
    case (1)
      ! Add to primal Neumann boundary
      if (rneumannBoundary%nregionsPrimal .eq. 0) then
        ! Create the structure if it does not exist.
        allocate (rneumannBoundary%p_rprimalBdTail)
        rneumannBoundary%p_rprimalBdHead => &
            rneumannBoundary%p_rprimalBdTail
      else
        ! Attach a new tail.
        allocate (rneumannBoundary%p_rprimalBdTail%p_nextBdRegion)
        rneumannBoundary%p_rprimalBdTail => &
            rneumannBoundary%p_rprimalBdTail%p_nextBdRegion
      end if
      
      ! Put the boundary region to there.
      rneumannBoundary%p_rprimalBdTail%rboundaryRegion = rboundaryRegion
      
      rneumannBoundary%nregionsPrimal = rneumannBoundary%nregionsPrimal + 1
      
    case (2)
      ! Add to dual Neumann boundary
      if (rneumannBoundary%nregionsDual .eq. 0) then
        ! Create the structure if it does not exist.
        allocate (rneumannBoundary%p_rdualBdTail)
        rneumannBoundary%p_rdualBdHead => &
            rneumannBoundary%p_rdualBdTail
      else
        ! Attach a new tail.
        allocate (rneumannBoundary%p_rdualBdTail%p_nextBdRegion)
        rneumannBoundary%p_rdualBdTail => &
            rneumannBoundary%p_rdualBdTail%p_nextBdRegion
      end if
      
      ! Put the boundary region to there.
      rneumannBoundary%p_rdualBdTail%rboundaryRegion = rboundaryRegion
      
      rneumannBoundary%nregionsDual = rneumannBoundary%nregionsDual + 1
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_releaseBoundaryList (rneumannBoundary)
  
!<description>
  ! Releases a rneumannBoundary structure that was allocated in
  ! sbc_assembleBDconditions.
!</desctiption>

!<inputoutput>
  ! Structure to release
  type(t_boundaryRegionList), intent(inout), target :: rneumannBoundary
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_bdRegionEntry), pointer :: p_rneumannRegion1,p_rneumannRegion2
    integer :: i
    
    ! Release the entries
    p_rneumannRegion1 => rneumannBoundary%p_rprimalBdHead
    do i=1,rneumannBoundary%nregionsPrimal
      p_rneumannRegion2 => p_rneumannRegion1
      p_rneumannRegion1 => p_rneumannRegion2%p_nextBDregion
      deallocate(p_rneumannRegion2)
    end do

    p_rneumannRegion1 => rneumannBoundary%p_rdualBdHead
    do i=1,rneumannBoundary%nregionsDual
      p_rneumannRegion2 => p_rneumannRegion1
      p_rneumannRegion1 => p_rneumannRegion2%p_nextBDregion
      deallocate(p_rneumannRegion2)
    end do
    
    ! Nothing inside anymore.
    rneumannBoundary%nregionsPrimal = 0
    rneumannBoundary%nregionsDual = 0

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine sbc_assembleBDconditions (roptcBDC,dtimePrimal,dtimeDual,&
      cbctype,rglobalData,casmFlags,rtimediscr,rspaceDiscr,&
      rdiscreteBC,rneumannBoundary,rdirichletControlBoundary,&
      rvectorDirichletBCCRHS)

!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
  !
  ! The bhasNeumannBoudary flag on every level is initialised according to
  ! whether there exist Neumann boundary components on the boundary or not.
!</description>
  
!<input>
  ! Boundary conditions in the problem.
  type(t_optcBDC), intent(in) :: roptcBDC

  ! Current simulation time for the primal equation.
  real(dp), intent(in) :: dtimePrimal

  ! Current simulation time for the dual equation.
  real(dp), intent(in) :: dtimeDual

  ! Type of boundary condition to assemble.
  ! CCSPACE_PRIMAL: BC for the primal space.
  ! CCSPACE_DUAL: BC for the dual space. If there are only 3 components
  !    in rspaceDiscr, the BC`s are assembled for a vector with 3 components
  !    which is assumed to be the dual space.
  ! CCSPACE_PRIMALDUAL: BC for primal and dual space.
  integer, intent(in) :: cbctype
  
  ! Global data, passed to callback routines
  type(t_globalData), intent(inout) :: rglobalData

  ! Assembly flags, determins what to assemble.
  ! SBC_BDC: Assembles all boundary conditions except Neumann BC's.
  ! SBC_NEUMANN: Assembles only Neumann boundary conditions. rneumannBoundary
  !     must be persent.
  ! SBC_ALL: Assembles both, standard and Neumann boundary conditions. Default.
  integer, intent(in) :: casmFlags
  
  ! A discretisation structure defining the underlying time discretisation.
  type(t_timeDiscretisation), intent(in) :: rtimeDiscr
  
  ! OPTIONAL: A discretisation structure defining the space discretisation
  ! of the current level.
  ! Can be omitted if only Neumann boundary is to be created.
  type(t_blockDiscretisation), intent(in) :: rspaceDiscr
  
  ! OPTIONAL: A t_discreteBC structure that receives a discretised version
  ! of the boundary boundary conditions. The structure should
  ! be empty; new BC's are simply added to the structure.
  ! Can be omitted if only Neumann boundary is to be created.
  type(t_discreteBC), intent(inout), optional :: rdiscreteBC
  
  ! OPTIONAL: A right-hand-side vector. If SBC_DIRICHLETBCC is specified,
  ! this RHS vector is modified according to the boundary conditions
  ! specified for the Dirichlet boundary control.
  type(t_vectorBlock), intent(inout), optional :: rvectorDirichletBCCRHS

!</input>

!<output>
  ! OPTIONAL: Returns a structure defining the Neumann boundary segments.
  ! Must be released by the caller using sbc_releaseBoundaryList!
  type(t_boundaryRegionList), intent(out), optional :: rneumannBoundary

  ! OPTIONAL: Returns a structure defining a list of boundary segments
  ! where boundary control is to be applied.
  ! Must be released by the caller using sbc_releaseBoundaryList!
  type(t_boundaryRegionList), intent(out), optional :: rdirichletControlBoundary
!</output>

!</subroutine>

    ! local variables
    integer :: i,ityp,ivalue,ibdComponent,isegment,iintervalEnds,iprimaldual
    integer :: ibctyp,icount,iexptyp,iid
    integer, dimension(2) :: IminIndex,imaxIndex
    integer, dimension(NDIM2D) :: IvelComp
    real(DP) :: dvalue,dpar1,dpar2
    character(LEN=PARLST_LENLINEBUF) :: cstr,cexpr,sbdex1,sbdex2,sbdex3,sbdex4
    character(LEN=PARLST_MLNAME) :: cname
    integer, dimension(NDIM2D) :: IvelEqns
    type(t_collection) :: rlocalCollection
    type(t_boundary), pointer :: p_rboundary
    type(t_linearForm) :: rlinformRhs
    
    ! A local collection we use for storing named parameters during the
    ! parsing process.
    type(t_collection) :: rcoll
    
    ! A collection structure passed to callback routines.
    type(t_collection), target :: rcallbackcollection
    
    ! Triangulation on currently highest level.
    type(t_triangulation), pointer :: p_rtriangulation

    ! A set of variables describing the analytic boundary conditions.
    type(t_boundaryRegion), target :: rboundaryRegion
    
    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection,p_rbdcond
    
    ! A compiled expression for evaluation at runtime
    type(t_fparser), target :: rparser
    
    select case (roptcBDC%p_rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
    
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
      call parlst_querysection(roptcBDC%p_rparamListBDC, &
          roptcBDC%ssectionBdExpressions, p_rsection)
      call parlst_querysection(roptcBDC%p_rparamListBDC, &
          roptcBDC%ssectionBdConditions, p_rbdcond)
      
      ! For intermediate storing of expression types, we use a local collection
      call collct_init (rcoll)
      
      ! Add a section to the collection that accepts the boundary expressions
      call collct_addsection (rcoll, roptcBDC%ssectionBdExpressions)
      
      ! Create a parser structure for as many expressions as configured
      call fparser_create (rparser,&
          parlst_querysubstrings (p_rsection, 'bdExpressions'))
      
      ! Add the parser to the collection
      call collct_setvalue_pars (rcoll, BDC_BDPARSER, rparser, &
                                  .true., 0, SEC_SBDEXPRESSIONS)
      
      ! Add the boundary expressions to the collection into the
      ! specified section.
      do i=1,parlst_querysubstrings (p_rsection, 'bdExpressions')
      
        call parlst_getvalue_string (p_rsection, 'bdExpressions', cstr, "", i)
        
        ! Get the type and decide on the identifier how to save the expression.
        read(cstr,*) cname,ityp
        
        select case (ityp)
        case (BDC_USERDEFID)
          ! Id of a hardcoded expression realised in the callback routines
          read(cstr,*) cname,ityp,ivalue
          call collct_setvalue_int (rcoll, cname, ivalue, .true., &
                                    0, SEC_SBDEXPRESSIONS)
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
                                     
        case default
          call output_line ('Expressions not implemented!', &
              OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
          call sys_halt()

        end select
        
        ! Put the type of the expression to the temporary collection section
        call collct_setvalue_int (rcoll, cname, ityp, .true.)
        
      end do
      
      ! Get the triangulation on the highest level
      p_rtriangulation => rspaceDiscr%p_rtriangulation
      p_rboundary => rspaceDiscr%p_rboundary
      
      ! Put some information to the quick access arrays for access
      ! in the callback routine.
      call collct_init(rcallbackcollection)

      rcoll%Iquickaccess(1) = 1 ! time-dependent
      rcoll%Dquickaccess(2) = rtimeDiscr%dtimeInit
      rcoll%Dquickaccess(3) = rtimeDiscr%dtimeMax
      
      rcoll%p_rnextCollection => rcallbackcollection
      
      ! Now to the actual boundary conditions.
      ! First primal (iprimaldual=1) then dual (iprimaldual=2)
      do iprimaldual = 1,2
        
        ! DquickAccess(4:) is reserved for BC specific information.
        !
        ! Put a link to a user defined collection into that local collection.
        ! That allows us to access it or to pass it to user defined callback
        ! functions.

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
                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                    ! Usually there's Neumann boundary in this region, but we can't be
                    ! sure. Check if, on the highest level, there's at least one edge
                    ! of the triangulation belonging to the boundary. If yes, we
                    ! have found Neumann boundary. If no, the segment is just too
                    ! small to be considered as Neumann boundary.
                    
                    !call bcasm_getEdgesInBCregion (p_rtriangulation,p_rboundary,&
                    !                              rboundaryRegion, &
                    !                              IminIndex,ImaxIndex,icount)
                    !if ((icount .gt. 0) .and. present(rneumannBoundary)) then
                    if (present(rneumannBoundary)) then
                      ! Add the bondary region to the Neumann boundary regions
                      ! if there is a structure present.
                      call sbc_addBoundaryRegion(&
                          rboundaryRegion,rneumannBoundary,iprimaldual)
                    end if
                  end if
                
                case (1)
                  if (iand(casmFlags,SBC_BDC) .ne. 0) then
                    ! Simple Dirichlet boundary.
                    ! Prescribed primal velocity; dual velocity os always zero.
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
                    
                    if (iprimaldual .eq. 1) then
                      ! Primal BC's
                      if ((cbctype .eq. CCSPACE_PRIMAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then
                      
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
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimePrimal
                          call user_initCollectForVecAssembly (&
                              rglobalData,iid,1,dtimePrimal,rcallbackcollection)
                          
                          ! Assemble the BC's.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)
                              
                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
                              
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
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimePrimal
                          call user_initCollectForVecAssembly (&
                              rglobalData,iid,2,dtimePrimal,rcallbackcollection)

                          ! Assemble the BC's.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)

                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        end if
                        
                      end if
                      
                    end if
                    
                    if (iprimaldual .eq. 2) then
                    
                      ! Dual BC's
                      if ((cbctype .eq. CCSPACE_DUAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then

                        ! Now the same thing again, this time separately for primal and dual
                        ! variables.
                        ! If a primal velocity is specified, Dirichlet-0-boundary conditions are
                        ! assumed for the corresponding dual.
                        if (sbdex1 .ne. '') then
                        
                          ! dual X-velocity if primal X-velocity exists
                          !
                          ! The 2nd element in IquickAccess saves the component number.
                          rcoll%IquickAccess(2) = 4
                          
                          ! Ditichlet-0 boundary
                          iexptyp = BDC_VALDOUBLE
                          rcoll%Dquickaccess(4) = 0.0_DP
                          rcoll%IquickAccess(3) = iexptyp
                          rcoll%SquickAccess(1) = ''
                          iid = 0
                          
                          rcoll%Dquickaccess(1) = dtimeDual
                          call user_initCollectForVecAssembly (rglobalData,iid,4,dtimeDual,rcallbackcollection)

                          ! Assemble the BC's.
                          ! If we only assemble dual BC's and there are 3 solution components,
                          ! we assume the vector to specify exactly the dual solution.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)
                              
                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
                              
                        end if
                        
                        if (sbdex2 .ne. '') then
                        
                          ! dual Y-velocity if primal Y-velocity exists
                          !
                          ! The 2nd element in IquickAccess saves the component number.
                          rcoll%IquickAccess(2) = 5
                          
                          ! Ditichlet-0 boundary
                          iexptyp = BDC_VALDOUBLE
                          rcoll%Dquickaccess(4) = 0.0_DP
                          rcoll%IquickAccess(3) = iexptyp
                          rcoll%SquickAccess(1) = ''
                          iid = 0
                          
                          rcoll%Dquickaccess(1) = dtimeDual
                          call user_initCollectForVecAssembly (rglobalData,iid,5,dtimeDual,rcallbackcollection)

                          ! Assemble the BC's.
                          ! If we only assemble dual BC's and there are 3 solution components,
                          ! we assume the vector to specify exactly the dual solution.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)

                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
                              
                        end if
                        
                      end if
                      
                    end if
                    
                  end if
                  
                case (2)
                
                  if (iand(casmFlags,SBC_BDC) .ne. 0) then
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
                            collct_getvalue_real (rcoll,sbdex1, 0, roptcBDC%ssectionBdExpressions)
                      case (BDC_EXPRESSION)
                        ! Expression. Write the identifier for the expression
                        ! as itag into the boundary condition structure.
                        rcoll%IquickAccess(4) = &
                            collct_getvalue_int (rcoll,sbdex1, 0, roptcBDC%ssectionBdExpressions)
                      end select
                    
                      IvelEqns = (/1,2/)
                      call bcasm_newPdropBConRealBd (&
                          rspaceDiscr,IvelEqns,rboundaryRegion,rdiscreteBC,&
                          cc_getBDconditionsNavSt2D,rcoll)
                    end if
                  
                  end if
                  
                case (3)
                  
                  if (iand(casmFlags,SBC_BDC) .ne. 0) then
                    ! Nonlinear slip boundary conditions.
                    IvelComp = (/1,2/)
                    call bcasm_newSlipBConRealBd (&
                        rspaceDiscr,IvelEqns(1:NDIM2D),rboundaryRegion,rdiscreteBC)
                  end if
                
                case (4)
                  if (iand(casmFlags,SBC_BDC) .ne. 0) then
                    ! Simple Dirichlet boundary, separate definition for all
                    ! solution components.
                    ! Read the line again, get the expressions for X- and Y-velocity
                    read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2,sbdex3,sbdex4
                    
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
                    
                    if (iprimaldual .eq. 1) then
                    
                      ! Primal BC's
                      if ((cbctype .eq. CCSPACE_PRIMAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then
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
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimePrimal
                          call user_initCollectForVecAssembly (&
                              rglobalData,iid,1,dtimePrimal,rcallbackcollection)

                          ! Assemble the BC's.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)

                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
                              
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
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimePrimal
                          call user_initCollectForVecAssembly (&
                              rglobalData,iid,2,dtimePrimal,rcallbackcollection)

                          ! Assemble the BC's.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)

                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        end if
                      end if
                      
                    end if
                    
                    if (iprimaldual .eq. 2) then
                    
                      ! Dual BC's
                      if ((cbctype .eq. CCSPACE_DUAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then
                      
                        ! Now the same thing again, this time separately for primal and dual
                        ! variables.
                        ! If a velocity is not specified, Dirichlet-0-boundary conditions are
                        ! assumed.
                        
                        if (sbdex3 .ne. '') then
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
                          select case (iexptyp)
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimeDual
                          call user_initCollectForVecAssembly (rglobalData,iid,4,dtimeDual,rcallbackcollection)

                          ! Assemble the BC's.
                          ! If we only assemble dual BC's and there are 3 solution components,
                          ! we assume the vector to specify exactly the dual solution.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)
                              
                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        end if
                        
                        if (sbdex4 .ne. '') then
                        
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
                          select case (iexptyp)
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex4, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimeDual
                          call user_initCollectForVecAssembly (rglobalData,iid,5,dtimeDual,rcallbackcollection)

                          ! Assemble the BC's.
                          ! If we only assemble dual BC's and there are 3 solution components,
                          ! we assume the vector to specify exactly the dual solution.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)
                              
                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        end if
                        
                      end if
                      
                    end if
                    
                  end if

                case (5)
                  ! Simple Dirichlet boundary for the primal equation.
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
                  
                  if (iand(casmFlags,SBC_BDC) .ne. 0) then
                    
                    if (iprimaldual .eq. 1) then
                    
                      ! Primal BC's
                      if ((cbctype .eq. CCSPACE_PRIMAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then
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
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimePrimal
                          call user_initCollectForVecAssembly (&
                              rglobalData,iid,1,dtimePrimal,rcallbackcollection)

                          ! Assemble the BC's.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)

                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
                              
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
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimePrimal
                          call user_initCollectForVecAssembly (&
                              rglobalData,iid,2,dtimePrimal,rcallbackcollection)

                          ! Assemble the BC's.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)

                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        end if
                      end if
                      
                    end if
                  
                  end if
                                 
                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                    ! Neumann boundary for the dual.
                    if (iprimaldual .eq. 2) then
                      !call bcasm_getEdgesInBCregion (p_rtriangulation,p_rboundary,&
                      !                              rboundaryRegion, &
                      !                              IminIndex,ImaxIndex,icount)
                      !if ((icount .gt. 0) .and. present(rneumannBoundary)) then
                      if (present(rneumannBoundary)) then
                        ! Add the bondary region to the Neumann boundary regions
                        ! if there is a structure present.
                        call sbc_addBoundaryRegion(&
                            rboundaryRegion,rneumannBoundary,iprimaldual)
                      end if
                    end if
                  end if
                                    
                case (6)
                  ! Simple Dirichlet boundary for the dual equation
                  ! Read the line again, get the expressions for X- and Y-velocity
                  read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex3,sbdex4
                  
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
                  
                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                    ! Neumann boundary for the primal.
                    if (iprimaldual .eq. 1) then
                      !call bcasm_getEdgesInBCregion (p_rtriangulation,rboundary,&
                      !                              rboundaryRegion, &
                      !                              IminIndex,ImaxIndex,icount)
                      !if ((icount .gt. 0) .and. present(rneumannBoundary)) then
                      if (present(rneumannBoundary)) then
                        ! Add the bondary region to the Neumann boundary regions
                        ! if there is a structure present.
                        call sbc_addBoundaryRegion(&
                            rboundaryRegion,rneumannBoundary,iprimaldual)
                      end if
                    end if
                  end if

                  if (iand(casmFlags,SBC_BDC) .ne. 0) then
                    if (iprimaldual .eq. 2) then
                    
                      ! Dual BC's
                      if ((cbctype .eq. CCSPACE_DUAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then
                      
                        ! Now the same thing again, this time separately for primal and dual
                        ! variables.
                        ! If a velocity is not specified, Dirichlet-0-boundary conditions are
                        ! assumed.
                        
                        if (sbdex3 .ne. '') then
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
                          select case (iexptyp)
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimeDual
                          call user_initCollectForVecAssembly (rglobalData,iid,4,dtimeDual,rcallbackcollection)

                          ! Assemble the BC's.
                          ! If we only assemble dual BC's and there are 3 solution components,
                          ! we assume the vector to specify exactly the dual solution.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)
                              
                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        end if
                        
                        if (sbdex4 .ne. '') then
                        
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
                          select case (iexptyp)
                          case (BDC_USERDEFID)
                            iid = collct_getvalue_int (rcoll, sbdex4, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                            ! Constant or parabolic profile
                            rcoll%Dquickaccess(4) = &
                                collct_getvalue_real (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
                          case (BDC_EXPRESSION)
                            ! Expression. Write the identifier for the expression
                            ! as itag into the boundary condition structure.
                            rcoll%IquickAccess(4) = &
                                collct_getvalue_int (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
                          end select
                        
                          rcoll%Dquickaccess(1) = dtimeDual
                          call user_initCollectForVecAssembly (rglobalData,iid,5,dtimeDual,rcallbackcollection)

                          ! Assemble the BC's.
                          ! If we only assemble dual BC's and there are 3 solution components,
                          ! we assume the vector to specify exactly the dual solution.
                          call bcasm_newDirichletBConRealBD (&
                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                              cc_getBDconditionsNavSt2D,rcoll)
                              
                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        end if
                        
                      end if
                      
                    end if
                    
                  end if
                  
                case (7)
                  
                  ! Dirichlet boundary control.
                  !
                  ! This is treated like Neumann boundary conditions here
                  ! but leads to a different implementation.
                  ! The corresponding region is added to the Neumann boundary
                  ! structure as well as to the Dirichlet boundary control
                  ! structure.
                  ! Adding it to the Neumann boundary conditions structure
                  ! is important, otherwise the components detect the problem
                  ! as pure Dirichlet and start filtering the pressure... 
!                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
!                    if (present(rneumannBoundary)) then
!                      call sbc_addBoundaryRegion(&
!                          rboundaryRegion,rneumannBoundary,iprimaldual)
!                    end if
!                  end if

                  ! For the primal equation, there is a special coupling
                  ! to the dual equation.
                  !
                  ! For the dual equation, we have Dirichlet-0 boundary conditions here.

                  if (iand(casmFlags,SBC_DIRICHLETBCC) .ne. 0) then
                  
                    if (iprimaldual .eq. 1) then
                      if (present(rdirichletControlBoundary)) then
                        ! Add the bondary region to the Dirichlet BCC boundary regions
                        ! if there is a structure present.
                        call sbc_addBoundaryRegion(&
                            rboundaryRegion,rdirichletControlBoundary,iprimaldual)
                      end if
                      
                      if (present(rvectorDirichletBCCRHS)) then
                        ! Modify the RHS according to the Dirichlet boundary control
                        ! conditions. Use a penalty approach to implement the Dirichlet
                        ! values. A linear form has to be applied to the RHS vector
                        ! in each component which incorporates the basic Dirichlet value.
                        ! The control then acts relative to the incoroprated value.
                        
                        ! Prepare the linear form.
                        rlinformRhs%itermCount = 1
                        rlinformRhs%Idescriptors(1) = DER_FUNC2D
                        rlinformRhs%Dcoefficients(1:rlinformRhs%itermCount)  = 1.0_DP
                        
                        ! Read the boundary data specification and prepare the callback routine
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
                        
                        ! The current boundary region is put to the collection
                        call collct_setvalue_bdreg (rcoll, "BDREG", rboundaryRegion, .true.)
                        
                        ! Primal BC's
                        if ((cbctype .eq. CCSPACE_PRIMAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then
                        
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
                            case (BDC_USERDEFID)
                              iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                            case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                              ! Constant or parabolic profile
                              rcoll%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%ddirichletBCPenalty * &
                                  collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                            case (BDC_EXPRESSION)
                              ! Expression. Write the identifier for the expression
                              ! as itag into the boundary condition structure.
                              rcoll%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%ddirichletBCPenalty
                              rcoll%IquickAccess(4) = &
                                  collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
                            end select
                          
                            rcoll%Dquickaccess(1) = dtimePrimal
                            call user_initCollectForVecAssembly (&
                                rglobalData,iid,1,dtimePrimal,rcallbackcollection)
                            
                            call linf_buildVectorScalarBdr2D (rlinformRhs, CUB_G4_1D, .false., &
                                rvectorDirichletBCCRHS%RvectorBlock(1),&
                                fcoeff_buildBCCRHS,rboundaryRegion, rcoll)
                                
                            call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
                                
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
                            case (BDC_USERDEFID)
                              iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
                            case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
                              ! Constant or parabolic profile
                              rcoll%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%ddirichletBCPenalty * &
                                  collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                            case (BDC_EXPRESSION)
                              ! Expression. Write the identifier for the expression
                              ! as itag into the boundary condition structure.
                              rcoll%IquickAccess(4) = &
                                  collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
                            end select
                          
                            rcoll%Dquickaccess(1) = dtimePrimal
                            call user_initCollectForVecAssembly (&
                                rglobalData,iid,2,dtimePrimal,rcallbackcollection)

                            ! Assemble the BC's.
                            call linf_buildVectorScalarBdr2D (rlinformRhs, CUB_G4_1D, .false., &
                                rvectorDirichletBCCRHS%RvectorBlock(1),&
                                fcoeff_buildBCCRHS,rboundaryRegion, rcoll)

                            call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                          end if
                          
                        end if
                        
                        ! Delete the boundary region again
                        call collct_deletevalue (rcoll,"BDREG")
                        
                      end if

                    end if
                    
                  end if
                  
                  if (iand(casmFlags,SBC_BDC) .ne. 0) then
                    if (iprimaldual .eq. 2) then
                    
                      ! Dual BC's
                      if ((cbctype .eq. CCSPACE_DUAL) .or. (cbctype .eq. CCSPACE_PRIMALDUAL)) then

                        ! dual X-velocity if primal X-velocity exists
                        
                        ! Ditichlet-0 boundary
                        iexptyp = BDC_VALDOUBLE
                        rcoll%Dquickaccess(4) = 0.0_DP
                        rcoll%IquickAccess(3) = iexptyp
                        rcoll%SquickAccess(1) = ''
                        iid = 0
                        
                        rcoll%Dquickaccess(1) = dtimeDual
                        call user_initCollectForVecAssembly (rglobalData,iid,4,dtimeDual,rcallbackcollection)

                        ! Assemble the BC's.
                        ! The 2nd element in IquickAccess saves the component number.
                        ! If we only assemble dual BC's and there are 3 solution components,
                        ! we assume the vector to specify exactly the dual solution.

                        rcoll%IquickAccess(2) = 4
                        call bcasm_newDirichletBConRealBD (&
                            rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                            cc_getBDconditionsNavSt2D,rcoll)

                        call user_doneCollectForAssembly (rglobalData,rcallbackcollection)

                        call user_initCollectForVecAssembly (rglobalData,iid,5,dtimeDual,rcallbackcollection)

                        rcoll%IquickAccess(2) = 5
                        call bcasm_newDirichletBConRealBD (&
                            rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
                            cc_getBDconditionsNavSt2D,rcoll)
                            
                        call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
                            
                      end if
                      
                    end if
                  
                  end if

                case default
                  call output_line ('Unknown boundary condition!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
                  call sys_halt()
                end select
                
                ! Move on to the next parameter value
                dpar1 = dpar2
                
              end if
                                                
            end do
          
          end if
          
        end do

      end do

      ! Release the parser object with all the expressions to be evaluated
      ! on the boundary.
      call fparser_release (rparser)

      ! Remove the boundary value parser from the collection
      call collct_deletevalue (rcoll, BDC_BDPARSER)
      
      ! Remove the boundary-expression section we added earlier,
      ! with all their content.
      call collct_deletesection(rcoll,roptcBDC%ssectionBdExpressions)

      ! Remove the temporary collection from memory.
      call collct_done(rcallbackcollection)
      call collct_done(rcoll)

    end select
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_buildBCCRHS (rdiscretisation, rform, &
                nelements, npointsPerElement, Dpoints, ibct, DpointPar, &
                IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
  
  use basicgeometry
  use collection
  use domainintegration
  use fsystem
  use scalarpde
  use spatialdiscretisation
  use triangulation
  
!<description>
  ! This subroutine is called during the vector assembly. It has to compute
  ! the coefficients in front of the terms of the linear form.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points and all the terms in the linear form
  ! the corresponding coefficients in front of the terms.
!</description>
  
!<input>
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  
  ! The linear form which is currently to be evaluated:
  type(t_linearForm), intent(in) :: rform
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  ! DIMENSION(dimension,npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! This is the number of the boundary component that contains the
  ! points in Dpoint. All points are on the same boundary component.
  integer, intent(in) :: ibct

  ! For every point under consideration, this specifies the parameter
  ! value of the point on the boundary component. The parameter value
  ! is calculated in LENGTH PARAMETRISATION!
  ! DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(in) :: DpointPar

  ! An array accepting the DOF`s on all elements in the test space.
  ! DIMENSION(\#local DOF`s in test space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! A list of all coefficients in front of all terms in the linear form -
  ! for all given points on all given elements.
  !   DIMENSION(itermCount,npointsPerElement,nelements)
  ! with itermCount the number of terms in the linear form.
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
  
!</subroutine>
  
    ! local variables
    integer, dimension(1) :: Icomponents
    real(DP), dimension(1) :: Dvalues
    integer :: i,j
    type(t_boundaryRegion), pointer :: p_rboundaryRegion
    real(DP) :: dpar, ddirichletBCPenalty

    ! Penalty parameter
    ddirichletBCPenalty = rcollection%Dquickaccess(4)

    ! Get the boundary region from the collection
    p_rboundaryRegion => collct_getvalue_bdreg (rcollection, "BDREG")

    ! Loop through all points and calculate the values.
    do i=1,nelements
      do j=1,npointsPerElement
        ! Current component
        Icomponents(1) = rcollection%IquickAccess(2)
        
        ! Calculate the parameter value in 0-1 parametrisation
        dpar = boundary_convertParameter(rdiscretisation%p_rboundary, ibct, &
            DpointPar(j,i), BDR_PAR_LENGTH, BDR_PAR_01)
        
        call cc_getBDconditionsNavSt2D(Icomponents,rdiscretisation,&
            p_rboundaryRegion,0,DISCBC_NEEDFUNC,ibct,dpar,&
            Dvalues,rcollection)
            
        ! Value must be mutiplied by the penalty parameter
        Dcoefficients(1,j,i) = ddirichletBCPenalty*Dvalues(1)
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getBDconditionsNavSt2D (Icomponents,rspaceDiscr,rboundaryRegion,ielement, &
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
  type(t_spatialDiscretisation), intent(IN)                   :: rspaceDiscr
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(IN)                                    :: ielement
  
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
  integer, intent(IN)                                    :: iwhere

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
  type(t_collection), intent(inout), optional      :: rcollection

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
      case (1,4,5,6,7)
        ! Simple Dirichlet BC's for either primal or primal/dual equation.
        ! Evaluate the expression iexprtyp.
        Dvalues(1) = evalBoundary (icomponent,rspaceDiscr, rboundaryRegion, &
            iexprtyp, rcollection%IquickAccess(4), rcollection%DquickAccess(4), &
            dwhere, rcollection%SquickAccess(1),dtime,&
            rcollection)
    
      case (2)
        ! Normal stress / pressure drop. Evaluate the  expression iexprtyp.
        Dvalues(1) = evalBoundary (icomponent,rspaceDiscr, rboundaryRegion, &
            iexprtyp, rcollection%IquickAccess(4), rcollection%DquickAccess(4), &
            dwhere, rcollection%SquickAccess(1),dtime,&
            rcollection)
            
      end select
      
    end select
  
  contains
  
    ! Auxiliary function: Evaluate a scalar expression on the boundary.
    
    real(DP) function evalBoundary (icomponent,rspaceDiscr, rboundaryRegion, &
                                    ityp, ivalue, dvalue, dpar, stag, dtime, rcollection)
    
    ! Solution component for which the expression is evaluated.
    ! 1 = X-velocity, 2 = y-velocity,...
    integer, intent(IN) :: icomponent
    
    ! Discretisation structure of the underlying discretisation
    type(t_spatialDiscretisation), intent(IN) :: rspaceDiscr
    
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
    type(t_collection), optional :: rcollection

      ! local variables
      real(DP) :: d,dx,dy
      character(LEN=PARLST_LENLINEBUF) :: sexpr
      real(DP), dimension(size(SEC_EXPRVARIABLES)) :: Rval
      
      select case (ityp)
      case (BDC_USERDEFID)
        ! Defined via analytic reference solution, identified by an id.
        ! The id is already specified in the collection.
        
        ! Call the user defined callback routine to evaluate the boudary value.
        ! No string tag is given which indicates that a reference function is to be
        ! evaluated.
        !
        ! As collection, we pass rcollection%p_rcollection here; this is a pointer
        ! to the application specific, global collection that may be of interest for
        ! callback routines. rcollection itself is actually a 'local' collection,
        ! a 'wrapper' for the rcollection of the application!
        call user_getBoundaryValues (&
            "",icomponent,rspaceDiscr,rboundaryRegion,&
            dpar, d, rcollection%p_rnextCollection)
        
        evalBoundary = d
        
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
        call user_getBoundaryValues (&
            stag,icomponent,rspaceDiscr,rboundaryRegion,&
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
        
        call boundary_getCoords(rspaceDiscr%p_rboundary, &
                                rboundaryRegion%iboundCompIdx, &
                                dpar, dx, dy)
        
        ! Get the local parameter value 0 <= d <= 1.
        ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
        ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
        ! although 0 <= dminpar <= max.par
        !      and 0 <= dmaxpar <= max.par!
        d = dpar
        if (d .lt. rboundaryRegion%dminParam) &
          d = d + boundary_dgetMaxParVal(rspaceDiscr%p_rboundary,&
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
        Rval(6) = boundary_convertParameter(rspaceDiscr%p_rboundary, &
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
          d = d + boundary_dgetMaxParVal(rspaceDiscr%p_rboundary,&
                                         rboundaryRegion%iboundCompIdx)
        d = d - rboundaryRegion%dminParam
    
        ! Normalise to 0..1 using the length of the parameter region.
        ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
        d = d / (rboundaryRegion%dmaxParam - rboundaryRegion%dminParam)
        
        evalBoundary = mprim_getParabolicProfile (d,1.0_DP,dvalue)
        
      case default
        call output_line ('Invalid boundary condition type!', &
            OU_CLASS_ERROR,OU_MODE_STD,'evalBoundary')
        call sys_halt()

      end select
    
    end function

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_assembleFBDconditions (dtime,rspaceDiscr,rtimeDiscr,cbctype,&
      rdiscreteFBC,rglobalData)

!<description>
  ! This parses the boundary conditions for fictitious boundary
  ! components in the problem and assembles them into rdiscreteFBC.
!</description>
  
!<input>
  ! A discretisation structure defining the discretisation of the current
  ! level.
  type(t_blockDiscretisation), intent(IN) :: rspaceDiscr

  ! A discretisation structure defining the underlying time discretisation.
  type(t_timeDiscretisation), intent(in) :: rtimeDiscr

  ! Global data, passed to callback routines
  type(t_globalData), intent(inout) :: rglobalData
!</input>

!<inputoutput>
  ! Current simulation time.
  real(dp), intent(in) :: dtime

  ! Type of boundary condition to assemble.
  ! CCSPACE_PRIMAL: BC for the primal space.
  ! CCSPACE_DUAL: BC for the dual space. If there are only 3 components
  !    in rspaceDiscr, the BC`s are assembled for a vector with 3 components
  !    which is assumed to be the dual space.
  ! CCSPACE_PRIMALDUAL: BC for primal and dual space.
  integer, intent(in) :: cbctype

  ! A t_discreteFBC structure that receives a discretised version
  ! of the fictitious boundary boundary conditions. The structure should
  ! be empty; new BC's are simply added to the structure.
  type(t_discreteFBC), intent(INOUT) :: rdiscreteFBC
!</inputoutput>

!</subroutine>

    ! An identifier array for the equations to be tackled when discretising
    ! boundary conditions on fictitious boundary components
    integer, dimension(2) :: Iequations
    
    Iequations = (/1,2/)    ! 1=x, 2=y-velocity
    ! CALL bcasm_newDirichletBConFBD (rspaceDiscr,Iequations,&
    !     user_getBoundaryValuesFBC,rcollection)

  end subroutine
  
end module
