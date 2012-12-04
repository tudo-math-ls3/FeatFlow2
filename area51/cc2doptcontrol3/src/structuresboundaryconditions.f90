!##############################################################################
!# ****************************************************************************
!# <name> structuresboundaryconditions </name>
!# ****************************************************************************
!#
!# <purpose>
!# Structures that encapsule the boundary conditions.
!# </purpose>
!##############################################################################

module structuresboundaryconditions

  use fsystem
  use storage
  use genoutput
  use paramlist
  use collection
  use fparser
  use mprimitives
  use spatialdiscretisation
  
  use boundary
  use discretebc
  use bcassemblybase
  
  use structuresdiscretisation
  use structuresoptcontrol
  
  implicit none
  
  private
  
!<constants>

!<constantblock description="Names of the sections used in the main collection">

  ! Section name of the section saving the boundary expressions, i.e. the
  ! expressions that are to be evaluated in each point on the boundary.
  character(LEN=COLLCT_MLSECTION), parameter, public :: SEC_SBDEXPRESSIONS = "BDEXPRESSIONS"

  ! Section name of the section saving the boundary conditiuon definitions
  ! for all boundary segments. Primal equation.
  character(LEN=COLLCT_MLSECTION), parameter, public :: SEC_SBDCONDITIONS_PRIM = "BDCONDITIONS_PRIMAL"

  ! Section name of the section saving the boundary conditiuon definitions
  ! for all boundary segments. Dual equation.
  character(LEN=COLLCT_MLSECTION), parameter, public :: SEC_SBDCONDITIONS_DUAL = "BDCONDITIONS_DUAL"

  ! Name of the parser object for boundary value expressions
  character(LEN=COLLCT_MLSECTION), parameter, public :: BDC_BDPARSER = "BDEXPRPARSER"
!</constantblock>

!<constantblock description="The different types of boundary conditions supported by this module">

  ! Automatic, analytical. Defined via analytical reference function.
  ! the associated "ivalue" is the id of the expression.
  integer, parameter, public :: BDC_USERDEFID = -3

  ! User defined expression, evaluated by calling the callback routine
  ! The associated "svalue" is the name of the expression.
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


!<constantblock description="Variables in expressions">

  ! Basic variables that are allowed in expressions.
  ! Variables that are not defined in the actual situation are set to 0.
  !
  ! X,Y,Z - coordinate of a point (z=0 in 2D case),
  ! L     - local parameter value in the range [0..1],
  ! R     - parameter value of a boundary point, 0-1 parametrisation,
  ! S     - parameter value of a boundary point, arc length parametrisation,
  ! TIME  - current simulation time
  !
  ! Depending on the situation, this list may be extended by situation
  ! specific variables or variables that are only available at runtime.
  character(LEN=10), dimension(9), parameter, public :: SEC_EXPRVARIABLES = &
    (/"X    ","Y    ","Z    ","L    ","R    ","S    ","TIME ", &
      "NX   ","NY   "/)

!</constantblock>

!<constantblock description="Assembly mode for boundary conditions">
  ! Assembles discrete boundary conditions.
  integer, parameter, public :: SBC_DISCRETEBC = 1

  ! Assembles strong Dirichlet boundary conditions
  integer, parameter, public :: SBC_DIRICHLETBC = 2

  ! Assembles Neumann boundary conditions.
  integer, parameter, public :: SBC_NEUMANN = 4

  ! Assembles Dirichlet boudary control boundary conditions.
  integer, parameter, public :: SBC_DIRICHLETBCC = 8
  
  ! Assembles all boundary conditions.
  integer, parameter, public :: SBC_ALL = SBC_DISCRETEBC + SBC_DIRICHLETBC + SBC_NEUMANN + SBC_DIRICHLETBCC
!</constantblock>

!</constants>
  
!<types>

!<typeblock>

  ! Type encapsuling the boundary conditions.
  type t_optcBDC
  
    ! Physics of the problem
    type(t_settings_physics), pointer :: p_rphysics => null()
    
    ! Reference to the optimal control parameters
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl => null()

    ! A parameter list that saves the boundary conditions.
    type(t_parlist), pointer :: p_rparList => null()

    ! A parser structure for compiled expression to be evaluated at runtime
    type(t_fparser) :: rparser
    
    ! A collection that saves expression types
    type(t_collection) :: rcollExprTypes

    ! Name of the section defining a set of expressions to be used
    ! for being evaluated on the boundary.
    character(len=SYS_STRLEN) :: ssectionBdExpr = ""

    ! Name of the section in rparamList describing the boundary conditions.
    ! Primal and dual equations.
    character(len=SYS_STRLEN) :: ssectionBdCondPrim = ""
    character(len=SYS_STRLEN) :: ssectionBdCondDual = ""

    ! Name of the section in rparamList describing the boundary conditions
    ! in the linearised equation
    character(len=SYS_STRLEN) :: ssectionBdCondPrimLin = ""
    character(len=SYS_STRLEN) :: ssectionBdCondDualLin = ""
    
    ! Name of the section in rparamList describing the boundary conditions
    ! in the Poincare-Steklov operator
    character(len=SYS_STRLEN) :: ssectionBdCondPCSteklov = ""

  end type

!</typeblock>

  public :: t_optcBDC

!<typeblock>

  ! Pointer to a t_optcBDC structure.
  type t_p_optcBDC
  
    type(t_optcBDC), pointer :: p_roptcBDC => null()
    
    ! Pointer to a boundary region, for inhomogeneous Neumann BC
    type(t_boundaryRegion), pointer :: p_rboundaryRegion => null()

  end type

!</typeblock>

  public :: t_p_optcBDC

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

  public :: t_bdRegionEntry

!<typeblock>

  ! Encapsules a set of boundary regions in the domain
  type t_boundaryRegionList
  
    ! Number of boundary regions
    integer :: nregions = 0

    ! Pointer to the head of a list of boundary regions
    type(t_bdRegionEntry), pointer :: p_rbdHead => null()

    ! Pointer to the tail of a list of boundary regions
    type(t_bdRegionEntry), pointer :: p_rbdTail => null()
  
  end type

!</typeblock>

  public :: t_boundaryRegionList

!<typeblock>

  ! Structure encapsuling discrete boundary conditions in space.
  type t_optcBDCSpace
  
    ! Discrete (Dirichlet-) boundary conditions, primal or dual space
    type(t_discreteBC) :: rdiscreteBC

    ! Boundary regions with do-nothing Neumann boundary conditions.
    type(t_boundaryRegionList) :: rneumannBoundary

    ! Boundary regions with Dirichlet boundary conditions.
    type(t_boundaryRegionList) :: rdirichletBoundary
  
    ! Boundary regions with L2 Dirichlet control boundary
    type(t_boundaryRegionList) :: rdirichletControlBoundaryL2

    ! Boundary regions with $H^{1/2}$ Dirichlet control boundary
    type(t_boundaryRegionList) :: rdirichletControlBoundaryH12
    
  end type
  
!</typeblock>

  public :: t_optcBDCSpace

!<typeblock>

  ! A hierarchy of discrete boundary condition structures in space.
  ! For every level of a hierarchy of FEM spaces, this structure contains
  ! a t_optcBDCSpace encapsuling the corresponding Boundary conditions in space.
  type t_optcBDCSpaceHierarchy
  
    ! Minimum level
    integer :: nlmin = 0
    
    ! Maximum level
    integer :: nlmax = 0
    
    ! Boundary conditions structures for all levels
    type(t_optcBDCSpace), dimension(:), pointer :: p_RoptcBDCspace => null()
  
  end type
  
!</typeblock>

  public :: t_optcBDCSpaceHierarchy

!</types>

  ! Initialises a boundary condition structure.
  public :: struc_initBDC

  ! Returns the type and information of a boundary expression.
  public :: struc_getBdExprInfo

  ! Cleans up information in the optimal control structure.
  public :: struc_doneBDC

  ! Adds boundary region to a list of boundary regions.
  public :: sbc_addBoundaryRegion
  
  ! Releases a rregionList structure that set up with
  ! sbc_addBoundaryRegion.
  public :: sbc_releaseBdRegionList 

  ! Determine all DOFs in a boundary region list.
  public :: sbc_getDofsInBdRegionList

  ! Obtains a pointer to the discrete boundary condition structure
  ! on level ilevel.
  public :: sbch_getDiscreteBC

  ! Obtain information about a bonudary condition definition
  public :: struc_getBDCinfo
  
  ! Obtain information about a boundary segment
  public :: struc_getSegmentInfo
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine struc_initBDC (roptcBDC,rparlist,rphysics,rsettingsOptControl,&
      ssectionBdExpr,ssectionBdCondPrim,ssectionBdCondDual,&
      ssectionBdCondPrimLin,ssectionBdCondDualLin,ssectionBoundaryCondPCSteklov)
  
!<description>
  ! Initialises a boundary condition structure.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in), target :: rparlist

  ! Structure defining the physics.
  type(t_settings_physics), target :: rphysics
  
  ! Optimal control parameter structure
  type(t_settings_optcontrol), target :: rsettingsOptControl

  ! Name of the section defining a set of expressions to be used
  ! for being evaluated on the boundary.
  character(len=*), intent(in) :: ssectionBdExpr

  ! Section defining the boundary conditions for the primal equations
  character(len=*), intent(in) :: ssectionBdCondPrim
  
  ! Section defining the boundary conditions for the dual equations
  character(len=*), intent(in) :: ssectionBdCondDual

  ! Section defining the boundary conditions for the linearised primal equations
  character(len=*), intent(in) :: ssectionBdCondPrimLin
  
  ! Section defining the boundary conditions for the linearised dual equations
  character(len=*), intent(in) :: ssectionBdCondDualLin

  ! Section defining the boundary conditions for the Poincare-Srteklov operator
  character(len=*), intent(in) :: ssectionBoundaryCondPCSteklov
!</input>

!<inputoutput>
  ! Boundary condition structure to be initialised
  type(t_optcBDC), intent(out) :: roptcBDC
!</inputoutput>

!</subroutine>

    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection
    integer :: ivalue, i, ityp
    real(DP) :: dvalue
    character(LEN=PARLST_LENLINEBUF) :: sexpr, sstr, sname

    roptcBDC%p_rphysics => rphysics
    roptcBDC%p_rparList => rparlist
    roptcBDC%p_rsettingsOptControl => rsettingsOptControl
    roptcBDC%ssectionBdExpr     = ssectionBdExpr
    roptcBDC%ssectionBdCondPrim = ssectionBdCondPrim
    roptcBDC%ssectionBdCondDual = ssectionBdCondDual
    roptcBDC%ssectionBdCondPrimLin = ssectionBdCondPrimLin
    roptcBDC%ssectionBdCondDualLin = ssectionBdCondDualLin
    roptcBDC%ssectionBdCondPCSteklov = ssectionBoundaryCondPCSteklov

    ! Create a parser structure for as many expressions as configured
    call parlst_querysection(rparlist, ssectionBdExpr, p_rsection)

    call fparser_create (roptcBDC%rparser,&
        parlst_querysubstrings (p_rsection, "bdExpressions"))
    
    ! For intermediate storing of expression types, we use a local collection
    call collct_init (roptcBDC%rcollExprTypes)
    
    ! Create a parser structure for as many expressions as configured
    call fparser_create (roptcBDC%rparser,&
        parlst_querysubstrings (p_rsection, "bdExpressions"))
    
    ! Add the boundary expressions to the collection into the
    ! specified section.
    do i=1,parlst_querysubstrings (p_rsection, "bdExpressions")
    
      call parlst_getvalue_string (p_rsection, "bdExpressions", sstr, "", i)
      
      ! Every expression has an associated ivalue and dvalue.
      ivalue = 0
      dvalue = 0.0_DP
      
      ! Get the type and decide on the identifier how to save the expression.
      read(sstr,*) sname,ityp
      
      select case (ityp)
      case (BDC_USERDEFID)
        ! Id of a hardcoded expression realised in the callback routines
        read(sstr,*) sname,ityp,ivalue
        
      case (BDC_USERDEF)
        ! Name of a hardcoded expression realised in the callback routines
        read(sstr,*) sname,ityp,sexpr
        
        ! The expression is saved to an associated "_sval" entry.
        call collct_setvalue_string (roptcBDC%rcollExprTypes, trim(sname)//"_sval",&
            sexpr, .true.)

      case (BDC_EXPRESSION)
        ! General expression; not implemented yet
        read(sstr,*) sname,ityp,sexpr
        
        ! Compile the expression; the expression gets number i
        call fparser_parseFunction (roptcBDC%rparser, i, sexpr, SEC_EXPRVARIABLES)
        
        ! i is the id of the expression, saved in the associated ivalue.
        ivalue = i
        
      case (BDC_VALDOUBLE)
        ! Real-value
        read(sstr,*) sname,ityp,dvalue
        
      case (BDC_VALINT)
        ! Integer-value
        read(sstr,*) sname,ityp,ivalue
                
      case (BDC_VALPARPROFILE)
        ! Parabolic profile with specified maximum velocity
        read(sstr,*) sname,ityp,dvalue
                                    
      case default
        call output_line ("Expressions not implemented!", &
            OU_CLASS_ERROR,OU_MODE_STD,"struc_initBDC")
        call sys_halt()

      end select
      
      ! Put the type of the expression to the main collection section
      call collct_setvalue_int (roptcBDC%rcollExprTypes, sname, ityp, .true.)
      
      ! Every expression has an associated ivalue and dvalue.
      call collct_setvalue_int (roptcBDC%rcollExprTypes, trim(sname)//"_ival", ivalue, .true.)
      call collct_setvalue_real (roptcBDC%rcollExprTypes, trim(sname)//"_dval", dvalue, .true.)
      
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine struc_getBdExprInfo (roptcBDC,sexpression,&
      ctype,iid,ivalue,dvalue,svalue)
  
!<description>
  ! Returns the type and information of a boundary expression.
!</description>

!<input>
  ! Boundary condition structure
  type(t_optcBDC), intent(inout) :: roptcBDC
  
  ! Name of the expression
  character(len=*), intent(in) :: sexpression
!</input>

!<output>
  ! Returns the type of the expression.
  integer, intent(out) :: ctype
  
  ! The associated id or =0, if no id is associated.
  integer, intent(out) :: iid
  
  ! The associated ivalue.
  integer, intent(out) :: ivalue
  
  ! The associated dvalue
  real(DP), intent(out) :: dvalue
  
  ! The associated svalue or "" if no svalue is associated.
  character(len=*), intent(out) :: svalue
!</output>

!</subroutine>

    logical :: bexists

    ! Get the type from the collection
    ctype = collct_getvalue_int (roptcBDC%rcollExprTypes, sexpression)
    
    ! Get the associated values.
    ivalue = 0
    ivalue = collct_getvalue_int (roptcBDC%rcollExprTypes, &
        trim(sexpression)//"_ival", bexists=bexists)
        
    dvalue = 0.0_DP
    dvalue = collct_getvalue_real (roptcBDC%rcollExprTypes, &
        trim(sexpression)//"_dval", bexists=bexists)
        
    svalue = ""
    call collct_getvalue_string (roptcBDC%rcollExprTypes, &
        trim(sexpression)//"_sval", svalue, bexists=bexists)
        
    ! If this is an expression, ivalue is the id.
    if (ctype .eq. BDC_EXPRESSION) then
      iid = ivalue
      ivalue = 0
    end if

  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine struc_getBDCinfo (roptcBDC,ssection,ibct,iindex,nsegments,bautomatic)

!<description>
  ! Returns information about the boundary condition definition on a boundary
  ! component.
!</description>
  
!<input>
  ! Boundary condition structure
  type(t_optcBDC), intent(inout) :: roptcBDC
  
  ! Name of the section with the boundary conditions.
  character(len=*), intent(in) :: ssection
  
  ! Number of the boundary component
  integer, intent(in) :: ibct
!</input>

!<output>
  ! Index of the boundary component. =0 if there is no definition for this 
  ! boundary component.
  integer, intent(out) :: iindex

  ! Number of segments that define the BCs on this bondary component.
  ! =0, if if there is no definition for this boundary component.
  integer, intent(out) :: nsegments
  
  ! Returns TRUE if the boundary component has 'automatic'
  ! bundary values. This is defined by a "-1" as the value
  ! of the boundary component identifier.
  logical, intent(out) :: bautomatic
!</output>

!</subroutine>

    ! local variables
    character(LEN=PARLST_MLDATA) :: cstr,cexpr
    type(t_parlstSection), pointer :: p_rsection
    integer :: i

    nsegments = 0
    bautomatic = .false.

    ! Parse the parameter "bdComponentX"
    write (cexpr,"(I10)") ibct
    cstr = "bdComponent" // adjustl(cexpr)
    
    ! Get the index of the parameter (in the parameter list)
    call parlst_querysection(roptcBDC%p_rparList, ssection, p_rsection)
    iindex = parlst_queryvalue (p_rsection, cstr)
    if (iindex .ne. 0) then
      ! Number of segments that define the BCs there.
      nsegments = parlst_querysubstrings (p_rsection, cstr)
      
      if (nsegments .eq. 0) then
        ! Ask the value. If it is =1, the entry if of type "automatic".
        call parlst_getvalue_int (p_rsection,iindex,i)
        bautomatic = i .eq. -1
      end if
    end if

  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine struc_getSegmentInfo (&
      roptcBDC,ssection,iindex,isegment,dpar,iintervalsEnd,cbcType,sparams)

!<description>
  ! Returns information about a boundary condition segment on a boundary
  ! component.
!</description>
  
!<input>
  ! Boundary condition structure
  type(t_optcBDC), intent(inout) :: roptcBDC
  
  ! Name of the section with the boundary conditions.
  character(len=*), intent(in) :: ssection

  ! Index of the boundary component, returned by cc_getBDCinfo.
  integer, intent(in) :: iindex

  ! Number of the segment where to return information for.
  integer, intent(in) :: isegment
!</input>

!<output>
  ! End parameter value of the segment
  real(DP), intent(out) :: dpar
  
  ! Definition of the endpoint inclusion.
  ! =0: NO endpoints included.
  ! =1: start point included
  ! =2: Ending point included
  ! =3: Start and ending point included into the segment.
  integer, intent(out) :: iintervalsEnd
  
  ! Type of the boundary conditions here.
  integer, intent(out) :: cbcType
  
  ! Additional parameters, boundary condition type dependent.
  ! Must be evaluated by the caller.
  character(LEN=*), intent(out) :: sparams
!</output>

!</subroutine>

    ! local variables
    character(LEN=PARLST_MLDATA) :: sstr,sstr2
    integer :: ntokens,istart
    type(t_parlstSection), pointer :: p_rsection

    ! Get the section
    call parlst_querysection(roptcBDC%p_rparlist, ssection, p_rsection)

    ! Get information about that segment.
    call parlst_getvalue_string (p_rsection, iindex, sstr, isubstring=isegment)
    
    call sys_counttokens (sstr,ntokens)
    if (ntokens .lt. 3) then
      call output_line ("Invalid definition of boundary conditions!", &
          OU_CLASS_ERROR,OU_MODE_STD,"sbc_getSegmentInfo")
      call output_line (""""//trim(sstr)//"""", &
          OU_CLASS_ERROR,OU_MODE_STD,"sbc_getSegmentInfo")
      call sys_halt()
    end if
    
    ! Read the segment parameters
    read(sstr,*) dpar,iintervalsEnd,cbcType

    ! Get the position of the parameters
    istart = 1
    call sys_getNextToken (sstr,sstr2,istart)
    call sys_getNextToken (sstr,sstr2,istart)
    call sys_getNextToken (sstr,sstr2,istart)
    
    ! Get the remaining parameters (if there are any)
    if (istart .ne. 0) then
      sparams = trim(sstr(istart:))
    else
      sparams = ""
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine struc_doneBDC (roptcBDC)
  
!<description>
  ! Cleans up information in the optimal control structure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up
  type(t_optcBDC), intent(inout) :: roptcBDC
!</inputoutput>

!</subroutine>

    ! Release the parser object.
    call fparser_release (roptcBDC%rparser)
    
    ! Release the collection
    call collct_done (roptcBDC%rcollExprTypes)
    
    ! Clean up the rest
    nullify(roptcBDC%p_rphysics)
    nullify(roptcBDC%p_rparList)
    roptcBDC%ssectionBdExpr     = ""
    roptcBDC%ssectionBdCondPrim = ""
    roptcBDC%ssectionBdCondDual = ""
    roptcBDC%ssectionBdCondPrimLin = ""
    roptcBDC%ssectionBdCondDualLin = ""

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_addBoundaryRegion (rboundaryRegion,rregionList)
  
!<description>
  ! Adds boundary region to a list of boundary regions.
!</desctiption>

!<input>
  ! Region to attach.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!</input>

!<inputoutput>
  ! List of boundary regions where to attach the region
  type(t_boundaryRegionList), intent(inout), target :: rregionList
!</inputoutput>

!</subroutine>

    if (rregionList%nregions .eq. 0) then
      ! Create the structure if it does not exist.
      allocate (rregionList%p_rbdTail)
      rregionList%p_rbdHead => rregionList%p_rbdTail
    else
      ! Attach a new tail.
      allocate (rregionList%p_rbdTail%p_nextBdRegion)
      rregionList%p_rbdTail => rregionList%p_rbdTail%p_nextBdRegion
    end if
    
    ! Put the boundary region to there.
    rregionList%p_rbdTail%rboundaryRegion = rboundaryRegion
    
    rregionList%nregions = rregionList%nregions + 1

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_releaseBdRegionList (rregionList)
  
!<description>
  ! Releases a rregionList structure that was allocated in
  ! sbc_assembleBDconditions.
!</desctiption>

!<inputoutput>
  ! Structure to release
  type(t_boundaryRegionList), intent(inout), target :: rregionList
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_bdRegionEntry), pointer :: p_rregion1,p_rregion2
    integer :: i
    
    ! Release the entries
    p_rregion1 => rregionList%p_rbdHead
    do i=1,rregionList%nregions
      p_rregion2 => p_rregion1
      p_rregion1 => p_rregion2%p_nextBDregion
      deallocate(p_rregion2)
    end do
    
    ! Nothing inside anymore.
    rregionList%nregions = 0

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine sbc_getDofsInBdRegionList (rspatialDiscr,rregionList,ndofs,IdofsArray)
  
!<description>
  ! Collects all DOFs in a boundary region list
!</desctiption>

!<input>
  ! The discretisation structure of the underlying discretisation.
  type(t_spatialDiscretisation), intent(in) :: rspatialDiscr

  ! List of boundary regions
  type(t_boundaryRegionList), intent(in) :: rregionList
!</input>

!<output>
  ! Number of DOFs
  integer, intent(out) :: ndofs
  
  ! OPTIONAL: Array where the DOFs are saved to. The array must be large enough.
  integer, dimension(:), intent(out), optional :: IdofsArray
!</output>

!</subroutine>

    ! local variables
    integer :: ndofsLocal
    type(t_bdRegionEntry), pointer :: p_rbdRegion

    if (rregionList%nregions .eq. 0) return
    
    if (.not. present(IdofsArray)) then
      ! Loop through all boundary regions and determine the number of DOFs.
      p_rbdRegion => rregionList%p_rbdHead
      
      ndofs = 0
      do while (associated(p_rbdRegion))
      
        call bcasm_getDOFsInBDRegion (rspatialDiscr, &
            p_rbdRegion%rboundaryRegion, ndofs=ndofsLocal)
            
        ndofs = ndofs + ndofsLocal
      
      end do
    else
      ! Loop through all boundary regions and save the DOFs.
      p_rbdRegion => rregionList%p_rbdHead
      
      ndofs = 0
      do while (associated(p_rbdRegion))
      
        call bcasm_getDOFsInBDRegion (rspatialDiscr, &
            p_rbdRegion%rboundaryRegion, ndofs=ndofsLocal,&
            IdofsArray=IdofsArray(ndofs+1:))
            
        ndofs = ndofs + ndofsLocal
      
      end do
    end if    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbch_getDiscreteBC (roptcBDCSpaceHierarchy,ilevel,p_rdiscreteBC)
  
!<description>
  ! Obtains a pointer to the discrete boundary condition structure
  ! on level ilevel.
!</desctiption>

!<input>
  ! Boundary condition hierarchy
  type(t_optcBDCSpaceHierarchy), intent(in), target :: roptcBDCSpaceHierarchy
  
  ! Level
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Pointer to the discrete boundary conditions
  type(t_discreteBC), pointer :: p_rdiscreteBC
!</inputoutput>

!</subroutine>

    p_rdiscreteBC => &
        roptcBDCSpaceHierarchy%p_RoptcBDCspace(ilevel)%rdiscretebc

  end subroutine
  
end module
