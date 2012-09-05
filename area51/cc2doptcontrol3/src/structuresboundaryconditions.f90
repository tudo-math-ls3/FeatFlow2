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
  character(LEN=10), dimension(7), parameter, public :: SEC_EXPRVARIABLES = &
    (/"X    ","Y    ","Z    ","L    ","R    ","S    ","TIME "/)

!</constantblock>

!<constantblock description="Assembly mode for boundary conditions">
  ! Assembles strong Dirichlet boundary conditions
  integer, parameter, public :: SBC_DIRICHLETBC = 1

  ! Assembles Neumann boundary conditions.
  integer, parameter, public :: SBC_NEUMANN = 2

  ! Assembles Dirichlet boudary control boundary conditions.
  integer, parameter, public :: SBC_DIRICHLETBCC = 4
  
  ! Assembles all boundary conditions.
  integer, parameter, public :: SBC_ALL = SBC_DIRICHLETBC + SBC_NEUMANN + SBC_DIRICHLETBCC
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

  end type

!</typeblock>

  public :: t_optcBDC

!<typeblock>

  ! Pointer to a t_optcBDC structure.
  type t_p_optcBDC
  
    type(t_optcBDC), pointer :: p_roptcBDC

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

  ! Encapsules all Neumann boundary parts in the domain
  type t_boundaryRegionList
  
    ! Number of Neumann boundary regions in the primal equation
    integer :: nregions = 0

    ! Pointer to the head of a list of Neumann boundary regions, primal equation
    type(t_bdRegionEntry), pointer :: p_rbdHead => null()

    ! Pointer to the tail of a list of Neumann boundary regions, primal equation
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
  
    ! Boundary regions with Dirichlet control boundary
    type(t_boundaryRegionList) :: rdirichletControlBoundary
    
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

  ! Evaluates an expression on the boundary.
  public :: struc_evalExpression

  ! Cleans up information in the optimal control structure.
  public :: struc_doneBDC

  ! Adds boundary region to a list of boundary regions.
  public :: sbc_addBoundaryRegion
  
  ! Releases a rregionList structure that was allocated in
  ! sbc_assembleBDconditions.
  public :: sbc_releaseBoundaryList 

  ! Obtains a pointer to the discrete boundary condition structure
  ! on level ilevel.
  public :: sbch_getDiscreteBC
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine struc_initBDC (roptcBDC,rparlist,rphysics,rsettingsOptControl,&
      ssectionBdExpr,ssectionBdCondPrim,ssectionBdCondDual,&
      ssectionBdCondPrimLin,ssectionBdCondDualLin)
  
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
      ctype,iid,ivalue,dvalue,svalue,bneedsParams)
  
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
  
  ! Defines whether or not a set of parameters is necessary for the
  ! evaluation of the expression.
  ! If this is set to FALSE, the expression can be evaluated without
  ! the need for parameters.
  ! If set to TRUE, the parameters must be specified.
  logical, intent(out) :: bneedsParams
!</output>

!</subroutine>

    logical :: bexists

    ! Get the type from the collection
    ctype = collct_getvalue_int (roptcBDC%rcollExprTypes, sexpression)
    
    ! Most types do not need parameters
    bneedsParams = .false.
    if ((ctype .eq. BDC_EXPRESSION) .or. (ctype .eq. BDC_VALPARPROFILE)) then
      bneedsParams = .true.
    end if
      
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

 ! ***************************************************************************

!<subroutine>

  real(DP) function struc_evalExpression (roptcBDC,cexprType,iexprid,&
      ivalue,dvalue,Dparams)
  
!<description>
  ! Evaluates an expression on the boundary.
  !
  ! This evaluates all expressions except for cexprtype=BDC_USERDEFID and
  ! BDC_USERDEF.
!</description>

!<input>
  ! Boundary condition structure
  type(t_optcBDC), intent(in) :: roptcBDC
  
  ! Type of the expression
  integer, intent(in) :: cexprType
  
  ! id of the expression
  integer, intent(in) :: iexprid
  
  ! Integer information value. Standard=0.
  ! For cexpr=BDC_EXPRESSION: Id of the expression.
  integer, intent(in) :: ivalue

  ! Double precision information value. Standard=0.0.
  ! For cexpr=BDC_VALDOUBLE/BDC_VALPARPROFILE: maximum velocity.
  real(DP), intent(in) :: dvalue
  
  ! An array with parameters corresponding to the possible expression
  ! variables defined by SEC_EXPRVARIABLES.
  ! This array only needs to be passed if bneedsparams was set to TRUE
  ! in the call to struc_getBdExprType.
  real(DP), dimension(:), intent(in), optional :: Dparams
!</input>

!</subroutine>

    select case (cexprType)
    case (BDC_USERDEFID,BDC_USERDEF)
      call output_line ("This expression must be evaluated manually!", &
          OU_CLASS_ERROR,OU_MODE_STD,"struc_evalExpression")
      call sys_halt()
    
    case (BDC_VALDOUBLE)
      ! A simple constant
      struc_evalExpression = dvalue

    case (BDC_EXPRESSION)
      ! A complex expression.

      ! Evaluate the expression. ivalue is the number of
      ! the expression to evaluate.
      if (.not. present(Dparams)) then
        call output_line ("Dparams not available!", &
            OU_CLASS_ERROR,OU_MODE_STD,"struc_evalExpression")
        call sys_halt()
      end if
      
      call fparser_evalFunction (roptcBDC%rparser, iexprid, Dparams, struc_evalExpression)
      
    case (BDC_VALPARPROFILE)
      ! A parabolic profile. dvalue expresses the
      ! maximum value of the profile.
      !
      ! Get the maximum of the profile.
      ! The evaluation point is L=Dparams(4) in [0,1], the
      ! local parameter value.
      struc_evalExpression = mprim_getParabolicProfile (Dparams(4),1.0_DP,dvalue)
      
    case default
      call output_line ("Invalid boundary condition type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"struc_evalExpression")
      call sys_halt()

    end select

  end function

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

  subroutine sbc_releaseBoundaryList (rregionList)
  
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
