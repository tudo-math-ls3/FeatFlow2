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
!#     -> Parses the definition of the BC`s from sections given by DAT files
!#        and sets up an discrete boundary condition description
!#
!# 2.) cc_assembleFBDconditions
!#     -> Assembles fictitious boundary boundary conditions.
!#
!# 3.) cc_assembleInhomNeumann
!#     -> Assembles inhomogeneous Neumann BCs into a RHS vector.
!#
!# Internal subroutines:
!#
!# 1.) cc_initAnalyticBC / cc_doneAnalyticBC
!#     -> initialise/clenaup a structure encapsuling the analytic BCs
!#
!# 2.) cc_getExprType
!#     -> Determine the type of an expression
!#
!# 3.) cc_getIvalue / cc_getDvalue / cc_getSvalue
!#     -> Determine information about an expression
!#
!# 4.) cc_getBDCinfo
!#     -> Obtain information about the boundary conditions on a boundary
!#        component
!#
!# 5.) cc_getSegmentInfo
!#     -> Determine the BCs on a boundary segment
!#
!# 6.) cc_initBCassembly / cc_doneBCassembly
!#     -> Initialise/Cleanup the BC assembly
!#
!# 7.) cc_getBCassemblyPointer
!#     -> For callback routines: Retrieve a BC assembly structure
!#        from a collection
!#
!# 8.) cc_getDirichletEdges
!#     -> Determine all edges with Dirichlet BCs
!#
!# 9.) cc_prepareExpressionEval
!#     -> Prepare the evaluation of an expression
!#
!# 10.) cc_evalBoundaryValue
!#      -> Evaluate an expression in a point on the boundary.
!#
!# 11.) cc_fcoeff_bdConditions
!#      -> Callback routine for Dirichlet BCs
!#
!# 12.) cc2d_fcoeff_inhomNeumann
!#      -> Callback routine for inhomogeneous Neumann BCs
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
  use scalarpde
  use derivatives
  use feevaluation
  
  use collection
  use convection
    
  use ccbasic
  use cccallback
  use ccmatvecassembly
  
  implicit none
  
  private

!<constants>

!<constantblock description="Names of the sections used in the main collection">

  ! Section name of the section saving the boundary expressions, i.e. the
  ! expressions that are to be evaluated in each point on the boundary.
  character(LEN=COLLCT_MLSECTION), parameter :: SEC_SBDEXPRESSIONS = "BDEXPRESSIONS"

  ! Name of the parser object for boundary value expressions
  character(LEN=COLLCT_MLSECTION), parameter :: BDC_BDPARSER = "BDEXPRPARSER"
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

!</constants>

!<types>

!<typeblock>

  ! COntains an analytic description of the boundary conditions
  ! based on expressions, including a parser object for parsing
  ! the expressions.
  type t_analyticBC
  
    ! A compiled expression for evaluation at runtime
    type(t_fparser) :: rparser
    
    ! A local collection we use for storing named parameters during the
    ! parsing process.
    type(t_collection) :: rbcCollection

    ! Section in the parameter list specifying the segment information.
    type(t_parlstSection), pointer :: p_rbdcond => null()

  end type

!</typeblock>

!<typeblock>

  ! A structure which is passed to the callback routine for assembly.
  type t_bcassemblyData
    
    ! Pointer to the problem data
    type(t_problem), pointer :: p_rproblem => null()
    
    ! Definition of the underlying boundary
    type(t_boundary), pointer :: p_rboundary => null()
    
    ! Discretisation structure currently under investigation.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! Current boundary region in process.
    type(t_boundaryRegion) :: rboundaryRegion
    
    ! The definition of the analytic boundary conditions.
    type(t_analyticBC) :: ranalyticBC
    
    ! Pointer to a user defined collection structure
    type(t_collection), pointer :: p_ruserCollection => null()
    
    ! Current time where to set up the BCs
    real(DP) :: dtime = 0.0_DP
    
    ! Minimum time
    real(DP) :: dtimeMin = 0.0_DP

    ! Maximum time
    real(DP) :: dtimeMax = 0.0_DP
    
    ! Weighting factor for the values on the boundary.
    real(DP) :: dweight = 1.0_DP
    
    !<!-- Boundary condition specific parameters -->
    
    ! Type of boundary condition
    integer :: cbcType = 0
    
    ! component under consideration (1=x-vel, 2=y-vel,...)
    integer :: icomponent = 0

    ! Whether or not the moving frame formulation is active.
    integer :: imovingFrame = 0
    
    !  X/Y-velocity of the moving frame
    real(DP), dimension(NDIM2D) :: DmframeVel = 0.0_DP
    
    !  X/Y-acceleration of the moving frame
    real(DP), dimension(NDIM2D) :: DmframeAcc = 0.0_DP
    
    !<!-- Parameters vor the evaluation of expressions -->
    
    ! expression type
    integer :: cexprType = 0
    
    ! Integer tag for the expression
    integer :: ivalue = 0
    
    ! Real value for the expression
    real(DP) :: dvalue = 0.0_DP
    
    ! Name of the expression
    character(len=PARLST_MLNAME) :: sexprName = ""
    
  end type

!</typeblock>

!<typeblock>
  ! Structure encapsuling a pointer to t_bcassemblyData, for being passed
  ! to callback routines.
  type t_p_bcassemblyData
    type(t_bcassemblyData), pointer :: p_rbcAssemblyData => null()
  end type
!</typeblock>

!</types>

  public :: cc_assembleBDconditions
  public :: cc_assembleFBDconditions
  public :: cc_assembleInhomNeumann

contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initAnalyticBC (rparlist,ranalyticBC)

!<description>
  ! Initialises a boundary condition structure for the evaluation
  ! of boundary conditions. Note: This is level independent
  ! and based on parser expressions.
!</description>
  
!<inputoutput>
  ! Parmeter list with parameters concerning boundary conditions.
  type(t_parlist), intent(inout), target :: rparlist
!</inputoutput>

!<output>
  ! Boundary condition object receiving an analytic definition
  ! of the boundary conditions.
  type(t_analyticBC), intent(out) :: ranalyticBC
!</output>

!</subroutine>

    ! local variables
    integer :: i,ityp,ivalue
    real(DP) :: dvalue
    character(LEN=PARLST_LENLINEBUF) :: cstr,cexpr
    character(LEN=PARLST_MLNAME) :: cname
    
    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection
    
    ! We need the analytic description of the boundary conditions.
    ! Initialise a structure for boundary conditions, which accepts this,
    ! on the heap.
    !
    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
    !
    ! Get the expression/bc sections from the boundary condition block
    call parlst_querysection(rparlist, "BDEXPRESSIONS", p_rsection)
    call parlst_querysection(rparlist, "BDCONDITIONS", ranalyticBC%p_rbdcond)
    
    ! For intermediate storing of expression types, we use a local collection
    call collct_init (ranalyticBC%rbcCollection)
    
    ! Add a section to the collection that accepts the boundary expressions
    call collct_addsection (ranalyticBC%rbcCollection, SEC_SBDEXPRESSIONS)
    
    ! Create a parser structure for as many expressions as configured
    call fparser_create (ranalyticBC%rparser,&
         parlst_querysubstrings (p_rsection, "bdExpressions"))
    
    ! Add the boundary expressions to the collection into the
    ! specified section.
    do i=1,parlst_querysubstrings (p_rsection, "bdExpressions")
    
      call parlst_getvalue_string (p_rsection, "bdExpressions", cstr, "", i)
      
      ! Get the type and decide on the identifier how to save the expression.
      read(cstr,*) cname,ityp
      
      select case (ityp)
      case (BDC_USERDEF)
        ! Name of a hardcoded expression realised in the callback routines
        read(cstr,*) cname,ityp,cexpr
        call collct_setvalue_string (ranalyticBC%rbcCollection, &
            cname, cexpr, .true., 0, SEC_SBDEXPRESSIONS)

      case (BDC_EXPRESSION)
        ! General expression; not implemented yet
        read(cstr,*) cname,ityp,cexpr
        
        ! Compile the expression; the expression gets number i
        call fparser_parseFunction (ranalyticBC%rparser, i, cexpr, EXPRVARIABLES)
        
        ! Add the number of the function in the parser object to
        ! the collection with the name of the expression
        call collct_setvalue_int (ranalyticBC%rbcCollection, &
            cname, i, .true., 0, SEC_SBDEXPRESSIONS)
        
      case (BDC_VALDOUBLE)
        ! Real-value
        read(cstr,*) cname,ityp,dvalue
        call collct_setvalue_real (ranalyticBC%rbcCollection, &
            cname, dvalue, .true., 0, SEC_SBDEXPRESSIONS)
        call fparser_defineConstant(cname, dvalue)

      case (BDC_VALINT)
        ! Integer-value
        read(cstr,*) cname,ityp,ivalue
        call collct_setvalue_int (ranalyticBC%rbcCollection, &
            cname, ivalue, .true., 0, SEC_SBDEXPRESSIONS)
        call fparser_defineConstant(cname, real(ivalue,DP))
                
      case (BDC_VALPARPROFILE)
        ! Parabolic profile with specified maximum velocity
        read(cstr,*) cname,ityp,dvalue
        call collct_setvalue_real (ranalyticBC%rbcCollection, &
            cname, dvalue, .true., 0, SEC_SBDEXPRESSIONS)
                                   
      case default
        call output_line ("Expressions not implemented!", &
            OU_CLASS_ERROR,OU_MODE_STD,"cc_initBDconditions")
        call sys_halt()

      end select
      
      ! Put the type of the expression to the temporary collection section
      call collct_setvalue_int (ranalyticBC%rbcCollection, cname, ityp, .true.)
      
    end do
    
  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine cc_getExprType (ranalyticBC,sexpression,iexprType)

!<description>
  ! Returns the type of an expression.
!</description>
  
!<input>
  ! Boundary condition object.
  type(t_analyticBC), intent(inout) :: ranalyticBC
  
  ! Name of the expression
  character(len=*), intent(in) :: sexpression
!</input>

!<output>
  ! Type of the expression
  integer, intent(out) :: iexprType
!</output>

!</subroutine>

    iexprType = collct_getvalue_int (ranalyticBC%rbcCollection, sexpression)

  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine cc_getIvalue (ranalyticBC,sexpression,ivalue)

!<description>
  ! Returns the ivalue for an expression.
!</description>
  
!<input>
  ! Boundary condition object.
  type(t_analyticBC), intent(inout) :: ranalyticBC
  
  ! Name of the expression
  character(len=*), intent(in) :: sexpression
!</input>

!<output>
  ! Integer value
  integer, intent(out) :: ivalue
!</output>

!</subroutine>

    ivalue = collct_getvalue_int (&
        ranalyticBC%rbcCollection, sexpression, 0, SEC_SBDEXPRESSIONS)

  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine cc_getDvalue (ranalyticBC,sexpression,dvalue)

!<description>
  ! Returns the dvalue for an expression.
!</description>
  
!<input>
  ! Boundary condition object.
  type(t_analyticBC), intent(inout) :: ranalyticBC
  
  ! Name of the expression
  character(len=*), intent(in) :: sexpression
!</input>

!<output>
  ! Real value
  real(DP), intent(out) :: dvalue
!</output>

!</subroutine>

    dvalue = collct_getvalue_real (ranalyticBC%rbcCollection, &
        sexpression, 0, SEC_SBDEXPRESSIONS)

  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine cc_getSvalue (ranalyticBC,sexpression,svalue)

!<description>
  ! Returns the svalue for an expression.
!</description>
  
!<input>
  ! Boundary condition object.
  type(t_analyticBC), intent(inout) :: ranalyticBC
  
  ! Name of the expression
  character(len=*), intent(in) :: sexpression
!</input>

!<output>
  ! Strinbg value
  character(len=*), intent(out) :: svalue
!</output>

!</subroutine>

    call collct_getvalue_string (ranalyticBC%rbcCollection, &
        sexpression, svalue, 0, SEC_SBDEXPRESSIONS)

  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine cc_getBDCinfo (ranalyticBC,ibct,iindex,nsegments)

!<description>
  ! Returns information about the boundary condition definition on a boundary
  ! component.
!</description>
  
!<input>
  ! Boundary condition object.
  type(t_analyticBC), intent(in) :: ranalyticBC
  
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
!</output>

!</subroutine>

    ! local variables
    character(LEN=PARLST_MLDATA) :: cstr,cexpr

    nsegments = 0

    ! Parse the parameter "bdComponentX"
    write (cexpr,"(I10)") ibct
    cstr = "bdComponent" // adjustl(cexpr)
    
    ! Get the index of the parameter (in the parameter list)
    iindex = parlst_queryvalue (ranalyticBC%p_rbdcond, cstr)
    if (iindex .ne. 0) then
      ! Number of segments that define the BCs there.
      nsegments = parlst_querysubstrings (ranalyticBC%p_rbdcond, cstr)
    end if

  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine cc_getSegmentInfo (&
      ranalyticBC,iindex,isegment,dpar,iintervalsEnd,cbcType,sparams)

!<description>
  ! Returns information about a boundary condition segment on a boundary
  ! component.
!</description>
  
!<input>
  ! Boundary condition object.
  type(t_analyticBC), intent(in) :: ranalyticBC
  
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

    ! Get information about that segment.
    call parlst_getvalue_string (ranalyticBC%p_rbdcond, iindex, sstr, isubstring=isegment)
    
    call sys_counttokens (sstr,ntokens)
    if (ntokens .lt. 3) then
      call output_line ("Invalid definition of boundary conditions!", &
          OU_CLASS_ERROR,OU_MODE_STD,"cc_getSegmentInfo")
      call output_line (""""//trim(sstr)//"""", &
          OU_CLASS_ERROR,OU_MODE_STD,"cc_getSegmentInfo")
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

  subroutine cc_doneAnalyticBC (ranalyticBC)

!<description>
  ! Cleans up an analytic BC structure.
!</description>
  
!<inputoutput>
  ! Structure to be cleaned up
  type(t_analyticBC), intent(inout) :: ranalyticBC
!</inputoutput>

!</subroutine>

    ! Release the content
    call fparser_release (ranalyticBC%rparser)
    call collct_done (ranalyticBC%rbcCollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initBCassembly (rproblem,rdiscretisation,ruserCollection,&
      rbcAssemblyData,rcollection)

!<description>
  ! Prepares the assembly of boundary conditions.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! User defined collection structure to be passed to callback routines
  type(t_collection), intent(in), target :: ruserCollection
  
  ! Discretisation structure currently under investigation
  type(t_blockDiscretisation), target :: rdiscretisation

!</input>
  
!<output>
  ! Boundary condition assembly data, to be set up
  type(t_bcassemblyData), intent(out), target :: rbcAssemblyData
!</output>

!<inputoutput>
  ! Collection structure where a pointer to rbcAssemblyData is saved to.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_p_bcassemblyData) :: rp_bcassemblyData
    
    ! Create rbcAssemblyData based on the given data.
    !
    ! At first, we need the analytic description of the boundary conditions.
    call cc_initAnalyticBC (rproblem%rparamlist,rbcAssemblyData%ranalyticBC)
    
    rbcAssemblyData%p_rproblem => rproblem
    rbcAssemblyData%p_rboundary => rproblem%rboundary
    rbcAssemblyData%p_ruserCollection => ruserCollection
    rbcAssemblyData%p_rdiscretisation => rdiscretisation
    rbcAssemblyData%dtime = rproblem%rtimedependence%dtime
    rbcAssemblyData%dtimeMin = rproblem%rtimedependence%dtimeInit
    rbcAssemblyData%dtimeMax = rproblem%rtimedependence%dtimeMax

    ! Create a pointer to this, encapsuled in a structure.
    rp_bcassemblyData%p_rbcAssemblyData => rbcAssemblyData
    
    ! Transfer this pointer into the IquickAccess array of the collection.
    ! We can get it from there later.
    rcollection%IquickAccess(:) = &
        transfer(rp_bcassemblyData,rcollection%IquickAccess(:),size(rcollection%IquickAccess(:)))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getBCassemblyPointer (rcollection,p_rbcAssemblyData)

!<description>
  ! Returns a pointer to the BC assembly data. The pointer is obtained
  ! from the given collection
!</description>

!<input>
  ! Collection structure. Must have been prepared with cc_initBCassembly.
  type(t_collection), intent(in) :: rcollection
!</input>

!<inputoutput>
  ! Pointer to be obtained.
  type(t_bcassemblyData), pointer :: p_rbcAssemblyData
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_p_bcassemblyData) :: rp_bcassemblyData

    ! Get back the pointer
    rp_bcassemblyData = transfer(rcollection%IquickAccess(:),rp_bcassemblyData)
    
    ! Return the pointer
    p_rbcAssemblyData => rp_bcassemblyData%p_rbcAssemblyData
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneBCassembly (rbcAssemblyData)

!<description>
  ! Cleans up the boundary condition assembly structures.
!</description>

!<inputoutput>
  ! Boundary condition assembly data, to be cleaned up
  type(t_bcassemblyData), intent(inout) :: rbcAssemblyData
!</inputoutput>

!</subroutine>

    ! Not much to do
    call cc_doneAnalyticBC (rbcAssemblyData%ranalyticBC)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleBDconditions (rproblem,rdiscretisation,rdynamicLevelInfo,&
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
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
  
  ! OPTIONAL: If this flag ist set to TRUE, the boundary conditions are
  ! assembled for postprocessing of a solution vector. When being set to FALSE
  ! or not present, the boundary condition are assembled for computation.
  logical, intent(in), optional :: bforPostprocessing
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! Collection structure to be passed to callback routines
  type(t_collection), intent(inout), target :: rcollection
  
  ! A t_dynamicLevelInfo structure that receives a discretised version
  ! of the boundary boundary conditions. The BC substructure should
  ! be empty; new BC`s are simply added to the structure.
  type(t_dynamicLevelInfo), intent(inout) :: rdynamicLevelInfo
!</inputoutput>

!</subroutine>

    ! local variables
    logical :: bNeumann
    integer :: i
    character(LEN=PARLST_MLDATA) :: sparams,sbdex1,sbdex2
    type(t_bcassemblyData) :: rbcAssemblyData
    type(t_collection) :: rlocalCollection
    real(DP) :: dpar1,dpar2
    integer, dimension(2) :: IminIndex,ImaxIndex
    integer :: iindex,nsegments,ibdComponent,isegment,iintervalEnds
    integer :: cbctype,icount
    integer, dimension(NDIM2D) :: IvelEqns
    integer(I32) :: casmComplexity
    
    ! Triangulation on currently highest level.
    type(t_triangulation), pointer :: p_rtriangulation

    ! Determine what to assemble
    casmComplexity = BCASM_DISCFORALL
    if (present(bforPostprocessing)) then
      ! Assemble only for the solution vector
      if (bforPostprocessing) casmComplexity = BCASM_DISCFORSOL
    endif
    
    ! Get the triangulation on the highest level
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For implementing boundary conditions, we use a `filter technique with
    ! discretised boundary conditions`. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, initialise the assembly.
    call cc_initBCassembly (rproblem,rdiscretisation,rcollection,&
        rbcAssemblyData,rlocalCollection)
    
    ! Is the moving-frame formulatino active?
    call parlst_getvalue_int (rproblem%rparamList,"CC-DISCRETISATION",&
        "imovingFrame",rbcAssemblyData%imovingFrame,0)
        
    if (rbcAssemblyData%imovingFrame .ne. 0) then
    
      ! Get the velocity and acceleration from the callback routine.
      ! Pass the user collection for this task.
      call getMovingFrameVelocity (&
          rbcAssemblyData%DmframeVel,rbcAssemblyData%DmframeAcc,rcollection)
          
    else
      rbcAssemblyData%DmframeVel(:) = 0.0_DP
      rbcAssemblyData%DmframeAcc(:) = 0.0_DP
    end if
    
    ! Now to the actual boundary conditions.
    !
    ! Log if there is Neumann boundary    
    bNeumann = .false.
    
    ! Loop through all boundary components we have.
    do ibdComponent = 1,boundary_igetNBoundComp(rbcAssemblyData%p_rboundary)
      
      ! Get information about that component: Number of segments,...
      call cc_getBDCinfo (rbcAssemblyData%ranalyticBC,ibdComponent,iindex,nsegments)
      
      ! Parse the segments.
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      if (iindex .eq. 0) then
        ! There is no parameter configuring the boundary condition on that
        ! component - so we have Neumann boundary there.
        bNeumann = .true.
      end if
      
      ! Parameter exists. Get the values in there.
      do isegment = 1,nsegments
        
        ! Get information about the next segment
        call cc_getSegmentInfo (rbcAssemblyData%ranalyticBC,&
            iindex,isegment,dpar2,iintervalEnds,cbctype,sparams)
            
        rbcAssemblyData%cbctype = cbctype
        
        ! Form a boundary condition segment that covers that boundary part
        if (dpar2 .ge. dpar1) then
          
          rbcAssemblyData%rboundaryRegion%dminParam = dpar1
          rbcAssemblyData%rboundaryRegion%dmaxParam = dpar2
          rbcAssemblyData%rboundaryRegion%iboundCompIdx = ibdComponent
          rbcAssemblyData%rboundaryRegion%dmaxParamBC = &
              boundary_dgetMaxParVal(rbcAssemblyData%p_rboundary, ibdComponent)
          rbcAssemblyData%rboundaryRegion%iproperties = iintervalEnds
          
          ! Now, which type of BC is to be created?
          select case (cbctype)
          
          ! -----------------------------------------------
          ! Homogeneous / Inhomogeneous Neumann BC
          ! -----------------------------------------------
          case (0,4)
            ! Usually there is Neumann boundary in this region, but we can not be
            ! sure. Check if, on the highest level, there is at least one edge
            ! of the triangulation belonging to the boundary. If yes, we
            ! have found Neumann boundary. If no, the segment is just too
            ! small to be considered as Neumann boundary.
            
            call bcasm_getElementsInBdRegion (&
                p_rtriangulation,rbcAssemblyData%rboundaryRegion, icount)
            if (icount .gt. 0) bNeumann = .true.
          
          ! -----------------------------------------------
          ! Dirichlet BC
          ! -----------------------------------------------
          case (1)
            ! Simple Dirichlet boundary
            ! Read the line again, get the expressions for X- and Y-velocity
            read(sparams,*) sbdex1,sbdex2
            
            ! For any string <> "", create the appropriate Dirichlet boundary
            ! condition and add it to the list of boundary conditions.
            !
            ! If the type is a double precision value, set the DquickAccess(4)
            ! to that value so it can quickly be accessed.
            if (sbdex1 .ne. "") then
              ! X-velocity
              !
              ! Component number, expression data
              rbcAssemblyData%icomponent = 1
              call cc_prepareExpressionEval (rbcAssemblyData,sbdex1)
            
              ! Assemble the BC`s.
              call bcasm_newDirichletBConRealBD (&
                  rdiscretisation,rbcAssemblyData%icomponent,&
                  rbcAssemblyData%rboundaryRegion,&
                  rdynamicLevelInfo%rdiscreteBC,&
                  cc_fcoeff_bdConditions,rlocalCollection,casmComplexity)
                  
            end if
            
            if (sbdex2 .ne. "") then
            
              ! Y-velocity
              !
              ! Component number, expression data
              rbcAssemblyData%icomponent = 2
              call cc_prepareExpressionEval (rbcAssemblyData,sbdex2)
            
              ! Assemble the BC`s.
              call bcasm_newDirichletBConRealBD (&
                  rdiscretisation,rbcAssemblyData%icomponent,&
                  rbcAssemblyData%rboundaryRegion,&
                  rdynamicLevelInfo%rdiscreteBC,&
                  cc_fcoeff_bdConditions,rlocalCollection,casmComplexity)

            end if
            
          ! -----------------------------------------------
          ! Pressure drop BC
          ! -----------------------------------------------
          case (2)
          
            ! Pressure drop boundary conditions.
            ! Read the line again to get the actual parameters
            read(sparams,*) sbdex1
            
            if (sbdex1 .ne. "") then
            
              ! PDrop BCs are a type of Neumann BCs since there is no Dirichlet value
              ! prescribed for the velocity. Check if there are edges in the region.
              ! If yes, create the BC and remember that we have Neumann BCs.
              
              call bcasm_getEdgesInBCregion (p_rtriangulation,rbcAssemblyData%p_rboundary,&
                  rbcAssemblyData%rboundaryRegion, IminIndex,ImaxIndex,icount)
              if (icount .gt. 0) then
              
                ! There is Neumann boundary
                bNeumann = .true.
                
                ! Component number, expression data
                rbcAssemblyData%icomponent = 0
                call cc_prepareExpressionEval (rbcAssemblyData,sbdex1)
                
                IvelEqns = (/1,2/)
                call bcasm_newPdropBConRealBd (&
                    rdiscretisation,IvelEqns,rbcAssemblyData%rboundaryRegion,&
                    rdynamicLevelInfo%rdiscreteBC,&
                    cc_fcoeff_bdConditions,rlocalCollection,casmComplexity)
                    
              end if

            end if
            
          ! -----------------------------------------------
          ! Slip BC
          ! -----------------------------------------------
          case (3)
            
            ! Nonlinear slip boundary conditions.
            IvelEqns = (/1,2/)
            call bcasm_newSlipBConRealBd (&
                rdiscretisation,IvelEqns(1:NDIM2D),rbcAssemblyData%rboundaryRegion,&
                rdynamicLevelInfo%rdiscreteBC,casmComplexity)
          
          case default
            call output_line ("Unknown boundary condition!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
            call sys_halt()
          end select
          
          ! Move on to the next parameter value
          dpar1 = dpar2
          
        end if
                                          
      end do
      
    end do

    ! The setting in the DAT file may overwrite our guess about Neumann boundaries.
    call parlst_getvalue_int (rbcAssemblyData%ranalyticBC%p_rbdcond, &
        "ineumannBoundary", i, -1)
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
    rdynamicLevelInfo%bhasNeumannBoundary = bNeumann
    
    ! Determine the edges with Dirichlet boundary conditions.
    ! These may be used for stabilisation or other operators that work
    ! on Dirichlet boundary parts differently than on Neumann boundary
    ! parts.
    if (rdynamicLevelInfo%hedgesDirichletBC .ne. ST_NOHANDLE) then
      call storage_free (rdynamicLevelInfo%hedgesDirichletBC)
    end if
    call cc_getDirichletEdges (rbcAssemblyData,p_rtriangulation,0,&
        rdynamicLevelInfo%hedgesDirichletBC,rdynamicLevelInfo%nedgesDirichletBC)
        
    ! Release the boundary condition assembly structure, finish.
    call cc_doneBCassembly (rbcAssemblyData)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_getDirichletEdges (rbcAssemblyData,rtriangulation,icomponent,hedges,ncount)

!<description>
  ! Calculates all edges in equation icomponent which are marked as Dirichlet
  ! boundary edges.
!</description>
  
!<input>
  ! Boudary condition assembly data
  type(t_bcassemblyData), intent(inout) :: rbcAssemblyData
  
  ! Underlying triangulation.
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! Number of the component where the edges should be computed for.
  ! =0: Any velocity component, =1: X-velocity, =2: Y-velocity
  integer, intent(in) :: icomponent
!</input>

!<output>
  ! Handle that identifies a memory block with all the edge numbers
  ! or =ST_NOHANDLE if no edges were found.
  integer, intent(out) :: hedges
  
  ! Number of edges with Dirichlet boundary attached.
  integer, intent(out) :: ncount
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Iedges, p_IedgesLocal
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer :: iindex, j, icount, hedgeslocal, isegment, nsegments
    integer :: ibdComponent,iintervalEnds,cbctype
    real(DP) :: dpar1,dpar2
    character(LEN=PARLST_MLDATA) :: sparams,sbdex1,sbdex2
    
    ! A set of variables describing the analytic boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    
    call storage_getbase_int2d(rtriangulation%h_IedgesAtElement,p_IedgesAtElement)

    ! Allocate memory for the edges. We have at most as many edges as
    ! there are on the boudary.
    call storage_new ("cc_getDirichletEdges","hedges",rtriangulation%NMBD,&
        ST_INT,hedges,ST_NEWBLOCK_NOINIT)
    ncount = 0
    call storage_getbase_int (hedges, p_Iedges)

    ! Allocate another array for local edge numbers.
    call storage_new ("cc_getDirichletEdges","hedgeslocal",rtriangulation%NMBD,&
        ST_INT,hedgeslocal,ST_NEWBLOCK_NOINIT)
    ncount = 0
    call storage_getbase_int (hedgeslocal, p_Iedgeslocal)
    
    ! Loop through all boundary components we have.
    do ibdComponent = 1,boundary_igetNBoundComp(rbcAssemblyData%p_rboundary)
      
      ! Get information about that component: Number of segments,...
      call cc_getBDCinfo (rbcAssemblyData%ranalyticBC,ibdComponent,iindex,nsegments)
      
      ! Parse the segments.
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      ! Parameter exists. Get the values in there.
      do isegment = 1,nsegments
        
        ! Get information about the next segment
        call cc_getSegmentInfo (&
            rbcAssemblyData%ranalyticBC,iindex,isegment,dpar2,iintervalEnds,cbctype,sparams)

        ! Form a boundary condition segment that covers that boundary part
        if (dpar2 .ge. dpar1) then
          
          ! Now, which type of BC is to be created?
          select case (cbctype)
          
          case (1)
            ! Simple Dirichlet boundary
            ! Read the line again, get the expressions for X- and Y-velocity
            read(sparams,*) sbdex1,sbdex2
            
            ! If the type is a double precision value, set the DquickAccess(4)
            ! to that value so it can quickly be accessed.
            if (((icomponent .eq. 0) .and. ((sbdex1 .ne. "") .or. (sbdex2 .ne. ""))) .or. &
                ((icomponent .eq. 1) .and. (sbdex1 .ne. "")) .or. &
                ((icomponent .eq. 2) .and. (sbdex2 .ne. ""))) then
            
              ! Add the edges.
              !
              ! Get the element numbers + the local edge numbers and compute the
              ! actual edge numbers.
              rboundaryRegion%dminParam = dpar1
              rboundaryRegion%dmaxParam = dpar2
              rboundaryRegion%iboundCompIdx = ibdComponent
              rboundaryRegion%dmaxParamBC = &
                boundary_dgetMaxParVal(rbcAssemblyData%p_rboundary, ibdComponent)
              rboundaryRegion%iproperties = iintervalEnds
              
              call bcasm_getElementsInBCregion (rtriangulation,rboundaryRegion, icount, &
                  IelList=p_Iedges(ncount+1:),IedgeLocal=p_Iedgeslocal)
              do j=1,icount
                p_Iedges(ncount+j) = p_IedgesAtElement(p_Iedgeslocal(j),p_Iedges(ncount+j))
              end do
              ncount = ncount + icount
              
            end if

          end select
          
          ! Move on to the next parameter value
          dpar1 = dpar2
          
        end if
                                          
      end do
    
    end do
    
    ! Release unused memory
    call storage_free (hedgeslocal)
    if (ncount .gt. 0) then
      call storage_realloc ("cc_getDirichletEdges", ncount, hedges, ST_NEWBLOCK_NOINIT, .true.)
    else
      call storage_free (hedges)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_fcoeff_bdConditions (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
      cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! "snapshot" of the (actually analytic) boundary conditions.
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
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(in)                                         :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                                         :: cinfoNeeded
  
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
  integer, intent(in)                                         :: iwhere

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
  real(DP), intent(in)                                        :: dwhere
    
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
  ! framework to ignore the node and treat it as "natural boundary condition"
  ! node.
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    type(t_bcassemblyData), pointer :: p_rbcAssemblyData

    ! Get a pointer to the boundary assembly data from the collection.
    call cc_getBCassemblyPointer (rcollection,p_rbcAssemblyData)

    ! Use boundary conditions from DAT files.
    select case (cinfoNeeded)
    
    case (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
          DISCBC_NEEDINTMEAN,DISCBC_NEEDNORMALSTRESS)
      
      ! Dirichlet boundary conditions
    
      ! -> 1=X-velocity, 2=Y-velocity.
      
      ! Return zero Dirichlet boundary values for all situations by default.
      Dvalues(1) = 0.0_DP

      ! Now, depending on the problem, calculate the return value.
      ! Which boundary condition do we have here?
      !
      ! Get the information from evalBoundary.
      ! Note: The information about the BC`s can be retrieved from the
      ! quick-access arrays in the collection as initialised above.
      select case (p_rbcAssemblyData%cbcType)
      case (1)
        ! Simple Dirichlet BCs.
        call cc_evalBoundaryValue (p_rbcAssemblyData,dwhere,Dvalues(1))
            
        if (p_rbcAssemblyData%imovingFrame .ne. 0) then
          ! Moving frame formulation active. Add the frame velocity to the
          ! Dirichlet boundary conditions.
          Dvalues(1) = Dvalues(1) + p_rbcAssemblyData%DmframeVel(p_rbcAssemblyData%icomponent)
        end if
    
      case (2)
        ! Normal stress / pressure drop. Evaluate the  expression iexprtyp.
        call cc_evalBoundaryValue (p_rbcAssemblyData,dwhere,Dvalues(1))
            
      end select
      
    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleInhomNeumann (rproblem,rcollection,rrhs,dweight)

!<description>
  ! Assembles inhomogeneous Neumann boundary conditions into the RHS vector
  ! of the system.
!</description>
  
!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! Collection structure to be passed to callback routines
  type(t_collection), intent(inout), target :: rcollection
  
  ! Weight for the inhomogeneous Neumann part.
  real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
  ! RHS vector to be modified.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    character(LEN=PARLST_MLDATA) :: sparams,sbdex1,sbdex2
    type(t_bcassemblyData) :: rbcAssemblyData
    type(t_collection) :: rlocalCollection
    real(DP) :: dpar1,dpar2
    integer :: iindex,nsegments,ibdComponent,isegment,iintervalEnds
    integer :: cbctype
    type(t_linearForm) :: rform
    
    ! For implementing boundary conditions, we use a `filter technique with
    ! discretised boundary conditions`. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, initialise the assembly.
    call cc_initBCassembly (rproblem,rrhs%p_rblockDiscr,&
        rcollection,rbcAssemblyData,rlocalCollection)
        
    ! Weight the result by the weighting factor.
    rbcAssemblyData%dweight = dweight
    
    ! Loop through all boundary components we have.
    do ibdComponent = 1,boundary_igetNBoundComp(rbcAssemblyData%p_rboundary)
      
      ! Get information about that component: Number of segments,...
      call cc_getBDCinfo (rbcAssemblyData%ranalyticBC,ibdComponent,iindex,nsegments)
      
      ! Parse the segments.
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      ! Parameter exists. Get the values in there.
      do isegment = 1,nsegments
        
        ! Get information about the next segment
        call cc_getSegmentInfo (rbcAssemblyData%ranalyticBC,&
            iindex,isegment,dpar2,iintervalEnds,cbctype,sparams)
        
        ! Form a boundary condition segment that covers that boundary part
        if (dpar2 .ge. dpar1) then
          
          rbcAssemblyData%rboundaryRegion%dminParam = dpar1
          rbcAssemblyData%rboundaryRegion%dmaxParam = dpar2
          rbcAssemblyData%rboundaryRegion%iboundCompIdx = ibdComponent
          rbcAssemblyData%rboundaryRegion%dmaxParamBC = &
              boundary_dgetMaxParVal(rbcAssemblyData%p_rboundary, ibdComponent)
          rbcAssemblyData%rboundaryRegion%iproperties = iintervalEnds
          
          ! Now, which type of BC is to be created?
          select case (cbctype)

          ! -----------------------------------------------
          ! Inhomogeneous Neumann BC
          ! -----------------------------------------------
          case (4)
            
            ! Invoke the assembly of the inhomogeneous Neumann BCs.
            ! This is just an integration on the boundary.
            
            ! Get the boundary expression
            read (sparams,*) sbdex1,sbdex2
            
            ! Use the 4x4 Gauss formula (hardcoded). Should be
            ! enough for most situations.
            rform%itermCount = 1
            rform%Dcoefficients(1) = 1.0_DP
            rform%Idescriptors(1) = DER_FUNC
            
            ! Assemble
            if (sbdex1 .ne. "") then
              rbcAssemblyData%icomponent = 1
              call cc_prepareExpressionEval (rbcAssemblyData,sbdex1)
              
              call linf_buildVectorScalarBdr2D (rform, CUB_G4_1D, .false., &
                  rrhs%RvectorBlock(1),cc2d_fcoeff_inhomNeumann,&
                  rbcAssemblyData%rboundaryRegion, rlocalCollection)
            end if

            if (sbdex2 .ne. "") then
              rbcAssemblyData%icomponent = 2
              call cc_prepareExpressionEval (rbcAssemblyData,sbdex2)
              
              call linf_buildVectorScalarBdr2D (rform, CUB_G4_1D, .false., &
                  rrhs%RvectorBlock(2),cc2d_fcoeff_inhomNeumann,&
                  rbcAssemblyData%rboundaryRegion, rlocalCollection)
            end if
          
          end select
          
          ! Move on to the next parameter value
          dpar1 = dpar2
          
        end if
                                          
      end do
      
    end do

    ! Release the boundary condition assembly structure, finish.
    call cc_doneBCassembly (rbcAssemblyData)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc2d_fcoeff_inhomNeumann (rdiscretisation, rform, &
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
    ! This subroutine is called during the calculation of boundary
    ! integrals. It calculates the inhomogeneous Neumann boundary conditions
    ! which are implemented into the RHS vector.
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
    integer :: ipt,iel
    real(DP) :: dpar
    type(t_bcassemblyData), pointer :: p_rbcAssemblyData

    ! Get a pointer to the boundary assembly data from the collection.
    call cc_getBCassemblyPointer (rcollection,p_rbcAssemblyData)

    do iel = 1,nelements
      do ipt = 1,npointsPerElement
      
        ! Get the parameter value in 0-1 parametrisation
        dpar = boundary_convertParameter(p_rbcAssemblyData%p_rboundary, ibct, &
            DpointPar(ipt,iel), BDR_PAR_LENGTH, BDR_PAR_01)
        
        ! Calculate the expression in the current point.
        call cc_evalBoundaryValue (p_rbcAssemblyData,dpar,Dcoefficients(1,ipt,iel))
      end do
    end do
            
  end subroutine

  ! **************************************************************************

!<subroutine>

  subroutine cc_prepareExpressionEval (rbcAssemblyData,sexpression)

!<description>
  ! Fetches information about an expression into the rbcAssemblyData structure
  ! such that the structure can be used for cc_evalBoundaryValue.
!</description>
  
!<input>
  ! Name of the expression
  character(len=*), intent(in) :: sexpression
!</input>

!<inputoutput>
  ! Boundary condition assembly data, to be prepared
  type(t_bcassemblyData), intent(inout) :: rbcAssemblyData
!</inputoutput>

!</subroutine>

    ! Expression name and type
    rbcAssemblyData%sexprName = sexpression
    
    call cc_getExprType (rbcAssemblyData%ranalyticBC,&
        rbcAssemblyData%sexprName,rbcAssemblyData%cexprType)
    
    ! Dquickaccess(4) / IquickAccess(4) saves information
    ! about the expression.
    select case (rbcAssemblyData%cexprType)
    case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
      ! Constant or parabolic profile
      call cc_getDvalue (rbcAssemblyData%ranalyticBC,&
          rbcAssemblyData%sexprName,rbcAssemblyData%dvalue)
          
    case (BDC_EXPRESSION)
      ! Expression. Write the identifier for the expression
      ! as ivalue into the boundary condition structure.
      call cc_getIvalue (rbcAssemblyData%ranalyticBC,&
          rbcAssemblyData%sexprName,rbcAssemblyData%ivalue)

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine cc_evalBoundaryValue (rbcAssemblyData,dpar,dresult)
  
!<description>
  ! Evaluates an expression in a point on the boundary.
!</description>
  
!<input>
  ! Boudary condition assembly data
  type(t_bcassemblyData), intent(inout) :: rbcAssemblyData

  ! Current parameter value of the point on the boundary.
  ! 0-1-parametrisation.
  real(DP), intent(in) :: dpar
!</input>

!<output>
  ! Value of the expression in the point.
  real(DP), intent(out) :: dresult
!</output>

!</function>

    ! local variables
    real(DP) :: d,dx,dy
    character(LEN=PARLST_MLDATA) :: stag
    real(DP), dimension(size(EXPRVARIABLES)) :: Rval
    integer :: cnormalMean
    
    dresult = 0.0_DP
    
    select case (rbcAssemblyData%cexprType)
    case (BDC_USERDEF)
      ! This is a hardcoded, user-defined identifier.
      ! Get the expression identifier from the collection, it is a string.
      call cc_getSvalue (rbcAssemblyData%ranalyticBC,rbcAssemblyData%sexprName,stag)
                              
      ! Call the user defined callback routine to evaluate the expression.
      ! Pass the user defined collection.
      call getBoundaryValues (stag,rbcAssemblyData%icomponent,&
          rbcAssemblyData%p_rdiscretisation%RspatialDiscr(rbcAssemblyData%icomponent),&
          rbcAssemblyData%rboundaryRegion,&
          dpar, dresult, rbcAssemblyData%p_ruserCollection)
      
    case (BDC_VALDOUBLE)
      ! A simple constant, given by dvalue
      dresult = rbcAssemblyData%dvalue

    case (BDC_EXPRESSION)
      ! A complex expression.
      !
      ! Set up an array with variables for evaluating the expression.
      ! Give the values in exactly the same order as specified
      ! by EXPRVARIABLES!
      Rval = 0.0_DP
      
      call boundary_getCoords(rbcAssemblyData%p_rboundary, &
          rbcAssemblyData%rboundaryRegion%iboundCompIdx, dpar, dx, dy)
      
      ! Get the local parameter value 0 <= d <= 1 in the boundary region.
      ! Note that if dpar < rboundaryRegion%dminParam, we have to add the maximum
      ! parameter value on the boundary to dpar as normally 0 <= dpar < max.par.
      ! although 0 <= dminpar <= max.par
      !      and 0 <= dmaxpar <= max.par!
      d = dpar
      if (d .lt. rbcAssemblyData%rboundaryRegion%dminParam) then
        d = d + boundary_dgetMaxParVal(rbcAssemblyData%p_rboundary,&
            rbcAssemblyData%rboundaryRegion%iboundCompIdx)
      end if
      d = d - rbcAssemblyData%rboundaryRegion%dminParam

      ! Normalise to 0..1 using the length of the parameter region.
      ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
      d = d / (rbcAssemblyData%rboundaryRegion%dmaxParam &
               - rbcAssemblyData%rboundaryRegion%dminParam)
      
      Rval(1) = dx
      Rval(2) = dy
      Rval(3) = 0.0_DP ! Not used.
      Rval(4) = d
      Rval(5) = dpar
      Rval(6) = boundary_convertParameter(rbcAssemblyData%p_rboundary, &
          rbcAssemblyData%rboundaryRegion%iboundCompIdx, dpar, &
          BDR_PAR_01, BDR_PAR_LENGTH)
      Rval(7) = rbcAssemblyData%dtime
      Rval(8:9) = rbcAssemblyData%DmframeVel(1:NDIM2D)
      Rval(10:11) = rbcAssemblyData%DmframeAcc(1:NDIM2D)
      
      ! The normal vector is a bit tricky. Normally, we take the "mean"
      ! setting.
      cnormalMean = BDR_NORMAL_MEAN
      
      ! Exception: If the parameter value corresponds to the beginning
      ! of the interval, we take the setting 'left'. On the right end,
      ! we take the setting 'right'.
      ! Reason: If this happens, the left/right point belongs to the interval,
      ! so it is more likely that the normal vector should be guided by
      ! the interval.
      if (dpar .eq. rbcAssemblyData%rboundaryRegion%dminParam) then
        cnormalMean = BDR_NORMAL_LEFT
      else if (dpar .eq. rbcAssemblyData%rboundaryRegion%dmaxParam) then
        cnormalMean = BDR_NORMAL_RIGHT
      end if
      
      ! Get the normal vector in that point
      call boundary_getNormalVec2D(rbcAssemblyData%p_rboundary, &
          rbcAssemblyData%rboundaryRegion%iboundCompIdx, dpar, Rval(12), Rval(13), &
          cnormalMean)
      
      ! Evaluate the expression. ivalue is the number of
      ! the expression to evaluate.
      call fparser_evalFunction (rbcAssemblyData%ranalyticBC%rparser, &
          rbcAssemblyData%ivalue, Rval, dresult)
      
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
      if (d .lt. rbcAssemblyData%rboundaryRegion%dminParam) then
        d = d + boundary_dgetMaxParVal(rbcAssemblyData%p_rboundary,&
            rbcAssemblyData%rboundaryRegion%iboundCompIdx)
      end if
      d = d - rbcAssemblyData%rboundaryRegion%dminParam
      
      ! Normalise to 0..1 using the length of the parameter region.
      ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
      d = d / (rbcAssemblyData%rboundaryRegion%dmaxParam &
               - rbcAssemblyData%rboundaryRegion%dminParam)
  
      dresult = mprim_getParabolicProfile (d,1.0_DP,rbcAssemblyData%dvalue)
    end select
    
    ! Weight the result by the weighting factor.
    dresult = dresult * rbcAssemblyData%dweight
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleFBDconditions (rproblem,rdiscretisation,rdynamicLevelInfo,rcollection)

!<description>
  ! This parses the boundary conditions for fictitious boundary
  ! components in the problem and assembles them into rdiscreteFBC.
!</description>
  
!<input>
  ! A discretisation structure defining the discretisation of the current
  ! level.
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! Collection structure to be passed to callback routines
  type(t_collection), intent(inout), target :: rcollection
  
  ! A t_dynamicLevelInfo structure that receives a discretised version
  ! of the boundary boundary conditions. The BC substructure should
  ! be empty; new BC`s are simply added to the structure.
  type(t_dynamicLevelInfo), intent(inout) :: rdynamicLevelInfo
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
    ! CALL bcasm_newDirichletBConFBD (rdiscretisation,Iequations,&
    !     rdynamicLevelInfo%rdiscreteFBC,getBoundaryValuesFBC,rcollection)

  end subroutine

end module