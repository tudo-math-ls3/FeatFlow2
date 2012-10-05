!##############################################################################
!# ****************************************************************************
!# <name> spatialbc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for parsing boundary conditions from a DAT
!# file
!#
!# The following files can be found here:
!#
!# 1.) sbc_assembleBDconditions
!#     -> Assembles the definition of the BCs from sections given by DAT files
!#        and sets up an analytical boundary condition description.
!#
!# 2.) sbc_releaseBdRegionList
!#     -> Release memory allocated in sbc_assembleBDconditions.
!#
!# 3.) sbc_implementDirichletBC
!#     -> Implement Dirichlet boundary conditions.
!# </purpose>
!##############################################################################

module spatialbc

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
  use vectorfilters
  use element

  use fespacehierarchybase
  use fespacehierarchy
  
  use collection
  use convection
    
  use constantsdiscretisation
  
  use structuresgeneral
  use structuresboundaryconditions
  use structuresoperatorasm
  use structuresoptcontrol
  use structuresdiscretisation
  
  use timediscretisation

  use user_callback
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! A structure which is passed to the callback routine for assembly.
  type t_bcassemblyData
    
    ! Boundary conditions in the problem.
    type(t_optcBDC), pointer :: p_roptcBDC => null()
    
    ! Definition of the underlying boundary
    type(t_boundary), pointer :: p_rboundary => null()
    
    ! Discretisation structure currently under investigation.
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr => null()

    ! Underlying time discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()

    ! Current boundary region in process.
    type(t_boundaryRegion) :: rboundaryRegion
    
    ! Pointer to a user defined collection structure
    type(t_collection), pointer :: p_ruserCollection => null()
    
    ! Current time where to set up the BCs
    real(DP) :: dtime = 0.0_DP
    
    ! Weighting factor for the values on the boundary.
    real(DP) :: dweight = 1.0_DP

    !<!-- Boundary condition specific parameters -->
    
    ! Type of boundary condition. A BDC_xxxx constant.
    integer :: cbcType = 0
    
    ! component under consideration (1=x-vel, 2=y-vel,...)
    integer :: icomponent = 0
    
    ! Type of boundary control currently in assembly.
    ! =0: No boundary control.
    ! =1: L2 Dirichlet boundary control.
    ! =2: $H^{1/2}$ Dirichlet boundary control.
    ! If the value is <> 0, the control must be specified in p_rcontrol.
    integer :: cboundaryControl = 0
    
    ! Pointer to a control vector on the boundary.
    type(t_vectorBlock), pointer :: p_rcontrol => null()

    !<!-- Parameters vor the evaluation of expressions -->
    
    ! Expression type. One of the BDC_xxxx constants.
    integer :: cexprType = 0
    
    ! expression id
    integer :: iid = 0
    
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

  ! Assemble boundary conditions
  public :: sbc_assembleBDconditions
  
  ! Resets the structure for discrete boundary conditions.
  public :: sbc_resetBCstructure

  ! Initialises a boundary condition hierarchy corresponding to a
  ! FE space hierarchy.
  public :: sbch_initBDCHierarchy
  
  ! Cleans up a boundary condition hierarchy.
  public :: sbch_doneBDCHierarchy

  ! Resets the hierarchy discrete boundary conditions.
  public  :: sbch_resetBCstructure
  
  ! Assemble the inhomogeneous Neumann boundary conditions into a RHS
  public :: sbc_assembleInhomNeumannRHS
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_initBCassembly (roptcBDC,rspaceDiscr,rtimeDiscr,dtime,&
      rbcAssemblyData,ruserCollection,rcollection)

!<description>
  ! Prepares the assembly of bonudary conditions.
!</description>

!<input>
  ! Boundary conditions in the problem.
  type(t_optcBDC), intent(inout), target :: roptcBDC
  
  ! User defined collection structure to be passed to callback routines
  type(t_collection), intent(in), target :: ruserCollection
  
  ! Underlying space discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
  
  ! Underlying time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
  
  ! Point im time, the BC should be assembled according to.
  real(DP), intent(in) :: dtime
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
    rbcAssemblyData%p_roptcBDC => roptcBDC
    rbcAssemblyData%p_rboundary => rspaceDiscr%p_rboundary
    rbcAssemblyData%p_ruserCollection => ruserCollection
    rbcAssemblyData%p_rspaceDiscr => rspaceDiscr
    rbcAssemblyData%p_rtimeDiscr => rtimeDiscr
    rbcAssemblyData%dtime = dtime

    ! Create a pointer to this, encapsuled in a structure.
    rp_bcassemblyData%p_rbcAssemblyData => rbcAssemblyData
    
    ! Transfer this pointer into the IquickAccess array of the collection.
    ! We can get it from there later.
    rcollection%IquickAccess(:) = &
        transfer(rp_bcassemblyData,rcollection%IquickAccess(:),size(rcollection%IquickAccess(:)))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_getBCassemblyPointer (rcollection,p_rbcAssemblyData)

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

  subroutine sbc_defineBdComponent (rbcAssemblyData,dpar1,dpar2,ibct,iintervalsEnd)

!<description>
  ! Initialises the boundary component structure in rbcAssemblyData.
!</description>

!<input>
  ! Start parameter value
  real(DP), intent(in) :: dpar1

  ! End parameter value
  real(DP), intent(in) :: dpar2

  ! Boundary component
  integer, intent(in) :: ibct
  
  ! Specifier that defines if the interval ends belong to the
  ! boundary region.
  ! =0: no
  ! =1: start belongs to the region
  ! =2: end belongs to the region
  ! =3: start and end belongs to the region
  integer, intent(in) :: iintervalsEnd

!</input>
  
!<inputoutput>
  ! Boundary condition assembly data, to be set up
  type(t_bcassemblyData), intent(inout) :: rbcAssemblyData
!</inputoutput>

!</subroutine>

    rbcAssemblyData%rboundaryRegion%dminParam = dpar1
    rbcAssemblyData%rboundaryRegion%dmaxParam = dpar2
    rbcAssemblyData%rboundaryRegion%iboundCompIdx = ibct
    rbcAssemblyData%rboundaryRegion%dmaxParamBC = &
        boundary_dgetMaxParVal(rbcAssemblyData%p_rboundary, ibct)
    rbcAssemblyData%rboundaryRegion%iproperties = iintervalsEnd

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_defineConstValue (rbcAssemblyData,dvalue)

!<description>
  ! Initialises an expression in rbcAssemblyData which returns a constant
  ! value dvalue.
!</description>

!<input>
  ! Value to be returned
  real(DP), intent(in) :: dvalue
!</input>
  
!<inputoutput>
  ! Boundary condition assembly data, to be set up
  type(t_bcassemblyData), intent(inout) :: rbcAssemblyData
!</inputoutput>

!</subroutine>

    rbcAssemblyData%cexprType = BDC_VALDOUBLE
    rbcAssemblyData%iid = 0
    rbcAssemblyData%ivalue = 0
    rbcAssemblyData%dvalue = dvalue
    rbcAssemblyData%sexprName = ""
    
    rbcAssemblyData%cboundaryControl = 0
    nullify(rbcAssemblyData%p_rcontrol)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_assembleBDconditions (roptcBDC,roptcBDCSpace,dtime,cequation,&
      copType,casmFlags,rspaceDiscr,rtimeDiscr,rglobalData,rvectorControl)

!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
  !
  ! The bhasNeumannBoudary flag on every level is initialised according to
  ! whether there exist Neumann boundary components on the boundary or not.
!</description>
  
!<input>
  ! Boundary conditions in the problem.
  type(t_optcBDC), intent(inout), target :: roptcBDC

  ! Current simulation time where to assemble the boundary conditions.
  real(dp), intent(in) :: dtime

  ! Underlying equation. One of the CEQ_xxxx constants.
  integer, intent(in) :: cequation

  ! Type of equation to be solved here. This can be 
  ! OPTP_PRIMAL, OPTP_DUAL, OPTP_PRIMALLIN or OPTP_DUALLIN,
  ! depending on which equation to solve.
  integer :: copType
  
  ! Assembly flags, determines what to assemble.
  ! One of the SBC_xxxx flags.
  integer, intent(in) :: casmFlags
  
  ! A discretisation structure defining the space discretisation
  ! of the current level.
  type(t_blockDiscretisation), intent(in) :: rspaceDiscr

  ! Time discretisation structure, dtime refers to.
  type(t_timeDiscretisation), intent(in) :: rtimeDiscr

  ! OPTIONAL; Global data, passed to callback routines
  ! Must be specified if casmFlags contains SBC_DISCRETEBC.
  type(t_globalData), intent(inout), optional :: rglobalData

  ! OPTIONAL: A control specifying the current control.
  ! Used for boundary control.
  ! Must be specified if casmFlags contains SBC_DISCRETEBC.
  type(t_vectorBlock), intent(in), target, optional :: rvectorControl
!</input>

!<inputoutput>
  ! Boundary condition structure that receives a definition of the boundary
  ! conditions. The assembled boundary conditions are appended to the
  ! elements in this structure.
  type(t_optcBDCSpace), intent(inout) :: roptcBDCSpace
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ibct, isegment, nsegments, nsegmentsParent, i, iintervalsEnd
    character(LEN=PARLST_LENLINEBUF) :: sbdex1,sbdex2,sbcsection,sparentSection
    character(LEN=SYS_STRLEN) :: sparams
    real(DP) :: dpar1, dpar2
    integer :: iindex,iindexParent
    logical :: bautomatic, bautomaticParent, bautomaticSegment
    type(t_collection), target :: rcollection, ruserCollection
    type(t_boundary), pointer :: p_rboundary
    type(t_bcassemblyData) :: rbcAssemblyData

    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection

    ! Fetch some parameters
    p_rboundary => rspaceDiscr%p_rboundary
    
    ! Get the expression/bc sections from the bondary condition block
    call parlst_querysection(roptcBDC%p_rparList, &
        roptcBDC%ssectionBdExpr, p_rsection)
        
    ! Get the section defining the primal or dual boundary conditions
    sbcsection = ""
    sparentSection = ""
    select case (copType)
    case (OPTP_PRIMAL)

      ! Primal boundary conditions
      sbcsection = roptcBDC%ssectionBdCondPrim

    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)

      ! Primal boundary conditions, linearised eqn.
      sbcsection = roptcBDC%ssectionBdCondPrimLin

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    case (OPTP_DUAL)
    
      ! Dual boundary conditions
      sbcsection = roptcBDC%ssectionBdCondDual

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    case (OPTP_DUALLIN,OPTP_DUALLIN_SIMPLE)
    
      ! Dual boundary conditions, linearised equation
      sbcsection = roptcBDC%ssectionBdCondDualLin

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    case (OPTP_PCSTEKLOV)
    
      ! Dual boundary conditions, linearised equation
      sbcsection = roptcBDC%ssectionBdCondPCSteklov

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    end select
        
    ! Initialise the user-defined assembly
    call collct_init (rusercollection)
    
    ! Basic initialisation of the BC assembly structure
    call sbc_initBCassembly (roptcBDC,rspaceDiscr,rtimeDiscr,dtime,&
        rbcAssemblyData,ruserCollection,rcollection)

    ! Loop through all boundary components we have.
    do ibct = 1,boundary_igetNBoundComp(p_rboundary)

      ! Get information about that component: Number of segments,...
      call struc_getBDCinfo (roptcBDC,sbcsection,&
          ibct,iindex,nsegments,bautomatic)
          
      if (bautomatic) then
        if (sparentSection .eq. "") then
          call output_line ("Automatic boundary conditions only allowed.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
          call sys_halt()
        end if
      end if

      if (sparentSection .ne. "") then
        ! Structure of the parent node        
        call struc_getBDCinfo (roptcBDC,sparentSection,&
            ibct,iindexParent,nsegmentsParent,bautomaticParent)

        ! Synchronise the number of segments.            
        if (nsegments .eq. 0) then
        
          nsegments = nsegmentsParent
        
        else if (bautomatic .and. (nsegments .ne. nsegmentsParent)) then

          call output_line ("Automatic boundary conditions invalid.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
          call sys_halt()

        end if
            
      end if

      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      ! Parameter exists. Get the values in there.
      do isegment = 1,nsegments
        
        ! Get information about the segment.
        !
        ! Is the component fully automatic? Is the segment automatic?
        bautomaticSegment = .false.
        
        if (bautomatic) then
          ! Get the info directly from the parent
          call struc_getSegmentInfo (&
                roptcBDC,sparentSection,iindexParent,isegment,dpar2,iintervalsEnd,&
                rbcAssemblyData%cbcType,sparams)
        else
          ! Get the information from the section
          call struc_getSegmentInfo (&
                roptcBDC,sbcSection,iindex,isegment,dpar2,iintervalsEnd,&
                rbcAssemblyData%cbcType,sparams)
          
          ! If the type of this segment is "automatic", get it from the parent
          if (rbcAssemblyData%cbcType .eq. -1) then
            bautomaticSegment = .true.
            call struc_getSegmentInfo (&
                  roptcBDC,sparentSection,iindexParent,isegment,dpar2,iintervalsEnd,&
                  rbcAssemblyData%cbcType,sparams)
          end if
        end if

        ! Form a boundary condition segment that covers that boundary part
        if (dpar2 .ge. dpar1) then

          ! Define the corresponding boundary region
          call sbc_defineBdComponent (rbcAssemblyData,dpar1,dpar2,ibct,iintervalsEnd)
          
          ! Type of equation?
          select case (cequation)
          
          ! ----------------------
          ! Stokes / Navier-Stokes
          ! ----------------------
          case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
          
            ! Do we have 'automatic' boundary conditions?
            if (.not. bautomatic .and. .not. bautomaticSegment) then
            
              ! -----------------------------------------
              ! Boundary conditions by definition
              ! -----------------------------------------

              ! Now, which type of BC is to be created?
              select case (rbcAssemblyData%cbcType)
              
              ! --------------------------------------------------
              ! Homogeneous/Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (0,9)
                if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                  ! Add the bondary region to the Neumann boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                end if
              
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)

                if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then

                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletBoundary)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                  
                    ! Simple Dirichlet boundary.
                    ! Get the definition of the boundary condition.
                    ! Read the line again, get the expressions for X- and Y-velocity
                    read(sparams,*) sbdex1,sbdex2
                  
                    ! For any string <> "", create the appropriate Dirichlet boundary
                    ! condition and add it to the list of boundary conditions.
                    
                    if (sbdex1 .ne. "") then
                      
                      ! X-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex1,&
                          rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                          rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                          
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = 1
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)
                      
                    end if
                        
                    if (sbdex2 .ne. "") then
                    
                      ! Y-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex2,&
                          rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                          rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                          
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = 2
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)
                      
                    end if
                    
                  end if ! SBC_DISCRETEBC

                end if

              ! --------------------------------------------------
              ! Dirichlet boundary control
              ! --------------------------------------------------
              case (7)

                if (iand(casmFlags,SBC_DIRICHLETBCC) .ne. 0) then

                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletControlBoundaryL2)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                    ! Init boundary control
                    rbcAssemblyData%cboundaryControl = 1
                    rbcAssemblyData%p_rcontrol => rvectorControl

                    ! Get the definition of the boundary condition.
                    ! Read the line again, get the expressions for X- and Y-velocity
                    read(sparams,*) sbdex1,sbdex2
                  
                    ! For any string <> "", create the appropriate Dirichlet boundary
                    ! condition and add it to the list of boundary conditions.
                    
                    if (sbdex1 .ne. "") then
                      
                      ! X-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex1,&
                          rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                          rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                          
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = 1
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)
                      
                    end if
                        
                    if (sbdex2 .ne. "") then
                    
                      ! Y-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex2,&
                          rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                          rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                          
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = 2
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)

                    end if
                    
                  end if ! SBC_DISCRETEBC

                end if

              ! --------------------------------------------------
              ! H^1/2 Dirichlet boundary control
              ! --------------------------------------------------
              case (8)

                if (iand(casmFlags,SBC_DIRICHLETBCC) .ne. 0) then

                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletControlBoundaryH12)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                    ! Init boundary control
                    rbcAssemblyData%cboundaryControl = 2
                    rbcAssemblyData%p_rcontrol => rvectorControl

                    ! Get the definition of the boundary condition.
                    ! Read the line again, get the expressions for X- and Y-velocity
                    read(sparams,*) sbdex1,sbdex2
                  
                    ! For any string <> "", create the appropriate Dirichlet boundary
                    ! condition and add it to the list of boundary conditions.

                    if (sbdex1 .ne. "") then
                      
                      ! X-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex1,&
                          rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                          rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                          
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = 1
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)
                      
                    end if
                        
                    if (sbdex2 .ne. "") then
                    
                      ! Y-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex2,&
                          rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                          rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                          
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = 2
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)

                    end if
                    
                  end if ! SBC_DISCRETEBC

                end if

              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
              end select ! cbctyp
              
            else

              ! -----------------------------------------
              ! Automatic boundary conditions
              ! -----------------------------------------

              ! The rule is:
              !  * Primal Neumann = Dual Neumann,
              !  * Primal Dirichlet = Dual Dirichlet-0

              select case (rbcAssemblyData%cbcType)
              
              ! --------------------------------------------------
              ! Neumann boundary conditions
              ! --------------------------------------------------
              case (0)
                if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                  ! Add the bondary region to the Neumann boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                end if
                
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)
                if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
                
                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletBoundary)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                    
                    ! Impose Dirichlet-0 boundary conditions in the dual equation
                    call sbc_defineConstValue (rbcAssemblyData,0.0_DP)
                  
                    ! X/Y-velocity
                    do i = 1,2
                    
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = i
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)

                    end do
                    
                  end if
                  
                end if
                
              ! --------------------------------------------------
              ! L2 Dirichlet boundary control
              ! --------------------------------------------------
              case (7)
              
                if (iand(casmFlags,SBC_DIRICHLETBCC) .ne. 0) then
                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletControlBoundaryL2)
                      
                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then

                    ! Impose Dirichlet-0 boundary conditions 
                    ! plus boundary control
                    call sbc_defineConstValue (rbcAssemblyData,0.0_DP)

                    if ((copType .eq. OPTP_PRIMALLIN) .or. (copType .eq. OPTP_PRIMALLIN_SIMPLE) &
                        .or. (copType .eq. OPTP_PCSTEKLOV)) then
                      ! Boundary control
                      rbcAssemblyData%cboundaryControl = 1
                      rbcAssemblyData%p_rcontrol => rvectorControl
                    end if

                    ! X/Y-velocity
                    do i = 1,2
                    
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = i
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)
                      
                    end do
                    
                  end if

                end if
                
              ! --------------------------------------------------
              ! H^1/2 Dirichlet boundary control
              ! --------------------------------------------------
              case (8)
              
                if (iand(casmFlags,SBC_DIRICHLETBCC) .ne. 0) then
                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletControlBoundaryH12)
                      
                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then

                    ! Impose Dirichlet-0 boundary conditions 
                    ! plus boundary control
                    call sbc_defineConstValue (rbcAssemblyData,0.0_DP)

                    if ((copType .eq. OPTP_PRIMALLIN) .or. (copType .eq. OPTP_PRIMALLIN_SIMPLE) &
                        .or. (copType .eq. OPTP_PCSTEKLOV)) then
                      ! Boundary control
                      rbcAssemblyData%cboundaryControl = 2
                      rbcAssemblyData%p_rcontrol => rvectorControl
                    end if

                    ! X/Y-velocity
                    do i = 1,2
                    
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = i
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)

                    end do
                  
                  end if
                
                end if

              case (9)
                ! Nothing to do here
                
              case default
                call output_line ("Cannot set up automatic boundary conditions", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
                
              end select
              
            end if
            
          ! ----------------------
          ! Heat equation
          ! ----------------------
          case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
          
            ! Do we have 'automatic' boundary conditions?
            if (.not. bautomatic .and. .not. bautomaticSegment) then
            
              ! -----------------------------------------
              ! Boundary conditions by definition
              ! -----------------------------------------
              
              ! Now, which type of BC is to be created?
              select case (rbcAssemblyData%cbcType)
              
              ! --------------------------------------------------
              ! Homogeneous/Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (0,9)
              
                if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                  ! Add the bondary region to the Neumann boundary regions
                  ! if there is a structure present.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                end if
              
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)

                if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then

                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletBoundary)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                  
                    ! Simple Dirichlet boundary.
                    ! Get the definition of the boundary condition.
                    ! Read the line again, get the expressions for the data field
                    read(sparams,*) sbdex1
                  
                    ! For any string <> "", create the appropriate Dirichlet boundary
                    ! condition and add it to the list of boundary conditions.
                    
                    if (sbdex1 .ne. "") then
                      
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex1,&
                          rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                          rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                          
                      ! Assemble the BCs.
                      rbcAssemblyData%icomponent = 1
                      call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                          rbcAssemblyData,rcollection,rglobalData)

                    end if
                    
                  end if ! SBC_DISCRETEBC
                      
                end if
                
              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
              end select ! cbctyp

            else

              ! -----------------------------------------
              ! Automatic boundary conditions
              ! -----------------------------------------

              ! The rule is:
              !  * Primal Neumann = Dual Neumann,
              !  * Primal Dirichlet = Dual Dirichlet-0

              select case (rbcAssemblyData%cbcType)
              
              ! --------------------------------------------------
              ! Neumann boundary conditions
              ! --------------------------------------------------
              case (0)
              
                if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                  ! Add the bondary region to the Neumann boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                end if
                
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)
              
                if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
                  
                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rbcAssemblyData%rboundaryRegion,roptcBDCSpace%rdirichletBoundary)
                      
                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                  
                    ! Impose Dirichlet-0 boundary conditions in the dual equation
                    call sbc_defineConstValue (rbcAssemblyData,0.0_DP)
                  
                    ! Assemble the BCs.
                    rbcAssemblyData%icomponent = 1
                    call sbc_assembleDirichletBC (roptcBDCSpace%rdiscreteBC,&
                        rbcAssemblyData,rcollection,rglobalData)
                    
                  end if
                
                end if
                
              case (9)
                ! Nothing to do here

              case default
                call output_line ("Cannot set up automatic boundary conditions", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
                
              end select

            end if            
          
          case default
            
            call output_line ("Unknown equation", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
            call sys_halt()
            
          end select ! equation
          
          ! Move on to the next parameter value
          dpar1 = dpar2
          
        end if
                                          
      end do ! isegment
      
    end do ! ibct

    call collct_done (rusercollection)
        
  end subroutine


  ! ***************************************************************************

!<subroutine>
  
  subroutine sbc_assembleDirichletBC (rdiscreteBC,rbcAssemblyData,rcollection,rglobalData)
  
!<description>
  ! Assembles Dirichlet boundary conditions
  ! using the settings in the rbcAssemblyData structure.
!</description>
  
!<input>
  ! Boundary condition assembly data.
  type(t_bcassemblyData), intent(inout), target :: rbcAssemblyData

  ! Collection structure prepared by sbc_initBCassembly.
  type(t_collection), intent(inout) :: rcollection

  ! OPTIONAL; Global data, passed to callback routines
  ! Must be specified if casmFlags contains SBC_DISCRETEBC.
  type(t_globalData), intent(inout), optional :: rglobalData
!</input>

!<output>
  ! Boundary condition structure that receives the discrete BC
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</output>

!</subroutine>

    ! Assemble the BCs.
    call user_initCollectForVecAssembly (&
        rglobalData,rbcAssemblyData%iid,rbcAssemblyData%icomponent,&
        rbcAssemblyData%dtime,rbcAssemblyData%p_ruserCollection)
    
    call bcasm_newDirichletBConRealBD (&
        rbcAssemblyData%p_rspaceDiscr,rbcAssemblyData%icomponent,&
        rbcAssemblyData%rboundaryRegion,&
        rdiscreteBC,cc_getDirBC,rcollection)
        
    call user_doneCollectForVecAssembly (rglobalData,rbcAssemblyData%p_ruserCollection)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine sbc_evalBoundaryValue (rbcAssemblyData,dpar,dresult)
  
!<description>
  ! Evaluates an expression in a point on the boundary.
!</description>
  
!<input>
  ! Boundary condition assembly data.
  type(t_bcassemblyData), intent(in), target :: rbcAssemblyData

  ! Current parameter value of the point on the boundary.
  ! 0-1-parametrisation.
  real(DP), intent(in) :: dpar
!</input>

!<output>
  ! Value of the expression in the point.
  real(DP), intent(out) :: dresult
!</output>

!</subroutine>

    ! local variables
    real(DP) :: d,dx,dy
    real(DP), dimension(size(SEC_EXPRVARIABLES)) :: Rval
    integer :: cnormalMean
    
    dresult = 0.0_DP
    
    select case (rbcAssemblyData%cexprType)
        
    case (BDC_USERDEFID)
      ! Defined via analytic reference solution, identified by an id.
      ! The id is already specified in the collection.
      
      ! Call the user defined callback routine to evaluate the boudary value.
      ! No string tag is given which indicates that a reference function is to be
      ! evaluated.
      call user_getBoundaryValues ("",rbcAssemblyData%icomponent,&
          rbcAssemblyData%p_rspaceDiscr%RspatialDiscr(rbcAssemblyData%icomponent),&
          rbcAssemblyData%rboundaryRegion,dpar,dresult,&
          rbcAssemblyData%p_ruserCollection)
      
    case (BDC_USERDEF)
      ! This is a hardcoded, user-defined identifier.
      ! Call the user defined callback routine to evaluate the expression.
      call user_getBoundaryValues (rbcAssemblyData%sexprName,rbcAssemblyData%icomponent,&
          rbcAssemblyData%p_rspaceDiscr%RspatialDiscr(rbcAssemblyData%icomponent),&
          rbcAssemblyData%rboundaryRegion,dpar,dresult,&
          rbcAssemblyData%p_ruserCollection)
      
    case (BDC_VALDOUBLE)
      ! A simple constant, given by dvalue
      dresult = rbcAssemblyData%dvalue

    case (BDC_EXPRESSION)
      ! A complex expression.
      !
      ! Set up an array with variables for evaluating the expression.
      ! Give the values in exactly the same order as specified
      ! by SEC_EXPRVARIABLES!
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
          rbcAssemblyData%rboundaryRegion%iboundCompIdx, dpar, Rval(8), Rval(9), &
          cnormalMean)
      
      ! Evaluate the expression. ivalue is the number of
      ! the expression to evaluate.
      call fparser_evalFunction (rbcAssemblyData%p_roptcBDC%rparser, &
          rbcAssemblyData%iid, Rval, dresult)
      
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

  subroutine cc_getDirBC (Icomponents,rspaceDiscr,rboundaryRegion,ielement, &
      cinfoNeeded, iwhere, dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates the Dirichlet boundary conditions and is
  ! used by the discretisation routines to generate a discrete
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
  ! framework to ignore the node and treat it as "natural boundary condition"
  ! node.
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    integer :: icontrolcomp
    type(t_settings_optcontrol), pointer :: p_roptControl
    type(t_settings_physics), pointer :: p_rphysics
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: iedge,iptpos,iv1,iv2
    real(DP) :: dpar1, dpar2, dpar, dval
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    type(t_vectorBlock), pointer :: p_rcontrol
    type(t_bcassemblyData), pointer :: p_rbcAssemblyData
    
    ! Get the assembly information structure.
    call sbc_getBCassemblyPointer (rcollection,p_rbcAssemblyData)
    
    p_roptControl => p_rbcAssemblyData%p_roptcBDC%p_rsettingsOptControl
    p_rphysics => p_rbcAssemblyData%p_roptcBDC%p_rphysics
    p_rcontrol => p_rbcAssemblyData%p_rcontrol
    
    ! Use boundary conditions from DAT files.
    select case (cinfoNeeded)
    
    case (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
          DISCBC_NEEDINTMEAN,DISCBC_NEEDNORMALSTRESS)
      
      ! Dirichlet boundary conditions.
      !
      ! Evaluate...
      call sbc_evalBoundaryValue (p_rbcAssemblyData,dwhere,Dvalues(1))
      
    end select
    
    select case (p_rbcAssemblyData%cboundaryControl)
    
    ! ---------------------------------
    ! L2 Dirichlet boudary control
    ! ---------------------------------
    case (1)
    
      if (p_roptControl%dalphaL2BdC .ge. 0.0_DP) then
        if (.not. associated(p_rcontrol)) then
          call output_line ("Control not specified.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBC")
          call sys_halt()
        end if
        
        select case (p_rphysics%cequation)
        case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
        
          ! Where is the boundary control?
          icontrolcomp = p_rbcAssemblyData%icomponent
          if (p_roptControl%dalphaDistC .ge. 0.0_DP) icontrolcomp = icontrolcomp + 2

          ! Quick and dirty implementation...
          
          ! Element type?
          select case (p_rcontrol%p_rblockDiscr%RspatialDiscr(icontrolcomp)%RelementDistr(1)%celement)
          case (EL_P2_1D)
            
            ! Type of info?
            select case (cinfoNeeded)
            case (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
                  DISCBC_NEEDINTMEAN,DISCBC_NEEDNORMALSTRESS)
              
              ! Search the boundary edge that contains the parameter value
              call tria_searchBoundaryEdgePar2D(rboundaryRegion%iboundCompIdx, dwhere, &
                  rspaceDiscr%p_rtriangulation, rspaceDiscr%p_rboundary, iedge,&
                  BDR_PAR_LENGTH)
              !print *,"Richtig?"
              
              ! Get the vector data
              call lsyssc_getbase_double (p_rcontrol%RvectorBlock(icontrolcomp),p_Ddata)
              
              ! Get the "relative" parameter value of dwhere.
              call storage_getbase_int (&
                  rspaceDiscr%p_rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
              call storage_getbase_double (&
                  rspaceDiscr%p_rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
                  
              iptpos = p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)
              
              ! First vertex number = edge number
              dpar1 = p_DvertexParameterValue(iptpos-1+iedge)
              iv1 = iedge
              
              ! Second vertex number = first vertex + 1 or overall first
              if (iedge .lt. p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx+1) &
                            -p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)) then
                dpar2 = p_DvertexParameterValue(iptpos-1+iedge+1)
                iv2 = iedge+1
              else
                dpar2 = boundary_dgetMaxParVal(rspaceDiscr%p_rboundary, &
                    rboundaryRegion%iboundCompIdx)
                iv2 = p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)
              end if
              
              ! Relative parameter value
              dpar = (boundary_convertParameter(rspaceDiscr%p_rboundary, &
                  rboundaryRegion%iboundCompIdx, dwhere, &
                  BDR_PAR_LENGTH, BDR_PAR_01)-dpar1)/(dpar2-dpar1)
                  
              ! Quadratic interpolation over [-1,1]
              call mprim_quadraticInterpolation (2*dpar-1.0_DP,&
                  p_Ddata(iv1),p_Ddata(iedge+rspaceDiscr%p_rtriangulation%nvbd),p_Ddata(iv2),dval)
              
              Dvalues(1) = Dvalues(1) + dval

            end select          
          
          case default
            call output_line ("Unsupported element.", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBC")
            call sys_halt()
          end select
        
        
        case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
          call output_line ("Equation not supported.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBC")
          call sys_halt()
        end select
        
      end if
    
    ! ---------------------------------
    ! H^1/2 Dirichlet boudary control
    ! ---------------------------------
    case (2)
    
      if (p_roptControl%dalphaH12BdC .ge. 0.0_DP) then
        if (.not. associated(p_rcontrol)) then
          call output_line ("Control not specified.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBC")
          call sys_halt()
        end if
        
        select case (p_rphysics%cequation)
        case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
        
          ! Where is the boundary control?
          icontrolcomp = p_rbcAssemblyData%icomponent
          if (p_roptControl%dalphaDistC .ge. 0.0_DP) icontrolcomp = icontrolcomp + 2
          if (p_roptControl%dalphaL2BdC .ge. 0.0_DP) icontrolcomp = icontrolcomp + 2

          ! Quick and dirty implementation...
          
          ! Element type?
          select case (p_rcontrol%p_rblockDiscr%RspatialDiscr(icontrolcomp)%RelementDistr(1)%celement)
          case (EL_P2_1D)
            
            ! Type of info?
            select case (cinfoNeeded)
            case (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
                  DISCBC_NEEDINTMEAN,DISCBC_NEEDNORMALSTRESS)
              
              ! Search the boundary edge that contains the parameter value
              call tria_searchBoundaryEdgePar2D(rboundaryRegion%iboundCompIdx, dwhere, &
                  rspaceDiscr%p_rtriangulation, rspaceDiscr%p_rboundary, iedge,&
                  BDR_PAR_LENGTH)
              !print *,"Richtig?"
              
              ! Get the vector data
              call lsyssc_getbase_double (p_rcontrol%RvectorBlock(icontrolcomp),p_Ddata)
              
              ! Get the "relative" parameter value of dwhere.
              call storage_getbase_int (&
                  rspaceDiscr%p_rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
              call storage_getbase_double (&
                  rspaceDiscr%p_rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
                  
              iptpos = p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)
              
              ! First vertex number = edge number
              dpar1 = p_DvertexParameterValue(iptpos-1+iedge)
              iv1 = iedge
              
              ! Second vertex number = first vertex + 1 or overall first
              if (iedge .lt. p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx+1) &
                            -p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)) then
                dpar2 = p_DvertexParameterValue(iptpos-1+iedge+1)
                iv2 = iedge+1
              else
                dpar2 = boundary_dgetMaxParVal(rspaceDiscr%p_rboundary, &
                    rboundaryRegion%iboundCompIdx)
                iv2 = p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)
              end if
              
              ! Relative parameter value
              dpar = (boundary_convertParameter(rspaceDiscr%p_rboundary, &
                  rboundaryRegion%iboundCompIdx, dwhere, &
                  BDR_PAR_LENGTH, BDR_PAR_01)-dpar1)/(dpar2-dpar1)
                  
              ! Quadratic interpolation over [-1,1]
              call mprim_quadraticInterpolation (2*dpar-1.0_DP,&
                  p_Ddata(iv1),p_Ddata(iedge+rspaceDiscr%p_rtriangulation%nvbd),p_Ddata(iv2),dval)
              
              Dvalues(1) = Dvalues(1) + dval

            end select          
          
          case default
            call output_line ("Unsupported element.", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBC")
            call sys_halt()
          end select
        
        
        case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
          call output_line ("Equation not supported.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBC")
          call sys_halt()
        end select
        
      end if
    
    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_resetBCstructure (roptcBDCSpace)
  
!<description>
  ! Resets the structure for discrete boundary conditions.
!</description>

!<inputoutput>
  ! Boundary condition structure to reset.
  type(t_optcBDCSpace), intent(inout) :: roptcBDCSpace
!</inputoutput>

!</subroutine>

    ! Reset the discrete boundary conditions
    call bcasm_clearDiscreteBC(roptcBDCSpace%rdiscreteBC)
    
    ! Release boundary region lists
    call sbc_releaseBdRegionList (roptcBDCSpace%rneumannBoundary)
    call sbc_releaseBdRegionList (roptcBDCSpace%rdirichletBoundary)
    call sbc_releaseBdRegionList (roptcBDCSpace%rdirichletControlBoundaryL2)
    call sbc_releaseBdRegionList (roptcBDCSpace%rdirichletControlBoundaryH12)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbch_initBDCHierarchy (roptcBDCHierarchy,nlmin,nlmax)
  
!<description>
  ! Initialises a boundary condition hierarchy corresponding to a
  ! FE space hierarchy.
!</desctiption>

!<input>
  ! Minimum level in the hierarchy.
  integer, intent(in) :: nlmin

  ! Maximum level in the hierarchy.
  integer, intent(in) :: nlmax
!</input>

!<output>
  ! Boundary condition hierarchy to initialise.
  type(t_optcBDCSpaceHierarchy), intent(out) :: roptcBDCHierarchy
!</output>

!</subroutine>

    integer :: i

    ! Set up the structure between level nlmin and nlmax
    roptcBDCHierarchy%nlmin = nlmin
    roptcBDCHierarchy%nlmax = nlmax
    
    allocate(roptcBDCHierarchy%p_RoptcBDCspace(roptcBDCHierarchy%nlmin:roptcBDCHierarchy%nlmax))

    ! Initialise the boundary condition structures.
    do i=roptcBDCHierarchy%nlmin,roptcBDCHierarchy%nlmax
      call bcasm_initDiscreteBC(roptcBDCHierarchy%p_RoptcBDCspace(i)%rdiscreteBC)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbch_doneBDCHierarchy (roptcBDCHierarchy)
  
!<description>
  ! Cleans up a boundary condition hierarchy.
!</desctiption>

!<inputoutput>
  ! Boundary condition hierarchy to clean up.
  type(t_optcBDCSpaceHierarchy), intent(inout) :: roptcBDCHierarchy
!</inputoutput>

!</subroutine>

    integer :: i

    ! Initialise the boundary condition structures.
    do i=roptcBDCHierarchy%nlmax,roptcBDCHierarchy%nlmin,-1
      call bcasm_releaseDiscreteBC(roptcBDCHierarchy%p_RoptcBDCspace(i)%rdiscreteBC)
    end do

    ! Release the memory
    deallocate(roptcBDCHierarchy%p_RoptcBDCspace)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbch_resetBCstructure (roptcBDCHierarchy)
  
!<description>
  ! Resets the hierarchy discrete boundary conditions.
!</description>

!<inputoutput>
  ! Boundary condition hierarchy to clean up.
  type(t_optcBDCSpaceHierarchy), intent(inout) :: roptcBDCHierarchy
!</inputoutput>

!</subroutine>

    integer :: i
    
    do i=roptcBDCHierarchy%nlmin,roptcBDCHierarchy%nlmax
      call sbc_resetBCstructure (roptcBDCHierarchy%p_RoptcBDCspace(i))
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getNeumannBdData (rdiscretisation, rform, &
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

    integer :: ipt,iel
    real(DP) :: dwhere
    type(t_bcassemblyData), pointer :: p_rbcAssemblyData
    
    ! Get the assembly information structure.
    call sbc_getBCassemblyPointer (rcollection,p_rbcAssemblyData)

    do iel = 1,nelements
      do ipt = 1,npointsPerElement

        ! Get the parameter value of the point in 0-1 parametrisation
        dwhere = boundary_convertParameter(p_rbcAssemblyData%p_rboundary, &
            p_rbcAssemblyData%rboundaryRegion%iboundCompIdx, DpointPar(ipt,iel), &
            BDR_PAR_LENGTH, BDR_PAR_01)
                  
        ! Evaluate the expression.
        call sbc_evalBoundaryValue (p_rbcAssemblyData,dwhere,Dcoefficients(1,ipt,iel))

      end do  
      
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_assembleInhomNeumannRHS (rrhs,roptcBDC,rspaceDiscr,rtimeDiscr,&
      dtime,cequation,copType,rglobalData)

!<description>
  ! Assembles inhomogeneous Neumann boundary data into a RHS vector.
!</description>
  
!<input>
  ! Boundary conditions in the problem.
  type(t_optcBDC), intent(inout), target :: roptcBDC

  ! Underlying space discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
  
  ! Underlying time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr

  ! Current simulation time where to assemble the boundary conditions.
  real(dp), intent(in) :: dtime

  ! Underlying equation. One of the CEQ_xxxx constants.
  integer, intent(in) :: cequation

  ! Type of equation to be solved here. This can be 
  ! OPTP_PRIMAL, OPTP_DUAL, OPTP_PRIMALLIN or OPTP_DUALLIN,
  ! depending on which equation to solve.
  integer :: copType
  
  ! OPTIONAL; Global data, passed to callback routines
  ! Must be specified if casmFlags contains SBC_DISCRETEBC.
  type(t_globalData), intent(inout), optional :: rglobalData
!</input>

!<inputoutput>
  ! RHS vector, the inhomogeneous Neumann boundary data should be implemented to.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ibct, isegment, nsegments, nsegmentsParent, iintervalsEnd
    character(LEN=PARLST_LENLINEBUF) :: sbdex1,sbdex2,sbcsection,sparentSection
    character(LEN=SYS_STRLEN) :: sparams
    real(DP) :: dpar1, dpar2
    integer :: iindex,iindexParent
    logical :: bautomatic, bautomaticParent, bautomaticSegment
    type(t_collection), target :: rcollection, ruserCollection
    type(t_boundary), pointer :: p_rboundary
    type(t_bcassemblyData) :: rbcAssemblyData
    type(t_linearform) :: rlinform

    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection

    ! Fetch some parameters
    p_rboundary => rspaceDiscr%p_rboundary
    
    ! Get the expression/bc sections from the bondary condition block
    call parlst_querysection(roptcBDC%p_rparList, &
        roptcBDC%ssectionBdExpr, p_rsection)
        
    ! Get the section defining the primal or dual boundary conditions
    sbcsection = ""
    sparentSection = ""
    select case (copType)
    case (OPTP_PRIMAL)

      ! Primal boundary conditions
      sbcsection = roptcBDC%ssectionBdCondPrim

    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)

      ! Primal boundary conditions, linearised eqn.
      sbcsection = roptcBDC%ssectionBdCondPrimLin

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    case (OPTP_DUAL)
    
      ! Dual boundary conditions
      sbcsection = roptcBDC%ssectionBdCondDual

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    case (OPTP_DUALLIN,OPTP_DUALLIN_SIMPLE)
    
      ! Dual boundary conditions, linearised equation
      sbcsection = roptcBDC%ssectionBdCondDualLin

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    case (OPTP_PCSTEKLOV)
    
      ! Dual boundary conditions, linearised equation
      sbcsection = roptcBDC%ssectionBdCondPCSteklov

      ! Primal boundary conditions
      sparentSection = roptcBDC%ssectionBdCondPrim

    end select
        
    ! Initialise the user-defined assembly
    call collct_init (rusercollection)
    
    ! Basic initialisation of the BC assembly structure
    call sbc_initBCassembly (roptcBDC,rspaceDiscr,rtimeDiscr,dtime,&
        rbcAssemblyData,ruserCollection,rcollection)

    ! Loop through all boundary components we have.
    do ibct = 1,boundary_igetNBoundComp(p_rboundary)

      ! Get information about that component: Number of segments,...
      call struc_getBDCinfo (rbcAssemblyData%p_roptcBDC,sbcsection,&
          ibct,iindex,nsegments,bautomatic)
          
      if (bautomatic) then
        if (sparentSection .eq. "") then
          call output_line (&
              "Automatic boundary conditions only allowed.", &
              OU_CLASS_ERROR,OU_MODE_STD,"sbc_assembleInhomNeumannRHS")
          call sys_halt()
        end if
      end if

      if (sparentSection .ne. "") then
        ! Structure of the parent node        
        call struc_getBDCinfo (rbcAssemblyData%p_roptcBDC,sparentSection,&
            ibct,iindexParent,nsegmentsParent,bautomaticParent)

        ! Synchronise the number of segments.            
        if (nsegments .eq. 0) then
        
          nsegments = nsegmentsParent
        
        else if (bautomatic .and. (nsegments .ne. nsegmentsParent)) then

          call output_line ("Automatic boundary conditions invalid.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
          call sys_halt()

        end if

      end if

      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      ! Parameter exists. Get the values in there.
      do isegment = 1,nsegments
        
        ! Get information about the segment.
        !
        ! Is the component fully automatic? Is the segment automatic?
        bautomaticSegment = .false.
        
        if (bautomatic) then
          ! Get the info directly from the parent
          call struc_getSegmentInfo (&
                roptcBDC,sparentSection,iindexParent,isegment,dpar2,iintervalsEnd,&
                rbcAssemblyData%cbcType,sparams)
        else
          ! Get the information from the section
          call struc_getSegmentInfo (&
                roptcBDC,sbcSection,iindex,isegment,dpar2,iintervalsEnd,&
                rbcAssemblyData%cbcType,sparams)
          
          ! If the type of this segment is "automatic", get it from the parent
          if (rbcAssemblyData%cbcType .eq. -1) then
            bautomaticSegment = .true.
            call struc_getSegmentInfo (&
                  roptcBDC,sparentSection,iindexParent,isegment,dpar2,iintervalsEnd,&
                  rbcAssemblyData%cbcType,sparams)
          end if
        end if
        
        ! Form a boundary condition segment that covers that boundary part
        if (dpar2 .ge. dpar1) then
          
          ! Define the corresponding boundary region
          call sbc_defineBdComponent (rbcAssemblyData,dpar1,dpar2,ibct,iintervalsEnd)
          
          ! Type of equation?
          select case (cequation)
          
          ! ----------------------
          ! Stokes / Navier-Stokes
          ! ----------------------
          case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
          
            ! Do we have 'automatic' boundary conditions?
            ! We only have to assemble something for non-automatic BCs!
            if (.not. bautomatic .and. .not. bautomaticSegment) then
            
              ! -----------------------------------------
              ! Boundary conditions by definition
              ! -----------------------------------------

              ! Now, which type of BC is to be created?
              select case (rbcAssemblyData%cbctype)
            
              ! --------------------------------------------------
              ! Other boundary conditions
              ! --------------------------------------------------
              case (0,1,7,8)
              
              ! --------------------------------------------------
              ! Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (9)

                ! Prepare the assembly
                !    
                ! Get the definition of the boundary condition.
                ! Read the line again, get the expressions for the data field
                read(sparams,*) sbdex1,sbdex2
              
                ! For any string <> "", create the appropriate Dirichlet boundary
                ! condition and add it to the list of boundary conditions.
                
                ! X-velocity
                if (sbdex1 .ne. "") then
                  
                  ! Get the expression information
                  call struc_getBdExprInfo (roptcBDC,sbdex1,&
                      rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                      rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                  
                  ! Assemble
                  rbcAssemblyData%icomponent = 1
                  call user_initCollectForVecAssembly (&
                      rglobalData,rbcAssemblyData%iid,rbcAssemblyData%icomponent,dtime,&
                      rusercollection)
                  
                  rlinform%itermCount = 1
                  rlinform%Idescriptors(1) = DER_FUNC
                  call linf_buildVectorScalarBdr2D (rlinform, CUB_G4_1D, .false., &
                      rrhs%RvectorBlock(rbcAssemblyData%icomponent),cc_getNeumannBdData,&
                      rbcAssemblyData%rboundaryRegion,rcollection)
                      
                  call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                
                end if

                ! Y-velocity
                if (sbdex1 .ne. "") then
                  
                  ! Get the expression information
                  call struc_getBdExprInfo (roptcBDC,sbdex1,&
                      rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                      rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                  
                  ! Assemble
                  rbcAssemblyData%icomponent = 2
                  call user_initCollectForVecAssembly (&
                      rglobalData,rbcAssemblyData%iid,rbcAssemblyData%icomponent,dtime,&
                      rusercollection)
                  
                  rlinform%itermCount = 1
                  rlinform%Idescriptors(1) = DER_FUNC
                  call linf_buildVectorScalarBdr2D (rlinform, CUB_G4_1D, .false., &
                      rrhs%RvectorBlock(rbcAssemblyData%icomponent),cc_getNeumannBdData,&
                      rbcAssemblyData%rboundaryRegion,rcollection)
                      
                  call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                
                end if

              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"sbc_assembleInhomNeumannRHS")
                call sys_halt()
              end select ! cbctyp
              
            end if ! bautomatic
          
          ! ----------------------
          ! Heat equation
          ! ----------------------
          case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
          
            ! Do we have 'automatic' boundary conditions?
            ! We only have to assemble something for non-automatic BCs!
            if (.not. bautomatic .and. .not. bautomaticSegment) then
            
              ! -----------------------------------------
              ! Boundary conditions by definition
              ! -----------------------------------------

              ! Now, which type of BC is to be created?
              select case (rbcAssemblyData%cbctype)
            
              ! --------------------------------------------------
              ! Other boundary conditions
              ! --------------------------------------------------
              case (0,1,7,8)
                
              ! --------------------------------------------------
              ! Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (9)

                ! Prepare the assembly
                !    
                ! Get the definition of the boundary condition.
                ! Read the line again, get the expressions for the data field
                read(sparams,*) sbdex1
              
                ! For any string <> "", create the appropriate Dirichlet boundary
                ! condition and add it to the list of boundary conditions.
                
                if (sbdex1 .ne. "") then
                  
                  ! Get the expression information
                  call struc_getBdExprInfo (roptcBDC,sbdex1,&
                      rbcAssemblyData%cexprType,rbcAssemblyData%iid,rbcAssemblyData%ivalue,&
                      rbcAssemblyData%dvalue,rbcAssemblyData%sexprName)
                  
                  ! Assemble
                  rbcAssemblyData%icomponent = 1
                  call user_initCollectForVecAssembly (&
                      rglobalData,rbcAssemblyData%iid,rbcAssemblyData%icomponent,dtime,&
                      rusercollection)
                  
                  rlinform%itermCount = 1
                  rlinform%Idescriptors(1) = DER_FUNC
                  call linf_buildVectorScalarBdr2D (rlinform, CUB_G4_1D, .false., &
                      rrhs%RvectorBlock(rbcAssemblyData%icomponent),cc_getNeumannBdData,&
                      rbcAssemblyData%rboundaryRegion,rcollection)
                      
                  call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                
                end if

              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"sbc_assembleInhomNeumannRHS")
                call sys_halt()
              end select ! cbctyp
              
            end if ! bautomatic

          case default
            
            call output_line ("Unknown equation", &
                OU_CLASS_ERROR,OU_MODE_STD,"sbc_assembleInhomNeumannRHS")
            call sys_halt()
            
          end select ! equation
          
          ! Move on to the next parameter value
          dpar1 = dpar2
          
        end if
                                          
      end do ! isegment
      
    end do ! ibct

    call collct_done (rusercollection)

  end subroutine

end module
