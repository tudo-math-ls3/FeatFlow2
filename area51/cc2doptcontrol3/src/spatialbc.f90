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

!</types>

  public :: sbc_assembleBDconditions
  public :: sbc_releaseBdRegionList
  !public :: sbc_implementDirichletBC
  
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

  subroutine sbc_assembleBDconditions (roptcBDC,roptcBDCSpace,dtime,cequation,&
      copType,casmFlags,rspaceDiscr,rglobalData,rvectorControl)

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
    
    integer :: ibdComponent, isegment, i, iprimal, iintervalEnds
    character(LEN=PARLST_LENLINEBUF) :: sstr,sexpr,svalue,sbdex1,sbdex2,sbdcomp
    real(DP) :: dpar1, dpar2
    integer :: ctype, ivalue, iid, ibctyp
    real(DP) :: dvalue
    logical :: bneedsParams
    type(t_collection), target :: rcollection, ruserCollection
    type(t_boundary), pointer :: p_rboundary
    type(t_p_optcBDC) :: r_t_p_optcBDC

    ! A boundary region defining where the boundary conditionis applied
    type(t_boundaryRegion), target :: rboundaryRegion

    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection,p_rbdcond,p_rbdcondPrimal

    ! Fetch some parameters
    p_rboundary => rspaceDiscr%p_rboundary
    
    ! Save a pointer to roptBDC ti r_t_p_optcBDC for later use.
    r_t_p_optcBDC%p_roptcBDC => roptcBDC

    ! Get the expression/bc sections from the bondary condition block
    call parlst_querysection(roptcBDC%p_rparList, &
        roptcBDC%ssectionBdExpr, p_rsection)
        
    ! Get the section defining the primal or dual boundary conditions
    select case (copType)
    case (OPTP_PRIMAL)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcond)

    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)

      ! Primal boundary conditions, linearised eqn.
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrimLin, p_rbdcond)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcondPrimal)

    case (OPTP_DUAL)
    
      ! Dual boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondDual, p_rbdcond)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcondPrimal)

    case (OPTP_DUALLIN,OPTP_DUALLIN_SIMPLE)
    
      ! Dual boundary conditions, linearised equation
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondDualLin, p_rbdcond)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcondPrimal)

    end select
        
    ! Initialise the user-defined assembly
    call collct_init (rusercollection)

    ! Loop through all boundary components we have.
    do ibdComponent = 1,boundary_igetNBoundComp(p_rboundary)

      ! Parse the parameter "bdComponentX"
      write (sexpr,"(I10)") ibdComponent
      sbdcomp = "bdComponent" // adjustl(sexpr)
      
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      i = parlst_queryvalue (p_rbdcond, sbdcomp)

      if (i .ne. 0) then
      
        ! get the corresponding index in the "primal" dection
        if (copType .ne. OPTP_PRIMAL) then
          iprimal = parlst_queryvalue (p_rbdcondPrimal, sbdcomp)
        end if
      
        ! Parameter exists. Get the values in there.
        do isegment = 1,parlst_querysubstrings (p_rbdcond, sbdcomp)
          
          call parlst_getvalue_string (p_rbdcond, i, sstr, isubstring=isegment)
          
          ! Read the segment parameters
          read(sstr,*) dpar2,iintervalEnds,ibctyp
          
          ! Form a boundary condition segment that covers that boundary part
          if (dpar2 .ge. dpar1) then
            
            rboundaryRegion%dminParam = dpar1
            rboundaryRegion%dmaxParam = dpar2
            rboundaryRegion%iboundCompIdx = ibdComponent
            rboundaryRegion%dmaxParamBC = &
              boundary_dgetMaxParVal(p_rboundary, ibdComponent)
            rboundaryRegion%iproperties = iintervalEnds
            
            ! Type of equation?
            select case (cequation)
            
            ! ----------------------
            ! Stokes / Navier-Stokes
            ! ----------------------
            case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
            
              ! Now, which type of BC is to be created?
              select case (ibctyp)
              
              ! ------------------------------------------------------------------
              ! Automatic boundary conditions. Dual and linearised equations only.
              ! ------------------------------------------------------------------
              case (-1)
                if (copType .eq. OPTP_PRIMAL) then
                  call output_line (&
                      "Automatic boundary conditions only allowed.", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                end if
                
                ! The rule is:
                !  * Primal Neumann = Dual Neumann,
                !  * Primal Dirichlet = Dual Dirichlet-0

                ! Get the boundary conditions of the primal equation
                call parlst_getvalue_string (&
                    p_rbdcondPrimal, iprimal, sstr, isubstring=isegment)
                
                ! Read the segment parameters
                read(sstr,*) dpar2,iintervalEnds,ibctyp
                
                select case (ibctyp)
                
                ! --------------------------------------------------
                ! Neumann boundary conditions
                ! --------------------------------------------------
                case (0)
                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                    ! Add the bondary region to the Neumann boundary regions.
                    call sbc_addBoundaryRegion(&
                        rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                  end if
                  
                ! --------------------------------------------------
                ! Dirichlet boundary conditions
                ! --------------------------------------------------
                case (1)
                  if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
                  
                    ! Add the bondary region to the Dirichlet boundary regions.
                    call sbc_addBoundaryRegion(&
                        rboundaryRegion,roptcBDCSpace%rdirichletBoundary)

                    if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                      ! Impose Dirichlet-0 boundary conditions in the dual equation
                    
                      rcollection%IquickAccess(2) = BDC_VALDOUBLE
                      rcollection%IquickAccess(3) = 0
                      rcollection%IquickAccess(4) = 0
                      rcollection%IquickAccess(5) = 0
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = 0.0_DP
                      rcollection%SquickAccess(1) = ""
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 0
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))

                      ! X-velocity
                      rcollection%IquickAccess(1) = 1
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,0,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)

                      ! Y-velocity
                      rcollection%IquickAccess(1) = 2
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,0,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                      
                    end if
                    
                  end if
                  
                ! --------------------------------------------------
                ! L2 Dirichlet boundary control
                ! --------------------------------------------------
                case (7)
                  if (iand(casmFlags,SBC_DIRICHLETBCC) .ne. 0) then
                    ! Add the bondary region to the Dirichlet boundary regions.
                    call sbc_addBoundaryRegion(&
                        rboundaryRegion,roptcBDCSpace%rdirichletControlBoundary)
                        
                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then

                      ! Impose Dirichlet-0 boundary conditions 
                      ! plus boundary control
                    
                      rcollection%IquickAccess(2) = BDC_VALDOUBLE
                      rcollection%IquickAccess(3) = 0
                      rcollection%IquickAccess(4) = 0
                      rcollection%IquickAccess(5) = 0
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = 0.0_DP
                      rcollection%SquickAccess(1) = ""
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 0
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))

                      if ((copType .eq. OPTP_PRIMALLIN) .or. (copType .eq. OPTP_PRIMALLIN_SIMPLE)) then
                        ! Prescribed boundary conditions
                        rcollection%IquickAccess(6) = 1
                        rcollection%p_rvectorQuickAccess1 => rvectorControl
                      end if

                      ! X-velocity
                      rcollection%IquickAccess(1) = 1
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,0,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)

                      ! Y-velocity
                      rcollection%IquickAccess(1) = 2
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,0,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                    
                    end if
                  
                  end if
                  
                case (8)
                  ! Nothing to do here
                  
                case default
                  call output_line ("Cannot set up automatic boundary conditions", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                  
                end select
                
              ! --------------------------------------------------
              ! Homogeneous/Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (0,8)
                if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                  ! Add the bondary region to the Neumann boundary regions.
                  call sbc_addBoundaryRegion(&
                      rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                end if
              
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)

                if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then

                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rboundaryRegion,roptcBDCSpace%rdirichletBoundary)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                  
                    ! Simple Dirichlet boundary.
                    ! Get the definition of the boundary condition.
                    ! Read the line again, get the expressions for X- and Y-velocity
                    read(sstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
                  
                    ! For any string <> "", create the appropriate Dirichlet boundary
                    ! condition and add it to the list of boundary conditions.
                    !
                    ! The IquickAccess array is set up as follows:
                    !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
                    !  IquickAccess(2) = expression type
                    !  IquickAccess(3) = iid
                    !  IquickAccess(4) = ivalue
                    !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
                    !  IquickAccess(6) = type of boundary control to impose.
                    !                    (=0: none, =1: Dirichlet L2)
                    !  IquickAccess(7:...) = The binary content of r_t_p_optcBDC.
                    !
                    ! The DquickAccess array is set up as follows:
                    !  DquickAccess(1) = dtime
                    !  DquickAccess(2) = dvalue
                    !
                    ! The SquickAccess array is set up as follows:
                    !  SquickAccess(1) = Name of the expression
                    
                    if (sbdex1 .ne. "") then
                      
                      ! X-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex1,&
                          ctype,iid,ivalue,dvalue,svalue,bneedsParams)
                          
                      rcollection%IquickAccess(1) = 1
                      rcollection%IquickAccess(2) = ctype
                      rcollection%IquickAccess(3) = iid
                      rcollection%IquickAccess(4) = ivalue
                      rcollection%IquickAccess(5) = 0
                      if (bneedsParams) rcollection%IquickAccess(5) = 1
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = dvalue
                      rcollection%SquickAccess(1) = svalue
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 0
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                      
                    end if
                        
                    if (sbdex2 .ne. "") then
                    
                      ! Y-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex2,&
                          ctype,iid,ivalue,dvalue,svalue,bneedsParams)
                          
                      rcollection%IquickAccess(1) = 2
                      rcollection%IquickAccess(2) = ctype
                      rcollection%IquickAccess(3) = iid
                      rcollection%IquickAccess(4) = ivalue
                      rcollection%IquickAccess(5) = 0
                      if (bneedsParams) rcollection%IquickAccess(5) = 1
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = dvalue
                      rcollection%SquickAccess(1) = svalue
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 0
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                      
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
                      rboundaryRegion,roptcBDCSpace%rdirichletControlBoundary)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                    ! Get the definition of the boundary condition.
                    ! Read the line again, get the expressions for X- and Y-velocity
                    read(sstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
                  
                    ! For any string <> "", create the appropriate Dirichlet boundary
                    ! condition and add it to the list of boundary conditions.
                    !
                    ! The IquickAccess array is set up as follows:
                    !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
                    !  IquickAccess(2) = expression type
                    !  IquickAccess(3) = iid
                    !  IquickAccess(4) = ivalue
                    !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
                    !  IquickAccess(6) = type of boundary control to impose.
                    !                    (=0: none, =1: Dirichlet L2)
                    !  IquickAccess(7:...) = The binary content of r_t_p_optcBDC.
                    !
                    ! The DquickAccess array is set up as follows:
                    !  DquickAccess(1) = dtime
                    !  DquickAccess(2) = dvalue
                    !
                    ! The SquickAccess array is set up as follows:
                    !  SquickAccess(1) = Name of the expression
                    !
                    ! The Quickvectors array is set up as follows:
                    !  p_rvectorQuickAccess1 => rvectorControl
                    
                    if (sbdex1 .ne. "") then
                      
                      ! X-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex1,&
                          ctype,iid,ivalue,dvalue,svalue,bneedsParams)
                          
                      rcollection%IquickAccess(1) = 1
                      rcollection%IquickAccess(2) = ctype
                      rcollection%IquickAccess(3) = iid
                      rcollection%IquickAccess(4) = ivalue
                      rcollection%IquickAccess(5) = 0
                      if (bneedsParams) rcollection%IquickAccess(5) = 1
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = dvalue
                      rcollection%SquickAccess(1) = svalue
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 1
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))
                      rcollection%p_rvectorQuickAccess1 => rvectorControl
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                      
                    end if
                        
                    if (sbdex2 .ne. "") then
                    
                      ! Y-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex2,&
                          ctype,iid,ivalue,dvalue,svalue,bneedsParams)
                          
                      rcollection%IquickAccess(1) = 2
                      rcollection%IquickAccess(2) = ctype
                      rcollection%IquickAccess(3) = iid
                      rcollection%IquickAccess(4) = ivalue
                      rcollection%IquickAccess(5) = 0
                      if (bneedsParams) rcollection%IquickAccess(5) = 1
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = dvalue
                      rcollection%SquickAccess(1) = svalue
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 1
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))
                      rcollection%p_rvectorQuickAccess1 => rvectorControl
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                      
                    end if
                    
                   end if ! SBC_DISCRETEBC

                end if

              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
              end select ! ibctyp

            ! ----------------------
            ! Heat equation
            ! ----------------------
            case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
            
              ! Now, which type of BC is to be created?
              select case (ibctyp)
              
              ! ------------------------------------------------------------------
              ! Automatic boundary conditions. Dual and linearised equations only.
              ! ------------------------------------------------------------------
              case (-1)
                if (copType .eq. OPTP_PRIMAL) then
                  call output_line (&
                      "Automatic boundary conditions only allowed.", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                end if
                
                ! The rule is:
                !  * Primal Neumann = Dual Neumann,
                !  * Primal Dirichlet = Dual Dirichlet-0

                ! Get the boundary conditions of the primal equation
                call parlst_getvalue_string (&
                    p_rbdcondPrimal, iprimal, sstr, isubstring=isegment)
                
                ! Read the segment parameters
                read(sstr,*) dpar2,iintervalEnds,ibctyp
                
                select case (ibctyp)
                
                ! --------------------------------------------------
                ! Neumann boundary conditions
                ! --------------------------------------------------
                case (0)
                
                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                    ! Add the bondary region to the Neumann boundary regions.
                    call sbc_addBoundaryRegion(&
                        rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                  end if
                  
                ! --------------------------------------------------
                ! Dirichlet boundary conditions
                ! --------------------------------------------------
                case (1)
                
                  if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
                    
                    ! Add the bondary region to the Dirichlet boundary regions.
                    call sbc_addBoundaryRegion(&
                        rboundaryRegion,roptcBDCSpace%rdirichletBoundary)
                        
                    if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                    
                      ! Impose Dirichlet-0 boundary conditions in the dual equation
                    
                      rcollection%IquickAccess(2) = BDC_VALDOUBLE
                      rcollection%IquickAccess(3) = 0
                      rcollection%IquickAccess(4) = 0
                      rcollection%IquickAccess(5) = 0
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = 0.0_DP
                      rcollection%SquickAccess(1) = ""
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 0
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))

                      ! Solution
                      rcollection%IquickAccess(1) = 1
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,0,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                      
                    end if
                  
                  end if
                  
                case (8)
                  ! Nothing to do here

                case default
                  call output_line ("Cannot set up automatic boundary conditions", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                  
                end select
                
              ! --------------------------------------------------
              ! Homogeneous/Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (0,8)
              
                if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                  ! Add the bondary region to the Neumann boundary regions
                  ! if there is a structure present.
                  call sbc_addBoundaryRegion(&
                      rboundaryRegion,roptcBDCSpace%rneumannBoundary)
                end if
              
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)

                if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then

                  ! Add the bondary region to the Dirichlet boundary regions.
                  call sbc_addBoundaryRegion(&
                      rboundaryRegion,roptcBDCSpace%rdirichletBoundary)

                  if (iand(casmFlags,SBC_DISCRETEBC) .ne. 0) then
                  
                    ! Simple Dirichlet boundary.
                    ! Get the definition of the boundary condition.
                    ! Read the line again, get the expressions for the data field
                    read(sstr,*) dvalue,iintervalEnds,ibctyp,sbdex1
                  
                    ! For any string <> "", create the appropriate Dirichlet boundary
                    ! condition and add it to the list of boundary conditions.
                    !
                    ! The IquickAccess array is set up as follows:
                    !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
                    !  IquickAccess(2) = expression type
                    !  IquickAccess(3) = iid
                    !  IquickAccess(4) = ivalue
                    !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
                    !  IquickAccess(7:...) = The binary content of r_t_p_optcBDC.
                    !
                    ! The DquickAccess array is set up as follows:
                    !  DquickAccess(1) = dtime
                    !  DquickAccess(2) = dvalue
                    !
                    ! The SquickAccess array is set up as follows:
                    !  SquickAccess(1) = Name of the expression
                    
                    if (sbdex1 .ne. "") then
                      
                      ! X-velocity
                      !
                      ! Get the expression information
                      call struc_getBdExprInfo (roptcBDC,sbdex1,&
                          ctype,iid,ivalue,dvalue,svalue,bneedsParams)
                          
                      rcollection%IquickAccess(1) = 1
                      rcollection%IquickAccess(2) = ctype
                      rcollection%IquickAccess(3) = iid
                      rcollection%IquickAccess(4) = ivalue
                      rcollection%IquickAccess(5) = 0
                      if (bneedsParams) rcollection%IquickAccess(5) = 1
                      rcollection%DquickAccess(1) = dtime
                      rcollection%DquickAccess(2) = dvalue
                      rcollection%SquickAccess(1) = svalue
                      rcollection%p_rnextCollection => ruserCollection
                      rcollection%IquickAccess(6) = 0
                      rcollection%IquickAccess(7:) = &
                          transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                        size(rcollection%IquickAccess(7:)))
                      
                      call user_initCollectForVecAssembly (&
                          rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                      
                      ! Assemble the BCs.
                      call bcasm_newDirichletBConRealBD (&
                          rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,&
                          roptcBDCSpace%rdiscreteBC,&
                          cc_getDirBCNavSt2D,rcollection)
                          
                      call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                      
                    end if
                    
                  end if ! SBC_DISCRETEBC
                      
                end if
                
              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
              end select ! ibctyp

            case default
              
              call output_line ("Unknown equation", &
                  OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
              call sys_halt()
              
            end select ! equation
            
            ! Move on to the next parameter value
            dpar1 = dpar2
            
          end if
                                            
        end do ! isegment
      
      end if
      
    end do ! ibdComponent

    call collct_done (rusercollection)


!    ! local variables
!    integer :: i,ityp,ivalue,ibdComponent,isegment,iintervalEnds,iprimaldual
!    integer :: ibctyp,icount,iexptyp,iid
!    integer, dimension(2) :: IminIndex,imaxIndex
!    integer, dimension(NDIM2D) :: IvelComp
!    real(DP) :: dvalue,dpar1,dpar2
!    character(LEN=PARLST_LENLINEBUF) :: cstr,cexpr,sbdex1,sbdex2,sbdex3,sbdex4
!    character(LEN=PARLST_MLNAME) :: cname
!    type(t_collection) :: rlocalCollection
!    type(t_boundary), pointer :: p_rboundary
!    type(t_linearForm) :: rlinformRhs
!    
!    ! A local collection we use for storing named parameters during the
!    ! parsing process.
!    type(t_collection) :: rcoll
!    
!    ! A collection structure passed to callback routines.
!    type(t_collection), target :: rcallbackcollection
!    
!    ! Triangulation on currently highest level.
!    type(t_triangulation), pointer :: p_rtriangulation
!
!    ! A set of variables describing the analytic boundary conditions.
!    type(t_boundaryRegion), target :: rboundaryRegion
!    
!    ! A pointer to the section with the expressions and the boundary conditions
!    type(t_parlstSection), pointer :: p_rsection,p_rbdcond
!    
!    ! A compiled expression for evaluation at runtime
!    type(t_fparser), target :: rparser
!    
!    select case (roptcBDC%p_rphysics%cequation)
!    case (0,1)
!      ! Stokes, Navier-Stokes, 2D
!    
!      ! For implementing boundary conditions, we use a "filter technique with
!      ! discretised boundary conditions". This means, we first have to calculate
!      ! a discrete version of the analytic BC, which we can implement into the
!      ! solution/RHS vectors using the corresponding filter.
!      !
!      ! At first, we need the analytic description of the boundary conditions.
!      ! Initialise a structure for boundary conditions, which accepts this,
!      ! on the heap.
!      !
!      ! We first set up the boundary conditions for the X-velocity, then those
!      ! of the Y-velocity.
!      !
!      ! Get the expression/bc sections from the bondary condition block
!      call parlst_querysection(roptcBDC%p_rparamListBDC, &
!          roptcBDC%ssectionBdExpressions, p_rsection)
!      call parlst_querysection(roptcBDC%p_rparamListBDC, &
!          roptcBDC%ssectionBdConditions, p_rbdcond)
!      
!      ! For intermediate storing of expression types, we use a local collection
!      call collct_init (rcoll)
!      
!      ! Add a section to the collection that accepts the boundary expressions
!      call collct_addsection (rcoll, roptcBDC%ssectionBdExpressions)
!      
!      ! Create a parser structure for as many expressions as configured
!      call fparser_create (rparser,&
!          parlst_querysubstrings (p_rsection, "bdExpressions"))
!      
!      ! Add the parser to the collection
!      call collct_setvalue_pars (rcoll, BDC_BDPARSER, rparser, &
!                                  .true., 0, SEC_SBDEXPRESSIONS)
!      
!      ! Add the boundary expressions to the collection into the
!      ! specified section.
!      do i=1,parlst_querysubstrings (p_rsection, "bdExpressions")
!      
!        call parlst_getvalue_string (p_rsection, "bdExpressions", cstr, "", i)
!        
!        ! Get the type and decide on the identifier how to save the expression.
!        read(cstr,*) cname,ityp
!        
!        select case (ityp)
!        case (BDC_USERDEFID)
!          ! Id of a hardcoded expression realised in the callback routines
!          read(cstr,*) cname,ityp,ivalue
!          call collct_setvalue_int (rcoll, cname, ivalue, .true., &
!                                    0, SEC_SBDEXPRESSIONS)
!        case (BDC_USERDEF)
!          ! Name of a hardcoded expression realised in the callback routines
!          read(cstr,*) cname,ityp,cexpr
!          call collct_setvalue_string (rcoll, cname, cexpr, .true., &
!                                      0, SEC_SBDEXPRESSIONS)
!
!        case (BDC_EXPRESSION)
!          ! General expression; not implemented yet
!          read(cstr,*) cname,ityp,cexpr
!          
!          ! Compile the expression; the expression gets number i
!          call fparser_parseFunction (rparser, i, cexpr, SEC_EXPRVARIABLES)
!          
!          ! Add the number of the function in the parser object to
!          ! the collection with the name of the expression
!          call collct_setvalue_int (rcoll, cname, i, .true., &
!                                    0, SEC_SBDEXPRESSIONS)
!          
!        case (BDC_VALDOUBLE)
!          ! Real-value
!          read(cstr,*) cname,ityp,dvalue
!          call collct_setvalue_real (rcoll, cname, dvalue, .true., &
!                                    0, SEC_SBDEXPRESSIONS)
!        case (BDC_VALINT)
!          ! Integer-value
!          read(cstr,*) cname,ityp,ivalue
!          call collct_setvalue_int (rcoll, cname, ivalue, .true., &
!                                    0, SEC_SBDEXPRESSIONS)
!                  
!        case (BDC_VALPARPROFILE)
!          ! Parabolic profile with specified maximum velocity
!          read(cstr,*) cname,ityp,dvalue
!          call collct_setvalue_real (rcoll, cname, dvalue, .true., &
!                                    0, SEC_SBDEXPRESSIONS)
!                                     
!        case default
!          call output_line ("Expressions not implemented!", &
!              OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
!          call sys_halt()
!
!        end select
!        
!        ! Put the type of the expression to the temporary collection section
!        call collct_setvalue_int (rcoll, cname, ityp, .true.)
!        
!      end do
!      
!      ! Get the triangulation on the highest level
!      p_rtriangulation => rspaceDiscr%p_rtriangulation
!      p_rboundary => rspaceDiscr%p_rboundary
!      
!      ! Put some information to the quick access arrays for access
!      ! in the callback routine.
!      call collct_init(rcallbackcollection)
!
!      rcoll%Iquickaccess(1) = 1 ! time-dependent
!      rcoll%Dquickaccess(2) = rtimeDiscr%dtimeInit
!      rcoll%Dquickaccess(3) = rtimeDiscr%dtimeMax
!      
!      rcoll%p_rnextCollection => rcallbackcollection
!      
!      ! Now to the actual boundary conditions.
!      ! First primal (iprimaldual=1) then dual (iprimaldual=2)
!      do iprimaldual = 1,2
!        
!        ! DquickAccess(4:) is reserved for BC specific information.
!        !
!        ! Put a link to a user defined collection into that local collection.
!        ! That allows us to access it or to pass it to user defined callback
!        ! functions.
!
!        ! Loop through all boundary components we have.
!        do ibdComponent = 1,boundary_igetNBoundComp(p_rboundary)
!
!          ! Parse the parameter "bdComponentX"
!          write (cexpr,"(I10)") ibdComponent
!          cstr = "bdComponent" // adjustl(cexpr)
!          
!          ! We start at parameter value 0.0.
!          dpar1 = 0.0_DP
!          
!          i = parlst_queryvalue (p_rbdcond, cstr)
!          if (i .ne. 0) then
!            ! Parameter exists. Get the values in there.
!            do isegment = 1,parlst_querysubstrings (p_rbdcond, cstr)
!              
!              call parlst_getvalue_string (p_rbdcond, i, cstr, isubstring=isegment)
!              ! Read the segment parameters
!              read(cstr,*) dpar2,iintervalEnds,ibctyp
!              
!              ! Form a boundary condition segment that covers that boundary part
!              if (dpar2 .ge. dpar1) then
!                
!                rboundaryRegion%dminParam = dpar1
!                rboundaryRegion%dmaxParam = dpar2
!                rboundaryRegion%iboundCompIdx = ibdComponent
!                rboundaryRegion%dmaxParamBC = &
!                  boundary_dgetMaxParVal(p_rboundary, ibdComponent)
!                rboundaryRegion%iproperties = iintervalEnds
!                
!                ! Now, which type of BC is to be created?
!                select case (ibctyp)
!                
!                case (0)
!                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
!                    if (present(rneumannBoundary)) then
!                      ! Add the bondary region to the Neumann boundary regions
!                      ! if there is a structure present.
!                      call sbc_addBoundaryRegion(&
!                          rboundaryRegion,rneumannBoundary,iprimaldual)
!                    end if
!                  end if
!                
!                case (1)
!                  if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
!                    ! Simple Dirichlet boundary.
!                    ! Prescribed primal velocity; dual velocity os always zero.
!                    ! Read the line again, get the expressions for X- and Y-velocity
!                    read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
!                    
!                    ! For any string <> "", create the appropriate Dirichlet boundary
!                    ! condition and add it to the list of boundary conditions.
!                    !
!                    ! The IquickAccess array is set up as follows:
!                    !  IquickAccess(1) = Type of boundary condition
!                    !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
!                    !  IquickAccess(3) = expression type
!                    !  IquickAccess(4) = expression identifier
!                    !
!                    ! The SquickAccess array is set up as follows:
!                    !  SquickAccess(1) = Name of the expression
!                    !
!                    rcoll%IquickAccess(1) = ibctyp
!                    
!                    if (iprimaldual .eq. 1) then
!                      ! Primal BCs
!                      if (cbctype .eq. CCSPACE_PRIMAL) then
!                      
!                        ! If the type is a double precision value, set the DquickAccess(4)
!                        ! to that value so it can quickly be accessed.
!                        if (sbdex1 .ne. "") then
!                          ! X-velocity
!                          !
!                          ! The 2nd element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 1
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex1)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex1
!                          
!                          ! Dquickaccess(3) / IquickAccess(3) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (&
!                              rglobalData,iid,1,dtime,rcallbackcollection)
!                          
!                          ! Assemble the BCs.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. "") then
!                        
!                          ! Y-velocity
!                          !
!                          ! The 1st element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 2
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex2)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex2
!                          
!                          ! Dquickaccess(4) / IquickAccess(4) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (&
!                              rglobalData,iid,2,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                      end if
!                      
!                    end if
!                    
!                    if (iprimaldual .eq. 2) then
!                    
!                      ! Dual BCs
!                      if (cbctype .eq. CCSPACE_DUAL) then
!
!                        ! Now the same thing again, this time separately for primal and dual
!                        ! variables.
!                        ! If a primal velocity is specified, Dirichlet-0-boundary conditions are
!                        ! assumed for the corresponding dual.
!                        if (sbdex1 .ne. "") then
!                        
!                          ! dual X-velocity if primal X-velocity exists
!                          !
!                          ! The 2nd element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 1
!                          
!                          ! Ditichlet-0 boundary
!                          iexptyp = BDC_VALDOUBLE
!                          rcoll%Dquickaccess(4) = 0.0_DP
!                          rcoll%IquickAccess(3) = iexptyp
!                          rcoll%SquickAccess(1) = ""
!                          iid = 0
!                          
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,4,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          ! If we only assemble dual BCs and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. "") then
!                        
!                          ! dual Y-velocity if primal Y-velocity exists
!                          !
!                          ! The 2nd element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 2
!                          
!                          ! Ditichlet-0 boundary
!                          iexptyp = BDC_VALDOUBLE
!                          rcoll%Dquickaccess(4) = 0.0_DP
!                          rcoll%IquickAccess(3) = iexptyp
!                          rcoll%SquickAccess(1) = ""
!                          iid = 0
!                          
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,5,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          ! If we only assemble dual BCs and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                      end if
!                      
!                    end if
!                    
!                  end if
!                  
!                case (4)
!                  if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
!                    ! Simple Dirichlet boundary, separate definition for all
!                    ! solution components.
!                    ! Read the line again, get the expressions for X- and Y-velocity
!                    read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2,sbdex3,sbdex4
!                    
!                    ! For any string <> "", create the appropriate Dirichlet boundary
!                    ! condition and add it to the list of boundary conditions.
!                    !
!                    ! The IquickAccess array is set up as follows:
!                    !  IquickAccess(1) = Type of boundary condition
!                    !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
!                    !  IquickAccess(3) = expression type
!                    !  IquickAccess(4) = expression identifier
!                    !
!                    ! The SquickAccess array is set up as follows:
!                    !  SquickAccess(1) = Name of the expression
!                    !
!                    rcoll%IquickAccess(1) = ibctyp
!                    
!                    if (iprimaldual .eq. 1) then
!                    
!                      ! Primal BCs
!                      if (cbctype .eq. CCSPACE_PRIMAL) then
!                        ! If the type is a double precision value, set the DquickAccess(4)
!                        ! to that value so it can quickly be accessed.
!                        if (sbdex1 .ne. "") then
!                          ! X-velocity
!                          !
!                          ! The 2nd element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 1
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex1)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex1
!                          
!                          ! Dquickaccess(3) / IquickAccess(3) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (&
!                              rglobalData,iid,1,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. "") then
!                        
!                          ! Y-velocity
!                          !
!                          ! The 1st element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 2
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex2)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex2
!                          
!                          ! Dquickaccess(4) / IquickAccess(4) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (&
!                              rglobalData,iid,2,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                      end if
!                      
!                    end if
!                    
!                    if (iprimaldual .eq. 2) then
!                    
!                      ! Dual BCs
!                      if (cbctype .eq. CCSPACE_DUAL) then
!                      
!                        ! Now the same thing again, this time separately for primal and dual
!                        ! variables.
!                        ! If a velocity is not specified, Dirichlet-0-boundary conditions are
!                        ! assumed.
!                        
!                        if (sbdex3 .ne. "") then
!                          ! X-velocity
!                          !
!                          ! The 2nd element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 1
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex3)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex3
!                          
!                          ! Dquickaccess(3) / IquickAccess(3) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,4,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          ! If we only assemble dual BCs and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                        if (sbdex4 .ne. "") then
!                        
!                          ! Y-velocity
!                          !
!                          ! The 1st element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 2
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex4)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex4
!                          
!                          ! Dquickaccess(4) / IquickAccess(4) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex4, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,5,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          ! If we only assemble dual BCs and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                      end if
!                      
!                    end if
!                    
!                  end if
!
!                case (5)
!                  ! Simple Dirichlet boundary for the primal equation.
!                  ! Read the line again, get the expressions for X- and Y-velocity
!                  read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
!                 
!                  ! For any string <> "", create the appropriate Dirichlet boundary
!                  ! condition and add it to the list of boundary conditions.
!                  !
!                  ! The IquickAccess array is set up as follows:
!                  !  IquickAccess(1) = Type of boundary condition
!                  !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
!                  !  IquickAccess(3) = expression type
!                  !  IquickAccess(4) = expression identifier
!                  !
!                  ! The SquickAccess array is set up as follows:
!                  !  SquickAccess(1) = Name of the expression
!                  !
!                  rcoll%IquickAccess(1) = ibctyp
!                  
!                  if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
!                    
!                    if (iprimaldual .eq. 1) then
!                    
!                      ! Primal BCs
!                      if (cbctype .eq. CCSPACE_PRIMAL) then
!                        ! If the type is a double precision value, set the DquickAccess(4)
!                        ! to that value so it can quickly be accessed.
!                        if (sbdex1 .ne. "") then
!                          ! X-velocity
!                          !
!                          ! The 2nd element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 1
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex1)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex1
!                          
!                          ! Dquickaccess(3) / IquickAccess(3) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (&
!                              rglobalData,iid,1,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. "") then
!                        
!                          ! Y-velocity
!                          !
!                          ! The 1st element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 2
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex2)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex2
!                          
!                          ! Dquickaccess(4) / IquickAccess(4) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (&
!                              rglobalData,iid,2,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                      end if
!                      
!                    end if
!                  
!                  end if
!                                 
!                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
!                    ! Neumann boundary for the dual.
!                    if (iprimaldual .eq. 2) then
!                      !call bcasm_getEdgesInBCregion (p_rtriangulation,p_rboundary,&
!                      !                              rboundaryRegion, &
!                      !                              IminIndex,ImaxIndex,icount)
!                      !if ((icount .gt. 0) .and. present(rneumannBoundary)) then
!                      if (present(rneumannBoundary)) then
!                        ! Add the bondary region to the Neumann boundary regions
!                        ! if there is a structure present.
!                        call sbc_addBoundaryRegion(&
!                            rboundaryRegion,rneumannBoundary,iprimaldual)
!                      end if
!                    end if
!                  end if
!                                    
!                case (6)
!                  ! Simple Dirichlet boundary for the dual equation
!                  ! Read the line again, get the expressions for X- and Y-velocity
!                  read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex3,sbdex4
!                  
!                  ! For any string <> "", create the appropriate Dirichlet boundary
!                  ! condition and add it to the list of boundary conditions.
!                  !
!                  ! The IquickAccess array is set up as follows:
!                  !  IquickAccess(1) = Type of boundary condition
!                  !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
!                  !  IquickAccess(3) = expression type
!                  !  IquickAccess(4) = expression identifier
!                  !
!                  ! The SquickAccess array is set up as follows:
!                  !  SquickAccess(1) = Name of the expression
!                  !
!                  rcoll%IquickAccess(1) = ibctyp
!                  
!                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
!                    ! Neumann boundary for the primal.
!                    if (iprimaldual .eq. 1) then
!                      !call bcasm_getEdgesInBCregion (p_rtriangulation,rboundary,&
!                      !                              rboundaryRegion, &
!                      !                              IminIndex,ImaxIndex,icount)
!                      !if ((icount .gt. 0) .and. present(rneumannBoundary)) then
!                      if (present(rneumannBoundary)) then
!                        ! Add the bondary region to the Neumann boundary regions
!                        ! if there is a structure present.
!                        call sbc_addBoundaryRegion(&
!                            rboundaryRegion,rneumannBoundary,iprimaldual)
!                      end if
!                    end if
!                  end if
!
!                  if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
!                    if (iprimaldual .eq. 2) then
!                    
!                      ! Dual BCs
!                      if (cbctype .eq. CCSPACE_DUAL) then
!                      
!                        ! Now the same thing again, this time separately for primal and dual
!                        ! variables.
!                        ! If a velocity is not specified, Dirichlet-0-boundary conditions are
!                        ! assumed.
!                        
!                        if (sbdex3 .ne. "") then
!                          ! X-velocity
!                          !
!                          ! The 2nd element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 1
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex3)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex3
!                          
!                          ! Dquickaccess(3) / IquickAccess(3) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll, sbdex3, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,4,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          ! If we only assemble dual BCs and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                        if (sbdex4 .ne. "") then
!                        
!                          ! Y-velocity
!                          !
!                          ! The 1st element in IquickAccess saves the component number.
!                          rcoll%IquickAccess(2) = 2
!                          
!                          ! IquickAccess(3) saves the type of the expression
!                          iexptyp = collct_getvalue_int (rcoll, sbdex4)
!                          rcoll%IquickAccess(3) = iexptyp
!                          
!                          ! The 1st element in the sting quick access array is
!                          ! the name of the expression to evaluate.
!                          rcoll%SquickAccess(1) = sbdex4
!                          
!                          ! Dquickaccess(4) / IquickAccess(4) saves information
!                          ! about the expression.
!                          select case (iexptyp)
!                          case (BDC_USERDEFID)
!                            iid = collct_getvalue_int (rcoll, sbdex4, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                            ! Constant or parabolic profile
!                            rcoll%Dquickaccess(4) = &
!                                collct_getvalue_real (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
!                          case (BDC_EXPRESSION)
!                            ! Expression. Write the identifier for the expression
!                            ! as itag into the boundary condition structure.
!                            rcoll%IquickAccess(4) = &
!                                collct_getvalue_int (rcoll,sbdex4, 0, roptcBDC%ssectionBdExpressions)
!                          end select
!                        
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,5,dtime,rcallbackcollection)
!
!                          ! Assemble the BCs.
!                          ! If we only assemble dual BCs and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                      end if
!                      
!                    end if
!                    
!                  end if
!                  
!                case (7)
!                  
!                  ! Dirichlet boundary control.
!                  !
!                  ! This is treated like Neumann boundary conditions here
!                  ! but leads to a different implementation.
!                  ! The corresponding region is added to the Neumann boundary
!                  ! structure as well as to the Dirichlet boundary control
!                  ! structure.
!                  ! Adding it to the Neumann boundary conditions structure
!                  ! is important, otherwise the components detect the problem
!                  ! as pure Dirichlet and start filtering the pressure... 
!!                  if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
!!                    if (present(rneumannBoundary)) then
!!                      call sbc_addBoundaryRegion(&
!!                          rboundaryRegion,rneumannBoundary,iprimaldual)
!!                    end if
!!                  end if
!
!                  ! For the primal equation, there is a special coupling
!                  ! to the dual equation.
!                  !
!                  ! For the dual equation, we have Dirichlet-0 boundary conditions here.
!
!                  if (iand(casmFlags,SBC_DIRICHLETBCC) .ne. 0) then
!                  
!                    if (iprimaldual .eq. 1) then
!                      if (present(rdirichletControlBoundary)) then
!                        ! Add the bondary region to the Dirichlet BCC boundary regions
!                        ! if there is a structure present.
!                        call sbc_addBoundaryRegion(&
!                            rboundaryRegion,rdirichletControlBoundary,iprimaldual)
!                      end if
!                      
!                      if (present(rvectorDirichletBCCRHS)) then
!                        ! Modify the RHS according to the Dirichlet boundary control
!                        ! conditions. Use a penalty approach to implement the Dirichlet
!                        ! values. A linear form has to be applied to the RHS vector
!                        ! in each component which incorporates the basic Dirichlet value.
!                        ! The control then acts relative to the incoroprated value.
!                        
!                        ! Prepare the linear form.
!                        rlinformRhs%itermCount = 1
!                        rlinformRhs%Idescriptors(1) = DER_FUNC2D
!                        rlinformRhs%Dcoefficients(1:rlinformRhs%itermCount)  = 1.0_DP
!                        
!                        ! Read the boundary data specification and prepare the callback routine
!                        read(cstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
!                        
!                        ! For any string <> "", create the appropriate Dirichlet boundary
!                        ! condition and add it to the list of boundary conditions.
!                        !
!                        ! The IquickAccess array is set up as follows:
!                        !  IquickAccess(1) = Type of boundary condition
!                        !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
!                        !  IquickAccess(3) = expression type
!                        !  IquickAccess(4) = expression identifier
!                        !
!                        ! The SquickAccess array is set up as follows:
!                        !  SquickAccess(1) = Name of the expression
!                        !
!                        rcoll%IquickAccess(1) = ibctyp
!                        
!                        ! The current boundary region is put to the collection
!                        call collct_setvalue_bdreg (rcoll, "BDREG", rboundaryRegion, .true.)
!                        
!                        ! Primal BCs
!                        if (cbctype .eq. CCSPACE_PRIMAL) then
!                        
!                          ! If the type is a double precision value, set the DquickAccess(4)
!                          ! to that value so it can quickly be accessed.
!                          if (sbdex1 .ne. "") then
!                            ! X-velocity
!                            !
!                            ! The 2nd element in IquickAccess saves the component number.
!                            rcoll%IquickAccess(2) = 1
!                            
!                            ! IquickAccess(3) saves the type of the expression
!                            iexptyp = collct_getvalue_int (rcoll, sbdex1)
!                            rcoll%IquickAccess(3) = iexptyp
!                            
!                            ! The 1st element in the sting quick access array is
!                            ! the name of the expression to evaluate.
!                            rcoll%SquickAccess(1) = sbdex1
!                            
!                            ! Dquickaccess(3) / IquickAccess(3) saves information
!                            ! about the expression.
!                            select case (iexptyp)
!                            case (BDC_USERDEFID)
!                              iid = collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                            case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                              ! Constant or parabolic profile
!                              rcoll%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%ddirichletBCPenalty * &
!                                  collct_getvalue_real (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                            case (BDC_EXPRESSION)
!                              ! Expression. Write the identifier for the expression
!                              ! as itag into the boundary condition structure.
!                              rcoll%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%ddirichletBCPenalty
!                              rcoll%IquickAccess(4) = &
!                                  collct_getvalue_int (rcoll, sbdex1, 0, roptcBDC%ssectionBdExpressions)
!                            end select
!                          
!                            rcoll%Dquickaccess(1) = dtime
!                            call user_initCollectForVecAssembly (&
!                                rglobalData,iid,1,dtime,rcallbackcollection)
!                            
!                            call linf_buildVectorScalarBdr2D (rlinformRhs, CUB_G4_1D, .false., &
!                                rvectorDirichletBCCRHS%RvectorBlock(1),&
!                                fcoeff_buildBCCRHS,rboundaryRegion, rcoll)
!                                
!                            call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!                                
!                          end if
!                          
!                          if (sbdex2 .ne. "") then
!                          
!                            ! Y-velocity
!                            !
!                            ! The 1st element in IquickAccess saves the component number.
!                            rcoll%IquickAccess(2) = 2
!                            
!                            ! IquickAccess(3) saves the type of the expression
!                            iexptyp = collct_getvalue_int (rcoll, sbdex2)
!                            rcoll%IquickAccess(3) = iexptyp
!                            
!                            ! The 1st element in the sting quick access array is
!                            ! the name of the expression to evaluate.
!                            rcoll%SquickAccess(1) = sbdex2
!                            
!                            ! Dquickaccess(4) / IquickAccess(4) saves information
!                            ! about the expression.
!                            select case (iexptyp)
!                            case (BDC_USERDEFID)
!                              iid = collct_getvalue_int (rcoll, sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                            case (BDC_VALDOUBLE,BDC_VALPARPROFILE)
!                              ! Constant or parabolic profile
!                              rcoll%Dquickaccess(4) = rglobalData%p_rsettingsOptControl%ddirichletBCPenalty * &
!                                  collct_getvalue_real (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                            case (BDC_EXPRESSION)
!                              ! Expression. Write the identifier for the expression
!                              ! as itag into the boundary condition structure.
!                              rcoll%IquickAccess(4) = &
!                                  collct_getvalue_int (rcoll,sbdex2, 0, roptcBDC%ssectionBdExpressions)
!                            end select
!                          
!                            rcoll%Dquickaccess(1) = dtime
!                            call user_initCollectForVecAssembly (&
!                                rglobalData,iid,2,dtime,rcallbackcollection)
!
!                            ! Assemble the BCs.
!                            call linf_buildVectorScalarBdr2D (rlinformRhs, CUB_G4_1D, .false., &
!                                rvectorDirichletBCCRHS%RvectorBlock(1),&
!                                fcoeff_buildBCCRHS,rboundaryRegion, rcoll)
!
!                            call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                          end if
!                          
!                        end if
!                        
!                        ! Delete the boundary region again
!                        call collct_deletevalue (rcoll,"BDREG")
!                        
!                      end if
!
!                    end if
!                    
!                  end if
!                  
!                  if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then
!                    if (iprimaldual .eq. 2) then
!                    
!                      ! Dual BCs
!                      if (cbctype .eq. CCSPACE_DUAL) then
!
!                        ! dual X-velocity if primal X-velocity exists
!                        
!                        ! Ditichlet-0 boundary
!                        iexptyp = BDC_VALDOUBLE
!                        rcoll%Dquickaccess(4) = 0.0_DP
!                        rcoll%IquickAccess(3) = iexptyp
!                        rcoll%SquickAccess(1) = ""
!                        iid = 0
!                        
!                        rcoll%Dquickaccess(1) = dtime
!                        call user_initCollectForVecAssembly (rglobalData,iid,4,dtime,rcallbackcollection)
!
!                        ! Assemble the BCs.
!                        ! The 2nd element in IquickAccess saves the component number.
!                        ! If we only assemble dual BCs and there are 3 solution components,
!                        ! we assume the vector to specify exactly the dual solution.
!
!                        rcoll%IquickAccess(2) = 1
!                        call bcasm_newDirichletBConRealBD (&
!                            rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                            cc_getBDconditionsNavSt2D,rcoll)
!
!                        call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!
!                        call user_initCollectForVecAssembly (rglobalData,iid,5,dtime,rcallbackcollection)
!
!                        rcoll%IquickAccess(2) = 2
!                        call bcasm_newDirichletBConRealBD (&
!                            rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                            cc_getBDconditionsNavSt2D,rcoll)
!                            
!                        call user_doneCollectForVecAssembly (rglobalData,rcallbackcollection)
!                            
!                      end if
!                      
!                    end if
!                  
!                  end if
!
!                case default
!                  call output_line ("Unknown boundary condition!", &
!                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
!                  call sys_halt()
!                end select
!                
!                ! Move on to the next parameter value
!                dpar1 = dpar2
!                
!              end if
!                                                
!            end do
!          
!          end if
!          
!        end do
!
!      end do
!
!      ! Release the parser object with all the expressions to be evaluated
!      ! on the boundary.
!      call fparser_release (rparser)
!
!      ! Remove the boundary value parser from the collection
!      call collct_deletevalue (rcoll, BDC_BDPARSER)
!      
!      ! Remove the boundary-expression section we added earlier,
!      ! with all their content.
!      call collct_deletesection(rcoll,roptcBDC%ssectionBdExpressions)
!
!      ! Remove the temporary collection from memory.
!      call collct_done(rcallbackcollection)
!      call collct_done(rcoll)
!
!    end select
        
  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine fcoeff_buildBCCRHS (rdiscretisation, rform, &
!                nelements, npointsPerElement, Dpoints, ibct, DpointPar, &
!                IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
!  
!  use basicgeometry
!  use collection
!  use domainintegration
!  use fsystem
!  use scalarpde
!  use spatialdiscretisation
!  use triangulation
!  
!!<description>
!  ! This subroutine is called during the vector assembly. It has to compute
!  ! the coefficients in front of the terms of the linear form.
!  !
!  ! The routine accepts a set of elements and a set of points on these
!  ! elements (cubature points) in in real coordinates.
!  ! According to the terms in the linear form, the routine has to compute
!  ! simultaneously for all these points and all the terms in the linear form
!  ! the corresponding coefficients in front of the terms.
!!</description>
!  
!!<input>
!  ! The discretisation structure that defines the basic shape of the
!  ! triangulation with references to the underlying triangulation,
!  ! analytic boundary boundary description etc.
!  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
!  
!  ! The linear form which is currently to be evaluated:
!  type(t_linearForm), intent(in) :: rform
!  
!  ! Number of elements, where the coefficients must be computed.
!  integer, intent(in) :: nelements
!  
!  ! Number of points per element, where the coefficients must be computed
!  integer, intent(in) :: npointsPerElement
!  
!  ! This is an array of all points on all the elements where coefficients
!  ! are needed.
!  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!  ! DIMENSION(dimension,npointsPerElement,nelements)
!  real(DP), dimension(:,:,:), intent(in) :: Dpoints
!
!  ! This is the number of the boundary component that contains the
!  ! points in Dpoint. All points are on the same boundary component.
!  integer, intent(in) :: ibct
!
!  ! For every point under consideration, this specifies the parameter
!  ! value of the point on the boundary component. The parameter value
!  ! is calculated in LENGTH PARAMETRISATION!
!  ! DIMENSION(npointsPerElement,nelements)
!  real(DP), dimension(:,:), intent(in) :: DpointPar
!
!  ! An array accepting the DOF`s on all elements in the test space.
!  ! DIMENSION(\#local DOF`s in test space,Number of elements)
!  integer, dimension(:,:), intent(in) :: IdofsTest
!
!  ! This is a t_domainIntSubset structure specifying more detailed information
!  ! about the element set that is currently being integrated.
!  ! It is usually used in more complex situations (e.g. nonlinear matrices).
!  type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!!</input>
!
!!<inputoutput>
!  ! Optional: A collection structure to provide additional
!  ! information to the coefficient routine.
!  type(t_collection), intent(inout), optional :: rcollection
!!</inputoutput>
!
!!<output>
!  ! A list of all coefficients in front of all terms in the linear form -
!  ! for all given points on all given elements.
!  !   DIMENSION(itermCount,npointsPerElement,nelements)
!  ! with itermCount the number of terms in the linear form.
!  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!!</output>
!  
!!</subroutine>
!  
!!    ! local variables
!!    integer, dimension(1) :: Icomponents
!!    real(DP), dimension(1) :: Dvalues
!!    integer :: i,j
!!    type(t_boundaryRegion), pointer :: p_rboundaryRegion
!!    real(DP) :: dpar, ddirichletBCPenalty
!!
!!    ! Penalty parameter
!!    ddirichletBCPenalty = rcollection%Dquickaccess(4)
!!
!!    ! Get the boundary region from the collection
!!    p_rboundaryRegion => collct_getvalue_bdreg (rcollection, "BDREG")
!!
!!    ! Loop through all points and calculate the values.
!!    do i=1,nelements
!!      do j=1,npointsPerElement
!!        ! Current component
!!        Icomponents(1) = rcollection%IquickAccess(2)
!!        
!!        ! Calculate the parameter value in 0-1 parametrisation
!!        dpar = boundary_convertParameter(rdiscretisation%p_rboundary, ibct, &
!!            DpointPar(j,i), BDR_PAR_LENGTH, BDR_PAR_01)
!!        
!!        call cc_getBDconditionsNavSt2D(Icomponents,rdiscretisation,&
!!            p_rboundaryRegion,0,DISCBC_NEEDFUNC,ibct,dpar,&
!!            Dvalues,rcollection)
!!            
!!        ! Value must be mutiplied by the penalty parameter
!!        Dcoefficients(1,j,i) = ddirichletBCPenalty*Dvalues(1)
!!      end do
!!    end do
!
!  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getDirBCNavSt2D (Icomponents,rspaceDiscr,rboundaryRegion,ielement, &
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

    integer :: icomponent,cexprtype,iid,ivalue,ccontrolType,icontrolcomp
    real(DP) :: dvalue
    logical :: bneedsParams
    real(DP), dimension(size(SEC_EXPRVARIABLES)) :: Dparams
    real(DP) :: d, dx, dy
    character(LEN=PARLST_LENLINEBUF) :: svalue
    type(t_p_optcBDC) :: r_t_p_optcBDC
    type(t_vectorBlock), pointer :: p_rcontrol
    type(t_settings_optcontrol), pointer :: p_roptControl
    type(t_settings_physics), pointer :: p_rphysics
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: iedge,iptpos,iv1,iv2
    real(DP) :: dpar1, dpar2, dpar, dval
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    
    real(DP) :: dtime

    ! The IquickAccess array is set up as follows:
    !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
    !  IquickAccess(2) = expression type
    !  IquickAccess(3) = iid
    !  IquickAccess(4) = ivalue
    !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
    !  IquickAccess(6) = Type of boundary control. =0: none, =1: Dirichlet L2
    !  IquickAccess(7) = Pointer to the t_optcBDC structure
    !
    ! The DquickAccess array is set up as follows:
    !  DquickAccess(1) = dvalue
    !  DquickAccess(2) = dtime
    !
    ! The SquickAccess array is set up as follows:
    !  SquickAccess(1) = Name of the expression
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    dtime        = rcollection%Dquickaccess(1)
    icomponent   = rcollection%Iquickaccess(1)
    cexprType    = rcollection%Iquickaccess(2)
    iid          = rcollection%Iquickaccess(3)
    ivalue       = rcollection%Iquickaccess(4)
    bneedsParams = rcollection%Iquickaccess(5) .ne. 0
    ccontrolType = rcollection%Iquickaccess(6)
    r_t_p_optcBDC   = transfer(rcollection%IquickAccess(7:),r_t_p_optcBDC)
    dvalue       = rcollection%Dquickaccess(2)
    p_rcontrol   => rcollection%p_rvectorQuickAccess1
    
    p_roptControl => r_t_p_optcBDC%p_roptcBDC%p_rsettingsOptControl
    p_rphysics => r_t_p_optcBDC%p_roptcBDC%p_rphysics
    
    if (bneedsParams) then
      ! Calculate the parameter array
      Dparams(:) = 0.0_DP
      
      call boundary_getCoords(rspaceDiscr%p_rboundary, &
                              rboundaryRegion%iboundCompIdx, &
                              dwhere, dx, dy)
      
      ! Get the local parameter value 0 <= d <= 1.
      ! Note that if dwhere < rboundaryRegion%dminParam, we have to add the maximum
      ! parameter value on the boundary to dpar as normally 0 <= dwhere < max.par.
      ! although 0 <= dminpar <= max.par
      !      and 0 <= dmaxpar <= max.par!
      d = dwhere
      if (d .lt. rboundaryRegion%dminParam) &
        d = d + boundary_dgetMaxParVal(rspaceDiscr%p_rboundary,&
                                        rboundaryRegion%iboundCompIdx)
      d = d - rboundaryRegion%dminParam
      
      ! Normalise to 0..1 using the length of the parameter region.
      ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
      d = d / (rboundaryRegion%dmaxParam - rboundaryRegion%dminParam)
      
      Dparams(1) = dx
      Dparams(2) = dy
      Dparams(3) = 0.0_DP ! does not exist, z-direction
      Dparams(4) = d
      Dparams(5) = dwhere
      Dparams(6) = boundary_convertParameter(rspaceDiscr%p_rboundary, &
                                          rboundaryRegion%iboundCompIdx, dwhere, &
                                          BDR_PAR_01, BDR_PAR_LENGTH)
      Dparams(7) = dtime
      
    end if

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
      icomponent = rcollection%IquickAccess(1)
      
      ! -> 1=X-velocity, 2=Y-velocity.
      
      ! Return zero Dirichlet boundary values for all situations by default.
      Dvalues(1) = 0.0_DP
      
      ! Type of expression?
      select case (cexprType)
      
      case (BDC_USERDEFID)
        ! Defined via analytic reference solution, identified by an id.
        ! The id is already specified in the collection.
        
        ! Call the user defined callback routine to evaluate the boudary value.
        ! No string tag is given which indicates that a reference function is to be
        ! evaluated.
        !
        ! As collection, we pass rcollection%p_rcollection here; this is a pointer
        ! to the application specific, global collection that may be of interest for
        ! callback routines. rcollection itself is actually a "local" collection,
        ! a "wrapper" for the rcollection of the application!
        call user_getBoundaryValues (&
            "",icomponent,rspaceDiscr,rboundaryRegion,&
            dwhere, Dvalues(1), rcollection%p_rnextCollection)
        
      case (BDC_USERDEF)
        ! This is a hardcoded, user-defined identifier.
        ! Call the user defined callback routine to evaluate the expression.
        !
        ! As collection, we pass rcollection%p_rcollection here; this is a pointer
        ! to the application specific, global collection that may be of interest for
        ! callback routines. rcollection itself is actually a "local" collection,
        ! a "wrapper" for the rcollection of the application!
        call user_getBoundaryValues (&
            svalue,icomponent,rspaceDiscr,rboundaryRegion,&
            dwhere, Dvalues(1), rcollection%p_rnextCollection)
      
      case default
        ! Default handling. Evaluate the expression.
        Dvalues(1) = struc_evalExpression (r_t_p_optcBDC%p_roptcBDC,cexprType,&
            iid,ivalue,dvalue,Dparams)
      
      end select
      
    end select
    
    select case (ccontrolType)
    
    ! ---------------------------------
    ! L2 Dirichlet boudary control
    ! ---------------------------------
    case (1)
    
      if (p_roptControl%dalphaL2BdC .ge. 0.0_DP) then
        if (.not. associated(p_rcontrol)) then
          call output_line ("Control not specified.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBCNavSt2D")
          call sys_halt()
        end if
        
        select case (p_rphysics%cequation)
        case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
        
          ! Where is the boundary control?
          icontrolcomp = icomponent
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
                  p_Ddata(iv1),p_Ddata(iedge+rspaceDiscr%p_rtriangulation%NVBD),p_Ddata(iv2),dval)
              
              Dvalues(1) = Dvalues(1) + dval

            end select          
          
          case default
            call output_line ("Unsupported element.", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBCNavSt2D")
            call sys_halt()
          end select
        
        
        case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
          call output_line ("Equation not supported.", &
              OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBCNavSt2D")
          call sys_halt()
        end select
        
      end if
    
    end select
  
  end subroutine

!  ! ***************************************************************************
!  
!!<subroutine>
!
!  subroutine sbc_implementDirichletBC (roptcBDC,dtime,cspace,&
!      rtimeDiscr,rspaceDiscr,rglobalData,&
!      rx,rb,rd,rdiscretebc)
!
!!<description>
!  ! Implements the boundary conditions at time dtime, for a solution
!  ! vector rx, a rhs vector rb and/or a defect vector rb.
!!</description>
!
!!<input>
!  ! Boundary conditions in the problem.
!  type(t_optcBDC), intent(inout) :: roptcBDC
!
!  ! Time where the BCs should be implemented.
!  real(DP), intent(IN) :: dtime
!
!  ! Type of solution space. CCSPACE_PRIMAL or CCSPACE_DUAL.
!  integer, intent(in) :: cspace
!
!  ! Time discretisation, dtime refers to.
!  type(t_timeDiscretisation), intent(in) :: rtimeDiscr
!
!  ! Space discretisation
!  type(t_blockDiscretisation), intent(in) :: rspaceDiscr
!
!  ! Global settings for callback routines.
!  type(t_globalData), intent(inout), target :: rglobalData
!!</input>
!
!!<inputoutput>
!  ! OPTIONAL: Solution vector
!  type(t_vectorBlock), intent(inout), optional :: rx
!
!  ! OPTIONAL: RHS vector
!  type(t_vectorBlock), intent(inout), optional :: rb
!
!  ! OPTIONAL: Defect vector
!  type(t_vectorBlock), intent(inout), optional :: rd
!
!  ! OPTIONAL: Boundary condition structure which receives the boudary
!  ! conditions. If present, it must have been initialised
!  ! by the caller. If present, it will be attached to the vectors.
!  type(t_discreteBC), intent(inout), optional :: rdiscreteBC
!!</inputoutput>
!
!!</subroutine>
!
!    ! Boundary condition structure which receives the boudary
!    ! conditions. 
!    type(t_discreteBC) :: rdiscreteBClocal
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Ddata
!  
!    ! DEBUG!!!
!    call lsysbl_getbase_double (rd,p_Ddata)
!
!    if (.not. present(rdiscreteBC)) then
!      ! Initialise the boundary conditions
!      call bcasm_initDiscreteBC(rdiscreteBClocal)
!
!      ! Assemble the BCs.
!      call sbc_assembleBDconditions (roptcBDC,dtime,cspace,&
!          rglobalData,SBC_DIRICHLETBC,&
!          rtimeDiscr,rspaceDiscr,rdiscreteBClocal)
!
!      ! Implement the boundary conditions into the vector(s).
!      if (present(rx)) then
!        call vecfil_discreteBCsol (rx,rdiscreteBClocal)
!      end if
!      
!      if (present(rb)) then
!        call vecfil_discreteBCrhs (rb,rdiscreteBClocal)
!      end if
!      
!      if (present(rd)) then
!        call vecfil_discreteBCdef (rd,rdiscreteBClocal)
!      end if
!
!      ! Release the BCs again.
!      call bcasm_releaseDiscreteBC(rdiscreteBClocal)
!
!    else
!
!      ! Assemble the BCs.
!      call sbc_assembleBDconditions (roptcBDC,dtime,cspace,&
!          rglobalData,SBC_DIRICHLETBC,&
!          rtimeDiscr,rspaceDiscr,rdiscreteBC)
!
!      ! Implement the boundary conditions into the vector(s).
!      ! Attach the structure.
!      if (present(rx)) then
!        call vecfil_discreteBCsol (rx,rdiscreteBC)
!        call lsysbl_assignDiscreteBC (rx,rdiscreteBC)
!      end if
!      
!      if (present(rb)) then
!        call vecfil_discreteBCrhs (rb,rdiscreteBC)
!        call lsysbl_assignDiscreteBC (rb,rdiscreteBC)
!      end if
!      
!      if (present(rd)) then
!        call vecfil_discreteBCdef (rd,rdiscreteBC)
!        call lsysbl_assignDiscreteBC (rd,rdiscreteBC)
!      end if
!
!    end if
!    
!  end subroutine

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
    call sbc_releaseBdRegionList (roptcBDCSpace%rdirichletControlBoundary)
    
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
    integer :: icomponent,cexprtype,iid,ivalue,ccontrolType,icontrolcomp
    real(DP) :: dvalue,dwhere
    logical :: bneedsParams
    real(DP), dimension(size(SEC_EXPRVARIABLES)) :: Dparams
    real(DP) :: d, dx, dy
    character(LEN=PARLST_LENLINEBUF) :: svalue
    type(t_p_optcBDC) :: r_t_p_optcBDC
    type(t_vectorBlock), pointer :: p_rcontrol
    type(t_settings_optcontrol), pointer :: p_roptControl
    type(t_settings_physics), pointer :: p_rphysics
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: iedge,iptpos,iv1,iv2
    real(DP) :: dpar1, dpar2, dpar, dval
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    type(t_boundaryRegion), pointer :: p_rboundaryRegion
    
    real(DP) :: dtime

    ! The IquickAccess array is set up as follows:
    !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
    !  IquickAccess(2) = expression type
    !  IquickAccess(3) = iid
    !  IquickAccess(4) = ivalue
    !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
    !  IquickAccess(6) = Type of boundary control. =0: none, =1: Dirichlet L2
    !  IquickAccess(7) = Pointer to the t_optcBDC structure
    !
    ! The DquickAccess array is set up as follows:
    !  DquickAccess(1) = dvalue
    !  DquickAccess(2) = dtime
    !
    ! The SquickAccess array is set up as follows:
    !  SquickAccess(1) = Name of the expression
    
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    dtime        = rcollection%Dquickaccess(1)
    icomponent   = rcollection%Iquickaccess(1)
    cexprType    = rcollection%Iquickaccess(2)
    iid          = rcollection%Iquickaccess(3)
    ivalue       = rcollection%Iquickaccess(4)
    bneedsParams = rcollection%Iquickaccess(5) .ne. 0
    ccontrolType = rcollection%Iquickaccess(6)
    r_t_p_optcBDC   = transfer(rcollection%IquickAccess(7:),r_t_p_optcBDC)
    dvalue       = rcollection%Dquickaccess(2)
    p_rcontrol   => rcollection%p_rvectorQuickAccess1
    p_rboundaryRegion => r_t_p_optcBDC%p_rboundaryRegion
    
    p_roptControl => r_t_p_optcBDC%p_roptcBDC%p_rsettingsOptControl
    p_rphysics => r_t_p_optcBDC%p_roptcBDC%p_rphysics
    
    do iel = 1,nelements
      do ipt = 1,npointsPerElement

        dwhere = DpointPar(ipt,iel)

        if (bneedsParams) then
          ! Calculate the parameter array
          Dparams(:) = 0.0_DP
          
          call boundary_getCoords(rdiscretisation%p_rboundary, ibct, &
                                  dwhere, dx, dy)
          
          ! Get the local parameter value 0 <= d <= 1.
          ! Note that if dwhere < rboundaryRegion%dminParam, we have to add the maximum
          ! parameter value on the boundary to dpar as normally 0 <= dwhere < max.par.
          ! although 0 <= dminpar <= max.par
          !      and 0 <= dmaxpar <= max.par!
          d = boundary_convertParameter(rdiscretisation%p_rboundary, &
                  ibct, dwhere, BDR_PAR_LENGTH, BDR_PAR_01)
          if (d .lt. p_rboundaryRegion%dminParam) &
            d = d + boundary_dgetMaxParVal(rdiscretisation%p_rboundary,ibct)
          d = d - p_rboundaryRegion%dminParam
          
          ! Normalise to 0..1 using the length of the parameter region.
          ! Necessary if a parabolic profile occurs in the inner of an edge e.g.
          d = d / (p_rboundaryRegion%dmaxParam - p_rboundaryRegion%dminParam)

          Dparams(1) = dx
          Dparams(2) = dy
          Dparams(3) = 0.0_DP ! does not exist, z-direction
          Dparams(4) = -1.0_DP
          Dparams(5) = boundary_convertParameter(rdiscretisation%p_rboundary, &
                         ibct, dwhere, BDR_PAR_LENGTH, BDR_PAR_01)
          Dparams(6) = dwhere
          Dparams(7) = dtime
          
        end if

        ! The IquickAccess array is set up as follows:
        !  IquickAccess(1) = Type of boundary condition
        !  IquickAccess(2) = component under consideration (1=x-vel, 2=y-vel,...)
        !  IquickAccess(3) = expression type
        !  IquickAccess(4) = expression identifier
        !
        ! Get from the current component of the PDE we are discretising:
        icomponent = rcollection%IquickAccess(1)
        
        ! -> 1=X-velocity, 2=Y-velocity.
        
        ! Return zero Dirichlet boundary values for all situations by default.
        Dcoefficients(1,ipt,iel) = 0.0_DP
        
        ! Type of expression?
        select case (cexprType)
        
        case (BDC_USERDEFID)
          ! Defined via analytic reference solution, identified by an id.
          ! The id is already specified in the collection.
          
          ! Call the user defined callback routine to evaluate the boudary value.
          ! No string tag is given which indicates that a reference function is to be
          ! evaluated.
          !
          ! As collection, we pass rcollection%p_rcollection here; this is a pointer
          ! to the application specific, global collection that may be of interest for
          ! callback routines. rcollection itself is actually a "local" collection,
          ! a "wrapper" for the rcollection of the application!
          call user_getBoundaryValues (&
              "",icomponent,rdiscretisation,p_rboundaryRegion,&
              dwhere, Dcoefficients(1,ipt,iel), rcollection%p_rnextCollection)
          
        case (BDC_USERDEF)
          ! This is a hardcoded, user-defined identifier.
          ! Call the user defined callback routine to evaluate the expression.
          !
          ! As collection, we pass rcollection%p_rcollection here; this is a pointer
          ! to the application specific, global collection that may be of interest for
          ! callback routines. rcollection itself is actually a "local" collection,
          ! a "wrapper" for the rcollection of the application!
          call user_getBoundaryValues (&
              svalue,icomponent,rdiscretisation,p_rboundaryRegion,&
              dwhere, Dcoefficients(1,ipt,iel), rcollection%p_rnextCollection)
        
        case default
          ! Default handling. Evaluate the expression.
          Dcoefficients(1,ipt,iel) = struc_evalExpression (r_t_p_optcBDC%p_roptcBDC,cexprType,&
              iid,ivalue,dvalue,Dparams)
        
        end select
        
!        select case (ccontrolType)
!        
!        ! ---------------------------------
!        ! L2 Dirichlet boudary control
!        ! ---------------------------------
!        case (1)
!        
!          if (p_roptControl%dalphaL2BdC .ge. 0.0_DP) then
!            if (.not. associated(p_rcontrol)) then
!              call output_line ("Control not specified.", &
!                  OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBCNavSt2D")
!              call sys_halt()
!            end if
!            
!            select case (p_rphysics%cequation)
!            case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
!            
!              ! Where is the boundary control?
!              icontrolcomp = icomponent
!              if (p_roptControl%dalphaDistC .ge. 0.0_DP) icontrolcomp = icontrolcomp + 2
!
!              ! Quick and dirty implementation...
!              
!              ! Element type?
!              select case (p_rcontrol%p_rblockDiscr%RspatialDiscr(icontrolcomp)%RelementDistr(1)%celement)
!              case (EL_P2_1D)
!                
!                ! Type of info?
!                select case (cinfoNeeded)
!                case (DISCBC_NEEDFUNC,DISCBC_NEEDFUNCMID,DISCBC_NEEDDERIV, &
!                      DISCBC_NEEDINTMEAN,DISCBC_NEEDNORMALSTRESS)
!                  
!                  ! Search the boundary edge that contains the parameter value
!                  call tria_searchBoundaryEdgePar2D(p_rboundaryRegion%iboundCompIdx, dwhere, &
!                      rspaceDiscr%p_rtriangulation, rspaceDiscr%p_rboundary, iedge,&
!                      BDR_PAR_LENGTH)
!                  !print *,"Richtig?"
!                  
!                  ! Get the vector data
!                  call lsyssc_getbase_double (p_rcontrol%RvectorBlock(icontrolcomp),p_Ddata)
!                  
!                  ! Get the "relative" parameter value of dwhere.
!                  call storage_getbase_int (&
!                      rspaceDiscr%p_rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
!                  call storage_getbase_double (&
!                      rspaceDiscr%p_rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
!                      
!                  iptpos = p_IboundaryCpIdx(p_rboundaryRegion%iboundCompIdx)
!                  
!                  ! First vertex number = edge number
!                  dpar1 = p_DvertexParameterValue(iptpos-1+iedge)
!                  iv1 = iedge
!                  
!                  ! Second vertex number = first vertex + 1 or overall first
!                  if (iedge .lt. p_IboundaryCpIdx(p_rboundaryRegion%iboundCompIdx+1) &
!                                -p_IboundaryCpIdx(p_rboundaryRegion%iboundCompIdx)) then
!                    dpar2 = p_DvertexParameterValue(iptpos-1+iedge+1)
!                    iv2 = iedge+1
!                  else
!                    dpar2 = boundary_dgetMaxParVal(rspaceDiscr%p_rboundary, &
!                        p_rboundaryRegion%iboundCompIdx)
!                    iv2 = p_IboundaryCpIdx(p_rboundaryRegion%iboundCompIdx)
!                  end if
!                  
!                  ! Relative parameter value
!                  dpar = (boundary_convertParameter(rspaceDiscr%p_rboundary, &
!                      p_rboundaryRegion%iboundCompIdx, dwhere, &
!                      BDR_PAR_LENGTH, BDR_PAR_01)-dpar1)/(dpar2-dpar1)
!                      
!                  ! Quadratic interpolation over [-1,1]
!                  call mprim_quadraticInterpolation (2*dpar-1.0_DP,&
!                      p_Ddata(iv1),p_Ddata(iedge+rspaceDiscr%p_rtriangulation%NVBD),p_Ddata(iv2),dval)
!                  
!                  Dvalues(1) = Dvalues(1) + dval
!
!                end select          
!              
!              case default
!                call output_line ("Unsupported element.", &
!                    OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBCNavSt2D")
!                call sys_halt()
!              end select
!            
!            
!            case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
!              call output_line ("Equation not supported.", &
!                  OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBCNavSt2D")
!              call sys_halt()
!            end select
!            
!          end if
!        
!        end select
        
      end do  
      
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_assembleInhomNeumannRHS (roptcBDC,rrhs,dtime,cequation,&
      copType,rglobalData)

!<description>
  ! Assembles inhomogeneous Neumann boundary data into a RHS vector.
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
    
    integer :: ibdComponent, isegment, i, iprimal, iintervalEnds
    character(LEN=PARLST_LENLINEBUF) :: sstr,sexpr,svalue,sbdex1,sbdex2,sbdcomp
    real(DP) :: dpar1, dpar2
    integer :: ctype, ivalue, iid, ibctyp
    real(DP) :: dvalue
    logical :: bneedsParams
    type(t_linearForm) :: rlinform
    type(t_collection), target :: rcollection, ruserCollection
    type(t_boundary), pointer :: p_rboundary
    type(t_p_optcBDC) :: r_t_p_optcBDC

    ! A boundary region defining where the boundary conditionis applied
    type(t_boundaryRegion), target :: rboundaryRegion

    ! A pointer to the section with the expressions and the boundary conditions
    type(t_parlstSection), pointer :: p_rsection,p_rbdcond,p_rbdcondPrimal

    ! Fetch some parameters
    p_rboundary => rrhs%p_rblockDiscr%p_rboundary
    
    ! Save a pointer to roptBDC ti r_t_p_optcBDC for later use.
    r_t_p_optcBDC%p_roptcBDC => roptcBDC
    r_t_p_optcBDC%p_rboundaryRegion => rboundaryRegion

    ! Get the expression/bc sections from the bondary condition block
    call parlst_querysection(roptcBDC%p_rparList, &
        roptcBDC%ssectionBdExpr, p_rsection)
        
    ! Get the section defining the primal or dual boundary conditions
    select case (copType)
    case (OPTP_PRIMAL)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcond)

    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)

      ! Primal boundary conditions, linearised eqn.
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrimLin, p_rbdcond)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcondPrimal)

    case (OPTP_DUAL)
    
      ! Dual boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondDual, p_rbdcond)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcondPrimal)

    case (OPTP_DUALLIN,OPTP_DUALLIN_SIMPLE)
    
      ! Dual boundary conditions, linearised equation
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondDualLin, p_rbdcond)

      ! Primal boundary conditions
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcondPrimal)

    end select
        
    ! Initialise the user-defined assembly
    call collct_init (rusercollection)

    ! Loop through all boundary components we have.
    do ibdComponent = 1,boundary_igetNBoundComp(p_rboundary)

      ! Parse the parameter "bdComponentX"
      write (sexpr,"(I10)") ibdComponent
      sbdcomp = "bdComponent" // adjustl(sexpr)
      
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      i = parlst_queryvalue (p_rbdcond, sbdcomp)

      if (i .ne. 0) then
      
        ! get the corresponding index in the "primal" dection
        if (copType .ne. OPTP_PRIMAL) then
          iprimal = parlst_queryvalue (p_rbdcondPrimal, sbdcomp)
        end if
      
        ! Parameter exists. Get the values in there.
        do isegment = 1,parlst_querysubstrings (p_rbdcond, sbdcomp)
          
          call parlst_getvalue_string (p_rbdcond, i, sstr, isubstring=isegment)
          
          ! Read the segment parameters
          read(sstr,*) dpar2,iintervalEnds,ibctyp
          
          ! Form a boundary condition segment that covers that boundary part
          if (dpar2 .ge. dpar1) then
            
            rboundaryRegion%dminParam = dpar1
            rboundaryRegion%dmaxParam = dpar2
            rboundaryRegion%iboundCompIdx = ibdComponent
            rboundaryRegion%dmaxParamBC = &
              boundary_dgetMaxParVal(p_rboundary, ibdComponent)
            rboundaryRegion%iproperties = iintervalEnds
            
            ! Type of equation?
            select case (cequation)
            
            ! ----------------------
            ! Stokes / Navier-Stokes
            ! ----------------------
            case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
            
              ! Now, which type of BC is to be created?
              select case (ibctyp)
              
              ! ------------------------------------------------------------------
              ! Automatic boundary conditions. Dual and linearised equations only.
              ! ------------------------------------------------------------------
              case (-1)
                if (copType .eq. OPTP_PRIMAL) then
                  call output_line (&
                      "Automatic boundary conditions only allowed.", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                end if
                
                ! The rule is:
                !  * Primal Neumann = Dual Neumann,
                !  * Primal Dirichlet = Dual Dirichlet-0

                ! Get the boundary conditions of the primal equation
                call parlst_getvalue_string (&
                    p_rbdcondPrimal, iprimal, sstr, isubstring=isegment)
                
                ! Read the segment parameters
                read(sstr,*) dpar2,iintervalEnds,ibctyp
                
                select case (ibctyp)
                
                ! --------------------------------------------------
                ! Neumann boundary conditions
                ! --------------------------------------------------
                case (0)
                  
                ! --------------------------------------------------
                ! Dirichlet boundary conditions
                ! --------------------------------------------------
                case (1)
                  
                ! --------------------------------------------------
                ! L2 Dirichlet boundary control
                ! --------------------------------------------------
                case (7)
                  
                ! --------------------------------------------------
                ! Inhomogeneous Neumann boundary conditions
                ! --------------------------------------------------
                case (8)
                  ! Nothing to do here
                  
                case default
                  call output_line ("Cannot set up automatic boundary conditions", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                  
                end select
                
              ! --------------------------------------------------
              ! Neumann boundary conditions
              ! --------------------------------------------------
              case (0)
              
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)

              ! --------------------------------------------------
              ! Dirichlet boundary control
              ! --------------------------------------------------
              case (7)

              ! --------------------------------------------------
              ! Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (8)
                ! Nothibg to do here

              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
              end select ! ibctyp

            ! ----------------------
            ! Heat equation
            ! ----------------------
            case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
            
              ! Now, which type of BC is to be created?
              select case (ibctyp)
              
              ! ------------------------------------------------------------------
              ! Automatic boundary conditions. Dual and linearised equations only.
              ! ------------------------------------------------------------------
              case (-1)
                if (copType .eq. OPTP_PRIMAL) then
                  call output_line (&
                      "Automatic boundary conditions only allowed.", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                end if
                
                ! The rule is:
                !  * Primal Neumann = Dual Neumann,
                !  * Primal Dirichlet = Dual Dirichlet-0

                ! Get the boundary conditions of the primal equation
                call parlst_getvalue_string (&
                    p_rbdcondPrimal, iprimal, sstr, isubstring=isegment)
                
                ! Read the segment parameters
                read(sstr,*) dpar2,iintervalEnds,ibctyp
                
                select case (ibctyp)
                
                ! --------------------------------------------------
                ! Neumann boundary conditions
                ! --------------------------------------------------
                case (0)
                
                ! --------------------------------------------------
                ! Dirichlet boundary conditions
                ! --------------------------------------------------
                case (1)
                
                ! --------------------------------------------------
                ! Dirichlet boundary control
                ! --------------------------------------------------
                case (7)

                ! --------------------------------------------------
                ! Inhomogeneous Neumann boundary conditions
                ! --------------------------------------------------
                case (8)
                  ! Nothibg to do here

                case default
                  call output_line ("Cannot set up automatic boundary conditions", &
                      OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                  call sys_halt()
                end select ! ibctyp

              ! --------------------------------------------------
              ! Neumann boundary conditions
              ! --------------------------------------------------
              case (0)
              
              ! --------------------------------------------------
              ! Dirichlet boundary conditions
              ! --------------------------------------------------
              case (1)
                
              ! --------------------------------------------------
              ! Inhomogeneous Neumann boundary conditions
              ! --------------------------------------------------
              case (8)

                ! Prepare the assembly
                !    
                ! Get the definition of the boundary condition.
                ! Read the line again, get the expressions for the data field
                read(sstr,*) dvalue,iintervalEnds,ibctyp,sbdex1
              
                ! For any string <> "", create the appropriate Dirichlet boundary
                ! condition and add it to the list of boundary conditions.
                !
                ! The IquickAccess array is set up as follows:
                !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
                !  IquickAccess(2) = expression type
                !  IquickAccess(3) = iid
                !  IquickAccess(4) = ivalue
                !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
                !  IquickAccess(7:...) = The binary content of r_t_p_optcBDC.
                !
                ! The DquickAccess array is set up as follows:
                !  DquickAccess(1) = dtime
                !  DquickAccess(2) = dvalue
                !
                ! The SquickAccess array is set up as follows:
                !  SquickAccess(1) = Name of the expression
                
                if (sbdex1 .ne. "") then
                  
                  ! X-velocity
                  !
                  ! Get the expression information
                  call struc_getBdExprInfo (roptcBDC,sbdex1,&
                      ctype,iid,ivalue,dvalue,svalue,bneedsParams)
                      
                  rcollection%IquickAccess(1) = 1
                  rcollection%IquickAccess(2) = ctype
                  rcollection%IquickAccess(3) = iid
                  rcollection%IquickAccess(4) = ivalue
                  rcollection%IquickAccess(5) = 0
                  if (bneedsParams) rcollection%IquickAccess(5) = 1
                  rcollection%DquickAccess(1) = dtime
                  rcollection%DquickAccess(2) = dvalue
                  rcollection%SquickAccess(1) = svalue
                  rcollection%p_rnextCollection => ruserCollection
                  rcollection%IquickAccess(6) = 0
                  rcollection%IquickAccess(7:) = &
                      transfer(r_t_p_optcBDC,rcollection%IquickAccess(7:),&
                                    size(rcollection%IquickAccess(7:)))
                  
                  call user_initCollectForVecAssembly (&
                      rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                  
                  ! Assemble the BCs.
                  rlinform%itermCount = 1
                  rlinform%Idescriptors(1) = DER_FUNC
                  call linf_buildVectorScalarBdr2D (rlinform, CUB_G4_1D, .false., &
                      rrhs%RvectorBlock(1),cc_getNeumannBdData,&
                      rboundaryRegion,rcollection)
                      
                  call user_doneCollectForVecAssembly (rglobalData,rusercollection)
                
                end if

              case default
                call output_line ("Unknown boundary condition!", &
                    OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
                call sys_halt()
              end select ! ibctyp

            case default
              
              call output_line ("Unknown equation", &
                  OU_CLASS_ERROR,OU_MODE_STD,"cc_parseBDconditions")
              call sys_halt()
              
            end select ! equation
            
            ! Move on to the next parameter value
            dpar1 = dpar2
            
          end if
                                            
        end do ! isegment
      
      end if
      
    end do ! ibdComponent

    call collct_done (rusercollection)

  end subroutine

end module
