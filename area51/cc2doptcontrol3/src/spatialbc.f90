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
!#     -> Assembles the definition of the BC's from sections given by DAT files
!#        and sets up an analytical boundary condition description.
!#
!# 2.) sbc_releaseBoundaryList
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
  
  use collection
  use convection
    
  use constantsdiscretisation
  
  use structuresgeneral
  use structuresboundaryconditions
  
  use timediscretisation

  use user_callback
  
  implicit none
  
  private
  
!<constants>

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
    integer :: nregions = 0

    ! Pointer to the head of a list of Neumann boundary regions, primal equation
    type(t_bdRegionEntry), pointer :: p_rbdHead => null()

    ! Pointer to the tail of a list of Neumann boundary regions, primal equation
    type(t_bdRegionEntry), pointer :: p_rbdTail => null()
  
  end type

!</typeblock>

!<types>

  public :: t_bdRegionEntry
  public :: t_boundaryRegionList
  public :: sbc_assembleBDconditions
  public :: sbc_releaseBoundaryList
  public :: sbc_implementDirichletBC

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sbc_addBoundaryRegion (rboundaryRegion,rregionList)
  
!<description>
  ! Adds bodunary region to a list of boundary regions.
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

  subroutine sbc_assembleBDconditions (roptcBDC,dtime,&
      cspace,rglobalData,casmFlags,rtimediscr,rspaceDiscr,&
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
  type(t_optcBDC), intent(inout), target :: roptcBDC

  ! Current simulation time where to assemble the boundary conditions.
  real(dp), intent(in) :: dtime

  ! Space for which to assemble the boundary conditions.
  ! CCSPACE_PRIMAL: BC for the primal space.
  ! CCSPACE_DUAL: BC for the dual space.
  integer, intent(in) :: cspace
  
  ! Global data, passed to callback routines
  type(t_globalData), intent(inout) :: rglobalData

  ! Assembly flags, determines what to assemble.
  ! One of the SBC_xxxx flags.
  integer, intent(in) :: casmFlags
  
  ! A discretisation structure defining the underlying time discretisation.
  type(t_timeDiscretisation), intent(in) :: rtimeDiscr
  
  ! OPTIONAL: A discretisation structure defining the space discretisation
  ! of the current level.
  ! Can be omitted if only Neumann bonudary is to be created.
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
    
    integer :: ibdComponent, isegment, i, iintervalEnds
    character(LEN=PARLST_LENLINEBUF) :: sstr,sexpr,svalue,sbdex1,sbdex2
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
    type(t_parlstSection), pointer :: p_rsection,p_rbdcond

    ! Fetch some parameters
    p_rboundary => rspaceDiscr%p_rboundary
    
    ! Save a pointer to roptBDC ti r_t_p_optcBDC for later use.
    r_t_p_optcBDC%p_roptcBDC => roptcBDC

    ! Get the expression/bc sections from the bondary condition block
    call parlst_querysection(roptcBDC%p_rparList, &
        roptcBDC%ssectionBdExpr, p_rsection)
        
    ! Get the section defining the primal or dual boundary conditions
    if (cspace .eq. CCSPACE_PRIMAL) then
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondPrim, p_rbdcond)
    else
      call parlst_querysection(roptcBDC%p_rparList, &
          roptcBDC%ssectionBdCondDual, p_rbdcond)
    end if
        
    ! Initialise the user-defined assembly
    call collct_init (rusercollection)

    ! Loop through all boundary components we have.
    do ibdComponent = 1,boundary_igetNBoundComp(p_rboundary)

      ! Parse the parameter 'bdComponentX'
      write (sexpr,'(I10)') ibdComponent
      sstr = 'bdComponent' // adjustl(sexpr)
      
      ! We start at parameter value 0.0.
      dpar1 = 0.0_DP
      
      i = parlst_queryvalue (p_rbdcond, sstr)
      if (i .ne. 0) then
        ! Parameter exists. Get the values in there.
        do isegment = 1,parlst_querysubstrings (p_rbdcond, sstr)
          
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
            
            ! Now, which type of BC is to be created?
            select case (ibctyp)
            
            case (0)
              if (iand(casmFlags,SBC_NEUMANN) .ne. 0) then
                if (present(rneumannBoundary)) then
                  ! Add the bondary region to the Neumann boundary regions
                  ! if there is a structure present.
                  call sbc_addBoundaryRegion(&
                      rboundaryRegion,rneumannBoundary)
                end if
              end if
            
            case (1)

              if (iand(casmFlags,SBC_DIRICHLETBC) .ne. 0) then

                ! Simple Dirichlet boundary.
                ! Get the definition of the boundary condition.
                ! Read the line again, get the expressions for X- and Y-velocity
                read(sstr,*) dvalue,iintervalEnds,ibctyp,sbdex1,sbdex2
              
                if ((cspace .eq. CCSPACE_DUAL) .and. (ibctyp .eq. -1)) then
                  
                  ! Automatic dual boundary conditions. Zero boundary conditions
                  ! where there are Dirichlet boundary conditions in the primal
                  ! equation.
                
                  rcollection%IquickAccess(2) = BDC_VALDOUBLE
                  rcollection%IquickAccess(3) = 0
                  rcollection%IquickAccess(4) = 0
                  rcollection%IquickAccess(5) = 0
                  rcollection%DquickAccess(1) = 0.0_DP
                  rcollection%DquickAccess(1) = dtime
                  rcollection%SquickAccess(1) = ""
                  rcollection%p_rnextCollection => ruserCollection
                  rcollection%IquickAccess(6:) = &
                      transfer(r_t_p_optcBDC,rcollection%IquickAccess(6:),&
                                     size(rcollection%IquickAccess(6:)))

                  ! X-velocity
                  rcollection%IquickAccess(1) = 1
                  
                  call user_initCollectForVecAssembly (&
                      rglobalData,0,rcollection%IquickAccess(1),dtime,rusercollection)
                  
                  ! Assemble the BC's.
                  call bcasm_newDirichletBConRealBD (&
                      rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,rdiscreteBC,&
                      cc_getDirBCNavSt2D,rcollection)
                      
                  call user_doneCollectForAssembly (rglobalData,rusercollection)

                  ! Y-velocity
                  rcollection%IquickAccess(1) = 2
                  
                  call user_initCollectForVecAssembly (&
                      rglobalData,0,rcollection%IquickAccess(1),dtime,rusercollection)
                  
                  ! Assemble the BC's.
                  call bcasm_newDirichletBConRealBD (&
                      rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,rdiscreteBC,&
                      cc_getDirBCNavSt2D,rcollection)
                      
                  call user_doneCollectForAssembly (rglobalData,rusercollection)
                
                else
                    
                  ! For any string <> '', create the appropriate Dirichlet boundary
                  ! condition and add it to the list of boundary conditions.
                  !
                  ! The IquickAccess array is set up as follows:
                  !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
                  !  IquickAccess(2) = expression type
                  !  IquickAccess(3) = iid
                  !  IquickAccess(4) = ivalue
                  !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
                  !  IquickAccess(6:...) = The binary content of r_t_p_optcBDC.
                  !
                  ! The DquickAccess array is set up as follows:
                  !  DquickAccess(1) = dvalue
                  !  DquickAccess(2) = dtime
                  !
                  ! The SquickAccess array is set up as follows:
                  !  SquickAccess(1) = Name of the expression
                  
                  if (sbdex1 .ne. '') then
                    
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
                    if (bneedsParams) rcollection%IquickAccess(4) = 1
                    rcollection%DquickAccess(1) = dvalue
                    rcollection%DquickAccess(1) = dtime
                    rcollection%SquickAccess(1) = svalue
                    rcollection%p_rnextCollection => ruserCollection
                    rcollection%IquickAccess(6:) = &
                        transfer(r_t_p_optcBDC,rcollection%IquickAccess(6:),&
                                       size(rcollection%IquickAccess(6:)))
                    
                    call user_initCollectForVecAssembly (&
                        rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                    
                    ! Assemble the BC's.
                    call bcasm_newDirichletBConRealBD (&
                        rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,rdiscreteBC,&
                        cc_getDirBCNavSt2D,rcollection)
                        
                    call user_doneCollectForAssembly (rglobalData,rusercollection)
                    
                  end if
                      
                  if (sbdex2 .ne. '') then
                  
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
                    if (bneedsParams) rcollection%IquickAccess(4) = 1
                    rcollection%DquickAccess(1) = dvalue
                    rcollection%DquickAccess(1) = dtime
                    rcollection%SquickAccess(1) = svalue
                    rcollection%p_rnextCollection => ruserCollection
                    rcollection%IquickAccess(6:) = &
                        transfer(r_t_p_optcBDC,rcollection%IquickAccess(6:),&
                                       size(rcollection%IquickAccess(6:)))
                    
                    call user_initCollectForVecAssembly (&
                        rglobalData,iid,rcollection%IquickAccess(1),dtime,rusercollection)
                    
                    ! Assemble the BC's.
                    call bcasm_newDirichletBConRealBD (&
                        rspaceDiscr,rcollection%IquickAccess(1),rboundaryRegion,rdiscreteBC,&
                        cc_getDirBCNavSt2D,rcollection)
                        
                    call user_doneCollectForAssembly (rglobalData,rusercollection)
                    
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
!      ! For implementing boundary conditions, we use a 'filter technique with
!      ! discretised boundary conditions'. This means, we first have to calculate
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
!          parlst_querysubstrings (p_rsection, 'bdExpressions'))
!      
!      ! Add the parser to the collection
!      call collct_setvalue_pars (rcoll, BDC_BDPARSER, rparser, &
!                                  .true., 0, SEC_SBDEXPRESSIONS)
!      
!      ! Add the boundary expressions to the collection into the
!      ! specified section.
!      do i=1,parlst_querysubstrings (p_rsection, 'bdExpressions')
!      
!        call parlst_getvalue_string (p_rsection, 'bdExpressions', cstr, "", i)
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
!          call output_line ('Expressions not implemented!', &
!              OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
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
!          ! Parse the parameter 'bdComponentX'
!          write (cexpr,'(I10)') ibdComponent
!          cstr = 'bdComponent' // adjustl(cexpr)
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
!                    ! For any string <> '', create the appropriate Dirichlet boundary
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
!                      ! Primal BC's
!                      if (cbctype .eq. CCSPACE_PRIMAL) then
!                      
!                        ! If the type is a double precision value, set the DquickAccess(4)
!                        ! to that value so it can quickly be accessed.
!                        if (sbdex1 .ne. '') then
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
!                          ! Assemble the BC's.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. '') then
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
!                          ! Assemble the BC's.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                      end if
!                      
!                    end if
!                    
!                    if (iprimaldual .eq. 2) then
!                    
!                      ! Dual BC's
!                      if (cbctype .eq. CCSPACE_DUAL) then
!
!                        ! Now the same thing again, this time separately for primal and dual
!                        ! variables.
!                        ! If a primal velocity is specified, Dirichlet-0-boundary conditions are
!                        ! assumed for the corresponding dual.
!                        if (sbdex1 .ne. '') then
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
!                          rcoll%SquickAccess(1) = ''
!                          iid = 0
!                          
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,4,dtime,rcallbackcollection)
!
!                          ! Assemble the BC's.
!                          ! If we only assemble dual BC's and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. '') then
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
!                          rcoll%SquickAccess(1) = ''
!                          iid = 0
!                          
!                          rcoll%Dquickaccess(1) = dtime
!                          call user_initCollectForVecAssembly (rglobalData,iid,5,dtime,rcallbackcollection)
!
!                          ! Assemble the BC's.
!                          ! If we only assemble dual BC's and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
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
!                    ! For any string <> '', create the appropriate Dirichlet boundary
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
!                      ! Primal BC's
!                      if (cbctype .eq. CCSPACE_PRIMAL) then
!                        ! If the type is a double precision value, set the DquickAccess(4)
!                        ! to that value so it can quickly be accessed.
!                        if (sbdex1 .ne. '') then
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
!                          ! Assemble the BC's.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. '') then
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
!                          ! Assemble the BC's.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                      end if
!                      
!                    end if
!                    
!                    if (iprimaldual .eq. 2) then
!                    
!                      ! Dual BC's
!                      if (cbctype .eq. CCSPACE_DUAL) then
!                      
!                        ! Now the same thing again, this time separately for primal and dual
!                        ! variables.
!                        ! If a velocity is not specified, Dirichlet-0-boundary conditions are
!                        ! assumed.
!                        
!                        if (sbdex3 .ne. '') then
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
!                          ! Assemble the BC's.
!                          ! If we only assemble dual BC's and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                        if (sbdex4 .ne. '') then
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
!                          ! Assemble the BC's.
!                          ! If we only assemble dual BC's and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
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
!                  ! For any string <> '', create the appropriate Dirichlet boundary
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
!                      ! Primal BC's
!                      if (cbctype .eq. CCSPACE_PRIMAL) then
!                        ! If the type is a double precision value, set the DquickAccess(4)
!                        ! to that value so it can quickly be accessed.
!                        if (sbdex1 .ne. '') then
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
!                          ! Assemble the BC's.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!                              
!                        end if
!                        
!                        if (sbdex2 .ne. '') then
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
!                          ! Assemble the BC's.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
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
!                  ! For any string <> '', create the appropriate Dirichlet boundary
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
!                      ! Dual BC's
!                      if (cbctype .eq. CCSPACE_DUAL) then
!                      
!                        ! Now the same thing again, this time separately for primal and dual
!                        ! variables.
!                        ! If a velocity is not specified, Dirichlet-0-boundary conditions are
!                        ! assumed.
!                        
!                        if (sbdex3 .ne. '') then
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
!                          ! Assemble the BC's.
!                          ! If we only assemble dual BC's and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!
!                        end if
!                        
!                        if (sbdex4 .ne. '') then
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
!                          ! Assemble the BC's.
!                          ! If we only assemble dual BC's and there are 3 solution components,
!                          ! we assume the vector to specify exactly the dual solution.
!                          call bcasm_newDirichletBConRealBD (&
!                              rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                              cc_getBDconditionsNavSt2D,rcoll)
!                              
!                          call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
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
!                        ! For any string <> '', create the appropriate Dirichlet boundary
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
!                        ! Primal BC's
!                        if (cbctype .eq. CCSPACE_PRIMAL) then
!                        
!                          ! If the type is a double precision value, set the DquickAccess(4)
!                          ! to that value so it can quickly be accessed.
!                          if (sbdex1 .ne. '') then
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
!                            call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!                                
!                          end if
!                          
!                          if (sbdex2 .ne. '') then
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
!                            ! Assemble the BC's.
!                            call linf_buildVectorScalarBdr2D (rlinformRhs, CUB_G4_1D, .false., &
!                                rvectorDirichletBCCRHS%RvectorBlock(1),&
!                                fcoeff_buildBCCRHS,rboundaryRegion, rcoll)
!
!                            call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
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
!                      ! Dual BC's
!                      if (cbctype .eq. CCSPACE_DUAL) then
!
!                        ! dual X-velocity if primal X-velocity exists
!                        
!                        ! Ditichlet-0 boundary
!                        iexptyp = BDC_VALDOUBLE
!                        rcoll%Dquickaccess(4) = 0.0_DP
!                        rcoll%IquickAccess(3) = iexptyp
!                        rcoll%SquickAccess(1) = ''
!                        iid = 0
!                        
!                        rcoll%Dquickaccess(1) = dtime
!                        call user_initCollectForVecAssembly (rglobalData,iid,4,dtime,rcallbackcollection)
!
!                        ! Assemble the BC's.
!                        ! The 2nd element in IquickAccess saves the component number.
!                        ! If we only assemble dual BC's and there are 3 solution components,
!                        ! we assume the vector to specify exactly the dual solution.
!
!                        rcoll%IquickAccess(2) = 1
!                        call bcasm_newDirichletBConRealBD (&
!                            rspaceDiscr,1+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                            cc_getBDconditionsNavSt2D,rcoll)
!
!                        call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!
!                        call user_initCollectForVecAssembly (rglobalData,iid,5,dtime,rcallbackcollection)
!
!                        rcoll%IquickAccess(2) = 2
!                        call bcasm_newDirichletBConRealBD (&
!                            rspaceDiscr,2+rspaceDiscr%ncomponents-3,rboundaryRegion,rdiscreteBC,&
!                            cc_getBDconditionsNavSt2D,rcoll)
!                            
!                        call user_doneCollectForAssembly (rglobalData,rcallbackcollection)
!                            
!                      end if
!                      
!                    end if
!                  
!                  end if
!
!                case default
!                  call output_line ('Unknown boundary condition!', &
!                      OU_CLASS_ERROR,OU_MODE_STD,'cc_parseBDconditions')
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
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    integer :: icomponent,cexprtype,iid,ivalue
    real(DP) :: dvalue
    logical :: bneedsParams
    real(DP), dimension(size(SEC_EXPRVARIABLES)) :: Dparams
    real(DP) :: d, dx, dy
    character(LEN=PARLST_LENLINEBUF) :: svalue
    type(t_p_optcBDC) :: r_t_p_optcBDC
    
    real(DP) :: dtime

    ! The IquickAccess array is set up as follows:
    !  IquickAccess(1) = component under consideration (1=x-vel, 2=y-vel,...)
    !  IquickAccess(2) = expression type
    !  IquickAccess(3) = iid
    !  IquickAccess(4) = ivalue
    !  IquickAccess(5) = 1, if parameters (x,y) are needed, =0 otherwise
    !  IquickAccess(6) = Pointer to the t_optcBDC structure
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
    r_t_p_optcBDC   = transfer(rcollection%IquickAccess(6:),r_t_p_optcBDC)
    dvalue       = rcollection%Dquickaccess(1)
    
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
        ! callback routines. rcollection itself is actually a 'local' collection,
        ! a 'wrapper' for the rcollection of the application!
        call user_getBoundaryValues (&
            "",icomponent,rspaceDiscr,rboundaryRegion,&
            dwhere, Dvalues(1), rcollection%p_rnextCollection)
        
      case (BDC_USERDEF)
        ! This is a hardcoded, user-defined identifier.
        ! Call the user defined callback routine to evaluate the expression.
        !
        ! As collection, we pass rcollection%p_rcollection here; this is a pointer
        ! to the application specific, global collection that may be of interest for
        ! callback routines. rcollection itself is actually a 'local' collection,
        ! a 'wrapper' for the rcollection of the application!
        call user_getBoundaryValues (&
            svalue,icomponent,rspaceDiscr,rboundaryRegion,&
            dwhere, Dvalues(1), rcollection%p_rnextCollection)
      
      case default
        ! Default handling. Evaluate the expression.
        Dvalues(1) = struc_evalExpression (r_t_p_optcBDC%p_roptcBDC,cexprType,&
            iid,ivalue,dvalue,Dparams)
      
      end select
      
    end select
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine sbc_implementDirichletBC (roptcBDC,dtime,cspace,&
      rtimeDiscr,rspaceDiscr,rglobalData,&
      rx,rb,rd,rdiscretebc)

!<description>
  ! Implements the boundary conditions at time dtime, for a solution
  ! vector rx, a rhs vector rb and/or a defect vector rb.
!</description>

!<input>
  ! Boundary conditions in the problem.
  type(t_optcBDC), intent(in) :: roptcBDC

  ! Time where the BC's should be implemented.
  real(DP), intent(IN) :: dtime

  ! Type of solution space. CCSPACE_PRIMAL or CCSPACE_DUAL.
  integer, intent(in) :: cspace

  ! Time discretisation, dtime refers to.
  type(t_timeDiscretisation), intent(in) :: rtimeDiscr

  ! Space discretisation
  type(t_blockDiscretisation), intent(in) :: rspaceDiscr

  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
!</input>

!<inputoutput>
  ! OPTIONAL: Solution vector
  type(t_vectorBlock), intent(inout), optional :: rx

  ! OPTIONAL: RHS vector
  type(t_vectorBlock), intent(inout), optional :: rb

  ! OPTIONAL: Defect vector
  type(t_vectorBlock), intent(inout), optional :: rd

  ! OPTIONAL: Boundary condition structure which receives the boudary
  ! conditions. If present, it must have been initialised
  ! by the caller. If present, it will be attached to the vectors.
  type(t_discreteBC), intent(inout), optional :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata
  
    ! DEBUG!!!
    call lsysbl_getbase_double (rd,p_Ddata)

    ! Boundary condition structure which receives the boudary
    ! conditions. 
    type(t_discreteBC) :: rdiscreteBClocal

    if (.not. present(rdiscreteBC)) then
      ! Initialise the boundary conditions
      call bcasm_initDiscreteBC(rdiscreteBClocal)

      ! Assemble the BC's.
      call sbc_assembleBDconditions (roptcBDC,dtime,cspace,&
          rglobalData,SBC_DIRICHLETBC,&
          rtimeDiscr,rspaceDiscr,rdiscreteBClocal)

      ! Implement the boundary conditions into the vector(s).
      if (present(rx)) then
        call vecfil_discreteBCsol (rx,rdiscreteBClocal)
      end if
      
      if (present(rb) then
        call vecfil_discreteBCrhs (rb,rdiscreteBClocal)
      end if
      
      if (present(rd)) then
        call vecfil_discreteBCdef (rd,rdiscreteBClocal)
      end if

      ! Release the BCs again.
      call bcasm_releaseDiscreteBC(rdiscreteBClocal)

    else

      ! Assemble the BC's.
      call sbc_assembleBDconditions (roptcBDC,dtime,cspace,&
          rglobalData,SBC_DIRICHLETBC,&
          rtimeDiscr,rspaceDiscr,rdiscreteBC)

      ! Implement the boundary conditions into the vector(s).
      ! Attach the structure.
      if (present(rx)) then
        call vecfil_discreteBCsol (rx,rdiscreteBC)
        call lsysbl_assignDiscreteBC (rx,rdiscreteBC)
      end if
      
      if (present(rb) then
        call vecfil_discreteBCrhs (rb,rdiscreteBC)
        call lsysbl_assignDiscreteBC (rb,rdiscreteBC)
      end if
      
      if (present(rd)) then
        call vecfil_discreteBCdef (rd,rdiscreteBC)
        call lsysbl_assignDiscreteBC (rd,rdiscreteBC)
      end if

    end if
    
  end subroutine

end module
