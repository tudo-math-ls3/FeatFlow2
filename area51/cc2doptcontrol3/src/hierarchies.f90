!##############################################################################
!# ****************************************************************************
!# <name> hierarchies </name>
!# ****************************************************************************
!#
!# <purpose>
!# 
!# This module contains a set of routines to maintain space and space-time
!# hierarchies of meshes and functions.
!#
!# 1.) init_initParamTria
!#     -> Initialise a mesh
!#
!# 2.) init_doneParamTria
!#     -> Release a mesh
!#
!# 3.) init_initTimeHierarchy
!#     -> Initialise a time hierarchy
!'
!# 4.) init_initSpaceDiscrHier
!#     -> Initialise a spatial hierarchy
!#
!# 5.) init_initSpaceTimeHierarchy
!#     -> Initialise a space-time hierarchy
!#
!# 6.) init_initSpacePrjHierarchy
!#     -> Initialise a projection hierarchy for prolongation and restriction
!#        in space
!#
!# 7.) init_initSpaceTimePrjHierarchy
!#     -> Initialise a projection hierarchy for prolongation and restriction
!#        in space-time
!# </purpose>
!##############################################################################

module hierarchies

  use fsystem
  use storage
  use genoutput
  use paramlist
  use basicgeometry
  use boundary
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  
  use spatialdiscretisation
  use timediscretisation
  use timescalehierarchy
  
  use collection
  
  use analyticsolution
  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimehierarchy
  
  use constantsdiscretisation
  use structuresdiscretisation
  use structuresgeneral
  use structuresoptcontrol
  
  use kktsystemspaces
  
  use spacediscretisation
  use spacetimeinterlevelprojection
  
  implicit none
  
  private

!<types>

!<typeblock>
  ! Type encapsuling the information used for creating a discretisation
  ! in space.
  type t_spaceDiscrParams
  
    ! Discretisation to create.
    integer :: cspace = CCSPACE_PRIMAL
    
    ! Pointer to physical parameters
    type(t_settings_physics), pointer :: p_rphysics => null()
    
    ! Pointer to structure defining the optimal control problem to calculate.
    ! can be NULL for cspace = CCSPACE_PRIMAL.
    type(t_settings_optcontrol), pointer :: p_roptControl => null()
    
    ! Pointer to discretisation settings
    type(t_settings_spacediscr), pointer :: p_rsettingsSpaceDiscr => null()
    
    ! Space discretisation hierarchy of the primal space or NULL if not available.
    type(t_feHierarchy), pointer :: p_rfeHierarchyPrimal => null()
    
  end type
  
!</typeblock>

  public :: t_spaceDiscrParams

!</types>
  
  ! Initialise a mesh
  public :: init_initParamTria
  
  ! Release a mesh
  public :: init_doneParamTria
  
  ! Initialise a time hierarchy
  public :: init_initTimeHierarchy
  
  ! Initialise a spatial hierarchy
  public :: init_initSpaceDiscrHier

  ! Releases a spatial hierarchy  
  public :: init_doneSpaceDiscrHier
  
  ! Initialise a space-time hierarchy
  public :: init_initSpaceTimeHierarchy
  
  ! Initialise a projection hierarchy for prolongation and restriction in space
  public :: init_initSpacePrjHierarchy
  
  ! Initialise a projection hierarchy for prolongation and restriction
  ! in space-time
  public :: init_initSpaceTimePrjHierarchy
  
  ! Releases a projection hierarchy.
  public :: init_doneSpaceTimePrjHierarchy
  
  ! Initialises a space hierarchy based on a coarse mesh and a refinement
  ! strategy specified in rrefinement.
  public :: init_initSpaceHierarchy
  
  ! Callback function which creates the discretisation on a level.
  public :: fgetDist1LvDiscr
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateSlicedQuadMesh (rtriangulation,ncellsX)
  
!<description>
  ! This routine generates a standard [0,1]^2 QUAD mesh with ncellsX cells in
  ! X-direction. By nature, these cells show an anisotropy of ncellsX:1.
!</description>
  
!<input>
  ! Number of cells in X-direction.
  integer, intent(IN) :: ncellsX
!</input>
    
!<output>
  ! Triangulation structure that receives the triangulation.
  type(t_triangulation), intent(OUT) :: rtriangulation
!</output>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:,:), pointer :: p_Idata2D
    integer, dimension(:), pointer :: p_Idata,p_IverticesAtBoundary,p_IboundaryCpIdx
    integer :: ivt, iel
    integer, dimension(2) :: Isize
    
    ! Initialise the basic mesh
    rtriangulation%ndim = NDIM2D
    rtriangulation%NEL = ncellsX
    rtriangulation%NVT = (ncellsX+1)*2
    rtriangulation%NMT = 0
    rtriangulation%NNVE = 4
    rtriangulation%NNEE = 4
    rtriangulation%NBCT = 1
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(TRIA_NVEQUAD2D) = rtriangulation%NEL
    
    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM2D, NVT)
    Isize = (/NDIM2D,rtriangulation%NVT/)
    call storage_new ('tria_read_tri2D', 'DCORVG', Isize, ST_DOUBLE, &
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
    
    ! Initialise the point coordinates.
    ! Odd vertices on the bottom, even vertices on top of the QUAD mesh,
    ! numbered from left to right.
    do ivt=0,ncellsX
      p_Ddata2D(1,2*ivt+1) = real(ivt,DP)/real(ncellsX,DP)
      p_Ddata2D(2,2*ivt+1) = 0.0_DP

      p_Ddata2D(1,2*ivt+2) = real(ivt,DP)/real(ncellsX,DP)
      p_Ddata2D(2,2*ivt+2) = 1.0_DP
    end do
    
    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriangulation%NNVE,rtriangulation%NEL/)
    call storage_new ('tria_read_tri2D', 'KVERT', Isize, ST_INT, &
        rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)

    ! Initialise the connectivity for the cells.
    do iel=0,ncellsX-1
      p_Idata2D(1,iel+1) = 2*iel+1
      p_Idata2D(2,iel+1) = 2*iel+3
      p_Idata2D(3,iel+1) = 2*iel+4
      p_Idata2D(4,iel+1) = 2*iel+2
    end do
    
    ! Allocate memory for InodalProperty
    call storage_new ('tria_read_tri2D', 'KNPR', &
        rtriangulation%NVT, ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_Idata)

    ! All vertices are on the boundary
    p_Idata(:) = 1
    
    ! Number of vertices on the boundary -- all of them
    rtriangulation%NVBD = rtriangulation%NVT
    
    ! Allocate memory for IverticesAtBoundary.
    call storage_new ('tria_generateBasicBoundary', &
        'KVBD', rtriangulation%NVBD, &
        ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    call storage_new ('tria_generateBasicBoundary', &
        'KBCT', rtriangulation%NBCT+1, &
        ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    call storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    call storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! The first element in p_IboundaryCpIdx is (as the head) always =1.
    p_IboundaryCpIdx(1) = 1
    p_IboundaryCpIdx(2) = 1 + rtriangulation%NVBD

    ! Initialise the numbers of the vertices on the boundary
    do ivt=0,ncellsX
      p_IverticesAtBoundary (ivt+1) = 2*ivt+1
      p_IverticesAtBoundary (2*(ncellsX+1)-ivt) = 2*ivt+2
    end do
    
    ! Allocate memory for  and DvertexParameterValue
    call storage_new ('tria_generateBasicBoundary', &
        'DVBDP', rtriangulation%NVBD, &
        ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
    
    ! Get the array where to store boundary parameter values.
    call storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
        
    ! Initialise the parameter values of the vertices. For the bottommost
    ! edge, they coincide wit the coordinate. For the topmost edge,
    ! that's 3 - x-coordinate.
    do ivt=1,ncellsX+1
      p_DvertexParameterValue (ivt) = p_Ddata2D(1,p_IverticesAtBoundary(ivt))
      p_DvertexParameterValue (ncellsX+1+ivt) = &
          3.0_DP - p_Ddata2D(1,p_IverticesAtBoundary(ncellsX+1+ivt))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initParamTria (rparlist,rboundary,rtriangulation,rdiscrTime,&
      ssectionSpace,ssectionTime)
  
!<description>
  ! Reads in / create the basic domain and coarse mesh in space and time.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the spatial mesh/domain can be found.
  character(len=*), intent(in) :: ssectionSpace

  ! Section where the parameters of the time mesh can be found.
  character(len=*), intent(in) :: ssectionTime
!</input>

!<output>
  ! Description of the boundary
  type(t_boundary), intent(out) :: rboundary
  
  ! Spatial coarse mesh. The result is a raw mesh.
  type(t_triangulation), intent(out) :: rtriangulation
  
  ! Time coarse mesh.
  ! The rdiscrTime%itag tag decides upon the type of one-step scheme, if
  ! a one-step scheme is chosen.
  ! =0: Standard, =1: Old (not respecting any minimisation problem)
  type(t_timeDiscretisation), intent(out) :: rdiscrTime
!</output>

!</subroutine>

  ! local variables
  integer :: imeshType,ncellsX
  
    ! Variable for a filename:
    character(LEN=60) :: sPRMFile, sTRIFile
    
    integer :: niterations
    real(dp) :: dtimeInit,dtimeMax,dtimeStepTheta
    integer :: ctimeStepScheme
    
    ! -----------------
    ! Space mesh/domain

    ! Get the .prm and the .tri file from the parameter list.
    ! note that parlst_getvalue_string returns us exactly what stands
    ! in the parameter file, so we have to apply READ to get rid of
    ! probable ""!
    call parlst_getvalue_string (rparlist,ssectionSpace,'sParametrisation',&
        sPRMFile,bdequote=.true.)
                              
    call parlst_getvalue_string (rparlist,ssectionSpace,'sMesh',&
        sTRIFile,bdequote=.true.)
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rboundary, sPrmFile)
        
    ! Now set up the basic triangulation. Which type of triangulation should
    ! be set up?
    call parlst_getvalue_int (rparlist,ssectionSpace,'imeshType',imeshType,0)
    select case (imeshType)
    case (0)
      ! Standard mesh, specified by a TRI file.
      call tria_readTriFile2D (rtriangulation, sTRIFile, rboundary)
    case (1)
      ! Sliced QUAD mesh with ncellsX cells on the coarse grid.
      call parlst_getvalue_int (rparlist,ssectionSpace,'ncellsX',ncellsX,1)
      call cc_generateSlicedQuadMesh(rtriangulation,ncellsX)
    case default
      call output_line ('Unknown mesh type!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initParamTria')
      call sys_halt()
    end select
    
    ! Create a standard mesh.
    call tria_initStandardMeshFromRaw(rtriangulation, rboundary)
    
    ! --------------------------
    ! Time mesh / discretisation
    
    ! Get the parameters
    call parlst_getvalue_int (rparlist,ssectionTime,'niterations', &
        niterations, 1000)
    call parlst_getvalue_double (rparlist,ssectionTime,'dtimestart', &
        dtimeInit, 0.0_DP)
    call parlst_getvalue_double (rparlist,ssectionTime,'dtimemax', &
        dtimeMax, 20.0_DP)
    call parlst_getvalue_int (rparlist, ssectionTime,'ctimeStepScheme', &
        ctimeStepScheme, 0)
    call parlst_getvalue_double (rparlist,ssectionTime,'dtimeStepTheta', &
        dtimeStepTheta, 1.0_DP)
    
    ! Initialise the coarse time discretisation
    select case (ctimeStepScheme)
    case (0)
      ! One-step scheme.
      call tdiscr_initOneStepTheta (rdiscrTime, &
          dtimeInit, dtimeMax, niterations, dtimeStepTheta)
          
      ! Default time stepping.
      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      rdiscrTime%itag = 1
    case (1)
      ! FS-Theta scheme.
      call tdiscr_initFSTheta (rdiscrTime, dtimeInit, dtimeMax, niterations)

    case (2)
      ! Old One-step scheme.
      call tdiscr_initOneStepTheta (rdiscrTime, &
          dtimeInit, dtimeMax, niterations, dtimeStepTheta)
          
      ! Old time stepping.
      ! itag=0/1 decides upon whether the new or old 1-step method is used.
      rdiscrTime%itag = 0

    case (3)
      ! dG(0)
      call tdiscr_initdG0 (rdiscrTime, dtimeInit, dtimeMax, niterations)
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_doneParamTria (rboundary,rtriangulation,rdiscrTime)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! Description of the boundary
  type(t_boundary), intent(inout) :: rboundary
  
  ! Coarse mesh. The result is a raw mesh.
  type(t_triangulation), intent(inout) :: rtriangulation

  ! Time coarse mesh.
  type(t_timeDiscretisation) :: rdiscrTime
!</inputoutput>

!</subroutine>

    ! Release the data.
    call tria_done (rtriangulation)
    call boundary_release (rboundary)
    call tdiscr_done(rdiscrTime)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceHierarchy (rboundary,rrefinement,&
      rtriaCoarse,rmeshHierarchy,ioutputLevel)
  
!<description>
  ! Initialises a space hierarchy based on a coarse mesh and a refinement
  ! strategy specified in rrefinement.
!</description>

!<input>
  ! Type of refinement
  type(t_settings_refinement), intent(in) :: rrefinement

  ! Description of the boundary
  type(t_boundary), intent(in) :: rboundary
  
  ! Spatial coarse mesh. The result is a raw mesh.
  type(t_triangulation), intent(in), target :: rtriaCoarse
  
  ! Output level during initialisation
  integer, intent(in) :: ioutputLevel
!</input>

!<output>
  ! Refinement specification.
  type(t_meshHierarchy), intent(out) :: rmeshHierarchy
!</output>

!</subroutine>

    ! Create the hierarchy.
    if (ioutputLevel .ge. 1) then
      call output_line ("Pre-refinement to level "//&
          trim(sys_siL(rrefinement%npreref+1,10))//".",bnolinebreak=.true.)
    end if

    call mshh_initHierarchy (rmeshHierarchy,rtriaCoarse,&
        rrefinement%npreref,rrefinement%nlevels,rboundary)
        
    if (ioutputLevel .ge. 1) then
      call output_line (" Done. Creating Hierarchy-Level: [1",bnolinebreak=.true.,&
          cdateTimeLogPolicy = OU_DTP_NONE)
    end if

    ! Refine the coarse mesh.
    call mshh_refineHierarchy2lv (rmeshHierarchy,rrefinement%nlevels,&
        rboundary=rboundary,bprint=ioutputLevel .ge. 1)

    if (ioutputLevel .ge. 1) then
      call output_line ("]",cdateTimeLogPolicy = OU_DTP_NONE)
    end if

    if (ioutputLevel .ge. 2) then
      call output_lbrk ()
      call output_line ('Mesh hierarchy statistics:')
      call output_line ("--------------------------")
      call mshh_printHierStatistics (rmeshHierarchy)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initTimeHierarchy (rrefinement,rtimeCoarse,rtimeHierarchy,ioutputLevel)
  
!<description>
  ! Initialises a space hierarchy based on a coarse mesh and a refinement
  ! strategy specified in rrefinement.
!</description>

!<input>
  ! Type of refinement
  type(t_settings_refinement), intent(in) :: rrefinement

  ! Spatial coarse mesh. The result is a raw mesh.
  type(t_timeDiscretisation), intent(in), target :: rtimeCoarse

  ! Output level during initialisation
  integer, intent(in) :: ioutputLevel
!</input>

!<output>
  ! Hierarchy of meshes in space.
  type(t_timescaleHierarchy), intent(out) :: rtimeHierarchy
!</output>

!</subroutine>

    if (ioutputLevel .ge. 1) then
      call output_line ("Creating Level: [1",bnolinebreak=.true.)
    end if

    ! Create the hierarchy.
    call tmsh_createHierarchy (rtimeCoarse,rtimeHierarchy,&
        rrefinement%npreref,rrefinement%nlevels,ioutputLevel .ge. 1)

    if (ioutputLevel .ge. 1) then
      call output_line ("]",cdateTimeLogPolicy = OU_DTP_NONE)
    end if

    if (ioutputLevel .ge. 2) then
      call output_lbrk ()
      call output_line ('Time hierarchy statistics:')
      call output_line ("--------------------------")
      call tmsh_printHierStatistics (rtimeHierarchy)
    end if

  end subroutine
  
  ! *****<**********************************************************************

!<subroutine>
  
  subroutine fgetDist1LvDiscr(ilevel,rtriangulation,rspaceDiscr,rboundary,rcollection)

!<description>
  ! Callback routine wrapper for spdsc_get1LevelDiscrNavSt2D. Allows to create
  ! a discretisation with the routines from fespacehierarchy.
  ! Expects the following information in rcollection:
  !   rcollection%IquickAccess(1) = ieltype
!</description>

  integer, intent(in) :: ilevel
  type(t_triangulation), intent(in) :: rtriangulation
  type(t_blockDiscretisation), intent(out) :: rspaceDiscr
  type(t_collection), intent(inout), optional :: rcollection
  type(t_boundary), intent(in), optional :: rboundary
    
!</subroutine>
    
    ! local variables
    type(t_spaceDiscrParams) :: rdiscrParams

    ! Recover the discretisation structure
    rdiscrParams = transfer(rcollection%DquickAccess(:),rdiscrParams)
    
    ! Which space is to be created?
    
    select case (rdiscrParams%cspace)
    
    ! --------------------
    ! Primal space
    ! --------------------
    case (CCSPACE_PRIMAL)
      call kktsp_initPrimalSpaceDiscr (rspaceDiscr,&
          rdiscrParams%p_rphysics,rdiscrParams%p_rsettingsSpaceDiscr,rtriangulation,rboundary)

    ! --------------------
    ! Dual space
    ! --------------------
    case (CCSPACE_DUAL)
      call kktsp_initDualSpaceDiscr (rspaceDiscr,&
          rdiscrParams%p_rfeHierarchyPrimal%p_rfeSpaces(ilevel)%p_rdiscretisation,&
          rdiscrParams%p_rphysics,rdiscrParams%p_roptControl)

    ! --------------------
    ! Control space
    ! --------------------
    case (CCSPACE_CONTROL)
      call kktsp_initControlSpaceDiscr (rspaceDiscr,&
          rdiscrParams%p_rfeHierarchyPrimal%p_rfeSpaces(ilevel)%p_rdiscretisation,&
          rdiscrParams%p_rsettingsSpaceDiscr,rdiscrParams%p_rphysics,&
          rdiscrParams%p_roptControl)
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceDiscrHier (&
      rfeHierarchyPrimal,rfeHierarchyDual,rfeHierarchyControl,&
      rphysics,roptControl,rsettingsDiscr,rmeshHierarchy,rboundary,&
      ioutputLevel)

!<description>
  ! Initialises the hierarchies for the space discretisation on all levels
  ! based on the parameters in the parameter list.
!</description>

!<input>
  ! Physics of the problem
  type(t_settings_physics), intent(in), target :: rphysics
  
  ! Structure defining the optimal control problem to calculate
  type(t_settings_optcontrol), intent(in), target :: roptControl

  ! Structure with space discretisation settings
  type(t_settings_spacediscr), intent(in), target :: rsettingsDiscr

  ! A mesh hierarchy with all available space meshes.
  type(t_meshHierarchy), intent(in), target :: rmeshHierarchy

  ! Description of the boundary
  type(t_boundary), intent(in), target :: rboundary

  ! Output level during initialisation
  integer, intent(in) :: ioutputLevel
!</input>

!<output>
  ! The FE space hierarch structure to create, primal space
  type(t_feHierarchy), intent(out), target :: rfeHierarchyPrimal

  ! The FE space hierarch structure to create, dual space
  type(t_feHierarchy), intent(out) :: rfeHierarchyDual

  ! The FE space hierarch structure to create, control space
  type(t_feHierarchy), intent(out) :: rfeHierarchyControl
!</output>

!</subroutine>
  
    ! local variables
    type(t_collection) :: rcollection
    type(t_spaceDiscrParams) :: rdiscrParams
    
    ! Read the parameters that define the underlying discretisation.
    ! We use fgetDist1LvDiscr to create the basic spaces.
    ! This routines expects the rcollection%IquickAccess array to be initialised
    ! as follows:
    !   rcollection%IquickAccess(1) = ieltype

    ! Initialise the structure encapsuling the parameters of the discretisation.
    ! The structure is placed in the collection

    rdiscrParams%p_rphysics => rphysics
    rdiscrParams%p_roptControl => roptControl
    rdiscrParams%p_rsettingsSpaceDiscr => rsettingsDiscr
    
    ! Create the hierarchy for the primal space
    rdiscrParams%cspace = CCSPACE_PRIMAL
    rdiscrParams%p_rfeHierarchyPrimal => null()
    rcollection%DquickAccess(:) = transfer(rdiscrParams,rcollection%DquickAccess(:))
    
    call fesph_createHierarchy (rfeHierarchyPrimal,&
        rmeshHierarchy%nlevels,rmeshHierarchy,&
        fgetDist1LvDiscr,rcollection,rboundary)

    ! Create the hierarchy for the dual space

    rdiscrParams%cspace = CCSPACE_DUAL
    rdiscrParams%p_rfeHierarchyPrimal => rfeHierarchyPrimal
    rcollection%DquickAccess(:) = transfer(rdiscrParams,rcollection%DquickAccess(:))

    call fesph_createHierarchy (rfeHierarchyDual,&
        rmeshHierarchy%nlevels,rmeshHierarchy,&
        fgetDist1LvDiscr,rcollection,rboundary)

    ! Create the hierarchy for the control space

    rdiscrParams%cspace = CCSPACE_CONTROL
    rdiscrParams%p_rfeHierarchyPrimal => rfeHierarchyPrimal
    rcollection%DquickAccess(:) = transfer(rdiscrParams,rcollection%DquickAccess(:))

    call fesph_createHierarchy (rfeHierarchyControl,&
        rmeshHierarchy%nlevels,rmeshHierarchy,&
        fgetDist1LvDiscr,rcollection,rboundary)
          
    if (ioutputLevel .ge. 2) then
      ! Print statistics about the discretisation
      call output_lbrk ()
      call output_line ("Space discretisation hierarchy statistics")
      call output_line ("-----------------------------------------")
      call output_line ("Primal space:")
      call output_lbrk ()
      call fesph_printHierStatistics (rfeHierarchyPrimal)
      call output_lbrk ()
      call output_line ("Dual space:")
      call output_lbrk ()
      call fesph_printHierStatistics (rfeHierarchyDual)
      call output_lbrk ()
      call output_line ("Control space:")
      call output_lbrk ()
      call fesph_printHierStatistics (rfeHierarchyControl)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_doneSpaceDiscrHier (rfeHierarchyPrimal,rfeHierarchyDual,rfeHierarchyControl)
  
!<description>
  ! Cleans up the discretisation hierarchies in rsettings.
!</description>

!<output>
  ! The FE space hierarch structure to clean up, primal space
  type(t_feHierarchy), intent(inout) :: rfeHierarchyPrimal

  ! The FE space hierarch structure to clean up, dual space
  type(t_feHierarchy), intent(inout) :: rfeHierarchyDual

  ! The FE space hierarch structure to clean up, control space
  type(t_feHierarchy), intent(inout) :: rfeHierarchyControl
!</output>

!</subroutine>

    ! Release all discretisation hierarchies.
    call fesph_releaseHierarchy(rfeHierarchyControl)
    call fesph_releaseHierarchy(rfeHierarchyDual)
    call fesph_releaseHierarchy(rfeHierarchyPrimal)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceTimeHierarchy (rparlist,ssection,&
      rrefinementSpace,rrefinementTime,&
      rfeHierPrimal,rfeHierDual,rfeHierControl,rtimeHierarchy,&
      rspaceTimeHierPrimal,rspaceTimeHierDual,rspaceTimeHierControl)
  
!<description>
  ! Creates a space-time hierarchy based on the parameters in the parameter list.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the parameters how to set up the
  ! space-time hierarchies.
  character(len=*), intent(in) :: ssection

  ! Settings that define the refinement in space
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! Settings that define the refinement in time
  type(t_settings_refinement), intent(in) :: rrefinementTime

  ! A mesh hierarchy with all available space meshes, primal space.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal

  ! A mesh hierarchy with all available space meshes, dual space.
  type(t_feHierarchy), intent(in) :: rfeHierDual

  ! A mesh hierarchy with all available space meshes, control space.
  type(t_feHierarchy), intent(in) :: rfeHierControl

  ! A hierarchy of time levels
  type(t_timescaleHierarchy), intent(in) :: rtimeHierarchy
!</input>

!<inputoutput>
  ! A space-time hierarchy based on the primal space
  type(t_spaceTimeHierarchy), intent(out) :: rspaceTimeHierPrimal
  
  ! A space-time hierarchy based on the dual space
  type(t_spaceTimeHierarchy), intent(out) :: rspaceTimeHierDual

  ! A space-time hierarchy based on the control space
  type(t_spaceTimeHierarchy), intent(out) :: rspaceTimeHierControl
!</inputoutput>

!</subroutine>

    integer :: ispacelevelcoupledtotimelevel,nmaxSimulRefLevel
    real(DP) :: dspacetimeRefFactor
    
    ! Get the parameters from the parameter list.
    !
    ! Should we couple space and time coarsening/refinement?
    call parlst_getvalue_int (rparlist,ssection,&
        'ispacelevelcoupledtotimelevel',ispacelevelcoupledtotimelevel,1)

    ! Parameters of the refinement?
    call parlst_getvalue_int (rparlist,ssection,&
        'nmaxSimulRefLevel',nmaxSimulRefLevel,0)
    
    if (nmaxSimulRefLevel .le. 0) then
      nmaxSimulRefLevel = rrefinementTime%nlevels+&
          nmaxSimulRefLevel-rrefinementTime%npreref
    end if
    nmaxSimulRefLevel = min(rrefinementTime%nlevels,max(1,nmaxSimulRefLevel))
        
    call parlst_getvalue_double (rparlist,ssection,&
        'dspacetimeRefFactor',dspacetimeRefFactor,1.0_DP)

    ! Create the hierarchies.
    call sth_initHierarchy (rspaceTimeHierPrimal,&
        rfeHierPrimal,rtimeHierarchy)

    call sth_initHierarchy (rspaceTimeHierDual,&
        rfeHierDual,rtimeHierarchy)

    call sth_initHierarchy (rspaceTimeHierControl,&
        rfeHierControl,rtimeHierarchy)
        
    select case (ispacelevelcoupledtotimelevel)
    case (0)
      ! Only in time, space level stays at max.
      dspacetimeRefFactor = SYS_INFINITY_DP
          
    case (1)
      ! Simultaneous refinement in space+time.

    case (2)
      ! Only in space, time level stays at max.
      dspacetimeRefFactor = 0.0_DP

    end select

    call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimal,&
        1,rrefinementSpace%nlevels,&
        1,rrefinementTime%nlevels,dspacetimeRefFactor)

    call sth_defineHierarchyByCoarsening (rspaceTimeHierDual,&
        1,rrefinementSpace%nlevels,&
        1,rrefinementTime%nlevels,dspacetimeRefFactor)

    call sth_defineHierarchyByCoarsening (rspaceTimeHierControl,&
        1,rrefinementSpace%nlevels,&
        1,rrefinementTime%nlevels,dspacetimeRefFactor)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpacePrjHierarchy (rprjHierarchy,rfehierarchy,&
      rphysics,rparlist,ssection)
  
!<description>
  ! Creates projection hierarchies for the interlevel projection in space.
!</description>
  
!<input>
  ! Underlying FEM hierarchy.
  type(t_feHierarchy), intent(in) :: rfehierarchy
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Parameter list with parameters about the projection.
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where parameters about the projection can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! The projection hierarchy to create.
  type(t_interlevelProjectionHier), intent(out) :: rprjHierarchy
!</output>

!</subroutine>

    ! local variables
    integer :: ilev

    ! Initialise a standard interlevel projection structure for every level
    call mlprj_initPrjHierarchy(rprjHierarchy,1,rfehierarchy%nlevels)
        
    do ilev=1,rfehierarchy%nlevels
    
      call mlprj_initPrjHierarchyLevel(rprjHierarchy,ilev,&
          rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation)

      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      if (ilev .gt. 1) then
        call getProlRest (rprjHierarchy%p_Rprojection(ilev), &
            rparlist, ssection)
      end if
          
    end do
    
    call mlprj_commitPrjHierarchy(rprjHierarchy)

  contains
  
    subroutine getProlRest (rprojection, rparamList, sname)
    
    ! Initialises an existing interlevel projection structure rprojection
    ! with parameters from the INI/DAT files. sname is the section in the
    ! parameter list containing parameters about prolongation restriction.

    ! Parameter list that contains the parameters from the INI/DAT file(s).
    type(t_parlist), intent(IN) :: rparamList
    
    ! Name of the section in the parameter list containing the parameters
    ! of the prolongation/restriction.
    character(LEN=*), intent(IN) :: sname

    ! An interlevel projection block structure containing an initial
    ! configuration of prolongation/restriction. The structure is modified
    ! according to the parameters in the INI/DAT file(s).
    type(t_interlevelProjectionBlock), intent(INOUT) :: rprojection

      ! local variables
      type(t_parlstSection), pointer :: p_rsection
      integer :: i1
      real(DP) :: d1

      ! Check that there is a section called sname - otherwise we
      ! cannot create anything!
      
      call parlst_querysection(rparamList, sname, p_rsection)

      if (.not. associated(p_rsection)) then
        ! We use the default configuration; stop here.
        return
      end if
      
      ! Now take a look which parameters appear in that section.

      select case (rphysics%cequation)
      
      ! ********************************************
      ! Stokes / Navier-Stokes
      ! ********************************************
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      
        ! Prolongation/restriction order for velocity components
        call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
        
        if (i1 .ne. -1) then
          ! Initialise order of prolongation/restriction for velocity components
          rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
          rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
          rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
        end if

        ! Prolongation/restriction order for pressure
        call parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
        
        if (i1 .ne. -1) then
          ! Initialise order of prolongation/restriction for pressure components
          rprojection%RscalarProjection(:,NDIM2D+1)%iprolongationOrder  = i1
          rprojection%RscalarProjection(:,NDIM2D+1)%irestrictionOrder   = i1
          rprojection%RscalarProjection(:,NDIM2D+1)%iinterpolationOrder = i1
        end if
        
        ! Prolongation/restriction variant for velocity components
        ! in case of Q1~ discretisation
        call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
        
        if (i1 .ne. -1) then
          rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
          rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1
        end if
        
        ! Aspect-ratio indicator in case of Q1~ discretisation
        ! with extended prolongation/restriction
        call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
        
        if (i1 .ne. 1) then
          rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
          rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
        end if

        ! Aspect-ratio bound for switching to constant prolongation/restriction
        ! in case of Q1~ discretisation with extended prolongation/restriction
        call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
        
        if (d1 .ne. 20.0_DP) then
          rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
          rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
        end if
        
      end select

    end subroutine

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceTimePrjHierarchy (rprjHierarchyBlock,rhierarchy,&
      rprojHierarchySpace,rphysics,roptcontrol,cspace,rparlist,ssection)
  
!<description>
  ! Creates projection hierarchies for the interlevel projection in space.
!</description>
  
!<input>
  ! Underlying space-time hierarchy.
  type(t_spaceTimeHierarchy), intent(in) :: rhierarchy
  
  ! Projection hierarchy in space
  type(t_interlevelProjectionHier), intent(in), target :: rprojHierarchySpace
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Optimal control parameters
  type(t_settings_optcontrol), intent(in) :: roptcontrol

  ! Type of space, this projection is set up for.
  ! =CCSPACE_PRIMAL  : Primal space, forward in time
  ! =CCSPACE_DUAL    : Dual space, backward in time
  ! =CCSPACE_CONTROL : Control space
  integer, intent(in) :: cspace

  ! Parameter list with parameters about the projection.
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where parameters about the projection can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! The projection hierarchy to create.
  type(t_sptiProjHierarchyBlock), intent(out) :: rprjHierarchyBlock
!</output>

!</subroutine>

    ! local variables
    integer :: ctypeProjection

    ! Type of prolongation/restriction in time
    call parlst_getvalue_int (rparlist, ssection, &
        'ctypeProjection', ctypeProjection, -1)
     
    call sptipr_initProjectionBlock (rprjHierarchyBlock,rhierarchy,&
        rprojHierarchySpace,rphysics,roptcontrol,cspace,ctypeProjection)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_doneSpaceTimePrjHierarchy (rprjHierarchyBlock)
  
!<description>
  ! Cleans up projection hierarchies for the interlevel projection in space.
!</description>
  
!<inputoutput>
  ! The projection hierarchy to create.
  type(t_sptiProjHierarchyBlock), intent(inout) :: rprjHierarchyBlock
!</inputoutput>

!</subroutine>

    call sptipr_doneProjectionBlock (rprjHierarchyBlock)
    
  end subroutine

end module
