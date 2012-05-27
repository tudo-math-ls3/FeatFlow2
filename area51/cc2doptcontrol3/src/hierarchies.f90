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

module initsolver

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
  
  use structuresdiscretisation
  use structuresgeneral
  
  use spacediscretisation
  use spacetimeinterlevelprojection
  
  implicit none
  
  private
  
  public :: init_initParamTria
  public :: init_doneParamTria
  public :: init_initTimeHierarchy
  public :: init_initSpaceDiscrHier
  public :: init_initSpacePrjHierarchy
  public :: init_initSpaceTimePrjHierarchy
  public :: init_initSpaceTimeHierarchy
  
contains

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
  type(t_timeDiscretisation) :: rdiscrTime
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
      call output_lbrk ()
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
      call output_lbrk ()
      call tmsh_printHierStatistics (rtimeHierarchy)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceDiscrHier (rparlist,rfeHierarchy,&
      rphysics,rsettingsSpaceDiscr,rboundary,&
      rmeshHierarchy,ioutputLevel)

!<description>
  ! Initialises the hierarchies for the space discretisation on all levels
  ! based on the parameters in the parameter list.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Physics of the problem
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Structure with space discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr

  ! Description of the boundary
  type(t_boundary), intent(in) :: rboundary

  ! A mesh hierarchy with all available space meshes.
  type(t_meshHierarchy), intent(in) :: rmeshHierarchy

  ! Output level during initialisation
  integer, intent(in) :: ioutputLevel
!</input>

!<output>
  ! The FE space hierarch structure to create
  type(t_feHierarchy), intent(out) :: rfeHierarchy
!</output>

!</subroutine>
  
    ! local variables
    type(t_collection) :: rcollection
    
    select case (rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
      
      ! Read the parameters that define the underlying discretisation.
      ! We use fget1LevelDiscretisation to create the basic spaces.
      ! This routines expects the rcollection%IquickAccess array to be initialised
      ! as follows:
      !   rcollection%IquickAccess(1) = ieltype

      ! Element type
      rcollection%IquickAccess(1) = rsettingsSpaceDiscr%ielementType

      ! Create an FE space hierarchy based on the existing mesh hierarchy.
      call fesph_createHierarchy (rfeHierarchy,&
          rmeshHierarchy%nlevels,rmeshHierarchy,&
          fgetDist1LvDiscrNavSt2D,rcollection,rboundary)
          
    end select
    
    if (ioutputLevel .ge. 2) then
      ! Print statistics about the discretisation
      call output_lbrk ()
      call output_line ("Space discretisation hierarchy statistics:")
      call output_lbrk ()
      call fesph_printHierStatistics (rfeHierarchy)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_doneSpaceDiscrHier (rfeHierarchy)
  
!<description>
  ! Cleans up the discretisation hierarchies in rsettings.
!</description>

!<output>
  ! The FE space hierarch structure to release
  type(t_feHierarchy), intent(out) :: rfeHierarchy
!</output>

!</subroutine>

    ! Release all discretisation hierarchies.
    call fesph_releaseHierarchy(rfeHierarchy)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceTimeHierarchy (rparlist,ssection,&
      rrefinementSpace,rrefinementTime,&
      rfeHierPrimal,rfeHierPrimalDual,rtimeHierarchy,&
      rspaceTimeHierPrimal,rspaceTimeHierPrimalDual)
  
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

  ! A mesh hierarchy with all available space meshes, only primal space.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal
  
  ! A mesh hierarchy with all available space meshes, primal + dual space.
  type(t_feHierarchy), intent(in) :: rfeHierPrimalDual

  ! A hierarchy of time levels
  type(t_timescaleHierarchy), intent(in) :: rtimeHierarchy
!</input>

!<inputoutput>
  ! A space-time hierarchy based on the primal/dual space
  type(t_spaceTimeHierarchy), intent(out) :: rspaceTimeHierPrimal
  
  ! A space-time hierarchy based on the primal+dual space
  type(t_spaceTimeHierarchy), intent(out) :: rspaceTimeHierPrimalDual
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

    call sth_initHierarchy (rspaceTimeHierPrimalDual,&
        rfeHierPrimalDual,rtimeHierarchy)
        
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

    call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimalDual,&
        1,rrefinementSpace%nlevels,&
        1,rrefinementTime%nlevels,dspacetimeRefFactor)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpacePrjHierarchy (rprjHierarchy,rfehierarchy,rparlist,ssection)
  
!<description>
  ! Creates projection hierarchies for the interlevel projection in space.
!</description>
  
!<input>
  ! Underlying FEM hierarchy.
  type(t_feHierarchy), intent(in) :: rfehierarchy
  
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

    end subroutine

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceTimePrjHierarchy (rprjHierarchy,rtimeDisr,rhierarchy,&
      rprojHierarchySpace,rphysics,rparlist,ssection)
  
!<description>
  ! Creates projection hierarchies for the interlevel projection in space.
!</description>
  
!<input>
  ! Underlying time (coarse) grid.
  type(t_timeDiscretisation), intent(in) :: rtimeDisr

  ! Underlying space-time hierarchy.
  type(t_spaceTimeHierarchy), intent(in) :: rhierarchy
  
  ! Projection hierarchy in space
  type(t_interlevelProjectionHier), intent(in), target :: rprojHierarchySpace
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Parameter list with parameters about the projection.
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where parameters about the projection can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! The projection hierarchy to create.
  type(t_sptiProjHierarchy), intent(out) :: rprjHierarchy
!</output>

!</subroutine>

    ! local variables
    integer :: ctypeProjection

    ! Type of prolongation/restriction in time
    call parlst_getvalue_int (rparlist, ssection, &
        'ctypeProjection', ctypeProjection, -1)
        
    call sptipr_initProjection (rprjHierarchy,rhierarchy,&
        rprojHierarchySpace,rphysics,ctypeProjection)
    
  end subroutine

end module