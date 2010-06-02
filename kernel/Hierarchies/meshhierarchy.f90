!##############################################################################
!# ****************************************************************************
!# <name> meshhierarchy </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module maintains a hierarchy of triangulations. 
!#
!# The following routines can be found here:
!#
!# 1.) mshh_initHierarchy
!#     -> Initialises a mesh hierarchy with a coarse mesh, read from a file,
!#        taken from the application or by extracting levels from another
!#        hierarchy.
!#
!# 2.) mshh_refineHierarchy2lv
!#     -> Refine the coarse mesh in a hierarchy and create a hierarchy
!#        by 2-level refinement
!#
!# 3.) mshh_releaseHierarchy
!#     -> Release a mesh hierarchy
!#
!# 4.) mshh_printHierStatistics
!#     -> Print statistics about the hierarchy.
!# </purpose>
!##############################################################################

module meshhierarchy

  use fsystem
  use genoutput
  use boundary
  use basicgeometry
  use triangulation
  
  implicit none
  
  private
  
  public :: t_meshHierarchy
  public :: mshh_initHierarchy
  public :: mshh_refineHierarchy2lv
  public :: mshh_releaseHierarchy
  public :: mshh_printHierStatistics
  
!<constants>

!<constantblock description = "Constants defining information in cflags about the hierarchy.">

  ! Hierarchy was created by 2-level refinement.
  integer(I32), parameter, public :: MSHH_FLG_2LVREF = 2_I32**0

  ! The coarse meshes share their point coordinates with the fine mesh.
  integer(I32), parameter, public :: MSHH_FLG_COMPRESSED = 2_I32**1
  
!</constantblock>

!<constantblock description = "Refinement flags">

  ! The coarse meshes share their point coordinates with the fine mesh.
  integer(I32), parameter, public :: MSHH_REF_SHAREDCOORDS = 2_I32**0

!</constantblock>


!</constants>


!<types>

!<typeblock>

  ! A mesh hierarchy structure that describes a hierarchy of triangulations.
  type t_meshHierarchy
  
    ! A flag that specifies which information in this structure.
    integer(I32) :: cflags = 0

    ! Reference to the underlying domain or NULL() if no domain is attached.
    type(t_boundary), pointer :: p_rboundary => null()

    ! Reference to the coarse mesh. 
    type(t_triangulation), pointer :: p_rcoarseMesh => null()

    ! Reference to the finest mesh. 
    type(t_triangulation), pointer :: p_rfineMesh => null()
    
    ! Current number of levels available in this structure.
    integer :: nlevels = 0

    ! Maximum number of levels available in this structure.
    integer :: nmaxlevels = 0
    
    ! Level information.
    type(t_triangulation), dimension(:), pointer :: p_Rtriangulations => null()
    
  end type

!</typeblock>

!</types>

  interface mshh_initHierarchy
    module procedure mshh_initHierarchyTria
    module procedure mshh_initHierarchyDup
    module procedure mshh_initHierarchyFromFile
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine mshh_refineHierarchy2lv (rmeshHierarchy,nlevels,&
      cflags,ctriaFlags,rboundary,bprint)

!<description>
  ! Refines the coarse mesh in a level hierarchy by 2-level refinement to
  ! create the fine meshes up to level nlevels.
  !
  ! The routine can also be applied to a pre-refined mesh hierarchy.
  ! In this case, it will only create the remaining meshes.
!</description>

!<input>
  ! Number of refinement levels.
  integer, intent(in) :: nlevels

  ! OPTIONAL: Flags that specify additional information about the refinement.
  ! One of the MSHH_REF_xxxx constants. 
  ! If not specified, MSHH_REF_SHAREDCOORDS is assumed.
  integer(I32), intent(in), optional :: cflags

  ! OPTIONAL: Bitfield of TRIA_R2LV_xxxx constants that allow to specify
  ! options for the refinement. If not specified, TRIA_R2LV_STANDARD
  ! is used as default.
  integer(I32), intent(in), optional :: ctriaFlags

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary
  
  ! OPTIONAL: Whether to print the current state of the assembly to
  ! the terminal. If set to TRUE, numbers "2 3 4..:" will be printed
  ! to the terminal for every level currently being processed.
  logical, intent(in), optional :: bprint
!</input>

!<inputoutput>
  ! The mesh space hierarchy structure to refine.
  type(t_meshHierarchy), intent(inout), target :: rmeshHierarchy
!</inputoutput>
  
!</subroutine>
  
    integer :: i
    logical :: boutput
    integer(I32) :: clocalflags,clocalTriaFlags
    type(t_boundary), pointer :: p_rboundary
    
    if ((rmeshHierarchy%nmaxLevels .eq. 0) .or. (rmeshHierarchy%nlevels .eq. 0)) then
      call output_line ('Hierarchy not initialised.', &
          OU_CLASS_ERROR,OU_MODE_STD,'mshh_initHierarchyFromFile')
      call sys_halt()
    end if    
    
    if (nlevels .le. 0) return
    
    boutput = .false.
    if (present(bprint)) boutput = bprint
    
    p_rboundary => rmeshHierarchy%p_rboundary
    if (present(rboundary)) p_rboundary => rboundary
    
    clocalflags = MSHH_REF_SHAREDCOORDS
    if (present(cflags)) clocalflags = cflags
    
    clocalTriaFlags = TRIA_R2LV_STANDARD
    if (present(ctriaFlags)) clocalTriaFlags = ctriaFlags
    
    ! Indicate the 2-level refinement if the mesh is not yet refined.
    ! Otherwise, we cannot guarantee that the previous levels are refined
    ! by 2-level refinement as well. But if that is the case, the flag
    ! is already set.
    if (rmeshHierarchy%nlevels .eq. 1) then
      rmeshHierarchy%cflags = ior(rmeshHierarchy%cflags,MSHH_FLG_2LVREF)
    
      if (iand(clocalflags,MSHH_REF_SHAREDCOORDS) .ne. 0) then
        ! Indicate the shared coordinates
        rmeshHierarchy%cflags = ior(rmeshHierarchy%cflags,MSHH_FLG_COMPRESSED)
      end if
    end if
    
    ! Initialise the FE spaces on all levels. Refine the mesh by 2-level 
    ! refinement to create all missing levels.
    do i=rmeshHierarchy%nlevels+1,max(rmeshHierarchy%nmaxlevels,nlevels)
    
      if (present(bprint)) then
        if (bprint) then
          ! Print current state.
          call output_line (" "//trim(sys_siL(i,10)),bnolinebreak=.true.,cdateTimeLogPolicy=OU_DTP_NONE)
        end if
      end if
    
      ! Refine.
      call tria_refine2LevelOrdering (rmeshHierarchy%p_Rtriangulations(i-1),&
          rmeshHierarchy%p_Rtriangulations(i),p_rboundary,clocalTriaFlags)
      call tria_initStandardMeshFromRaw (&
          rmeshHierarchy%p_Rtriangulations(i),p_rboundary)
      
      ! Probably share the coordinates of the coarse mesh with the fine mesh.
      if (iand(clocalflags,MSHH_REF_SHAREDCOORDS) .ne. 0) then
        call tria_compress2LevelOrdHierarchy (&
            rmeshHierarchy%p_Rtriangulations(i),&
            rmeshHierarchy%p_Rtriangulations(i-1))
      end if
    end do
    
    ! Now, we have nlevels levels available.
    rmeshHierarchy%nlevels = max(rmeshHierarchy%nmaxlevels,nlevels)

    rmeshHierarchy%p_rfineMesh => rmeshHierarchy%p_Rtriangulations(rmeshHierarchy%nlevels)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mshh_initHierarchyTria (rmeshHierarchy,rtriangulation,npreref,nmaxLevels,&
      rboundary,cdupFlag)

!<description>
  ! Initialises a mesh hierarchy based on rtriangulation being
  ! the coarse mesh.
  ! The hierarchy is set up to hold up to nmaxLevels levels. No refinement is done.
!</description>

!<input>
  ! A triangulation structure that defines the coarse mesh.
  type(t_triangulation), intent(in), target :: rtriangulation

  ! Number of pre-refinements to be applied to rtriangulation to get the coarse mesh.
  integer, intent(in) :: npreref
  
  ! Maximum number of refinement levels in the structure
  integer, intent(in) :: nmaxLevels

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary

  ! OPTIONAL: Duplication flag. Defines how data is copied from 
  ! rtriangulation to the coarse mesh in rmeshHierarchy. One of the TR_SHARE_xxxx
  ! constants. If not specified, TR_SHARE_ALL is assumed, i.e. the coarse
  ! mesh shares all information with rtriangulation.
  integer(I32), intent(in), optional :: cdupFlag
!</input>

!<output>
  ! The mesh space hierarchy structure to create
  type(t_meshHierarchy), intent(out) :: rmeshHierarchy
!</output>
  
!</subroutine>
  
    if (nmaxLevels .le. 0) return
    
    ! Initialise the structure
    rmeshHierarchy%cflags = 0
    rmeshHierarchy%nmaxLevels = nmaxLevels
    if (present(rboundary)) then
      rmeshHierarchy%p_rboundary => rboundary
    end if
    allocate(rmeshHierarchy%p_Rtriangulations(rmeshHierarchy%nmaxLevels))
    
    ! The coarse mesh pointer points to the first level.
    rmeshHierarchy%nlevels = 1
    rmeshHierarchy%p_rcoarseMesh => rmeshHierarchy%p_Rtriangulations(1)
    rmeshHierarchy%p_rfineMesh => rmeshHierarchy%p_Rtriangulations(1)
    
    ! Create the first level. Duplicate the coarse mesh, share all information.
    if (present(cdupFlag)) then
      call tria_duplicate (rtriangulation, &
          rmeshHierarchy%p_Rtriangulations(1),cdupFlag,.false.)
    else
      call tria_duplicate (rtriangulation, &
          rmeshHierarchy%p_Rtriangulations(1),TR_SHARE_ALL,.false.)
    end if
    
    ! Pre-refine.
    if (npreref .gt. 0) then
      call tria_quickRefine2LevelOrdering(npreref,&
          rmeshHierarchy%p_Rtriangulations(1),rmeshHierarchy%p_rboundary)
      call tria_initStandardMeshFromRaw (&
          rmeshHierarchy%p_Rtriangulations(1),rmeshHierarchy%p_rboundary)
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mshh_initHierarchyDup (rmeshHierarchySource,rmeshHierarchyDest,&
      nminLevel, nmaxLevel, nmaxLevels, cdupFlag)

!<description>
  ! Creates a mesh hierarchy by copying another mesh hierarchy.
!</description>

!<input>
  ! A source mesh hierarchy.
  type(t_meshHierarchy), intent(in) :: rmeshHierarchySource
  
  ! OPTIONAL: The minimum level from rmeshHierarchySource that should be
  ! used as coarse mesh in rmeshHierarchyDest. 
  ! If not specified, this defaults to 1.
  integer, intent(in), optional :: nminLevel

  ! OPTIONAL: The maximum level from rmeshHierarchySource that should be
  ! used as fine mesh in rmeshHierarchyDest.
  ! If not specified, this defaults to the number of levels in rmeshHierarchySource.
  integer, intent(in), optional :: nmaxLevel

  ! OPTIONAL: Total maximum number of levels, rmeshHierarchyDest should support.
  ! If not specified, this defaults to the number of meshes to copy from
  ! rmeshHierarchySource.
  integer, intent(in), optional :: nmaxLevels
  
  ! OPTIONAL: Duplication flag. Defines how data is copied from 
  ! rmeshHierarchySource to rmeshHierarchyDest. One of the TR_SHARE_xxxx
  ! constants. If not specified, TR_SHARE_ALL is assumed.
  integer(I32), intent(in), optional :: cdupFlag
!</input>

!<output>
  ! The FE space hierarchy structure to create
  type(t_meshHierarchy), intent(out), target :: rmeshHierarchyDest
!</output>
  
!</subroutine>

    integer :: ntotalLevels,i,nmin,nmax
    integer(I32) :: cdup
  
    cdup = TR_SHARE_ALL
    if (present(cdupFlag)) cdup = cdupFlag
  
    ! Copy data.
    rmeshHierarchyDest = rmeshHierarchySource
  
    nmin = 1
    if (present(nminLevel)) nmin = max(nmin,nminLevel)
    nmin = min(nmin,rmeshHierarchySource%nlevels)
    
    nmax = rmeshHierarchySource%nlevels
    if (present(nmaxLevel)) nmax = min(nmax,nmaxLevel)
    nmax = max(nmax,1)
    
    ! How much levels should we allocate?
    ntotalLevels = max(1,nmax-nmin+1)
    if (present(nmaxLevels)) ntotalLevels = max(nmaxLevels,ntotalLevels)
    
    ! Allocate local data
    rmeshHierarchyDest%nmaxLevels = ntotalLevels
    rmeshHierarchyDest%nlevels = nmax-nmin+1
    allocate(rmeshHierarchyDest%p_Rtriangulations(rmeshHierarchyDest%nmaxLevels))
    
    ! Copy mesh data of the highest mesh as defined by cdupFlag.
    call tria_duplicate (&
        rmeshHierarchySource%p_Rtriangulations(nmax),&
        rmeshHierarchyDest%p_Rtriangulations(rmeshHierarchyDest%nlevels),cdup)
        
    ! Copy mesh data of the lower levels. 
    do i=rmeshHierarchyDest%nlevels-1,1,-1
    
      if (iand(rmeshHierarchySource%cflags,MSHH_FLG_2LVREF+MSHH_FLG_COMPRESSED) &
          .eq. MSHH_FLG_2LVREF+MSHH_FLG_COMPRESSED) then
      
        ! If we have a compressed 2-level hierarchy, don't duplicate the
        ! coordinates in any case. They are replaced anyway.
        call tria_duplicate (&
            rmeshHierarchySource%p_Rtriangulations(nmax-(rmeshHierarchyDest%nlevels-i)),&
            rmeshHierarchyDest%p_Rtriangulations(i),&
            ior(cdup,TR_SHARE_DVERTEXCOORDS))
            
        ! Re-compress the hierarchy.
        call tria_compress2LevelOrdHierarchy (&
            rmeshHierarchyDest%p_Rtriangulations(i+1),&
            rmeshHierarchyDest%p_Rtriangulations(i))
      else
      
        ! Just duplicate.
        call tria_duplicate (&
            rmeshHierarchySource%p_Rtriangulations(nmax-(rmeshHierarchyDest%nlevels-i)),&
            rmeshHierarchyDest%p_Rtriangulations(i),&
            cdup)

      end if
      
    end do
    
    ! The coarse mesh pointer points to the first level.
    rmeshHierarchyDest%p_rcoarseMesh => rmeshHierarchyDest%p_Rtriangulations(1)
    rmeshHierarchyDest%p_rfineMesh => &
        rmeshHierarchyDest%p_Rtriangulations(rmeshHierarchyDest%nlevels)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mshh_initHierarchyFromFile (rmeshHierarchy,nmaxLevels,sTRIfile,ndim,&
      npreref,rboundary,ctriaFlags)

!<description>
  ! Initialises a mesh hierarchy by reading in a .TRI file and use it
  ! as the coarse mesh. The hierarchy is set up to hold up to nmaxLevels
  ! levels. No refinement is done.
!</description>

!<input>
  ! Maximum number of refinement levels in the hierarchy.
  integer, intent(in) :: nmaxLevels
  
  ! Name/path of a triangulation file that specifies the coarse mesh.
  character(len=*), intent(in) :: sTRIFile

  ! Dimension of the triangulation. NDIM2D for 2D, NDIM3D for 3D.
  integer, intent(in) :: ndim

  ! OPTIONAL: Number of prerefinements to get the coarse mesh. Standard is =0.
  integer, intent(in), optional :: npreref

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary

  ! OPTIONAL: Bitfield of TRIA_R2LV_xxxx constants that allow to specify
  ! options for the refinement. If not specified, TRIA_R2LV_STANDARD
  ! is used as default.
  integer(I32), intent(in), optional :: ctriaFlags
!</input>

!<output>
  ! The FE space hierarchy structure to create
  type(t_meshHierarchy), intent(out), target :: rmeshHierarchy
!</output>
  
!</subroutine>
  
    if (nmaxLevels .le. 0) return
    
    ! Initialise the structure
    rmeshHierarchy%cflags = 0
    rmeshHierarchy%nlevels = 1
    rmeshHierarchy%nmaxLevels = nmaxLevels
    if (present(rboundary)) then
      rmeshHierarchy%p_rboundary => rboundary
    end if
    allocate(rmeshHierarchy%p_Rtriangulations(rmeshHierarchy%nmaxLevels))

    ! The coarse mesh pointer points to the first level.
    rmeshHierarchy%p_rcoarseMesh => rmeshHierarchy%p_Rtriangulations(1)
    rmeshHierarchy%p_rfineMesh => rmeshHierarchy%p_Rtriangulations(1)

    ! Read the file and pre-refine.
    select case (ndim)
    case (NDIM2D)
      call tria_readTriFile2D (rmeshHierarchy%p_rcoarseMesh, sTRIFile, rboundary)
      
    case default
      call output_line ('Dimension not supported.', &
          OU_CLASS_ERROR,OU_MODE_STD,'mshh_initHierarchyFromFile')
      call sys_halt()
    end select
    
    call tria_quickRefine2LevelOrdering(npreref,rmeshHierarchy%p_rcoarseMesh,&
        rboundary,ctriaflags)
    call tria_initStandardMeshFromRaw(rmeshHierarchy%p_rcoarseMesh, rboundary)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mshh_releaseHierarchy (rmeshHierarchy)

!<description>
  ! Releases a FE hierarchy.
!</description>

!<inputoutput>
  ! The FE space hierarchy structure to release.
  type(t_meshHierarchy), intent(inout) :: rmeshHierarchy
!</inputoutput>
  
!</subroutine>
  
    integer :: i
    
    if (rmeshHierarchy%nmaxLevels .eq. 0)then
      ! Hierarchy not initialised.
      return
    end if    

    ! Release all levels
    do i=rmeshHierarchy%nlevels,1,-1
      call tria_done (rmeshHierarchy%p_Rtriangulations(i))
    end do
    
    ! Clean up the rest
    deallocate (rmeshHierarchy%p_Rtriangulations)
    
    nullify(rmeshHierarchy%p_rcoarseMesh)
    nullify(rmeshHierarchy%p_rfineMesh)
    nullify(rmeshHierarchy%p_rboundary)
    
    rmeshHierarchy%nlevels = 0
    rmeshHierarchy%nmaxLevels = 0
    rmeshHierarchy%cflags = 0
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mshh_printHierStatistics (rmeshHierarchy,iminlevel)

!<description>
  ! Writes statistics about the mesh hierarchy to the terminal.
!</description>

!<input>
  ! The FE space hierarchy.
  type(t_meshHierarchy), intent(in) :: rmeshHierarchy
  
  ! Optional: Number of the minimum level in the hierarchy.
  integer, intent(in), optional :: iminlevel
!</input>
  
!</subroutine>
  
    integer :: i,iofs
    
    iofs = 0
    if (present(iminlevel)) iofs = iminlevel-1
  
    do i=1,rmeshHierarchy%nlevels
      call tria_infoStatistics (rmeshHierarchy%p_Rtriangulations(i),i .eq. 1,i+iofs)
    end do
    
  end subroutine

end module
