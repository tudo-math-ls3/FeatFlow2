!##############################################################################
!# ****************************************************************************
!# <name> fespacehierarchy </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module maintains a hierarchy of FE spaces.
!#
!# A FE space is realised as a combination of triangulation and
!# discretisation. The problem is that an application may create a
!# discretisation/triangulation and might want to set up a FE space
!# without allocating new memory. The routines in this module provide the
!# functionality to automatically create and refine meshes and discretisation
!# structures when needed and recycling existing information as much as
!# possible:
!#  - If a FE space is created by a trangulation, the triangulation is reused
!#    and only refined if necessary. A new discretisation is created.
!#  - If a FE space is created by an existing discretisation, the 
!#    discretisation is reused and only refined if necessary.
!#  - If a FE space is created based on a set of existing discretisations/
!#    triangulations, the routines try to reuse the existing structures and
!#    only refine the mesh/create a new discretisation if they don't find
!#    a proper existing one.
!#
!# All creation routines accept a callback routine fgetdiscr which is called
!# in case a new discretisation must be created based on a mesh. The caller
!# has to provide the functionality to create the discretisation as it is
!# needed.
!#
!# The following routines can be found here. These add functionality
!# to the module fespacehierarchybase.f90:
!#
!# 1.) fesph_createFEspace
!#     -> Creates a FE space based on a triangulation, discretisation, TRI file,
!#        FE space hierarchy or a mesh hierarchy.
!#
!# 2.) fesph_releaseFEspace
!#     -> Releases a FE space
!#
!# 3.) fesph_createHierarchy
!#     -> Creates a FE space hierarchy based on a triangulation, 
!#        discretisation, TRI-file or FE space
!#
!# 4.) fesph_releaseHierarchy
!#     -> Releases a FE space hierarchy
!#
!# 5.) fesph_concatFeSpaces
!#     -> Concatenates two FE spaces
!#
!# 6.) fesph_deriveFeSpaces
!#     -> Extracts a set of FE spaces from an FE space
!#
!# 7.) fesph_concatFeHierarchies
!#     -> Concatenates two FE spaces on every level of a hierarchy
!#
!# 8.) fesph_deriveFeHierarchy
!#     -> Extracts a set of FE spaces from an FE space on every level of 
!#        a hierarchy
!#
!# 9.) fesph_infoStatistics
!#     -> Print statistics about an FE space
!#
!# 10.) fesph_printHierStatistics
!#     -> Print statistics about an FE space hierarchy
!# </purpose>
!##############################################################################

module fespacehierarchy

  use fsystem
  use genoutput
  use boundary
  use basicgeometry
  use triangulation
  use meshhierarchy
  use spatialdiscretisation
  use collection
  use dofmapping
  
  use fespacehierarchybase
  
  implicit none
  
  private
  
  public :: fesph_createFEspace
  public :: fesph_releaseFEspace
  public :: fesph_concatFeSpaces
  public :: fesph_deriveFeSpaces
  public :: fesph_createHierarchy
  public :: fesph_releaseHierarchy
  public :: fesph_concatFeHierarchies
  public :: fesph_deriveFeHierarchy
  public :: fesph_infoStatistics
  public :: fesph_printHierStatistics
  
  interface fesph_createFEspace
    module procedure fesph_createFEspaceTria
    module procedure fesph_createFEspaceDiscr
    module procedure fesph_createFEspaceMeshH
    module procedure fesph_createFEspaceHier
    module procedure fesph_createFEspaceRef
    module procedure fesph_createFEspaceFromFile
  end interface
    
  interface fesph_createHierarchy
    module procedure fesph_createHierarchyTria
    module procedure fesph_createHierarchyDiscr
    module procedure fesph_createHierarchyMeshH
    module procedure fesph_createHierarchyDup
  end interface

contains

! ***************************************************************************

!<subroutine>

  subroutine fesph_createFEspaceHier (rfeSpace,ilevel,rfeHierarchy,&
      fgetDiscr,rcollection,rboundary)

!<description>
  ! Creates a FE space.
  ! ilevel describes the refinement level of the underlying mesh.
  ! If possible, the solution vector uses a discretisation specified in
  ! the FE hierarchy rfeHierarchy.
  ! If ilevel specifies a level that does not exist in the hierarchy,
  ! the underlying mesh is automatically refined and a suitable 
  ! discretisation structure is automatically created.
!</description>
 
!<input>
  ! An existing FE hierarchy.
  type(t_feHierarchy), intent(inout) :: rfeHierarchy
  
  ! Refinement level. Level 1 corresponds to the coarse mesh in Rdiscr(1).
  ! A value <= 0 specifies a level relative to the maximum, resulting in
  ! level SIZE(Rdiscr)+ilevel.
  integer, intent(in) :: ilevel
  
  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary
  
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpace
!</output>
  
!</subroutine>
  
    integer :: iactlevel
    
    iactlevel = ilevel
    if (ilevel .le. 0) iactlevel = rfeHierarchy%nlevels-ilevel
    
    if (iactlevel .lt. 1) then
      call output_line ('Refinement level < 1.', &
          OU_CLASS_ERROR,OU_MODE_STD,'fesph_createFEspaceHier')
      call sys_halt()
    end if
  
    ! Is the level in the hierarchy?
    if (iactlevel .le. rfeHierarchy%nlevels) then
      ! Take the existing hierarchy and remember that we are a copy.
      rfeSpace = rfeHierarchy%p_rfeSpaces(iactlevel)
      rfeSpace%cflags = FESPH_SHAREDTRIA+FESPH_SHAREDDISCR
      rfeSpace%p_rboundary => rfeHierarchy%p_rboundary
    else
      ! Refine the maximum level and create a new FE space.
      call fesph_createFEspaceRef (rfeSpace,iactlevel,&
          rfeHierarchy%p_rfeSpaces(rfeHierarchy%nlevels),&
          rfeHierarchy%nlevels,fgetDiscr,rcollection)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createFEspaceMeshH (rfeSpace,ilevel,rmeshHierarchy,&
      fgetDiscr,rcollection,rboundary)
  
!<description>
  ! Creates a FE space based on a mesh hierarchy.
  ! ilevel describes the refinement level of the underlying mesh.
  ! If possible, the solution vector uses a triangulation specified in
  ! rmeshHerarchy. If ilevel specifies a level which cannot be found in 
  ! rmeshHerarchy, the underlying mesh is automatically refined and a suitable 
  ! discretisation structure is automatically created.
!</description>

!<input>  
  ! A mesh hierarchy
  type(t_meshHierarchy), intent(in) :: rmeshHierarchy
  
  ! Refinement level. Level 1 corresponds to the coarse mesh in rmeshHierarchy.
  ! A value <= 0 specifies a level relative to the maximum, resulting in
  ! level t_meshHierarchy%nlevels+ilevel.
  integer, intent(in) :: ilevel

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface

  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpace
!</output>
  
!</subroutine>
  
    integer :: iactlevel
    
    iactlevel = ilevel
    if (ilevel .le. 0) iactlevel = rmeshHierarchy%nlevels-ilevel
  
    if (iactlevel .lt. 1) then
      call output_line ('Refinement level < 1.', &
          OU_CLASS_ERROR,OU_MODE_STD,'fesph_createFEspaceMeshH')
      call sys_halt()
    end if
  
    if (iactlevel .le. rmeshHierarchy%nlevels) then
      ! If we are in the range of Rdiscr, simply take the existing discretisation
      ! Remember discretisation and trriangulation.
      rfeSpace%cflags = FESPH_SHAREDTRIA
      rfeSpace%p_rtriangulation => rmeshHierarchy%p_Rtriangulations(iactlevel)
    else
      ! Refine the given mesh and create a new discretisation based on the new mesh.
      ! This information is released later on together with the structure.
      rfeSpace%cflags = 0
      allocate(rfeSpace%p_rtriangulation)
      
      call tria_duplicate (rmeshHierarchy%p_rfineMesh, &
          rfeSpace%p_rtriangulation,TR_SHARE_ALL,.false.)
      
      call tria_quickRefine2LevelOrdering(iactlevel-rmeshHierarchy%nlevels,&
          rfeSpace%p_rtriangulation,rboundary)
      call tria_initStandardMeshFromRaw(rfeSpace%p_rtriangulation, rboundary)
    end if

    if (present(rboundary)) then
      rfeSpace%p_rboundary => rboundary
    end if

    allocate(rfeSpace%p_rdiscretisation)
    call fgetDiscr(rfeSpace%p_rtriangulation,&
        rfeSpace%p_rdiscretisation,rboundary,rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createFEspaceTria (rfeSpace,ilevelTri,&
      rtriangulation,ilevel,fgetDiscr,rcollection,rboundary)
  
!<description>
  ! Creates a FE space based on a triangulation.
  ! rtriangulation is a given triangulation corresponding to level ilevelTri.
  ! ilevel>=ilevelTri describes the refinement level of the underlying mesh.
  ! If possible, the solution vector uses a triangulation specified in
  ! rtriangulation. If ilevel specifies a level > ilevelTri
  ! the underlying mesh is automatically refined and a suitable 
  ! discretisation structure is automatically created.
!</description>

!<input>  
  ! A basic triangulation corresponding to level ilevelTri.
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! Refinement level. Must be >= ilevelTri.
  ! If ilevel > ilevelTri, the mesh is automatically refined up to this level.
  integer, intent(in) :: ilevel

  ! A level identifier that specifies the refinement level of rtriangulation.
  integer, intent(in) :: ilevelTri
  
  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface

  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpace
!</output>
  
!</subroutine>
  
    if (ilevel .lt. ilevelTri) then
      call output_line ('Refinement level < level of the triangulation.', &
          OU_CLASS_ERROR,OU_MODE_STD,'fesph_createFEspaceTria')
      call sys_halt()
    end if
  
    if (ilevel .eq. ilevelTri) then
      ! If the level corresponds to the existing triangulation, use it.
      rfeSpace%cflags = FESPH_SHAREDTRIA
      rfeSpace%p_rtriangulation => rtriangulation
    else
      ! Refine the given mesh and create a new discretisation based on the new mesh.
      ! This information is released later on together with the structure.
      rfeSpace%cflags = 0
      allocate(rfeSpace%p_rtriangulation)
      
      call tria_duplicate (rtriangulation, &
          rfeSpace%p_rtriangulation,TR_SHARE_ALL,.false.)
      
      call tria_quickRefine2LevelOrdering(ilevel-ilevelTri,&
          rfeSpace%p_rtriangulation,rboundary)
      call tria_initStandardMeshFromRaw(rfeSpace%p_rtriangulation, rboundary)
    end if

    if (present(rboundary)) then
      rfeSpace%p_rboundary => rboundary
    end if

    ! Create a new discretisation based on the new mesh.
    allocate(rfeSpace%p_rdiscretisation)
    call fgetDiscr(rfeSpace%p_rtriangulation,&
        rfeSpace%p_rdiscretisation,rboundary,rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createFEspaceDiscr (rfeSpace,ilevel,&
      rdiscretisation,ilevelDiscr,fgetDiscr,rcollection)
  
!<description>
  ! Creates a FE space based on a discretisation.
  ! rtriangulation is a given triangulation corresponding to level ilevelTri.
  ! ilevel>=ilevelTri describes the refinement level of the underlying mesh.
  ! If possible, the solution vector uses a triangulation specified in
  ! rtriangulation. If ilevel specifies a level > ilevelTri
  ! the underlying mesh is automatically refined and a suitable 
  ! discretisation structure is automatically created.
!</description>

!<input>  
  ! A basic discretisation corresponding to level ilevelDiscr.
  type(t_blockDiscretisation), intent(in), target :: rdiscretisation
  
  ! A level identifier that specifies the refinement level of rdiscretisation.
  integer, intent(in) :: ilevelDiscr
  
  ! Refinement level. Must be >= ilevelDiscr.
  ! If ilevel > ilevelDiscr, the mesh is automatically refined up to this level.
  integer, intent(in) :: ilevel

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface

  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpace
!</output>
  
!</subroutine>
  
    if (ilevel .lt. ilevelDiscr) then
      call output_line ('Refinement level < level of the triangulation.', &
          OU_CLASS_ERROR,OU_MODE_STD,'fesph_createFEspaceTria')
      call sys_halt()
    end if
  
    rfeSpace%p_rboundary => rdiscretisation%p_rboundary
    
    if (ilevel .eq. ilevelDiscr) then
      ! If the level corresponds to the existing triangulation, use it.
      rfeSpace%cflags = FESPH_SHAREDTRIA+FESPH_SHAREDDISCR
      rfeSpace%p_rtriangulation => rdiscretisation%p_rtriangulation
      rfeSpace%p_rdiscretisation => rdiscretisation
    else
      ! Refine the given mesh and create a new discretisation based on the new mesh.
      ! This information is released later on together with the structure.
      rfeSpace%cflags = 0
      allocate(rfeSpace%p_rtriangulation)
      
      call tria_duplicate (rdiscretisation%p_rtriangulation, &
          rfeSpace%p_rtriangulation,TR_SHARE_ALL,.false.)
      
      call tria_quickRefine2LevelOrdering(ilevel-ilevelDiscr,&
          rfeSpace%p_rtriangulation,rdiscretisation%p_rboundary)
      call tria_initStandardMeshFromRaw(rfeSpace%p_rtriangulation, &
          rdiscretisation%p_rboundary)

      ! Create a new discretisation based on the new mesh.
      allocate(rfeSpace%p_rdiscretisation)
      if (associated(rdiscretisation%p_rboundary)) then
        call fgetDiscr(rfeSpace%p_rtriangulation,&
            rfeSpace%p_rdiscretisation,rdiscretisation%p_rboundary,rcollection)
      else
        call fgetDiscr(rfeSpace%p_rtriangulation,&
            rfeSpace%p_rdiscretisation,rcollection=rcollection)
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createFEspaceRef (rfeSpace,ilevel,rfeSpaceTemplate,&
      ilevelTemplate,fgetDiscr,rcollection)
  
!<description>
  ! Creates a FE space by refinment of an existing FE space.
  ! rfeSpaceTemplate is a template FE space which is refined up to level ilevel,
  ! where ilevel=1 corresponds to the unrefined coarse mesh.
  !
  ! For ilevel=ilevelTemplate, the source structure is copied.
!</description>

!<input>  
  ! A template FE space structure used for creating rfeSpace.
  type(t_feSpaceLevel), intent(in) :: rfeSpaceTemplate

  ! Refinement level of the template FE-space rfeSpaceTemplate
  integer, intent(in) :: ilevelTemplate

  ! Refinement level. Level 1 corresponds to the coarse mesh in rmeshHierarchy.
  ! A value <= 0 specifies a level relative to the maximum, resulting in
  ! level rmeshHierarchy%nlevels%nlevels+ilevel.
  integer, intent(in) :: ilevel

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface

  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpace
!</output>
  
!</subroutine>
  
    if (ilevel .lt. ilevelTemplate) then
      call output_line ('Refinement level < 1.', &
          OU_CLASS_ERROR,OU_MODE_STD,'fesph_createFEspaceRef')
      call sys_halt()
    end if

    rfeSpace%p_rboundary => rfeSpaceTemplate%p_rboundary
    
    if (ilevel .eq. ilevelTemplate) then
      ! If the level corresponds to the existing one, use it.
      rfeSpace%cflags = FESPH_SHAREDDISCR+FESPH_SHAREDTRIA
      rfeSpace%p_rtriangulation => rfeSpaceTemplate%p_rtriangulation
      rfeSpace%p_rdiscretisation => rfeSpaceTemplate%p_rdiscretisation
    else
      ! Refine the given mesh and create a new discretisation based on the new mesh.
      ! This information is released later on together with the structure.
      rfeSpace%cflags = 0
      allocate(rfeSpace%p_rtriangulation)
      allocate(rfeSpace%p_rdiscretisation)
      
      call tria_duplicate (rfeSpaceTemplate%p_rtriangulation, &
          rfeSpace%p_rtriangulation,TR_SHARE_ALL,.false.)
      
      if (associated(rfeSpace%p_rboundary)) then
        call tria_quickRefine2LevelOrdering(ilevel-ilevelTemplate,&
            rfeSpace%p_rtriangulation,rfeSpace%p_rboundary)
        call tria_initStandardMeshFromRaw(rfeSpace%p_rtriangulation, &
            rfeSpace%p_rboundary)

        call fgetDiscr(rfeSpace%p_rtriangulation,&
            rfeSpace%p_rdiscretisation,rfeSpace%p_rboundary,rcollection)
      else
        call tria_quickRefine2LevelOrdering(ilevel-ilevelTemplate,&
            rfeSpace%p_rtriangulation)
        call tria_initStandardMeshFromRaw(rfeSpace%p_rtriangulation)

        call fgetDiscr(rfeSpace%p_rtriangulation,&
            rfeSpace%p_rdiscretisation,rcollection=rcollection)
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createFEspaceFromFile (rfeSpace,ilevel,sTRIFile,ndim,&
      fgetDiscr,rcollection,rboundary)

!<description>
  ! Creates a FE space based on a triangulation read from hard disc.
  ! ilevel describes the refinement level of the underlying mesh.
  ! sTRIFile is the name of the mesh the discretisation type.
!</description>
  
!<input>
  ! Name/path of a triangulation file that specifies the coarse mesh.
  character(len=*), intent(in) :: sTRIFile
  
  ! Dimension of the triangulation. NDIM2D for 2D, NDIM3D for 3D.
  integer, intent(in) :: ndim
  
  ! Refinement level. Level 1 corresponds to the coarse mesh specified by sTRIFile.
  ! Mus be >= 1.
  integer, intent(in) :: ilevel

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpace
!</output>
  
!</subroutine>
  
    if (ilevel .lt. 1) then
      call output_line ('Refinement level < 1.', &
          OU_CLASS_ERROR,OU_MODE_STD,'fesph_createFEspaceFromFile')
      call sys_halt()
    end if

    ! Read the file and refine/set up discretisation by hand.
    rfeSpace%cflags = 0
    allocate(rfeSpace%p_rtriangulation)
    allocate(rfeSpace%p_rdiscretisation)
    
    if (present(rboundary)) then
      rfeSpace%p_rboundary => rboundary
    end if
    
    select case (ndim)
    case (NDIM2D)
      call tria_readTriFile2D (rfeSpace%p_rtriangulation, sTRIFile, rboundary)
    case default
      call output_line ('Dimension not supported.', &
          OU_CLASS_ERROR,OU_MODE_STD,'fesph_createFEspaceFromFile')
      call sys_halt()
    end select
    
    call tria_quickRefine2LevelOrdering(ilevel-1,rfeSpace%p_rtriangulation,rboundary)
    call tria_initStandardMeshFromRaw(rfeSpace%p_rtriangulation, rboundary)
    
    call fgetDiscr(rfeSpace%p_rtriangulation,&
        rfeSpace%p_rdiscretisation,rboundary,rcollection)
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_concatFeSpaces (rfeSpace1,rfeSpace2,rfeSpace,rtriangulation)
  
!<description>
  ! Concatenates two FE space to a common FE space (cross product).
!</description>

!<input>  
  ! First FE space
  type(t_feSpaceLevel), intent(in) :: rfeSpace1

  ! Second FE space
  type(t_feSpaceLevel), intent(in) :: rfeSpace2

  ! OPTIONAL: Reference to the triangulation to use.
  ! If not specified, the triangulation in rfeSpace1 is used.
  type(t_triangulation), intent(in), target, optional :: rtriangulation
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpace
!</output>
  
!</subroutine>
  
    ! Copy the main data from the first space
    rfeSpace = rfeSpace1
    
    ! Probably use rtriangulation
    if (present(rtriangulation)) then
      rfeSpace%p_rtriangulation => rtriangulation
    end if
    
    ! Discretisation is new. Triangulation is old -- in any case.
    rfeSpace%cflags = iand(rfeSpace%cflags,not(FESPH_SHAREDDISCR))
    rfeSpace%cflags = ior(rfeSpace%cflags,FESPH_SHAREDTRIA)
    
    ! Create a new discretisation
    allocate(rfeSpace%p_rdiscretisation)
    call spdiscr_concatBlockDiscr(&
        rfeSpace1%p_rdiscretisation,rfeSpace2%p_rdiscretisation,&
        rfeSpace%p_rdiscretisation,rtriangulation)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_deriveFeSpaces (rfeSpaceSource,rfeSpaceDest,ifirstBlock,ilastBlock,&
      rtriangulation)
  
!<description>
  ! Creates an FE space as a subset of another FE space.
!</description>

!<input>  
  ! Source FE space
  type(t_feSpaceLevel), intent(in) :: rfeSpaceSource

  ! Number of the block in rfeSpaceSource that should be
  ! used as first block in rfeSpaceDest.
  integer, intent(in) :: ifirstBlock

  ! Number of the last block in rfeSpaceSource that should be
  ! used as last block in rfeSpaceDest.
  integer, intent(in) :: ilastBlock
  
  ! OPTIONAL: Reference to the triangulation to use.
  ! If not specified, the triangulation in rfeSpaceSource is used.
  type(t_triangulation), intent(in), target, optional :: rtriangulation
!</input>

!<output>
  ! A new FE space structure.
  type(t_feSpaceLevel), intent(out) :: rfeSpaceDest
!</output>
  
!</subroutine>
  
    ! Copy the main data from the first space
    rfeSpaceDest = rfeSpaceSource
    
    ! Probably use rtriangulation
    if (present(rtriangulation)) then
      rfeSpaceDest%p_rtriangulation => rtriangulation
    end if
    
    ! Discretisation is new. Triangulation is old -- in any case.
    rfeSpaceDest%cflags = iand(rfeSpaceDest%cflags,not(FESPH_SHAREDDISCR))
    rfeSpaceDest%cflags = ior(rfeSpaceDest%cflags,FESPH_SHAREDTRIA)
    
    ! Create a new discretisation
    allocate(rfeSpaceDest%p_rdiscretisation)
    call spdiscr_deriveBlockDiscr(rfeSpaceSource%p_rdiscretisation,&
        rfeSpaceDest%p_rdiscretisation,ifirstBlock,ilastBlock)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_releaseFEspace (rfeSpace)

!<description>
  ! Releases a FE space.
!</description>

!<inputoutput>
  ! The FE space structure to be released.
  type(t_feSpaceLevel), intent(inout) :: rfeSpace
!</inputoutput>
  
!</subroutine>
  
    ! Release all information
    if (iand(rfeSpace%cflags,FESPH_SHAREDDISCR) .eq. 0) then
      call spdiscr_releaseBlockDiscr (rfeSpace%p_rdiscretisation)
      deallocate(rfeSpace%p_rdiscretisation)
    else
      nullify(rfeSpace%p_rdiscretisation)
    end if

    if (iand(rfeSpace%cflags,FESPH_SHAREDTRIA) .eq. 0) then
      call tria_done (rfeSpace%p_rtriangulation)
      deallocate(rfeSpace%p_rtriangulation)
    else
      nullify(rfeSpace%p_rtriangulation)
    end if
    
    nullify(rfeSpace%p_rboundary)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createHierarchyTria (rfeHierarchy,nlevels,rtriangulation,&
      fgetDiscr,rcollection,rboundary,nmaxLevels)

!<description>
  ! Creates a FE space hierarchy based on a triangulation by refinement.
!</description>

!<input>
  ! A triangulation structure that defines the coarse mesh.
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! Number of refinement levels.
  integer, intent(in) :: nlevels

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary
  
  ! OPTIONAL: Maximum number of levels, the structure should support.
  integer, intent(in), optional :: nmaxLevels
!</input>

!<output>
  ! The FE space hierarch structure to create
  type(t_feHierarchy), intent(out) :: rfeHierarchy
!</output>
  
!</subroutine>
  
    integer :: i
    
    if (nlevels .le. 0) return
    
    ! Initialise the structure
    rfeHierarchy%cflags = 0
    rfeHierarchy%nlevels = nlevels
    if (present(nmaxLevels)) then
      rfeHierarchy%nmaxLevels = max(nlevels,nmaxLevels)
    else
      rfeHierarchy%nmaxLevels = nlevels
    end if
    if (present(rboundary)) then
      rfeHierarchy%p_rboundary => rboundary
    end if
    allocate(rfeHierarchy%p_rfeSpaces(rfeHierarchy%nmaxLevels))
    
    ! Initialise a mesh hierarchy corresponding to the FE hierarchy
    ! we are creating here.
    call mshh_initHierarchy(rfeHierarchy%rmeshHierarchy,&
        rtriangulation,0,rfeHierarchy%nmaxLevels,rboundary,TR_SHARE_ALL)
        
    ! Refine to get all spatial meshes.
    call mshh_refineHierarchy2lv (rfeHierarchy%rmeshHierarchy,nlevels)
    
    ! Initialise all the FE space levels.
    do i=1,nlevels
      call fesph_createFEspace (rfeHierarchy%p_rfeSpaces(i),i,&
          rfeHierarchy%rmeshHierarchy%p_Rtriangulations(i),i,&
          fgetDiscr,rcollection,rboundary)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createHierarchyDiscr (rfeHierarchy,nlevels,rdiscretisation,&
      fgetDiscr,rcollection,rboundary,nmaxLevels)

!<description>
  ! Creates a FE space hierarchy based on a discretisation by refinement.
!</description>

!<input>
  ! A discretisation structure that defines the coarse discretisation.
  type(t_blockDiscretisation), intent(in), target :: rdiscretisation
  
  ! Number of refinement levels.
  integer, intent(in) :: nlevels

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary

  ! OPTIONAL: Maximum number of levels, the structure should support.
  integer, intent(in), optional :: nmaxLevels
!</input>

!<output>
  ! The FE space hierarch structure to create
  type(t_feHierarchy), intent(out) :: rfeHierarchy
!</output>
  
!</subroutine>
  
    integer :: i
    
    if (nlevels .le. 0) return
    
    ! Initialise the structure
    rfeHierarchy%cflags = 0
    rfeHierarchy%nlevels = nlevels
    if (present(nmaxLevels)) then
      rfeHierarchy%nmaxLevels = max(nlevels,nmaxLevels)
    else
      rfeHierarchy%nmaxLevels = nlevels
    end if
    if (present(rboundary)) then
      rfeHierarchy%p_rboundary => rboundary
    end if
    allocate(rfeHierarchy%p_rfeSpaces(rfeHierarchy%nmaxLevels))
    
    ! Initialise a mesh hierarchy corresponding to the FE hierarchy
    ! we are creating here.
    call mshh_initHierarchy(rfeHierarchy%rmeshHierarchy,&
        rdiscretisation%p_rtriangulation,0,rfeHierarchy%nmaxLevels,&
        rboundary,TR_SHARE_ALL)
        
    ! Refine to get all spatial meshes.
    call mshh_refineHierarchy2lv (rfeHierarchy%rmeshHierarchy,nlevels)
    
    ! Initialise all the FE space levels.
    ! The first level can be a pointer to the existing discretisation.
    call fesph_createFEspaceDiscr (rfeHierarchy%p_rfeSpaces(i),1,&
      rdiscretisation,1,fgetDiscr,rcollection)
    do i=2,nlevels
      call fesph_createFEspace (rfeHierarchy%p_rfeSpaces(i),i,&
          rfeHierarchy%rmeshHierarchy%p_Rtriangulations(i),i,&
          fgetDiscr,rcollection,rboundary)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createHierarchyMeshH (rfeHierarchy,nlevels,rmeshHierarchy,&
      fgetDiscr,rcollection,rboundary,nmaxLevels)

!<description>
  ! Creates a FE space hierarchy based on a mesh hierarchy.
  ! If necessary, new levels are automatically created.
!</description>

!<input>
  ! An underlying mesh hierarchy.
  type(t_meshHierarchy), intent(in) :: rmeshHierarchy
  
  ! Number of refinement levels in the hierarchy. Must be >= 1.
  integer, intent(in) :: nlevels

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), target, optional :: rboundary

  ! OPTIONAL: Maximum number of levels, the structure should support.
  integer, intent(in), optional :: nmaxLevels
!</input>

!<output>
  ! The FE space hierarch structure to create
  type(t_feHierarchy), intent(out) :: rfeHierarchy
!</output>
  
!</subroutine>
  
    integer :: i
    
    if (nlevels .le. 0) return
    
    ! Initialise the structure
    rfeHierarchy%cflags = 0
    rfeHierarchy%nlevels = nlevels
    if (present(nmaxLevels)) then
      rfeHierarchy%nmaxLevels = max(nlevels,nmaxLevels)
    else
      rfeHierarchy%nmaxLevels = nlevels
    end if
    if (present(rboundary)) then
      rfeHierarchy%p_rboundary => rboundary
    end if
    allocate(rfeHierarchy%p_rfeSpaces(rfeHierarchy%nmaxLevels))

    ! Duplicate the mesh hierarchy.
    call mshh_initHierarchy(rmeshHierarchy,rfeHierarchy%rmeshHierarchy,&
        nmaxLevels=rfeHierarchy%nmaxLevels)
        
    ! Refine up to level nlevels.
    call mshh_refineHierarchy2lv (rfeHierarchy%rmeshHierarchy,nlevels,&
        rboundary=rboundary)
    
    ! Initialise all the FE space levels.
    do i=1,nlevels
      call fesph_createFEspace (rfeHierarchy%p_rfeSpaces(i),i,&
          rfeHierarchy%rmeshHierarchy%p_Rtriangulations(i),i,&
          fgetDiscr,rcollection,rboundary)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createHierarchyDup (rfeHierarchySource,rfeHierarchyDest,&
      nminLevel, nmaxLevel, nmaxLevels, fgetDiscr,rcollection)

!<description>
  ! Creates a FE space hierarchy based on am existing FE hierarchy.
  ! If necessary, new levels are automatically created.
!</description>

!<input>
  ! A source FE hierarchy.
  type(t_feHierarchy), intent(in) :: rfeHierarchySource
  
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
  
  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
  ! The FE space hierarch structure to create
  type(t_feHierarchy), intent(out) :: rfeHierarchyDest
!</output>
  
!</subroutine>
  
    integer :: i
    integer :: ntotalLevels,nmin,nmax
    
    nmin = 1
    if (present(nminLevel)) nmin = max(nmin,nminLevel)
    nmin = min(nmin,rfeHierarchySource%nlevels)
    
    nmax = rfeHierarchySource%nlevels
    if (present(nmaxLevel)) nmax = min(nmax,nmaxLevel)
    nmax = max(nmax,1)
    
    ! How much levels should we allocate?
    ntotalLevels = max(1,nmax-nmin+1)
    if (present(nmaxLevels)) ntotalLevels = max(nmaxLevels,ntotalLevels)

    ! Initialise the structure
    rfeHierarchyDest%cflags = 0
    rfeHierarchyDest%nmaxLevels = ntotalLevels
    rfeHierarchyDest%nlevels = nmax-nmin+1
    rfeHierarchyDest%p_rboundary => rfeHierarchySource%p_rboundary
    allocate(rfeHierarchyDest%p_rfeSpaces(rfeHierarchyDest%nmaxLevels))
    
    ! Copy data of the levels
    do i=1,rfeHierarchyDest%nlevels
    
      call fesph_createFEspaceRef (rfeHierarchyDest%p_rfeSpaces(i),i,&
          rfeHierarchySource%p_rfeSpaces(nmax-(rfeHierarchyDest%nlevels-i)),i,&
          fgetDiscr,rcollection)

    end do
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_createHierarchyFromFile (rfeHierarchy,nlevels,sTRIfile,ndim,&
      fgetDiscr,rcollection,npreref,nmaxLevels,rboundary)

!<description>
  ! Creates a FE space hierarchy based on a triangulation by refinement.
!</description>

!<input>
  ! Number of refinement levels.
  integer, intent(in) :: nlevels
  
  ! Name/path of a triangulation file that specifies the coarse mesh.
  character(len=*), intent(in) :: sTRIFile

  ! Dimension of the triangulation. NDIM2D for 2D, NDIM3D for 3D.
  integer, intent(in) :: ndim

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: Number of prerefinements to get the coarse mesh. Standard is =0.
  integer, intent(in), optional :: npreref

  ! OPTIONAL: Maximum number of levels, the structure should support.
  integer, intent(in), optional :: nmaxLevels

  ! OPTIONAL: A boundary object that defines the domain.
  type(t_boundary), intent(in), optional :: rboundary
!</input>

!<output>
  ! The FE space hierarch structure to create
  type(t_feHierarchy), intent(out) :: rfeHierarchy
!</output>
  
!</subroutine>
  
    integer :: i
    
    if (nlevels .le. 0) return
    
    ! Initialise the structure
    rfeHierarchy%nlevels = nlevels
    allocate(rfeHierarchy%p_rfeSpaces(nlevels))

    ! Initialise the mesh hierarchy.
    call mshh_initHierarchy(rfeHierarchy%rmeshHierarchy,nmaxLevels,sTRIfile,ndim,&
      npreref,rboundary)
        
    ! Refine up to level nlevels.
    call mshh_refineHierarchy2lv (rfeHierarchy%rmeshHierarchy,nlevels,&
        rboundary=rboundary)
    
    ! Initialise all the FE space levels.
    do i=1,nlevels
      call fesph_createFEspace (rfeHierarchy%p_rfeSpaces(i),i,&
          rfeHierarchy%rmeshHierarchy%p_Rtriangulations(i),i,&
          fgetDiscr,rcollection,rboundary)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_releaseHierarchy (rfeHierarchy)

!<description>
  ! Releases a FE hierarchy.
!</description>

!<inputoutput>
  ! The FE space hierarch structure to release.
  type(t_feHierarchy), intent(inout) :: rfeHierarchy
!</inputoutput>
  
!</subroutine>
  
    integer :: i
    
    ! Release all levels
    do i=rfeHierarchy%nlevels,1,-1
      call fesph_releaseFEspace (rfeHierarchy%p_rfeSpaces(i))
    end do
    
    ! Clean up the rest
    deallocate (rfeHierarchy%p_rfeSpaces)
    
    ! Release the mesh hierarchy if initialised.
    call mshh_releaseHierarchy(rfeHierarchy%rmeshHierarchy)
    
    rfeHierarchy%nlevels = 0
    rfeHierarchy%nmaxLevels = 0
    rfeHierarchy%cflags = 0
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_concatFeHierarchies (rfeHier1,rfeHier2,rfeHier)
  
!<description>
  ! Concatenates two FE hierarchies to a common FE hierarchy (cross product
  ! on every level).
!</description>

!<input>  
  ! First FE hierarchy
  type(t_feHierarchy), intent(in) :: rfeHier1

  ! Second FE hierarchy
  type(t_feHierarchy), intent(in) :: rfeHier2
!</input>

!<output>
  ! A new FE hierarchy structure.
  type(t_feHierarchy), intent(out) :: rfeHier
!</output>
  
!</subroutine>

    integer :: i

    ! As a basis, copy the first hierarchy
    rfeHier = rfeHier1
    
    ! The mesh hierarchy must be copied.
    call mshh_initHierarchy (rfeHier1%rmeshHierarchy,rfeHier%rmeshHierarchy)

    ! Allocate new memory for the levels and create them
    allocate(rfeHier%p_rfeSpaces(rfeHier%nlevels))
    
    do i=1,rfeHier%nlevels
      call fesph_concatFeSpaces (rfeHier1%p_rfeSpaces(i),rfeHier2%p_rfeSpaces(i),&
          rfeHier%p_rfeSpaces(i),rfeHier2%rmeshHierarchy%p_Rtriangulations(i))
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_deriveFeHierarchy (rfeHierSource,rfeHierDest,ifirstBlock,ilastBlock)
  
!<description>
  ! Creates an FE hierarchy as a subset of another FE hierarchy.
!</description>

!<input>  
  ! Source FE hierarchy
  type(t_feHierarchy), intent(in) :: rfeHierSource

  ! Number of the block on every level in rfeHierSource that should be
  ! used as first block on every level in rfeHierDest.
  integer, intent(in) :: ifirstBlock

  ! Number of the last block on every level in rfeHierSource that should be
  ! used as last block on every level in rfeHierDest.
  integer, intent(in) :: ilastBlock
!</input>

!<output>
  ! A new FE space structure.
  type(t_feHierarchy), intent(out) :: rfeHierDest
!</output>
  
!</subroutine>

    integer :: i

    ! As a basis, copy the source hierarchy
    rfeHierDest = rfeHierSource  

    ! The mesh hierarchy must be copied.
    call mshh_initHierarchy (rfeHierSource%rmeshHierarchy,rfeHierDest%rmeshHierarchy)

    ! Allocate new memory for the levels and create them
    allocate(rfeHierDest%p_rfeSpaces(rfeHierDest%nlevels))
    
    do i=1,rfeHierDest%nlevels
      call fesph_deriveFeSpaces (rfeHierSource%p_rfeSpaces(i),rfeHierDest%p_rfeSpaces(i),&
          ifirstBlock,ilastBlock,rfeHierDest%rmeshHierarchy%p_Rtriangulations(i))
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_infoStatistics (rfeSpaceLevel,bheadline,ilevel)
  
!<description>
  ! Prints out statistical information of the given FE space to the
  ! terminal. The output is formatted as a table with an optional headline
  ! and an optional level identifier.
!</description>

!<input>
  ! Triangulation structure.
  type(t_feSpaceLevel), intent(in) :: rfeSpaceLevel
  
  ! OPTIONAL: Print out a headline above the statistical data.
  ! =FALSE: do not print = standard.
  ! =TRUE: print a headline.
  logical, intent(in), optional :: bheadline
  
  ! OPTIONAL: Level identifier.
  ! If specified, an additional column 'Level' is added to the front of
  ! the statistics table. ilevel is printed to this column.
  integer, intent(in), optional :: ilevel
!</input>

!</subroutine>

    integer :: j

    if (bheadline) then
      ! Print a headline
      if (present(ilevel)) then
        call output_line("Lv.",bnoLineBreak=.true.)
      else
        call output_line ("",bnoLineBreak=.true.)
      end if
      do j=1,rfeSpaceLevel%p_rdiscretisation%ncomponents
        call output_line("    #dof(x"//trim(sys_si(j,1))//")",bnoLineBreak=.true.,&
            cdateTimeLogPolicy=OU_DTP_NONE)
      end do
      call output_line(" #dof(total)",cdateTimeLogPolicy=OU_DTP_NONE)
      call output_line("---",bnoLineBreak=.true.)
      do j=1,rfeSpaceLevel%p_rdiscretisation%ncomponents
        call output_line("------------",bnoLineBreak=.true.,cdateTimeLogPolicy=OU_DTP_NONE)
      end do
      call output_line("------------",cdateTimeLogPolicy=OU_DTP_NONE)
    end if
    
    ! Print statistics about that level
    if (present(ilevel)) then
      call output_line (trim(sys_si(ilevel,3)),bnoLineBreak=.true.)
    else
      call output_line ("",bnoLineBreak=.true.)
    end if
    do j=1,rfeSpaceLevel%p_rdiscretisation%ncomponents
      call output_line(&
          trim(sys_si(dof_igetNDofGlob(&
              rfeSpaceLevel%p_rdiscretisation%RspatialDiscr(j)),12)),&
              bnoLineBreak=.true.,cdateTimeLogPolicy=OU_DTP_NONE)
    end do

    call output_line(trim(sys_si(dof_igetNDofGlobBlock(&
        rfeSpaceLevel%p_rdiscretisation),12)),cdateTimeLogPolicy=OU_DTP_NONE)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fesph_printHierStatistics (rhierarchy)

!<description>
  ! Writes statistics about the mesh hierarchy to the terminal.
!</description>

!<inputoutput>
  ! FE hierarchy
  type(t_feHierarchy), intent(in) :: rhierarchy
!</inputoutput>
  
!</subroutine>
  
    integer :: i
  
    do i=1,rhierarchy%nlevels
      call fesph_infoStatistics (rhierarchy%p_rfeSpaces(i),(i .eq. 1),i)
    end do
    
  end subroutine

end module
