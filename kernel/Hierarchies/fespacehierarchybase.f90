!##############################################################################
!# ****************************************************************************
!# <name> fespacehierarchybase </name>
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
!# </purpose>
!##############################################################################

module fespacehierarchybase

!$use omp_lib
  use fsystem
  use boundary
  use basicgeometry
  use triangulation
  use meshhierarchy
  use spatialdiscretisation
  
  implicit none
  
  private
  
  public :: t_feSpaceLevel
  public :: t_feHierarchy
  
!<constants>

!<constantblock description = "Constants defining shared information">

  ! Triangulation is shared.
  integer(I32), parameter, public :: FESPH_SHAREDTRIA = 2_I32**0

  ! Discretisation is shared.
  integer(I32), parameter, public :: FESPH_SHAREDDISCR = 2_I32**1

!</constantblock>

!</constants>


!<types>

!<typeblock>

  ! A structure that specifies a FE space which can be used to create
  ! FE vectors/matrices.
  type t_feSpaceLevel
  
    ! A shared flag that specifies which information in this structure
    ! is shared with information from outside and which information
    ! belongs to this structure.
    integer(I32) :: cflags = 0
    
    ! Reference to the underlying domain or NULL() if no domain is attached.
    type(t_boundary), pointer :: p_rboundary => null()
    
    ! Reference to the underlying triangulation.
    ! This is either a direct reference or refers to allocated memory
    ! on the heap if the triangulation had been created by refinement
    ! in the create routines.
    type(t_triangulation), pointer :: p_rtriangulation => null()
    
    ! Reference to the underlying discretisation.
    ! This is either a direct reference or refers to allocated memory
    ! on the heap if the discretisation had been created by
    ! the create routines.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()
    
  end type

!</typeblock>

!<typeblock>

  ! A FE hierarchy that describes a hierarchy of FE spaces.
  type t_feHierarchy
  
    ! A shared flag that specifies which information in this structure
    ! is shared with information from outside and which information
    ! belongs to this structure.
    integer(I32) :: cflags = 0

    ! Reference to the underlying domain or NULL() if no domain is attached.
    type(t_boundary), pointer :: p_rboundary => null()

    ! An underlying mesh hierarchy.
    type(t_meshHierarchy) :: rmeshHierarchy
    
    ! Number of levels available in this structure.
    integer :: nlevels = 0
    
    ! Maximum number of available levels available in this structure.
    integer :: nmaxLevels = 0
    
    ! Level information.
    type(t_feSpaceLevel), dimension(:), pointer :: p_rfeSpaces => null()
    
  end type

!</typeblock>

!</types>

end module
