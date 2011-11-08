!##############################################################################
!# ****************************************************************************
!# <name> ccdomaincontrol </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all domain-dependent routines, such as reading
!# the .TRI filename of the coarse mesh, performing mesh correction after a
!# mesh has been refined or calulculating mesh regions describing the domain's
!# boundary region.
!#
!# If a new domain is to be added, the routines in this module have to be
!# modified to work with the new domain.
!#
!# 1.) ccdc_init
!#     -> Initialises the domain control.
!#
!# 2.) ccdc_done
!#     -> Releases the domain control.
!#
!# 3.) ccdc_getNumRegions
!#     -> Returns the number of regions of the selected domain.
!#
!# 4.) ccdc_calcBoundaryMeshRegion
!#     -> Calculates a mesh region containing all desired boundary vertices,
!#        edges and faces which belong to specific domain boundary regions.
!#
!# 5.) ccdc_correctMesh
!#     -> Corrects a mesh that has been refined, e.g. projects vertices
!#        which belong to the boundary onto the domain's boundary.
!#
!# 6.) ccdc_calcParProfile
!#     -> Calculates the parabolic profile for a given domain region.
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-
!# = How to add a new domain =
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-
!# - To Do -
!#
!# </purpose>
!##############################################################################

module ccdomaincontrol

  use fsystem
  use paramlist
  use collection
  use triangulation
  use meshregion
    
  use ccbasic
  
  ! Domain-Dependend includes
  use dom3d_cube
  use dom3d_c3d0
  use dom3d_c3d1
  
  implicit none

!<constants>

!<constantblock description="Domain Identifiers">
  
  ! Domain Identifier for CUBE-domain
  integer, parameter :: CCDC_ID_DOM3D_CUBE = 1
  
  ! Domain Identifier for C3D0-domain
  integer, parameter :: CCDC_ID_DOM3D_C3D0 = 2
  
  ! Domain Identifier for C3D1-domain
  integer, parameter :: CCDC_ID_DOM3D_C3D1 = 3

!</constantblock>

!</constants>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine ccdc_init (rproblem)
  
!<description>
  ! This routine initialises the domain control.
  ! If the selected domain needs some additional initialisation as parsing
  ! something from the parameter list or adding information into the
  ! collection, this can be performed in this routine.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout)                 :: rproblem
!</inputoutput>

!</subroutine>

  integer :: idomain
  
    ! Read in the domain identifier from the parameter list
    call parlst_getvalue_int (rproblem%rparamList,'DOMAIN','iDomain',&
                              idomain,-1)

    ! Store the domain id to the problem structure
    rproblem%idomain = idomain
    
    ! And store it to the collection structure
    call collct_setvalue_int (rproblem%rcollection,'IDOMAIN',idomain,.true.)
    
    ! Place domain-dependent initialisation here...

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine ccdc_done (rproblem)
  
!<description>
  ! This routine releases the domain control.
  ! If the selected domain has allocated any addition information on the heap,
  ! in the collection or elsewhere, it has to release the memory in this
  ! routine.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout)                 :: rproblem
!</inputoutput>

!</subroutine>

    ! Place domain-dependent destruction here...
    ! INTEGER :: idomain
    ! CALL parlst_getvalue_int (rproblem%rparamList,'DOMAIN','iDomain',&
    !                           idomain,-1)

    ! Now we can remove the domain identifier from the collection
    call collct_deletevalue(rproblem%rcollection, 'IDOMAIN')

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine ccdc_getNumRegions (rproblem, inumRegions)
  
!<description>
  ! This routine returns the number of boundary regions of the domain that
  ! was selected in the domain.dat file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout)                 :: rproblem
!</inputoutput>

!<output>
  ! The number of boundary regions of the domain
  integer, intent(out)                           :: inumRegions
!</output>

!</subroutine>

    ! What domain do we have here?
    select case(rproblem%idomain)
    case (CCDC_ID_DOM3D_CUBE)
      inumRegions = 6
      
    case (CCDC_ID_DOM3D_C3D0)
      inumRegions = 7
    
    case (CCDC_ID_DOM3D_C3D1)
      inumRegions = 7
    
    case DEFAULT
      print *, 'ERROR: ccdc_getNumRegions: Unknown Domain identifier'
      call sys_halt()
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine ccdc_calcBoundaryMeshRegion (rproblem, rmeshRegion, &
                                          rtriangulation, Iregions)
  
!<description>
  ! This routine calculates a mesh region containing all vertices, edges and
  ! faces that belong to a set of specified domain boundary regions.
!</description>

!<input>
  ! The triangulation on which the mesh region is to be defined.
  type(t_triangulation), target, intent(in)      :: rtriangulation
  
  ! An array holding the indices of all domain regions which should be added
  ! into the mesh region.
  integer, dimension(:), intent(in)              :: Iregions
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout)                 :: rproblem
!</inputoutput>

!<output>
  ! The mesh region
  type(t_meshRegion), intent(out)                :: rmeshRegion
!</output>

!</subroutine>

    ! Call the corresponding routine
    select case(rproblem%idomain)
    case (CCDC_ID_DOM3D_CUBE)
      call dom3d_cube_calcMeshRegion(rmeshRegion,rtriangulation,Iregions)
    
    case (CCDC_ID_DOM3D_C3D0)
      call dom3d_c3d0_calcMeshRegion(rmeshRegion,rtriangulation,Iregions)
    
    case (CCDC_ID_DOM3D_C3D1)
      call dom3d_c3d1_calcMeshRegion(rmeshRegion,rtriangulation,Iregions)

    case DEFAULT
      print *, 'ERROR: ccdc_calcBoundaryMeshRegion: Unknown Domain identifier'
      call sys_halt()
      
    end select

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine ccdc_correctMesh (rproblem, rtriangulation)
  
!<description>
  ! This routine corrects a mesh after it has been read in from a .TRI file
  ! or after it has been refined.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout)                 :: rproblem

  ! The triangulation which is to be corrected.
  type(t_triangulation), intent(inout)           :: rtriangulation
!</inputoutput>

!</subroutine>

    ! Call the corresponding routine
    select case(rproblem%iDomain)
    case (CCDC_ID_DOM3D_C3D0)
      call dom3d_c3d0_correctMesh(rtriangulation)
    
    ! We will not abort if the domain does not support mesh correction
    ! as most of the simple domains do not require any mesh corrections.
    
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine ccdc_calcParProfile (dvalue, iregion, dx, dy, dz, rcollection)
  
!<description>
  ! This routine calculates the parabolic profile for a specified domain
  ! boundary region.
!</description>

!<input>
  ! Specifies the index of the boundary region of which the parabolic profile
  ! is to be calculated.
  integer, intent(in)                            :: iregion
  
  ! The coordinates of the point in which the parabolic profile is to be
  ! evaluated.
  real(DP), intent(in)                           :: dx, dy, dz
!</input>

!<inputoutput>
  ! A collection structure holding the domain identifier.
  type(t_collection), intent(inout)              :: rcollection
!</inputoutput>

!<output>
  ! The result of the evaluation of the parabolic profile.
  real(DP), intent(out)                          :: dvalue
!</output>

!</subroutine>

  ! The domain id of the currently selection domain.
  integer(I32) :: idomain
    
    ! Get the index of the domain from the collection
    idomain = collct_getvalue_int(rcollection, 'IDOMAIN')

    ! And call the corresponding routine
    select case(idomain)
    case (CCDC_ID_DOM3D_CUBE)
      call dom3d_cube_calcParProfile(dvalue, iregion, dx, dy, dz)
    
    case (CCDC_ID_DOM3D_C3D0)
      call dom3d_c3d0_calcParProfile(dvalue, iregion, dx, dy, dz)
    
    case (CCDC_ID_DOM3D_C3D1)
      call dom3d_c3d1_calcParProfile(dvalue, iregion, dx, dy, dz)
    
    case DEFAULT
      print *, 'ERROR: ccdc_calcParProfile: Unknown Domain identifier'
      call sys_halt()
    
    end select
    
  end subroutine

end module
