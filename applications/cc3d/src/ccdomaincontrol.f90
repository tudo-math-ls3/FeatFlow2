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

MODULE ccdomaincontrol

  USE fsystem
  USE paramlist
  USE collection
  USE triangulation
  USE meshregion
    
  USE ccbasic
  
  ! Domain-Dependend includes
  USE dom3d_cube
  USE dom3d_c3d0
  USE dom3d_c3d1
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Domain Identifiers">
  
  ! Domain Identifier for CUBE-domain
  INTEGER, PARAMETER :: CCDC_ID_DOM3D_CUBE = 1
  
  ! Domain Identifier for C3D0-domain
  INTEGER, PARAMETER :: CCDC_ID_DOM3D_C3D0 = 2
  
  ! Domain Identifier for C3D1-domain
  INTEGER, PARAMETER :: CCDC_ID_DOM3D_C3D1 = 3

!</constantblock>

!</constants>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ccdc_init (rproblem)
  
!<description>
  ! This routine initialises the domain control.
  ! If the selected domain needs some additional initialisation as parsing
  ! something from the parameter list or adding information into the
  ! collection, this can be performed in this routine.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT)                 :: rproblem
!</inputoutput>

!</subroutine>

  INTEGER :: idomain
  
    ! Read in the domain identifier from the parameter list
    CALL parlst_getvalue_int (rproblem%rparamList,'DOMAIN','iDomain',&
                              idomain,-1)

    ! Store the domain id to the problem structure
    rproblem%idomain = idomain
    
    ! And store it to the collection structure
    CALL collct_setvalue_int (rproblem%rcollection,'IDOMAIN',idomain,.TRUE.)
    
    ! Place domain-dependent initialisation here...

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ccdc_done (rproblem)
  
!<description>
  ! This routine releases the domain control.
  ! If the selected domain has allocated any addition information on the heap,
  ! in the collection or elsewhere, it has to release the memory in this
  ! routine.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT)                 :: rproblem
!</inputoutput>

!</subroutine>

    ! Place domain-dependent destruction here...
    ! INTEGER :: idomain
    ! CALL parlst_getvalue_int (rproblem%rparamList,'DOMAIN','iDomain',&
    !                           idomain,-1)

    ! Now we can remove the domain identifier from the collection
    CALL collct_deletevalue(rproblem%rcollection, 'IDOMAIN')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE ccdc_getNumRegions (rproblem, inumRegions)
  
!<description>
  ! This routine returns the number of boundary regions of the domain that
  ! was selected in the domain.dat file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT)                 :: rproblem
!</inputoutput>

!<output>
  ! The number of boundary regions of the domain
  INTEGER, INTENT(OUT)                           :: inumRegions
!</output>

!</subroutine>

    ! What domain do we have here?
    SELECT CASE(rproblem%idomain)
    CASE (CCDC_ID_DOM3D_CUBE)
      inumRegions = 6
      
    CASE (CCDC_ID_DOM3D_C3D0)
      inumRegions = 7
    
    CASE (CCDC_ID_DOM3D_C3D1)
      inumRegions = 7
    
    CASE DEFAULT
      PRINT *, 'ERROR: ccdc_getNumRegions: Unknown Domain identifier'
      CALL sys_halt()
      
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ccdc_calcBoundaryMeshRegion (rproblem, rmeshRegion, &
                                          rtriangulation, Iregions)
  
!<description>
  ! This routine calculates a mesh region containing all vertices, edges and
  ! faces that belong to a set of specified domain boundary regions.
!</description>

!<input>
  ! The triangulation on which the mesh region is to be defined.
  TYPE(t_triangulation), TARGET, INTENT(IN)      :: rtriangulation
  
  ! An array holding the indices of all domain regions which should be added
  ! into the mesh region.
  INTEGER, DIMENSION(:), INTENT(IN)              :: Iregions
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT)                 :: rproblem
!</inputoutput>

!<output>
  ! The mesh region
  TYPE(t_meshRegion), INTENT(OUT)                :: rmeshRegion
!</output>

!</subroutine>

    ! Call the corresponding routine
    SELECT CASE(rproblem%idomain)
    CASE (CCDC_ID_DOM3D_CUBE)
      CALL dom3d_cube_calcMeshRegion(rmeshRegion,rtriangulation,Iregions)
    
    CASE (CCDC_ID_DOM3D_C3D0)
      CALL dom3d_c3d0_calcMeshRegion(rmeshRegion,rtriangulation,Iregions)
    
    CASE (CCDC_ID_DOM3D_C3D1)
      CALL dom3d_c3d1_calcMeshRegion(rmeshRegion,rtriangulation,Iregions)

    CASE DEFAULT
      PRINT *, 'ERROR: ccdc_calcBoundaryMeshRegion: Unknown Domain identifier'
      CALL sys_halt()
      
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE ccdc_correctMesh (rproblem, rtriangulation)
  
!<description>
  ! This routine corrects a mesh after it has been read in from a .TRI file
  ! or after it has been refined.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT)                 :: rproblem

  ! The triangulation which is to be corrected.
  TYPE(t_triangulation), INTENT(INOUT)           :: rtriangulation
!</inputoutput>

!</subroutine>

    ! Call the corresponding routine
    SELECT CASE(rproblem%iDomain)
    CASE (CCDC_ID_DOM3D_C3D0)
      CALL dom3d_c3d0_correctMesh(rtriangulation)
    
    ! We will not abort if the domain does not support mesh correction
    ! as most of the simple domains do not require any mesh corrections.
    
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE ccdc_calcParProfile (dvalue, iregion, dx, dy, dz, rcollection)
  
!<description>
  ! This routine calculates the parabolic profile for a specified domain
  ! boundary region.
!</description>

!<input>
  ! Specifies the index of the boundary region of which the parabolic profile
  ! is to be calculated.
  INTEGER, INTENT(IN)                            :: iregion
  
  ! The coordinates of the point in which the parabolic profile is to be
  ! evaluated.
  REAL(DP), INTENT(IN)                           :: dx, dy, dz
!</input>

!<inputoutput>
  ! A collection structure holding the domain identifier.
  TYPE(t_collection), INTENT(INOUT)              :: rcollection
!</inputoutput>

!<output>
  ! The result of the evaluation of the parabolic profile.
  REAL(DP), INTENT(OUT)                          :: dvalue
!</output>

!</subroutine>

  ! The domain id of the currently selection domain.
  INTEGER(I32) :: idomain
    
    ! Get the index of the domain from the collection
    idomain = collct_getvalue_int(rcollection, 'IDOMAIN')

    ! And call the corresponding routine
    SELECT CASE(idomain)
    CASE (CCDC_ID_DOM3D_CUBE)
      CALL dom3d_cube_calcParProfile(dvalue, iregion, dx, dy, dz)
    
    CASE (CCDC_ID_DOM3D_C3D0)
      CALL dom3d_c3d0_calcParProfile(dvalue, iregion, dx, dy, dz)
    
    CASE (CCDC_ID_DOM3D_C3D1)
      CALL dom3d_c3d1_calcParProfile(dvalue, iregion, dx, dy, dz)
    
    CASE DEFAULT
      PRINT *, 'ERROR: ccdc_calcParProfile: Unknown Domain identifier'
      CALL sys_halt()
    
    END SELECT
    
  END SUBROUTINE

END MODULE