!##############################################################################
!# ****************************************************************************
!# <name> meshmodification </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a collection of different (more or less) simple mesh
!# modification routines.
!# The following routines can be found here:
!#
!# 1.) meshmod_disturbMesh
!#     -> Applies stochastical grid disturbance
!#
!# </purpose>
!##############################################################################

MODULE meshmodification

  USE triangulation

  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE meshmod_disturbMesh (rtriangulation,damount,InodalProperty)
  
!<description>
  ! Applies stochastical grid disturbance to a given mesh. damount is a
  ! value in the range 0..1 and specifies the amount of disturbance
  ! that should be applied to the mesh rtriangulation.
  !
  ! Boundary vertices are not disturbed.
!</description>

!<input>
  ! Amount of stochastical grid disturbance to be applied to rtriangulation.
  ! Range 0..1; e.g. 0.2 = 20%.
  REAL(DP), INTENT(IN) :: damount
  
  ! OPTIONAL: Nodal property of the points, specifying which points should
  ! be disturbed and which not.
  ! A point i with InodalProperty(i)=0 is disturbed. All points with
  ! InodalProperty<>0 are not disturbed.
  ! If not specified, the nodal property array from rtriangulation is used,
  ! thus all points in the inner are disturbed while all points on the
  ! boundary are kept as they are.
  INTEGER(I32), DIMENSION(:), INTENT(IN), OPTIONAL, TARGET :: InodalProperty
!</input>

!<inputoutput>
  ! Mesh whose grid points should be disturbed.
  ! THe grid point coordinates are modified by this routine.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(I32), DIMENSION(:), POINTER :: p_InodalProperty
    REAL(DP) :: dh,dhdist
    INTEGER(PREC_VERTEXIDX) :: ivt

    ! Mesh width
    dh = 1.0_DP/(SQRT(REAL(rtriangulation%NVT,DP))-1.0_DP)
    
    ! Amount of distortion
    dhdist = damount * dh
    
    ! Get arrays for coordinates / boundary definition
    CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    IF (PRESENT(InodalProperty)) THEN
      p_InodalProperty => InodalProperty
    ELSE
      CALL storage_getbase_int (rtriangulation%h_InodalProperty,&
          p_InodalProperty)
    END IF
    
    ! Loop over the vertices and disturb them
    DO ivt = 1,rtriangulation%NVT
    
      ! Only modify inner points.
      IF (p_InodalProperty(ivt) .EQ. 0) THEN
        p_DvertexCoords(:,ivt) = &
          p_DvertexCoords(:,ivt) + REAL((-1)**MOD(ivt,17),DP)*dhdist
      END IF
    
    END DO

  END SUBROUTINE

END MODULE
