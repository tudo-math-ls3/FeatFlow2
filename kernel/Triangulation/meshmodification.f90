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

module meshmodification

  use triangulation

  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine meshmod_disturbMesh (rtriangulation,damount,InodalProperty)
  
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
  real(DP), intent(IN) :: damount
  
  ! OPTIONAL: Nodal property of the points, specifying which points should
  ! be disturbed and which not.
  ! A point i with InodalProperty(i)=0 is disturbed. All points with
  ! InodalProperty<>0 are not disturbed.
  ! If not specified, the nodal property array from rtriangulation is used,
  ! thus all points in the inner are disturbed while all points on the
  ! boundary are kept as they are.
  integer(I32), dimension(:), intent(IN), optional, target :: InodalProperty
!</input>

!<inputoutput>
  ! Mesh whose grid points should be disturbed.
  ! THe grid point coordinates are modified by this routine.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer(I32), dimension(:), pointer :: p_InodalProperty
    real(DP) :: dh,dhdist
    integer(PREC_VERTEXIDX) :: ivt

    ! Mesh width
    dh = 1.0_DP/(sqrt(real(rtriangulation%NVT,DP))-1.0_DP)
    
    ! Amount of distortion
    dhdist = damount * dh
    
    ! Get arrays for coordinates / boundary definition
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    if (present(InodalProperty)) then
      p_InodalProperty => InodalProperty
    else
      call storage_getbase_int (rtriangulation%h_InodalProperty,&
          p_InodalProperty)
    end if
    
    ! Loop over the vertices and disturb them
    do ivt = 1,rtriangulation%NVT
    
      ! Only modify inner points.
      if (p_InodalProperty(ivt) .eq. 0) then
        p_DvertexCoords(:,ivt) = &
          p_DvertexCoords(:,ivt) + real((-1)**mod(ivt,17),DP)*dhdist
      end if
    
    end do

  end subroutine

end module
