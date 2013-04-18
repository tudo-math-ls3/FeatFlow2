!##############################################################################
!# Tutorial 007a: Identify points in the mesh
!##############################################################################

module tutorial007a

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  use meshregion

  implicit none
  private
  
  public :: start_tutorial007a

contains

  ! ***************************************************************************

  subroutine start_tutorial007a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_meshRegion) :: rmeshRegion
    integer :: i
    integer, dimension(:), pointer :: p_Idata
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 007a")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Identify all points on the boundary.
    ! =================================
    
    call mshreg_createFromNodalProp(rmeshRegion, rtriangulation,MSHREG_IDX_ALL)

    ! =================================
    ! Output
    ! =================================

    ! ---------------------------------
    call output_line ("Points on the boundary:")
    
    ! Get the list of points.
    call storage_getbase_int (rmeshRegion%h_IvertexIdx,p_Idata)

    ! Print the points in the mesh region
    do i=1,rmeshRegion%NVT
      call output_line (trim(sys_siL(p_Idata(i),10)))
    end do

    ! ---------------------------------
    call output_lbrk()
    call output_line ("Edges on the boundary:")
    
    ! Get the list of points.
    call storage_getbase_int (rmeshRegion%h_IedgeIdx,p_Idata)

    ! Print the points in the mesh region
    do i=1,rmeshRegion%NMT
      call output_line (trim(sys_siL(p_Idata(i),10)))
    end do

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
