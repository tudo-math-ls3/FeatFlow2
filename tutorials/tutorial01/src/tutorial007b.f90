!##############################################################################
!# Tutorial 007b: Identify points in the mesh using an expression
!##############################################################################

module tutorial007b

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  use meshregion
  use fparser

  implicit none
  private
  
  public :: start_tutorial007b

contains

  ! ***************************************************************************

  subroutine start_tutorial007b

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_meshRegion) :: rmeshRegion
    integer :: i
    integer, dimension(:), pointer :: p_Idata
    character(len=SYS_STRLEN) :: scondition
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 007b")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Initialise the parser.
    ! =================================

    call fparser_init()

    ! =================================
    ! Identify all points on the boundary.
    ! =================================
    
    ! The condition fetches all points on the right/top.
    scondition = "(X=1) | (Y=1)"
    
    call mshreg_createFromExpression(rmeshRegion, rtriangulation, &
        MSHREG_IDX_ALL, .true., scondition)

    ! =================================
    ! Output
    ! =================================

    ! ---------------------------------
    call output_line ("Points on the right/top boundary:")
    
    ! Get the list of points.
    call storage_getbase_int (rmeshRegion%h_IvertexIdx,p_Idata)

    ! Print the points in the mesh region
    do i=1,rmeshRegion%NVT
      call output_line (trim(sys_siL(p_Idata(i),10)))
    end do

    ! ---------------------------------
    call output_lbrk()
    call output_line ("Edges on the right/top boundary:")
    
    ! Get the list of points.
    call storage_getbase_int (rmeshRegion%h_IedgeIdx,p_Idata)

    ! Print the edges in the mesh region
    do i=1,rmeshRegion%NMT
      call output_line (trim(sys_siL(p_Idata(i),10)))
    end do

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the mesh region
    call mshreg_done(rmeshRegion)

    ! Clean up the parser.
    call fparser_done()
    
    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
