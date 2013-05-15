!##############################################################################
!# Tutorial 014a: Analysis of a 2D brick mesh
!##############################################################################

module tutorial014a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  
  use triangulation
  use meshgeneration

  implicit none
  private
  
  public :: start_tutorial014a

contains

  ! ***************************************************************************

  subroutine start_tutorial014a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:,:), pointer :: h_IelementsAtEdge
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DelementVolume
    integer, dimension(:), pointer :: p_Iidx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    
    real(DP) :: dx,dy,dvol
    integer :: i, j, j1, j2, j3, j4

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 014a")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format to work with it.
    ! First create a 3x3-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 2, 2)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Print basic information about the mesh.
    ! =================================
    
    call output_line ("Mesh analysis: ")
    call output_line ("--------------")
    
    ! ---------------------------------
    ! Print basic information
    call output_line ("Mesh dimension           : " // trim(sys_siL( rtriangulation%ndim, 10)) )
    call output_lbrk()
    call output_line ("Number of elements       : " // trim(sys_siL( rtriangulation%NEL , 10)) )
    call output_line ("Number of vertices       : " // trim(sys_siL( rtriangulation%NVT , 10)) )
    call output_line ("Number of edges          : " // trim(sys_siL( rtriangulation%NMT , 10)) )
    call output_lbrk()
    call output_line ("Number of boundary comp. : " // trim(sys_siL( rtriangulation%NBCT, 10)) )
    call output_line ("Number of vert. on bd.   : " // trim(sys_siL( rtriangulation%NVBD, 10)) )
    call output_line ("Number of edges on bd.   : " // trim(sys_siL( rtriangulation%NMBD, 10)) )
    call output_lbrk()
    call output_line ("Bounding box             : " // &
        trim(sys_sdL( rtriangulation%DboundingBoxMin(1), 2)) // "/"   // &
        trim(sys_sdL( rtriangulation%DboundingBoxMin(2), 2)) // " - " // &
        trim(sys_sdL( rtriangulation%DboundingBoxMax(1), 2)) // "/"   // &
        trim(sys_sdL( rtriangulation%DboundingBoxMax(2), 2)) )
    call output_lbrk()
    
    ! ---------------------------------
    ! Access the coordinates of vertex 5.
    ! Get the pointer p_DvertexCoords associated to the handle h_DvertexCoords.
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    dx = p_DvertexCoords(1,5)
    dy = p_DvertexCoords(2,5)
    
    call output_line ("Coords of vertex 5       : " // trim(sys_sdL( dx,10 )) // " " &
                                                    // trim(sys_sdL( dy,10 )) )

    ! ---------------------------------
    ! Access boundary status of vertex 1, 2 and 5
    call storage_getbase_int (rtriangulation%h_InodalProperty,p_InodalProperty)
    j1 = p_InodalProperty(1)
    j2 = p_InodalProperty(2)
    j3 = p_InodalProperty(5)
    
    call output_line ("Boundary comp. vertex 1  : " // trim(sys_siL( j1,10 )) )
    call output_line ("Boundary comp. vertex 2  : " // trim(sys_siL( j2,10 )) )
    call output_line ("Boundary comp. vertex 5  : " // trim(sys_siL( j3,10 )) // " (0=inner node)")

    ! ---------------------------------            
    ! Access vertices on element 1
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    j1 = p_IverticesAtElement(1,1)
    j2 = p_IverticesAtElement(2,1)
    j3 = p_IverticesAtElement(3,1)
    j4 = p_IverticesAtElement(4,1)
    
    call output_line ("Vertices on element 1    : " // trim(sys_siL( j1,10 )) // " " &
                                                    // trim(sys_siL( j2,10 )) // " " &
                                                    // trim(sys_siL( j3,10 )) // " " &
                                                    // trim(sys_siL( j4,10 )) )
    
    ! ---------------------------------
    ! Access edges on element 1
    call storage_getbase_int2d (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    j1 = p_IedgesAtElement(1,1)
    j2 = p_IedgesAtElement(2,1)
    j3 = p_IedgesAtElement(3,1)
    j4 = p_IedgesAtElement(4,1)
    
    call output_line ("Edges on element 1       : " // trim(sys_siL( j1,10 )) // " " &
                                                    // trim(sys_siL( j2,10 )) // " " &
                                                    // trim(sys_siL( j3,10 )) // " " &
                                                    // trim(sys_siL( j4,10 )) )
    
    ! ---------------------------------
    ! Access neighbouring elements of element 1
    call storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
    j1 = p_IneighboursAtElement(1,1)
    j2 = p_IneighboursAtElement(2,1)
    j3 = p_IneighboursAtElement(3,1)
    j4 = p_IneighboursAtElement(4,1)
    
    call output_line ("Neighbours of element 1  : " // trim(sys_siL( j1,10 )) // " " &
                                                    // trim(sys_siL( j2,10 )) // " " &
                                                    // trim(sys_siL( j3,10 )) // " " &
                                                    // trim(sys_siL( j4,10 )) )
    
    ! ---------------------------------
    ! Access vertices on edge 3
    call storage_getbase_int2d (rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    j1 = p_IverticesAtEdge(1,1)
    j2 = p_IverticesAtEdge(2,1)
    
    call output_line ("Vertices on edge 3       : " // trim(sys_siL( j1,10 )) // " " &
                                                    // trim(sys_siL( j2,10 )) )

    ! ---------------------------------
    ! Access elements left and right of edge 3
    call storage_getbase_int2d (rtriangulation%h_IelementsAtEdge,h_IelementsAtEdge)
    j1 = h_IelementsAtEdge(1,1)
    j2 = h_IelementsAtEdge(2,1)
    
    call output_line ("Elements on edge 3       : " // trim(sys_siL( j1,10 )) // " " &
                                                    // trim(sys_siL( j2,10 )) )

    ! ---------------------------------
    ! Volume of element 1 and of the mesh
    call storage_getbase_double (rtriangulation%h_DelementVolume,p_DelementVolume)
    
    dvol = p_DelementVolume(1)
    call output_line ("Volume of element 1      : " // trim(sys_sdL( dvol,10 )) )

    dvol = p_DelementVolume(rtriangulation%NEL + 1)
    call output_line ("Volume of the mesh       : " // trim(sys_sdL( dvol,10 )) )

    ! ---------------------------------
    ! Elements adjacent to vertex 5 and 6
    call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,p_Iidx)
    call storage_getbase_int (rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)

    call output_line ("Elements on vertex 5     :" , bnolinebreak=.true.)
    
    do i = p_Iidx(5), p_Iidx(6)-1
      j = p_IelementsAtVertex(i)
      call output_line ( " " // trim(sys_siL( j, 10 )) , bnolinebreak = (i .ne. p_Iidx(6)-1) )
    end do

    call output_line ("Elements on vertex 6     :" , bnolinebreak=.true.)
    
    do i = p_Iidx(6), p_Iidx(7)-1
      j = p_IelementsAtVertex(i)
      call output_line ( " " // trim(sys_siL( j, 10 )) , bnolinebreak = (i .ne. p_Iidx(7)-1) )
    end do

    ! =================================
    ! Cleanup
    ! =================================
    
    call tria_done (rtriangulation)
    
  end subroutine

end module
