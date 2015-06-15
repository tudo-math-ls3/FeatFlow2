!##############################################################################
!# Tutorial 014b: Analysis of a 2D mesh with boundary
!##############################################################################

module tutorial014b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  use boundary
  use triangulation
  use meshgeneration

  implicit none
  private

  public :: start_tutorial014b

contains

  ! ***************************************************************************

  subroutine start_tutorial014b

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    character(LEN=SYS_STRLEN) :: spredir

    integer, dimension(:), pointer :: p_InodalProperty
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:), pointer :: p_Iidx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:), pointer :: p_DedgeParameterValue

    real(DP) :: dx,dy,dlambda
    integer :: i,j,j1,j2,j3

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 014b")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Read the underlying domain
    ! and the mesh
    ! =================================

    if (sys_getenv_string("PREDIR",spredir)) then
      call boundary_read_prm(rboundary, trim(spredir)//"/bench1.prm")
    else
      call boundary_read_prm(rboundary, "pre/bench1.prm")
    end if

    ! The mesh must always be in "standard" format to work with it.
    ! First read, then convert to standard, based on rboundary.
    if (sys_getenv_string("PREDIR",spredir)) then
      call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
      call tria_readTriFile2D (rtriangulation, trim(spredir)//"pre/bench1.tri", rboundary)
    else
      call tria_readTriFile2D (rtriangulation, "pre/bench1.tri", rboundary)
    end if
    call tria_initStandardMeshFromRaw (rtriangulation, rboundary)

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
    ! Access the coordinates of vertices 42, 1, 13
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)

    dx = p_DvertexCoords(1,42)
    dy = p_DvertexCoords(2,42)
    call output_line ("Coords of vertex 42      : " // trim(sys_sdL( dx,10 )) // " " &
                                                    // trim(sys_sdL( dy,10 )) )

    dx = p_DvertexCoords(1,1)
    dy = p_DvertexCoords(2,1)
    call output_line ("Coords of vertex 1       : " // trim(sys_sdL( dx,10 )) // " " &
                                                    // trim(sys_sdL( dy,10 )) )

    dx = p_DvertexCoords(1,13)
    dy = p_DvertexCoords(2,13)
    call output_line ("Coords of vertex 13      : " // trim(sys_sdL( dx,10 )) // " " &
                                                    // trim(sys_sdL( dy,10 )) )

    ! ---------------------------------
    ! Access boundary status of the vertices 42, 1, 13
    call storage_getbase_int (rtriangulation%h_InodalProperty,p_InodalProperty)
    j1 = p_InodalProperty(42)
    j2 = p_InodalProperty(1)
    j3 = p_InodalProperty(13)

    call output_line ("Boundary comp. vertex 1  : " // trim(sys_siL( j1,10 )) // " (1=outer box)")
    call output_line ("Boundary comp. vertex 2  : " // trim(sys_siL( j2,10 )) // " (2=circle)")
    call output_line ("Boundary comp. vertex 5  : " // trim(sys_siL( j3,10 )) // " (0=inner node)")

    ! ---------------------------------
    ! Number of vertices on the boundary components
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_Iidx)

    j = p_Iidx(2)-p_Iidx(1)
    call output_line ("#Vertices, bd. comp. 1   : " // trim(sys_siL( j,10 )) )

    j = p_Iidx(3)-p_Iidx(2)
    call output_line ("#Vertices, bd. comp. 2   : " // trim(sys_siL( j,10 )) )

    ! ---------------------------------
    ! Vertices on the circle
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)

    call output_line ("Vertices on the circle   :" , bnolinebreak=.true.)

    do i = p_Iidx(2), p_Iidx(3)-1
      j = p_IverticesAtBoundary(i)
      call output_line ( " " // trim(sys_siL( j, 10 )) , bnolinebreak = (i .ne. p_Iidx(3)-1) )
    end do

    ! ---------------------------------
    ! Edges on the circle
    call storage_getbase_int (rtriangulation%h_IelementsAtVertex,p_IedgesAtBoundary)

    call output_line ("Edges on the circle      :" , bnolinebreak=.true.)

    do i = p_Iidx(2), p_Iidx(3)-1
      j = p_IedgesAtBoundary(i)
      call output_line ( " " // trim(sys_siL( j, 10 )) , bnolinebreak = (i .ne. p_Iidx(3)-1) )
    end do

    ! ---------------------------------
    ! Elements on the circle
    call storage_getbase_int (rtriangulation%h_IelementsAtBoundary,p_IelementsAtBoundary)

    call output_line ("Elements on the circle   :" , bnolinebreak=.true.)

    do i = p_Iidx(2), p_Iidx(3)-1
      j = p_IelementsAtBoundary(i)
      call output_line ( " " // trim(sys_siL( j, 10 )) , bnolinebreak = (i .ne. p_Iidx(3)-1) )
    end do

    ! ---------------------------------
    ! Parameter values of the vertices on the circle
    call storage_getbase_double (rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)

    call output_line ("Par. val. of circle verts:" , bnolinebreak=.true.)

    do i = p_Iidx(2), p_Iidx(3)-1
      dlambda = p_DvertexParameterValue(i)
      call output_line ( " " // trim(sys_sdL( dlambda, 10 )) )
    end do

    ! ---------------------------------
    ! Parameter values of the edges (=edge midpoints) on the circle
    call storage_getbase_double (rtriangulation%h_DedgeParameterValue,p_DedgeParameterValue)

    call output_line ("Par. val. of circle edges:" , bnolinebreak=.true.)

    do i = p_Iidx(2), p_Iidx(3)-1
      dlambda = p_DedgeParameterValue(i)
      call output_line ( " " // trim(sys_sdL( dlambda, 10 )) )
    end do

    ! ---------------------------------
    ! Find vertices 11 and 13 on the boundary.
    ! Vertex 11 is on the circle, vertex 13 is an inner vertex.
    call tria_searchBoundaryVertex(13, rtriangulation, i)
    call output_line ("Boundary-idx of vert. 13 : " // trim(sys_siL( i, 10 )) // " (0=inner vertex)")

    call tria_searchBoundaryVertex(11, rtriangulation, i)
    call output_line ("Boundary-idx of vert. 11 : " // trim(sys_siL( i, 10 )) // " (0=inner vertex)")

    ! With the index i, we can print the parameter value of vertex 11
    dlambda = p_DvertexParameterValue(i)
    call output_line ("Corresp. parameter value : " // trim(sys_sdL( dlambda, 10 )) )

    ! ---------------------------------
    ! Find edges 1 and 3 on the boundary.
    call tria_searchBoundaryEdge(1, rtriangulation, i)
    call output_line ("Boundary-idx of edge 1   : " // trim(sys_siL( i, 10 )) // " (0=inner edge)")

    dlambda = p_DedgeParameterValue(i)
    call output_line ("Corresp. parameter value : " // trim(sys_sdL( dlambda, 10 )) )

    call tria_searchBoundaryEdge(3, rtriangulation, i)
    call output_line ("Boundary-idx of edge 3   : " // trim(sys_siL( i, 10 )) // " (0=inner edge)")

    ! =================================
    ! Cleanup
    ! =================================

    ! Release the triangulation and the boundary definition
    call tria_done (rtriangulation)

    call boundary_release(rboundary)

  end subroutine

end module
