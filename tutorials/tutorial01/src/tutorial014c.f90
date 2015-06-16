!##############################################################################
!# Tutorial 014c: Analysis of a mixed 2d tri/quad mesh
!##############################################################################

module tutorial014c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  use boundary
  use triangulation
  use meshgeneration

  use ucd

  implicit none
  private

  public :: start_tutorial014c

contains

  ! ***************************************************************************

  subroutine start_tutorial014c

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    character(LEN=SYS_STRLEN) :: spredir,spostdir
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: i
    type(t_ucdExport) :: rexport

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 014c")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Read the underlying domain
    ! and the mesh
    ! =================================

    if (sys_getenv_string("PREDIR",spredir)) then
      call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
    else
      call boundary_read_prm(rboundary, "pre/QUAD.prm")
    end if

    ! Read a mixed tri/quad mesh that uses the above parametrisation
    if (sys_getenv_string("PREDIR",spredir)) then
      call tria_readTriFile2D (rtriangulation, trim(spredir)//"/QUAD_TRIA.tri", rboundary)
    else
      call tria_readTriFile2D (rtriangulation, "pre/QUAD_TRIA.tri", rboundary)
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

    i = rtriangulation%InelOfType(TRIA_NVETRI2D)
    call output_line ("Number of triangles      : " // trim(sys_siL( i, 10)) )

    i = rtriangulation%InelOfType(TRIA_NVEQUAD2D)
    call output_line ("Number of quads          : " // trim(sys_siL( i, 10)) )

    call output_lbrk()

    ! ---------------------------------
    ! For elements 1 and 5, print the number of vertices
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

    ! p_IverticesAtElement has always NNVE entries (NNVE=maximum number of
    ! vertices per element). For triangles in a quad mesh, there are 3 items
    ! used in p_IverticesAtElement, the array is filled with zero to NNVE items.
    ! Loop starting from NNVE and count down; reaching the first nonzero entry
    ! gives the number of vertices -- 3 for triangle, 4 for quad.

    ! Element 1 - a triangle.
    do i=rtriangulation%NNVE, 1, -1
      if (p_IverticesAtElement(i, 1) .ne. 0) exit
    end do

    call output_line ("Number of verts, elem. 1 : " // trim(sys_siL( i, 10)) )

    ! Element 5 - a quad
    do i=rtriangulation%NNVE, 1, -1
      if (p_IverticesAtElement(i, 5) .ne. 0) exit
    end do

    call output_line ("Number of verts, elem. 5 : " // trim(sys_siL( i, 10)) )

    ! =================================
    ! Create a VTK file for this mesh
    ! =================================

    call output_lbrk()
    call output_line ("Writing postprocessing files...")

    ! Open / write / close; write the solution to a VTK file.
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(spostdir)//"/tutorial014c.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial014c.vtk")
    end if
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================

    ! Release the triangulation and the boundary definition
    call tria_done (rtriangulation)

    call boundary_release(rboundary)

  end subroutine

end module
