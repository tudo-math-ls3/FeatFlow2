!##############################################################################
!# Tutorial 012d: Create a mesh hierarchy
!##############################################################################

module tutorial012d

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary
  use triangulation
  use ucd
  
  use meshhierarchy

  implicit none
  private
  
  public :: start_tutorial012d

contains

  ! ***************************************************************************

  subroutine start_tutorial012d

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    type(t_ucdExport) :: rexport
    integer :: i
    
    type(t_meshHierarchy) :: rmeshHierarchy

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 012d")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Read the underlying domain
    ! and the mesh
    ! =================================

    call boundary_read_prm(rboundary, "pre/QUAD.prm")
    
    ! The mesh must always be in "standard" format to work with it.
    ! First read, then convert to standard.
    call tria_readTriFile2D (rtriangulation, "pre/QUAD.tri", rboundary)
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! =================================
    ! Create a hierarchy of 5 meshes with
    ! rtriangulation as coarse mesh.
    ! =================================
    
    ! Create a hierarchy with 5 levels,
    ! initialise the coarse mesh with rtriangulation
    call mshh_initHierarchy (rmeshHierarchy,rtriangulation,0,5,rboundary)
    
    ! Refine the coarse mesh 5x to get the hierarchy.
    call mshh_refineHierarchy2lv (rmeshHierarchy,5,rboundary=rboundary)
    
    ! =================================
    ! Write the meshes
    ! =================================
    
    do i=1,rmeshHierarchy%nlevels
      call output_line ("Writing file 'post/tutorial012d_level"//trim(sys_siL(i,10))//".vtk'.")

      ! Open / write / close
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
          "post/tutorial012d_level"//trim(sys_siL(i,10))//".vtk")
      call ucd_write (rexport)
      call ucd_release (rexport)
    end do
    
    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the mesh hierarchy
    call mshh_releaseHierarchy (rmeshHierarchy)

    ! Release the triangulation    
    call tria_done (rtriangulation)
    
    ! Release the boundary definition
    call boundary_release(rboundary)
  end subroutine

end module
