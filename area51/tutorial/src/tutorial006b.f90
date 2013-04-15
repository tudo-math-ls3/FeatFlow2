!##############################################################################
!# Tutorial 006b: Read a basic 2D mesh and write a VTK file.
!##############################################################################

module tutorial006b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary
  use triangulation
  use ucd

  implicit none
  private
  
  public :: start_tutorial006b

contains

  ! ***************************************************************************

  subroutine start_tutorial006b

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    type(t_ucdExport) :: rexport

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006b")
    call output_separator (OU_SEP_MINUS)
    
    call output_line ("Writing file 'gmv/tutorial006b.vtk'.")

    ! =================================
    ! Read the underlying domain
    ! and the mesh
    ! =================================

    call boundary_read_prm(rboundary, "pre/bench1.prm")
    
    ! The mesh must always be in "standard" format. First read,
    ! then convert to standard.
    call tria_readTriFile2D (rtriangulation, "pre/bench1.tri", rboundary)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Write a VTK file with the mesh
    ! =================================

    ! Open / write / close
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "gmv/tutorial006b.vtk")
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    call tria_done (rtriangulation)
    call boundary_release (rboundary)
    
  end subroutine

end module
