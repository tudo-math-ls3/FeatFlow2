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
    character(LEN=SYS_STRLEN) :: spredir,spostdir

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006b")
    call output_separator (OU_SEP_MINUS)
    
        if (sys_getenv_string("POSTDIR",spostdir)) then
      call output_line ("Writing file "//trim(spostdir)//"'tutorial006b.vtk'.")
    else
      call output_line ("Writing file './post/tutorial006b.vtk'.")
    end if

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
      call tria_readTriFile2D (rtriangulation, trim(spredir)//"/bench1.tri", rboundary)
    else
      call tria_readTriFile2D (rtriangulation, "pre/bench1.tri", rboundary)
    end if
    call tria_initStandardMeshFromRaw (rtriangulation, rboundary)

    ! =================================
    ! Write a VTK file with the mesh
    ! =================================

    ! Open / write / close
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial006b.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "post/tutorial006b.vtk")
    end if
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    call tria_done (rtriangulation)
    call boundary_release (rboundary)
    
  end subroutine

end module
