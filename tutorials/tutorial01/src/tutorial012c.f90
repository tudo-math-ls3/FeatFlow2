!##############################################################################
!# Tutorial 012c: Fast refinement
!##############################################################################

module tutorial012c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary
  use triangulation
  use ucd

  implicit none
  private
  
  public :: start_tutorial012c

contains

  ! ***************************************************************************

  subroutine start_tutorial012c

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    type(t_ucdExport) :: rexport
    character(LEN=SYS_STRLEN) :: spredir,spostdir

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 012c")
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
    ! First read, then convert to standard.
    if (sys_getenv_string("PREDIR",spredir)) then
      call tria_readTriFile2D (rtriangulation, trim(spredir)//"/bench1.tri", rboundary)
    else
      call tria_readTriFile2D (rtriangulation, "pre/bench1.tri", rboundary)
    end if
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! =================================
    ! Write a VTK file with the coarse mesh
    ! =================================
    call output_line ("Writing file 'post/tutorial012c_level1.vtk'.")

    ! Open / write / close
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial012c_level1.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "post/tutorial012c_level1.vtk")
    end if
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! =================================
    ! 4x quick refinement by 2-level ordering
    ! including rboundary.
    ! =================================
    
    call tria_quickRefine2LevelOrdering(4, rtriangulation, rboundary)
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! =================================
    ! Write the refined mesh
    ! =================================
    
    call output_line ("Writing file 'post/tutorial012c_level5.vtk'.")

    ! Open / write / close
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial012c_level5.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "post/tutorial012c_level5.vtk")
    end if
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! =================================
    ! Cleanup
    ! =================================

    ! Release the triangulation    
    call tria_done (rtriangulation)
    
    ! Release the boundary definition
    call boundary_release(rboundary)
  end subroutine

end module
