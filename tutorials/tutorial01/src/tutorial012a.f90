!##############################################################################
!# Tutorial 012a: Refine a brick triangulation
!##############################################################################

module tutorial012a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use triangulation
  use meshgeneration
  use ucd

  implicit none
  private

  public :: start_tutorial012a

contains

  ! ***************************************************************************

  subroutine start_tutorial012a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_ucdExport) :: rexport
    character(LEN=SYS_STRLEN) :: spostdir

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 012a")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format to work with it.
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Write a VTK file with the mesh
    ! =================================
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call output_line ("Writing file '"//trim(spostdir)//"/tutorial012a_level1.vtk'.")
    else
      call output_line ("Writing file './post/tutorial012a_level1.vtk'.")
    end if

    ! Open / write / close
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial012a_level1.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "./post/tutorial012a_level1.vtk")
    end if
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Refine by 2-level ordering
    ! =================================

    call tria_refine2LevelOrdering(rtriangulation)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Write a VTK file with the mesh
    ! =================================
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call output_line ("Writing file '"//trim(spostdir)//"/tutorial012a_level2.vtk'.")
    else
      call output_line ("Writing file './post/tutorial012a_level2.vtk'.")
    end if

    ! Open / write / close
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial012a_level2.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "./post/tutorial012a_level2.vtk")
    end if
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Refine by 2-level ordering
    ! =================================

    call tria_refine2LevelOrdering(rtriangulation)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Write a VTK file with the mesh
    ! =================================
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call output_line ("Writing file '"//trim(spostdir)//"/tutorial012a_level3.vtk'.")
    else
      call output_line ("Writing file './post/tutorial012a_level3.vtk'.")
    end if

    ! Open / write / close
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial012a_level3.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "./post/tutorial012a_level3.vtk")
    end if
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================

    call tria_done (rtriangulation)

  end subroutine

end module
