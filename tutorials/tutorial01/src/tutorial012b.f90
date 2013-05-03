!##############################################################################
!# Tutorial 012b: Refine a complex mesh without and with boundary
!##############################################################################

module tutorial012b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary
  use triangulation
  use ucd

  implicit none
  private
  
  public :: start_tutorial012b

contains

  ! ***************************************************************************

  subroutine start_tutorial012b

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation, rtria1, rtria2
    type(t_ucdExport) :: rexport

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 012b")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Read the underlying domain
    ! and the mesh
    ! =================================

    call boundary_read_prm(rboundary, "pre/bench1.prm")
    
    ! The mesh must always be in "standard" format to work with it.
    ! First read, then convert to standard, based on rboundary.
    call tria_readTriFile2D (rtriangulation, "pre/bench1.tri", rboundary)
    call tria_initStandardMeshFromRaw (rtriangulation, rboundary)

    ! =================================
    ! Write a VTK file with the coarse mesh
    ! =================================
    call output_line ("Writing file 'post/tutorial012b_level1.vtk'.")

    ! Open / write / close
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial012b_level1.vtk")
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! =================================
    ! Refine twice by 2-level ordering
    ! without rboundary.
    ! Note that the circle is not
    ! represented appropriately this way.
    ! =================================
    
    ! Refine rtriangulation -> rtria1
    call tria_refine2LevelOrdering(rtriangulation,rtria1)
    call tria_initStandardMeshFromRaw (rtria1)
    
    ! Refine rtria1 again
    call tria_refine2LevelOrdering(rtria1)
    call tria_initStandardMeshFromRaw (rtria1)
    
    ! =================================
    ! Refine twice by 2-level ordering
    ! with specifying rboundary.
    ! This will give a reconstruction
    ! of the circle!
    ! =================================
    
    ! Refine rtriangulation -> rtria2
    call tria_refine2LevelOrdering(rtriangulation,rtria2,rboundary)
    call tria_initStandardMeshFromRaw (rtria2,rboundary)
    
    ! Refine rtria2 again
    call tria_refine2LevelOrdering(rtria2,rboundary=rboundary)
    call tria_initStandardMeshFromRaw (rtria2,rboundary)

    ! =================================
    ! Write the refines meshes
    ! =================================
    
    ! ----- 1ST MESH -----
    call output_line ("Writing file 'post/tutorial012b_level3a.vtk'.")

    ! Open / write / close
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtria1,&
                       "post/tutorial012b_level3a.vtk")
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! ----- 2ND MESH -----
    call output_line ("Writing file 'post/tutorial012b_level3v.vtk'.")

    ! Open / write / close
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtria2,&
                       "post/tutorial012b_level3b.vtk")
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Write the refined meshes
    ! =================================
    call output_line ("Writing file 'post/tutorial012b_level1.vtk'.")

    ! Open / write / close
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial012b_level2.vtk")
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================

    ! Release the triangulations
    call tria_done (rtria2)
    call tria_done (rtria1)
    call tria_done (rtriangulation)
    
    ! Release the boundary defintiion
    call boundary_release(rboundary)
  end subroutine

end module
