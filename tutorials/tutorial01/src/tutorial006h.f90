!##############################################################################
!# Tutorial 006h: Discretise with Q1, manually visualise a function
!##############################################################################

module tutorial006h

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use bilinearformevaluation
  
  use matrixio
  use vectorio
  use ucd

  implicit none
  private
  
  public :: start_tutorial006h

contains

  ! ***************************************************************************

  subroutine start_tutorial006h

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_vectorScalar) :: rx
    type(t_ucdExport) :: rexport

    integer :: ivt
    real(DP) :: dx,dy
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006h")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rdiscretisation
    ! which describes the discretisation.
    ! =================================

    call spdiscr_initDiscr_simple (rdiscretisation,EL_Q1_2D,rtriangulation)

    ! =================================
    ! Create a scalar vector.
    ! =================================

    ! Create a vector.
    call lsyssc_createVector (rdiscretisation,rx)
    
    ! =================================
    ! Fill the vector with data.
    ! =================================
    
    ! Get a pointer to the data.
    call lsyssc_getbase_double (rx,p_Ddata)
    
    ! Get a pointer to the point coordinates.
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    ! Set the entries of the vector according to the function
    !    u(x,y) = x^2 * y^2
    do ivt=1,rx%NEQ
      dx = p_DvertexCoords(1,ivt)
      dy = p_DvertexCoords(2,ivt)
      p_Ddata(ivt) = dx**2 * dy**2
    end do

    ! =================================
    ! Write a VTK file with the mesh
    ! and this solution vector.
    ! =================================
    call output_line ("Writing file 'post/tutorial006h.vtk'.")

    ! Open
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial006h.vtk")
                       
    ! Pass the vector as solution.
    call ucd_addVectorByVertex (rexport, "x", UCD_VAR_STANDARD, rx)
          
    ! Write / close             
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsyssc_releaseVector (rx)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
