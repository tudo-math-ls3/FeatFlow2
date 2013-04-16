!##############################################################################
!# Tutorial 006h: Discretise with Q1, manually visualise a function
!##############################################################################

module tutorial006h

  ! Include basic Feat-2 modules
  use fsystem
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

    integer :: i,j
    real(DP) :: dx,dy
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
    call meshgen_rectangular2DQuadMesh (rtriangulation, 1.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rdiscretisation
    ! which describes the discretisation.
    ! =================================

    call spdiscr_initDiscr_simple (rdiscretisation,EL_Q1,rtriangulation)

    ! =================================
    ! Create a scalar vector.
    ! =================================

    ! Create a vector.
    call lsyssc_createVector (rdiscretisation,rx)
    
    ! =================================
    ! Fill the vector with data.
    ! =================================
    
    ! Get a pointer to the data
    call lsyssc_getbase_double (rx,p_Ddata)
    
    ! Set the entries of the vector according to the function
    !    f(x,y) = x^2 * y^2
    do i=0,4
      do j=0,4
        dx = real(j,DP) / 4.0_DP
        dy = real(i,DP) / 4.0_DP
        p_Ddata(i*5+j+1) = dx**2 * dy**2
      end do
    end do

    ! =================================
    ! Write a VTK file with the mesh
    ! and this solution vector.
    ! =================================
    call output_line ("Writing file 'gmv/tutorial006h.vtk'.")

    ! Open
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "gmv/tutorial006h.vtk")
                       
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
    
    ! Remease the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
