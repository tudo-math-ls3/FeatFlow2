!##############################################################################
!# Tutorial 006f: Discretise with Q1, create a vector
!##############################################################################

module tutorial006f

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use spatialdiscretisation
  use linearsystemscalar
  
  use vectorio

  implicit none
  private
  
  public :: start_tutorial006f

contains

  ! ***************************************************************************

  subroutine start_tutorial006f

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_vectorScalar) :: rx

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006f")
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
    ! Create a scalar vector
    ! =================================

    call lsyssc_createVector (rdiscretisation,rx)
    
    ! Fill the vector with "1.0".
    call lsyssc_clearVector (rx,1.0_DP)
    
    ! =================================
    ! Output to files. 
    ! =================================
    call output_line ("Writing vector to text files...")

    ! Write the vector to a text file.
    call vecio_writeVectorHR (rx, "vector", .true., 0, &
        "data/tutorial006f_vector.txt", "(E11.2)")

    ! Write the vector to a MATLAB file.
    call vecio_spyVector(&
        "data/tutorial006f_vector","vector",rx,.true.)
    
    ! Write the vector to a MAPLE file
    call vecio_writeVectorMaple (rx, "vector", .true., 0, &
        "data/tutorial006f_vector.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the matrix
    call lsyssc_releaseVector (rx)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
