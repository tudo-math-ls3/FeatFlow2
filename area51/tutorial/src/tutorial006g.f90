!##############################################################################
!# Tutorial 006g: Discretise with Q1, access entries in a matrix and a vector.
!##############################################################################

module tutorial006g

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

  implicit none
  private
  
  public :: start_tutorial006g

contains

  ! ***************************************************************************

  subroutine start_tutorial006g

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_vectorScalar) :: rx
    type(t_matrixScalar) :: rmatrix
    
    integer :: i
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Kcol, p_Kld

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006g")
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

    call spdiscr_initDiscr_simple (rdiscretisation,EL_Q1_2D,rtriangulation)

    ! =================================
    ! Create a scalar vector.
    ! Create a CSR matrix corresponding
    ! to the Q1 space.
    ! =================================

    ! Create a vector.
    call lsyssc_createVector (rdiscretisation,rx)
    
    ! Create the matrix structure.
    call bilf_createMatrixStructure (rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    ! Allocate memory for the matrix entries.
    call lsyssc_allocEmptyMatrix (rmatrix)
    
    ! =================================
    ! Fill the vector with data.
    ! =================================
    
    ! Get a pointer to the data
    call lsyssc_getbase_double (rx,p_Ddata)
    
    ! Initialise the data: 1,2,3,...
    do i=1,rx%NEQ
      p_Ddata(i) = real(i,DP)
    end do

    ! =================================
    ! Fill the matrix with data.
    ! =================================
    
    ! Get a pointer to the matrix data, the column numbers
    ! and the row indices that are usesd in the CSR format.
    call lsyssc_getbase_double (rmatrix,p_Ddata)
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)

    ! Initialise the data: Column number.
    do i=1,rmatrix%NA
      p_Ddata(i) = real( p_Kcol(i) ,DP )
    end do
    
    ! Fill the first row with 1.0.
    do i=p_Kld(1), p_Kld(2)-1
      p_Ddata(i) = 1.0_DP
    end do

    ! Fill the last row with "NEQ".
    do i=p_Kld(rmatrix%NEQ), p_Kld(rmatrix%NEQ+1)-1
      p_Ddata(i) = real( rmatrix%NEQ, DP )
    end do

    ! =================================
    ! Output to files. 
    ! =================================
    call output_line ("Writing to text files...")

    ! Write the vector to a text file.
    call vecio_writeVectorHR (rx, "vector", .true., 0, &
        "post/tutorial006g_vector.txt", "(E11.2)")

    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial006g_matrix.txt", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsyssc_releaseVector (rx)
    
    ! Release the matrix
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
