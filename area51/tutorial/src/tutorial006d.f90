!##############################################################################
!# Tutorial 006d: Discretise with Q1, create a Laplace matrix
!##############################################################################

module tutorial006d

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use bilinearformevaluation
  use derivatives
  use stdoperators
  
  use matrixio

  implicit none
  private
  
  public :: start_tutorial006d

contains

  ! ***************************************************************************

  subroutine start_tutorial006d

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_matrixScalar) :: rmatrix

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006d")
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
    ! Create a CSR-matrix for the above
    ! discretisation.
    ! =================================

    ! Create the matrix structure
    call bilf_createMatrixStructure (rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    ! Allocate memory for the content
    call lsyssc_allocEmptyMatrix (rmatrix)
    
    ! =================================
    ! Discretise Laplace and Mass matrix,
    ! output to text files.
    ! =================================
    call output_line ("Writing matrix to a text file...")

    ! -------------------------------------------
    ! Discretise the Laplace matrix.
    call lsyssc_clearMatrix (rmatrix)
    
    call stdop_assembleLaplaceMatrix (rmatrix)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial006d_laplace.txt", "(E15.5)")

    ! -------------------------------------------
    ! Discretise the Mass matrix.
    call lsyssc_clearMatrix (rmatrix)
    
    call stdop_assembleSimpleMatrix (rmatrix,DER_FUNC2D,DER_FUNC2D)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial006d_mass.txt", "(E15.5)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the matrix
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
