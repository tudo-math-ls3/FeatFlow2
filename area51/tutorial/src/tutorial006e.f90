!##############################################################################
!# Tutorial 006e: Create Laplace with different cubature formulas.
!##############################################################################

module tutorial006e

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
  
  public :: start_tutorial006e

contains

  ! ***************************************************************************

  subroutine start_tutorial006e

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_matrixScalar) :: rmatrix
    type(t_scalarCubatureInfo), target :: rcubatureInfo_TRZ
    type(t_scalarCubatureInfo), target :: rcubatureInfo_G1X1
    type(t_scalarCubatureInfo), target :: rcubatureInfo_G2X2
    type(t_scalarCubatureInfo), target :: rcubatureInfo_G3X3
    type(t_scalarCubatureInfo), target :: rcubatureInfo_G4X4
    type(t_scalarCubatureInfo), target :: rcubatureInfo_G5X5

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006c")
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
    ! Create a CSR-matrix for the above
    ! discretisation.
    ! =================================

    ! Create the matrix structure
    call bilf_createMatrixStructure (rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    ! Allocate memory for the content.
    call lsyssc_allocEmptyMatrix (rmatrix)
    
    ! =================================
    ! Create a set of cubature structures
    ! that represent different cubature 
    ! formulas.
    ! =================================
    
    ! Trapezoidal rule and Gauss formulas up to Gauss 5x5.
    call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo_TRZ,CUB_GEN_AUTO_TRZ)
    call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo_G1 ,CUB_GEN_AUTO_G1)
    call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo_G2 ,CUB_GEN_AUTO_G2)
    call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo_G3 ,CUB_GEN_AUTO_G3)
    call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo_G4 ,CUB_GEN_AUTO_G4)
    call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo_G5 ,CUB_GEN_AUTO_G5)
    
    ! =================================
    ! Discretise Laplace matrix, output to text files.
    ! =================================
    call output_line ("Writing matrices to text files...")

    ! -------------------------------------------
    ! Discretise the Laplace matrix. Trapezoidal rule.
    call lsyssc_clearMatrix (rmatrix)
    call stdop_assembleLaplaceMatrix (rmatrix,rcubatureInfo=rcubatureInfo_TRZ)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "data/tutorial006e_trz.txt", "(E15.5)")

    ! -------------------------------------------
    ! Discretise the Laplace matrix. Gauss 1x1.
    call lsyssc_clearMatrix (rmatrix)
    call stdop_assembleLaplaceMatrix (rmatrix,rcubatureInfo=rcubatureInfo_G1)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "data/tutorial006e_g1.txt", "(E15.5)")

    ! -------------------------------------------
    ! Discretise the Laplace matrix. Gauss 2x2
    call lsyssc_clearMatrix (rmatrix)
    call stdop_assembleLaplaceMatrix (rmatrix,rcubatureInfo=rcubatureInfo_G2)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "data/tutorial006e_g2.txt", "(E15.5)")

    ! -------------------------------------------
    ! Discretise the Laplace matrix. Gauss 3x3.
    call lsyssc_clearMatrix (rmatrix)
    call stdop_assembleLaplaceMatrix (rmatrix,rcubatureInfo=rcubatureInfo_G3)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "data/tutorial006e_g3.txt", "(E15.5)")

    ! -------------------------------------------
    ! Discretise the Laplace matrix. Gauss 5x5.
    call lsyssc_clearMatrix (rmatrix)
    call stdop_assembleLaplaceMatrix (rmatrix,rcubatureInfo=rcubatureInfo_G4)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "data/tutorial006e_g4.txt", "(E15.5)")

    ! -------------------------------------------
    ! Discretise the Laplace matrix. Gauss 5x5.
    call lsyssc_clearMatrix (rmatrix)
    call stdop_assembleLaplaceMatrix (rmatrix,rcubatureInfo=rcubatureInfo_G5)
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "data/tutorial006e_g5.txt", "(E15.5)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the matrix
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release the cubature formula structures.
    call spdiscr_releaseCubStructure (rcubatureInfo_TRZ)
    call spdiscr_releaseCubStructure (rcubatureInfo_G1)
    call spdiscr_releaseCubStructure (rcubatureInfo_G2)
    call spdiscr_releaseCubStructure (rcubatureInfo_G3)
    call spdiscr_releaseCubStructure (rcubatureInfo_G4)
    call spdiscr_releaseCubStructure (rcubatureInfo_G5)
    
    ! Remease the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
