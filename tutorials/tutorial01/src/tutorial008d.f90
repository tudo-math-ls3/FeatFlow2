!##############################################################################
!# Tutorial 008d: Create a 2x2 block system with Mass and Laplace on the diag.
!##############################################################################

module tutorial008d

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use derivatives
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use stdoperators
  
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use collection
  
  use matrixio

  implicit none
  private
  
  public :: start_tutorial008d

contains

  !****************************************************************************

!<subroutine>

  subroutine fassembleLocalMatrices(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Callback routine. Calculates the local matrices of an operator.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! local variables
    type(t_collection) :: rlocalCollect

    ! Put the mass matrix to (1,1).
    rlocalCollect%DquickAccess(1) = 1.0_DP    ! Multiplier
    rlocalCollect%IquickAccess(1) = 1         ! x-position
    rlocalCollect%IquickAccess(2) = 1         ! y-position
    
    ! Call the assembly routine
    call bma_fcalc_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,revalVectors,rlocalCollect)
    
    ! Put the Laplace matrix to (2,2)
    rlocalCollect%DquickAccess(1) = 1.0_DP    ! Multiplier
    rlocalCollect%IquickAccess(1) = 2         ! x-position
    rlocalCollect%IquickAccess(2) = 2         ! y-position
    
    ! Call the assembly routine
    call bma_fcalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,revalVectors,rlocalCollect)

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial008d

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_matrixScalar) :: rtemplateMatrix
    type(t_matrixBlock) :: rmatrix
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 008d")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Create a block discretisation that encapsules
    ! two equations, both discretised with Q1.
    ! =================================

    ! Set up a Q1 discretisation
    call spdiscr_initDiscr_simple (rspatialDiscr,EL_Q1_2D,rtriangulation)

    ! Set up a block discretisation for two equations.
    call spdiscr_initBlockDiscr (rblockDiscr,2,rtriangulation)

    ! Both equations are Q1. Initialise with the Q1 discretisation.
    call spdiscr_duplicateDiscrSc (rspatialDiscr,rblockDiscr%RspatialDiscr(1))
    call spdiscr_duplicateDiscrSc (rspatialDiscr,rblockDiscr%RspatialDiscr(2))
    
    ! =================================
    ! Create a 2x2 block matrix with
    ! entries in (1,1) and (2,2)
    ! =================================
    
    ! Create a template matrix for the FEM space. CSR structure.
    call bilf_createMatrixStructure (rspatialDiscr,LSYSSC_MATRIX9,rtemplateMatrix)

    ! Use the block discretisation to create a basic system matrix.
    call lsysbl_createMatBlockByDiscr (rblockDiscr,rmatrix)

    ! Create (1,1) and (2,2) using the template matrix.
    ! Structure is "shared" (no new memory is allocated), content is allocated.
    call lsyssc_duplicateMatrix (rtemplateMatrix,rmatrix%RmatrixBlock(1,1),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    call lsyssc_duplicateMatrix (rtemplateMatrix,rmatrix%RmatrixBlock(2,2),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    ! =================================
    ! Create Mass and Laplace.
    ! Use block assembly routines and a
    ! callback routine.
    ! =================================
    
    ! Clear the matrix
    call lsysbl_clearMatrix (rmatrix)
    
    ! Assemble the matrix using our callback routine above.
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,fassembleLocalMatrices)

    ! =================================
    ! Output of the matrix structure
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial008d_matrix.txt", "(E11.2)")

    ! Write the matrix to a MATLAB file.
    call matio_spyBlockMatrix(&
        "post/tutorial008d_matrix","matrix",rmatrix,.true.)
    
    ! Write the matrix to a MAPLE file
    call matio_writeBlockMatrixMaple (rmatrix, "matrix", 0, &
        "post/tutorial008d_matrix.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the matrix
    call lsysbl_releaseMatrix (rmatrix)
    
    ! ... and the FEM template matrix
    call lsyssc_releaseMatrix (rtemplateMatrix)
    
    ! Release the block discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
