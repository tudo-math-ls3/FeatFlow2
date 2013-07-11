!##############################################################################
!# Tutorial 010a: Create a 2x2 block system with Mass and Laplace on the diag.
!##############################################################################

module tutorial010a

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use derivatives
  use element
  use cubature
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
  
  public :: start_tutorial010a

contains

  ! ***************************************************************************

  subroutine start_tutorial010a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_matrixScalar) :: rtemplateMatrix
    type(t_matrixBlock) :: rmatrix
    
    type(t_collection) :: rcollection

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 010a")
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

    ! Set up a block discretisation for two equations. Both equations are Q1.
    call spdiscr_initBlockDiscr (rblockDiscr,rtriangulation)
    call spdiscr_appendBlockComponent (rblockDiscr,rspatialDiscr)
    call spdiscr_appendBlockComponent (rblockDiscr,rspatialDiscr)
    call spdiscr_commitBlockDiscr (rblockDiscr)
    
    ! =================================
    ! Use a 3-point Gauss Formula for the assembly
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Create a 2x2 block matrix with
    ! entries in (1,1) and (2,2)
    ! =================================
    
    ! Create a template matrix for the FEM space. CSR structure.
    call bilf_createMatrixStructure (rspatialDiscr,LSYSSC_MATRIX9,rtemplateMatrix)

    ! Use the block discretisation to create a basic system matrix.
    call lsysbl_createMatrix (rblockDiscr,rmatrix)

    ! Create (1,1) and (2,2) using the template matrix.
    ! Structure is "shared" (no new memory is allocated), content is allocated.
    call lsysbl_duplicateMatrix (rtemplateMatrix,rmatrix,1,1,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    call lsysbl_duplicateMatrix (rtemplateMatrix,rmatrix,2,2,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    ! =================================
    ! Create Mass and Laplace.
    ! Use default block assembly routines.
    ! =================================
    
    ! Clear the matrix
    call lsysbl_clearMatrix (rmatrix)
    
    ! Put the mass matrix to (1,1).
    rcollection%DquickAccess(1) = 1.0_DP    ! Multiplier
    rcollection%IquickAccess(1) = 1         ! y-position
    rcollection%IquickAccess(2) = 1         ! x-position
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,bma_fcalc_mass,rcollection,&
        rcubatureInfo=rcubatureInfo)
    
    ! Put the Laplace matrix to (2,2)
    rcollection%DquickAccess(1) = 1.0_DP    ! Multiplier
    rcollection%IquickAccess(1) = 2         ! y-position
    rcollection%IquickAccess(2) = 2         ! x-position
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,bma_fcalc_laplace,rcollection,&
        rcubatureInfo=rcubatureInfo)

    ! =================================
    ! Output of the matrix structure
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial010a_matrix.txt", "(E11.2)")

    ! Write the matrix to a MATLAB file.
    call matio_spyBlockMatrix(&
        "post/tutorial010a_matrix.m","matrix",rmatrix,.true.)
    
    ! Write the matrix to a MAPLE file
    call matio_writeBlockMatrixMaple (rmatrix, "matrix", 0, &
        "post/tutorial010a_matrix.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo)

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
