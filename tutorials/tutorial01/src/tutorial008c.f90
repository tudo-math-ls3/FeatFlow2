!##############################################################################
!# Tutorial 008c: Create a 2x2 block system with diffusion/convection/reaction
!##############################################################################

module tutorial008c

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

  use scalarpde
  use derivatives
  use bilinearformevaluation  
  
  use matrixio

  implicit none
  private
  
  public :: start_tutorial008c

contains

  ! ***************************************************************************

  subroutine start_tutorial008c

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_matrixScalar) :: rtemplateMatrix
    type(t_matrixBlock) :: rmatrix

    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_bilinearForm) :: rform

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 008c")
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
    ! Clear the matrix
    ! =================================

    call lsysbl_clearMatrix (rmatrix)

    ! =================================
    ! Use a 3-point Gauss Formula for the assembly
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Block (1,1): -Laplace(u) + beta grad(u)
    ! =================================

    ! Set up a bilinear for for Laplace + convection. There are 4 terms...    
    rform%itermCount = 4
    
    ! We have constant coefficients, given in rform%Dcoefficients.
    rform%ballCoeffConstant = .true.

    ! Term 1/2: -Laplace  =>  (1 phi_x psi_x  +  1 phi_y psi_y)
    rform%Dcoefficients(1)  = 1.0_DP
    rform%Idescriptors(1,1) = DER_DERIV2D_X
    rform%Idescriptors(2,1) = DER_DERIV2D_X
    
    rform%Dcoefficients(2)  = 1.0_DP
    rform%Idescriptors(1,2) = DER_DERIV2D_Y
    rform%Idescriptors(2,2) = DER_DERIV2D_Y
    
    ! Term 3/4: Convection  beta*grad  =>  beta_1 phi_x psi  +  beta_2 phi_y psi
    ! We take beta=(2,0.5).
    rform%Dcoefficients(3)  = 2.0_DP
    rform%Idescriptors(1,3) = DER_DERIV2D_X
    rform%Idescriptors(2,3) = DER_FUNC2D
    
    rform%Dcoefficients(4)  = 0.5_DP
    rform%Idescriptors(1,4) = DER_DERIV2D_Y
    rform%Idescriptors(2,4) = DER_FUNC2D

    ! Build the matrix block (1,1)
    call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(1,1),rcubatureInfo)

    ! =================================
    ! Block (2,2): -Laplace(u) + u
    ! =================================

    ! Set up a bilinear for for Laplace + convection. There are 4 terms...    
    rform%itermCount = 3
    
    ! We have constant coefficients:
    rform%ballCoeffConstant = .true.

    ! Term 1/2: -Laplace  =>  (1 phi_x psi_x  +  1 phi_y psi_y)
    rform%Dcoefficients(1)  = 1.0_DP
    rform%Idescriptors(1,1) = DER_DERIV2D_X
    rform%Idescriptors(2,1) = DER_DERIV2D_X
    
    rform%Dcoefficients(2)  = 1.0_DP
    rform%Idescriptors(1,2) = DER_DERIV2D_Y
    rform%Idescriptors(2,2) = DER_DERIV2D_Y
    
    ! Term 3/4: Reaction.  =>  (1 phi psi)
    rform%Dcoefficients(3)  = 1.0_DP
    rform%Idescriptors(1,3) = DER_FUNC2D
    rform%Idescriptors(2,3) = DER_FUNC2D
    
    ! Build the matrix block (2,2)
    call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(2,2),rcubatureInfo)

    ! =================================
    ! Output of the matrix
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial008c_matrix.txt", "(E11.2)")

    ! Write the matrix to a MATLAB file.
    call matio_spyBlockMatrix(&
        "post/tutorial008c_matrix","matrix",rmatrix,.true.)
    
    ! Write the matrix to a MAPLE file
    call matio_writeBlockMatrixMaple (rmatrix, "matrix", 0, &
        "post/tutorial008c_matrix.maple", "(E11.2)")

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
