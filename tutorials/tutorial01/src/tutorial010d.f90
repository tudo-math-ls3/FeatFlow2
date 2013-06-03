!##############################################################################
!# Tutorial 010d: Create a 2x2 block system, nonconstant coefficients.
!##############################################################################

module tutorial010d

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
  
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use collection
  
  use matrixio

  implicit none
  private
  
  public :: start_tutorial010d

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

    ! Local variables
    real(DP) :: dbasI, dbasJ
    real(DP) :: dbasIx, dbasJx, dbasIy, dbasJy
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix22
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData11
    real(DP), dimension(:,:,:), pointer :: p_DpointCoords

    real(DP) :: df, dg, dx, dy

    ! Get a pointer to cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get pointers to local data
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial
    p_DbasTest => RmatrixData(1,1)%p_DbasTest
    
    ! Get the coordinates of the cubature points
    p_DpointCoords => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the matrix data.
    p_rmatrixData11 => RmatrixData(1,1)
    
    ! Get the matrix entries.
    p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
    p_DlocalMatrix22 => RmatrixData(2,2)%p_Dentry

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Coordinates of the current cubature point
        dx = p_DpointCoords(1,icubp,iel)
        dy = p_DpointCoords(2,icubp,iel)
        
        ! We set up two matrices:
        !
        !    A11 = f * Mass
        !    A22 = g * (-Laplace)
        !
        ! with
        !
        !    f = 1 + x^2
        !    g = 1 + y^2
        !
        ! Define f and g in the current cubature point
        df = 1.0_DP + dx**2
        dg = 1.0_DP + dy**2

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Psi_i:
        do idofe=1,p_rmatrixData11%ndofTest

          ! Fetch the contributions of the (test) basis functions Psi_i
          ! into dbasI
          dbasI  = p_DbasTest(idofe,DER_FUNC2D,icubp,iel)
          dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the (trial) basis function Phi_j:
          do jdofe=1,p_rmatrixData11%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ  = p_DbasTrial(jdofe,DER_FUNC2D,icubp,iel)
            dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
            dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.

            ! Mass to the block (1,1)
            p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * df * dbasJ*dbasI

            ! Laplace to the block (2,2)
            p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * dg * ( dbasJx*dbasIx + dbasJy*dbasIy )

          end do ! jdofe

        end do ! idofe

      end do ! icubp

    end do ! iel
        
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial010d

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_matrixScalar) :: rtemplateMatrix
    type(t_matrixBlock) :: rmatrix
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 010d")
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
    ! Use a 3-point Gauss Formula for the assembly
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Create Mass and Laplace.
    ! Use block assembly routines and a
    ! callback routine which applies
    ! the complete assembly in one step.
    ! =================================
    
    ! Clear the matrix
    call lsysbl_clearMatrix (rmatrix)
    
    ! Assemble the matrix using our callback routine above.
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,fassembleLocalMatrices,&
        rcubatureInfo=rcubatureInfo)

    ! =================================
    ! Output of the matrix structure
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial010d_matrix.txt", "(E11.2)")

    ! Write the matrix to a MATLAB file.
    call matio_spyBlockMatrix(&
        "post/tutorial010d_matrix.m","matrix",rmatrix,.true.)
    
    ! Write the matrix to a MAPLE file
    call matio_writeBlockMatrixMaple (rmatrix, "matrix", 0, &
        "post/tutorial010d_matrix.maple", "(E11.2)")

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
