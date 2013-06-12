!##############################################################################
!# Tutorial 008d: Create a 2x2 block system with Mass/Laplace, nonlinear coeff.
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
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock

  use scalarpde
  use derivatives
  use bilinearformevaluation  
  use domainintegration
  use collection
  
  use matrixio

  implicit none
  private
  
  public :: start_tutorial008d

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_Matrix (rdiscretisationTrial,&
                rdiscretisationTest, rform, nelements, npointsPerElement,&
                Dpoints, IdofsTrial, IdofsTest, rdomainIntSubset,&
                Dcoefficients, rcollection)

  use basicgeometry
  use collection
  use domainintegration
  use scalarpde
  use spatialdiscretisation
  use triangulation
  use fsystem

!<description>
  ! This subroutine is called during the matrix assembly. It has to compute
  ! the coefficients in front of the terms of the bilinear form.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in real coordinates.
  ! According to the terms in the bilinear form, the routine has to compute
  ! simultaneously for all these points and all the terms in the bilinear form
  ! the corresponding coefficients in front of the terms.
!</description>

!<input>
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.; trial space.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.; test space.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

  ! The bilinear form which is currently being evaluated:
  type(t_bilinearForm), intent(in) :: rform

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement

  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  ! DIMENSION(dimension,npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTrial

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in test space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! A list of all coefficients in front of all terms in the bilinear form -
  ! for all given points on all given elements.
  !   DIMENSION(itermCount,npointsPerElement,nelements)
  ! with itermCount the number of terms in the bilinear form.
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    integer :: imatrixType
    real(DP) :: dx,dy
    integer :: ivt,iel
    
    ! The collection tells us what to assemble.
    imatrixType = rcollection%IquickAccess(1)

    select case (imatrixType)
    
    ! Mass matrix
    case (1)
    
      ! Loop over all points and elements. Return the function in each point
      do iel=1,nelements
        do ivt=1,npointsPerElement
        
          ! Point coordinates
          dx = Dpoints(1,ivt,iel)
          dy = Dpoints(2,ivt,iel)
          
          ! Initialise the coefficients
          Dcoefficients(1,ivt,iel) = 1.0_DP + 16.0_DP * dx * (1.0_DP-dx) * dy * (1.0_DP-dy)
        
        end do
      end do

    ! Laplace matrix
    case (2)
    
      ! Loop over all points and elements. Return the function in each point
      do iel=1,nelements
        do ivt=1,npointsPerElement
        
          ! Point coordinates
          dx = Dpoints(1,ivt,iel)
          dy = Dpoints(2,ivt,iel)
        
          ! Initialise the coefficients. Both the same.
          Dcoefficients(1,ivt,iel) = 1.0_DP + SIN( SYS_PI * dx ) * SIN( SYS_PI * dy )
          Dcoefficients(2,ivt,iel) = Dcoefficients(1,ivt,iel)

        end do
      end do
      
    end select

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial008d

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_matrixScalar) :: rtemplateMatrix
    type(t_matrixBlock) :: rmatrix

    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_bilinearForm) :: rform
    type(t_collection) :: rcollection

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
    call lsysbl_duplicateMatrix (rtemplateMatrix,rmatrix,1,1,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    call lsysbl_duplicateMatrix (rtemplateMatrix,rmatrix,2,2,&
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
    ! Block (1,1): f Mass
    ! =================================

    ! Set up a bilinear for for Mass. There is 1 term...    
    rform%itermCount = 1
    
    ! We have nonconstant coefficients
    rform%ballCoeffConstant = .false.

    ! Term 1/2: Mass  =>  (f phi, psi)
    rform%Idescriptors(1,1) = DER_FUNC2D
    rform%Idescriptors(2,1) = DER_FUNC2D
    
    ! Via the collection tell the callback routine to assemble Mass.
    rcollection%IquickAccess(1) = 1

    ! Build the matrix block (1,1)
    call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(1,1),rcubatureInfo,&
        fcoeff_Matrix, rcollection)

    ! =================================
    ! Block (2,2): g (-Laplace)
    ! =================================

    ! Set up a bilinear for for Laplace...    
    rform%itermCount = 2
    
    ! We have constant coefficients:
    rform%ballCoeffConstant = .false.

    ! Term 1/2: -Laplace  =>  (g phi_x, psi_x)  +  (g phi_y, psi_y)
    rform%Idescriptors(1,1) = DER_DERIV2D_X
    rform%Idescriptors(2,1) = DER_DERIV2D_X
    
    rform%Idescriptors(1,2) = DER_DERIV2D_Y
    rform%Idescriptors(2,2) = DER_DERIV2D_Y
    
    ! Via the collection tell the callback routine to assemble Laplace.
    rcollection%IquickAccess(1) = 2
    
    ! Build the matrix block (2,2)
    call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(2,2),rcubatureInfo,&
        fcoeff_Matrix, rcollection)

    ! =================================
    ! Output of the matrix
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial008d_matrix.txt", "(E11.2)")

    ! Write the matrix to a MATLAB file.
    call matio_spyBlockMatrix(&
        "post/tutorial008d_matrix.m","matrix",rmatrix,.true.)
    
    ! Write the matrix to a MAPLE file
    call matio_writeBlockMatrixMaple (rmatrix, "matrix", 0, &
        "post/tutorial008d_matrix.maple", "(E11.2)")

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
