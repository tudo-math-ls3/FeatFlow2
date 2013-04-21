!##############################################################################
!# Tutorial 008e: Create a 2x2 block system with Mass/Laplace, nonlinear coeff.
!##############################################################################

module tutorial008e

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
  use feevaluation
  use collection
  
  use matrixio

  implicit none
  private
  
  public :: start_tutorial008e

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_Matrix (rdiscretisationTrial,&
                rdiscretisationTest, rform, nelements, npointsPerElement,&
                Dpoints, IdofsTrial, IdofsTest, rdomainIntSubset,&
                Dcoefficients, rcollection)

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
    integer :: ivt,iel
    type(t_vectorBlock), pointer :: p_rcoeffVector
    real(DP), dimension(:,:,:), pointer :: p_Dcoeff
    
    ! The collection tells us what to assemble.
    imatrixType = rcollection%IquickAccess(1)
    
    ! Get the coefficient vector from the collection
    p_rcoeffVector => rcollection%p_rvectorQuickAccess1
    
    ! Get the temp memory provided by bilf_buildMatrixScalar
    ! for the coefficients. They were reserved due to ntemp=2.
    p_Dcoeff => rdomainIntSubset%p_DtempArrays 
    
    ! --------------------------------------------------
    ! Calculate the coefficients in the cubature points.
    ! Use the two temporary arrays provided by bilf_buildMatrixScalar.
    !
    ! 1st component
    call fevl_evaluate_sim (p_rcoeffVector%RvectorBlock(1),rdomainIntSubset, &
        DER_FUNC2D, p_Dcoeff(:,:,1))
    
    ! 2nd component
    call fevl_evaluate_sim (p_rcoeffVector%RvectorBlock(2),rdomainIntSubset, &
        DER_FUNC2D, p_Dcoeff(:,:,2))

    ! Assemble
    select case (imatrixType)
    
    ! Mass matrix
    case (1)
    
      ! Loop over all points and elements. Return the function in each point
      do iel=1,nelements
        do ivt=1,npointsPerElement
        
          ! Initialise the coefficients by the 1st coordinate:
          Dcoefficients(1,ivt,iel) = p_Dcoeff(ivt,iel,1)
        
        end do
      end do

    ! Laplace matrix
    case (2)
    
      ! Loop over all points and elements. Return the function in each point
      do iel=1,nelements
        do ivt=1,npointsPerElement
        
          ! Initialise the coefficients by the 2nd coordinate.
          ! Both the same.
          Dcoefficients(1,ivt,iel) = p_Dcoeff(ivt,iel,2)   ! for (g phi_x,psi_x)
          Dcoefficients(2,ivt,iel) = p_Dcoeff(ivt,iel,2)   ! for (g phi_y,psi_y)

        end do
      end do
      
    end select

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial008e

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_matrixScalar) :: rtemplateMatrix
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock), target :: rcoeffVector

    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_bilinearForm) :: rform
    type(t_collection) :: rcollection
    
    integer :: i,ntemp
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 008e")
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
    ! Clear the matrix
    ! =================================

    call lsysbl_clearMatrix (rmatrix)

    ! =================================
    ! Use a 3-point Gauss Formula for the assembly
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Create a coefficient block vector rcoeffVector = ( 1+x^2, 1+y^2 )
    ! =================================
    
    ! Create a vector
    call lsysbl_createVector (rblockDiscr,rcoeffVector)

    ! Get vertex positions
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    ! Get pointers to the 1st and 2nd subvector
    call lsyssc_getbase_double (rcoeffVector%RvectorBlock(1), p_Ddata1)
    call lsyssc_getbase_double (rcoeffVector%RvectorBlock(2), p_Ddata2)
    
    ! Initialise the vector
    do i=1,rtriangulation%NVT
      p_Ddata1(i) = 1.0_DP + p_DvertexCoords(1,i)**2
      p_Ddata2(i) = 1.0_DP + p_DvertexCoords(2,i)**2
    end do

    ! =================================
    ! Block (1,1): f Mass
    ! =================================

    ! Set up a bilinear for for Mass. There is 1 term...    
    rform%itermCount = 1
    
    ! We have nonconstant coefficients
    rform%ballCoeffConstant = .false.

    ! Term 1/2: Mass  =>  (f phi, psi)
    rform%Dcoefficients(1)  = 1.0_DP
    rform%Idescriptors(1,1) = DER_FUNC2D
    rform%Idescriptors(2,1) = DER_FUNC2D
    
    ! Via the collection tell the callback routine to assemble Mass.
    ! Pass rcoeffVector via the collection as nonconstant coefficient. 
    rcollection%IquickAccess(1) = 1
    rcollection%p_rvectorQuickAccess1 => rcoeffVector

    ! Build the matrix block (1,1). Reserve two temporary arrays for the coefficients.
    ntemp = 2
    call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(1,1),rcubatureInfo,&
        fcoeff_Matrix, rcollection, ntemp)

    ! =================================
    ! Block (2,2): g (-Laplace)
    ! =================================

    ! Set up a bilinear for for Laplace...    
    rform%itermCount = 2
    
    ! We have constant coefficients:
    rform%ballCoeffConstant = .false.

    ! Term 1/2: -Laplace  =>  (g phi_x, psi_x)  +  (g phi_y, psi_y)
    rform%Dcoefficients(1)  = 1.0_DP
    rform%Idescriptors(1,1) = DER_DERIV2D_X
    rform%Idescriptors(2,1) = DER_DERIV2D_X
    
    rform%Dcoefficients(2)  = 1.0_DP
    rform%Idescriptors(1,2) = DER_DERIV2D_Y
    rform%Idescriptors(2,2) = DER_DERIV2D_Y
    
    ! Via the collection tell the callback routine to assemble Laplace.
    ! Pass rvector via the collection as nonconstant coefficient. 
    rcollection%IquickAccess(1) = 2
    rcollection%p_rvectorQuickAccess1 => rcoeffVector
    
    ! Build the matrix block (2,2). Reserve two temporary arrays for the coefficients.
    ntemp = 2
    call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(2,2),rcubatureInfo,&
        fcoeff_Matrix, rcollection, ntemp)

    ! =================================
    ! Output of the matrix
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial008e_matrix.txt", "(E11.2)")

    ! Write the matrix to a MATLAB file.
    call matio_spyBlockMatrix(&
        "post/tutorial008e_matrix","matrix",rmatrix,.true.)
    
    ! Write the matrix to a MAPLE file
    call matio_writeBlockMatrixMaple (rmatrix, "matrix", 0, &
        "post/tutorial008e_matrix.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the coefficient vector
    call lsysbl_releaseVector (rcoeffVector)

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
