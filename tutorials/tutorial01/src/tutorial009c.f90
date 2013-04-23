!##############################################################################
!# Tutorial 009c: Assemble a block vector. RHS given as FEM function.
!##############################################################################

module tutorial009c

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
  use linearformevaluation
  use domainintegration
  use feevaluation
  
  use collection
  use vectorio

  implicit none
  private
  
  public :: start_tutorial009c

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_RHS (rdiscretisation, rform, &
                nelements, npointsPerElement, Dpoints, &
                IdofsTest, rdomainIntSubset, &
                Dcoefficients, rcollection)

  use basicgeometry
  use collection
  use domainintegration
  use fsystem
  use scalarpde
  use spatialdiscretisation
  use triangulation

!<description>
  ! This subroutine is called during the vector assembly. It returns
  ! the values of a RHS function f in a set of (cubature) points
  ! on a set of elements.
!</description>

!<input>
  ! Underlying discretisation structure.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! The linear form which is currently to be evaluated:
  type(t_linearForm), intent(in) :: rform

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement

  ! Array of all points (x/y coords) on all elements where the values are needed.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Array with degrees of freedom of the test function.
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! A t_domainIntSubset structure specifying more detailed assembly information.
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
  ! Optional: A collection structure to provide additional information.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! A list of all coefficients in front of all terms in the linear form -
  ! for all given points on all given elements.
  !   DIMENSION(itermCount,npointsPerElement,nelements)
  ! with itermCount the number of terms in the linear form.
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    integer :: ielement, ipoint, iblock
    type(t_vectorBlock), pointer :: p_rcoeffVector
    real(DP), dimension(:,:,:), pointer :: p_Dcoeff
    
    ! Get the coefficient vector (f_1, f_2) from the collection
    p_rcoeffVector => rcollection%p_rvectorQuickAccess1

    ! Get the temp memory provided by bilf_buildVectorScalar
    ! for the coefficients. This is reserved due to ntemp=1.
    p_Dcoeff => rdomainIntSubset%p_DtempArrays 

    ! Get the ID of the RHS from the assembly.
    iblock = rcollection%IquickAccess(1)
    
    select case (iblock)
    
    ! -----------------------------------------------------
    ! First block. f(x,y) = f_1
    case (1)

      ! Evaluafe f_1 into p_Dcoeff(:,:,1)
      call fevl_evaluate_sim (p_rcoeffVector%RvectorBlock(1),rdomainIntSubset, &
          DER_FUNC2D, p_Dcoeff(:,:,1))

      do ielement = 1,nelements
        do ipoint = 1,npointsPerElement
        
          ! Write f(x,y)=f_1 into Dcoefficients(1,:,:)
          Dcoefficients(1,ipoint,ielement) = p_Dcoeff(ipoint,ielement,1)
        
        end do
      end do

    ! -----------------------------------------------------
    ! 2nd block. f(x,y) = f_2
    case (2)
    
      ! Evaluafe f_2 into p_Dcoeff(:,:,1)
      call fevl_evaluate_sim (p_rcoeffVector%RvectorBlock(2),rdomainIntSubset, &
          DER_FUNC2D, p_Dcoeff(:,:,1))

      do ielement = 1,nelements
        do ipoint = 1,npointsPerElement
        
          ! Write f(x,y)=f_2 into Dcoefficients(1,:,:)
          Dcoefficients(1,ipoint,ielement) = p_Dcoeff(ipoint,ielement,1)
        
        end do
      end do
      
    end select
    
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial009c

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_linearForm) :: rlinform
    type(t_collection) :: rcollection
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_vectorBlock) :: rrhs
    type(t_vectorBlock), target :: rcoeffVector

    integer :: i, ntemp
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 009c")
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
    ! Create a vector with 2 blocks
    ! =================================
    
    ! Use the block discretisation to create a block vector.
    call lsysbl_createVector (rblockDiscr,rrhs)

    ! =================================
    ! Clear the vector
    ! =================================
    
    call lsysbl_clearVector (rrhs)

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
    ! Use a 3-point Gauss Formula for the RHS
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Assemble a RHS vector.
    ! =================================
    
    ! Prepare a linear form structure for one term in the RHS: (f, phi).
    ! The test function is just phi without derivative.
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC2D
    
    ! Pass rcoeffVector via the collection as nonconstant function f. 
    rcollection%p_rvectorQuickAccess1 => rcoeffVector

    ! Build the RHS using fcoeff_RHS above. Pass the block to assemble in the collection.
    ! rcollection%DquickAccess/IquickAccess/... can arbitrarily be used to pass values.
    
    ! First block. Reserve one temporary array for the coefficients.
    ntemp = 1
    rcollection%IquickAccess(1) = 1
    call linf_buildVectorScalar (rlinform,.false.,&
        rrhs%RvectorBlock(1),rcubatureInfo,fcoeff_RHS,rcollection,ntemp)

    ! 2nd block. Reserve one temporary array for the coefficients.
    ntemp = 1
    rcollection%IquickAccess(1) = 2
    call linf_buildVectorScalar (rlinform,.false.,&
        rrhs%RvectorBlock(2),rcubatureInfo,fcoeff_RHS,rcollection,ntemp)

    ! =================================
    ! Output of the vector
    ! =================================
    call output_line ("Writing vectors to text files...")
    
    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rrhs, "vector", .true., 0, &
        "post/tutorial009c_rhs.txt", "(E11.2)")

    ! Write the vector to a MATLAB file.
    call vecio_spyBlockVector(&
        "post/tutorial009c_vector","vector",rrhs,.true.)
    
    ! Write the vector to a MAPLE file
    call vecio_writeBlockVectorMaple (rrhs, "vector", .true., 0,&
        "post/tutorial009c_vector.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the vector
    call lsysbl_releaseVector (rrhs)
    
    ! Release the coefficient vector
    call lsysbl_releaseVector (rcoeffVector)

    ! Release the block discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
