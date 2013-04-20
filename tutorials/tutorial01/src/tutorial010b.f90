!##############################################################################
!# Tutorial 010b: Assemble a block vector.
!##############################################################################

module tutorial010b

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
  
  use collection
  use vectorio

  implicit none
  private
  
  public :: start_tutorial010b

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
    real(DP) :: dx, dy
    
    ! Get the ID of the RHS from the assembly.
    iblock = rcollection%IquickAccess(1)
    
    select case (iblock)
    
    ! -----------------------------------------------------
    ! First block. f(x,y) = 1
    case (1)
    
      do ielement = 1,nelements
        do ipoint = 1,npointsPerElement
        
          ! x/y coordinate
          dx = Dpoints(1,ipoint,ielement)
          dy = Dpoints(2,ipoint,ielement)
          
          ! Write f(x,y)=d*16*x*(1-x)*y*(1-y) into Dcoefficients(1,:,:)
          Dcoefficients(1,ipoint,ielement) = 1.0_DP
        
        end do
      end do

    ! -----------------------------------------------------
    ! 2nd block. f(x,y) = 16x(1-x)y(1-y)
    case (2)
    
      do ielement = 1,nelements
        do ipoint = 1,npointsPerElement
        
          ! x/y coordinate
          dx = Dpoints(1,ipoint,ielement)
          dy = Dpoints(2,ipoint,ielement)
          
          ! Write f(x,y)=d*16*x*(1-x)*y*(1-y) into Dcoefficients(1,:,:)
          Dcoefficients(1,ipoint,ielement) = 16.0_DP * dx * (1.0_DP - dx) * dy * (1.0_DP - dy)
        
        end do
      end do
      
    end select
    
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial010b

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_linearForm) :: rlinform
    type(t_collection) :: rcollection
    type(t_scalarCubatureInfo), target :: rcubatureInfo_G3
    type(t_vectorBlock) :: rrhs

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 010b")
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
    ! Use a 3-point Gauss Formula for the RHS
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo_G3,CUB_GEN_AUTO_G3)

    ! =================================
    ! Assemble a RHS vector.
    ! =================================
    
    ! Prepare a linear form structure for one term in the RHS: (f, phi).
    ! The test function is just phi without derivative.
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC2D
    
    ! Build the RHS using fcoeff_RHS above. Pass the block to assemble in the collection.
    ! rcollection%DquickAccess/IquickAccess/... can arbitrarily be used to pass values.
    
    ! First block
    rcollection%IquickAccess(1) = 1
    call linf_buildVectorScalar (rlinform,.true.,&
        rrhs%RvectorBlock(1),rcubatureInfo_G3,fcoeff_RHS,rcollection)

    ! 2nd block
    rcollection%IquickAccess(1) = 2
    call linf_buildVectorScalar (rlinform,.true.,&
        rrhs%RvectorBlock(2),rcubatureInfo_G3,fcoeff_RHS,rcollection)

    ! =================================
    ! Output of the vector
    ! =================================
    call output_line ("Writing vector to a text file...")
    
    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rrhs, "vector", .true., 0, &
        "post/tutorial010b_rhs.txt", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo_G3)

    ! Release the vector
    call lsysbl_releaseVector (rrhs)
    
    ! Release the block discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
