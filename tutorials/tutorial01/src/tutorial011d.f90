!##############################################################################
!# Tutorial 011d: Assemble a vector with block meth., FE-function/deriv. as coeff.
!##############################################################################

module tutorial011d

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
  
  use vectorio

  implicit none
  private
  
  public :: start_tutorial011d

contains

  !****************************************************************************

!<subroutine>

  subroutine fcoeff_rhs(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the RHS.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
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
    real(DP) :: dbasI, df,dfx, dfy, dx, dy
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_DpointCoords
    real(DP), dimension(:,:,:), pointer :: p_Df
    real(DP), dimension(:,:), pointer :: p_Dcoeff1, p_Dcoeff2
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_DpointCoords => rassemblyData%revalElementSet%p_DpointsReal

    ! We set up two vectors:
    !
    !    rhs_1 = f + sin(dx*dy*PI)
    !    rhs_2 = g + sin(dx*dy*PI)
    !
    ! with
    !
    !    f = 1 + x*y         ( given in rcoeffVectors%p_RvectorData(1) )
    !    g = 1 + f_x * f_y
    !
    ! Get the nonconstant coefficients for f. They are
    ! automatically evaluated for us, the data can be found
    ! in rcoeffVectors%p_RvectorData(:)%p_Ddata.
    !
    ! Due to ideriv=1 below, f, f_x and f_y are available here as well.
    p_Df => revalVectors%p_RvectorData(1)%p_Ddata
    
    ! In rcoeffVectors%p_RvectorData(2), we have temporary memory available.
    ! Get it, we precompute the coefficients to this.
    p_Dcoeff1 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,1)
    p_Dcoeff2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,2)
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Coordinates of the current cubature point
        dx = p_DpointCoords(1,icubp,iel)
        dy = p_DpointCoords(2,icubp,iel)

        ! Get the values of f, f_x, f_y in the cubature point.
        df  = p_Df(icubp,iel,DER_FUNC2D)
        dfx = p_Df(icubp,iel,DER_DERIV2D_X)
        dfy = p_Df(icubp,iel,DER_DERIV2D_Y)
        
        ! Precompute the coefficients.
        p_Dcoeff1(icubp,iel) = df + sin(dx*dy*SYS_PI)
        p_Dcoeff2(icubp,iel) = 1.0_DP + dfx * dfy + sin(dx*dy*SYS_PI)

      end do
      
    end do

    ! Get the data arrays of the subvectors
    p_rvectorData => RvectorData(1)
    p_DbasTest => RvectorData(1)%p_DbasTest

    ! Get pointers to the local vectors for block 1 + 2
    p_DlocalVector1 => RvectorData(1)%p_Dentry
    p_DlocalVector2 => RvectorData(2)%p_Dentry

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Psi_i:
        do idofe=1,p_rvectorData%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Psi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
          
          ! Multiply the values of the basis functions
          ! (1st derivatives) by the cubature weight and sum up
          ! into the local vectors.
          p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
              p_DcubWeight(icubp,iel) * p_Dcoeff1(icubp,iel) * dbasI
          
          p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
              p_DcubWeight(icubp,iel) * p_Dcoeff2(icubp,iel) * dbasI

        end do ! idofe
          
      end do ! icubp
    
    end do ! iel
      
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial011d

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_fev2Vectors) :: rcoeffVectors
    type(t_vectorScalar) :: rcoeffVector
    type(t_vectorBlock) :: rrhs
    
    integer :: i, ideriv
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DdataCoeff

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 011d")
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
    ! Use a 3-point Gauss Formula for the assembly
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Create a coefficient vector rcoeffVector = ( 1+x*y )
    ! =================================
    
    ! Create a scalar coefficient vector.
    call lsyssc_createVector (rspatialDiscr,rcoeffVector)

    ! Get vertex positions
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    ! Get pointers to the 1st and 2nd subvector
    call lsyssc_getbase_double (rcoeffVector, p_DdataCoeff)
    
    ! Initialise the vector 1+x*y
    do i=1,rtriangulation%NVT
      p_DdataCoeff(i) = 1.0_DP + p_DvertexCoords(1,i)*p_DvertexCoords(2,i)
    end do

    ! =================================
    ! Create a RHS vector.
    ! =================================
    
    ! Use the block discretisation to create a vector.
    call lsysbl_createVector (rblockDiscr,rrhs)
    
    ! =================================
    ! Assemble rrhs.
    ! =================================
    
    ! Clear the vector.
    call lsysbl_clearVector (rrhs)

    ! Set up a vector evaluation structure that automatically evaluates rcoeffVector
    ! in the cubature points during the assembly. 
    ! Evaluate function values + 1st derivatives, as we need the latter.
    ideriv = 1
    call fev2_addVectorToEvalList(rcoeffVectors,rcoeffVector,ideriv)

    ! Add two dummy vectors for temporary calculations.
    call fev2_addDummyVectorToEvalList(rcoeffVectors,2)

    ! Calculate the RHS.
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,fcoeff_rhs,&
        revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)

    ! Release the evaluation structure.
    call fev2_releaseVectorList(rcoeffVectors)

    ! =================================
    ! Output of the matrix structure
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rrhs, "vector", .true., 0, &
        "post/tutorial011d_vector.txt", "(E11.2)")

    ! Write the vector to a MATLAB file.
    call vecio_spyBlockVector(&
        "post/tutorial011d_vector","vector",rrhs,.true.)
    
    ! Write the vector to a MAPLE file
    call vecio_writeBlockVectorMaple (rrhs, "vector", .true., 0,&
        "post/tutorial011d_vector.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the coefficient vector
    call lsyssc_releaseVector (rcoeffVector)

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
