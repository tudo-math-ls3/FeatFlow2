!##############################################################################
!# Tutorial 023d: Calculate an L2 and H1 error with the block integration
!##############################################################################

module tutorial023d

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use derivatives
  use collection
  
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  implicit none
  private
  
  public :: start_tutorial023d

contains

  !****************************************************************************

!<subroutine>

  subroutine fcalc_L2error(Dintvalues,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the L2 error of the function passed by revalVectors
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

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

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalues
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dx,dy,dval
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Initial integral value
    Dintvalues(1) = 0.0_DP

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the coordinates of the point
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Get the value of the error
        dval = p_Dfunc(icubp,iel) - 1.0_DP/exp(dx)
        
        ! Multiply the values by the cubature weight and sum up
        ! into the integral value
        Dintvalues(1) = Dintvalues(1) + p_DcubWeight(icubp,iel) * dval**2
          
      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fcalc_H1error(Dintvalues,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the H1 (semi-)error of the function passed by revalVectors
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

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

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalues
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dx,dy,dderivX,dderivY
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Initial integral value
    Dintvalues(1) = 0.0_DP

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_DderivX => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_X)
    p_DderivY => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_Y)

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the coordinates of the point
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)
          
        ! Get the value of the error
        dderivX = p_DderivX(icubp,iel) - (-1.0_DP/exp(dx))
        dderivY = p_DderivY(icubp,iel) - (0)
        
        ! Multiply the values by the cubature weight and sum up
        ! into the integral value
        Dintvalues(1) = Dintvalues(1) + p_DcubWeight(icubp,iel) * &
            ( dderivX**2 + dderivY**2 )
              
      end do ! icubp
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial023d

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_vectorBlock) :: rx
    type(t_fev2Vectors) :: rcoeffVectors

    integer :: ivt
    real(DP) :: dx, dintvalue
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 023d")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 9x9-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 8, 8)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rspatialDiscr
    ! which describes the discretisation.
    ! We generate a 1x1 system with Q1.
    ! =================================

    ! Create a spatial discretisation with Q1
    call spdiscr_initDiscr_simple (rspatialDiscr,EL_Q1_2D,rtriangulation)
    
    ! Create a block discretisation with 1 block Q1.
    call spdiscr_initBlockDiscr (rblockDiscr,rtriangulation)
    call spdiscr_appendBlockComponent (rblockDiscr,rspatialDiscr)
    call spdiscr_commitBlockDiscr (rblockDiscr)

    ! =================================
    ! Create a scalar vector.
    ! =================================

    ! Create a vector.
    call lsysbl_createVector (rblockDiscr,rx)
    
    ! =================================
    ! Fill the vector with data.
    ! =================================
    
    ! Get a pointer to the data.
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Ddata)
    
    ! Get a pointer to the point coordinates.
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    ! Set the entries of the vector according to the function
    !    u(x,y) = 1/exp(x)
    do ivt=1,rx%NEQ
      dx = p_DvertexCoords(1,ivt)
      p_Ddata(ivt) = 1.0_DP / exp(dx)
    end do

    ! =================================
    ! Define cubature formula
    ! =================================

    ! Use a Gauss 3x3 formula for the discretisation.
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Calculate the L2-error of u
    ! =================================
    
    ! Pass u via rcoeffVectors to bma_buildIntegral
    call fev2_addVectorToEvalList (rcoeffVectors,rx%RvectorBlock(1),0)

    ! Calculate ||u_h-u||_L2^2
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
        fcalc_L2error,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)

    call fev2_releaseVectorList(rcoeffVectors)
        
    ! Take the square root to get ||u_h-u||
    dintvalue = sqrt(dintvalue)
        
    call output_line ("L2-error = "//trim(sys_sdEL(dintvalue,10)))

    ! =================================
    ! Calculate the H1-(semi)error of u
    ! =================================
    
    ! Pass u and Du via rcoeffVectors to bma_buildIntegral
    call fev2_addVectorToEvalList (rcoeffVectors,rx%RvectorBlock(1),1)

    ! Calculate |u_h-u|_H1^2 = ||D(u_h-u)||_L2^2
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
        fcalc_H1error,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)

    call fev2_releaseVectorList(rcoeffVectors)
        
    ! Take the square root to get |u_h-u|_H1
    dintvalue = sqrt(dintvalue)
        
    call output_line ("H1-error = "//trim(sys_sdEL(dintvalue,10)))

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsysbl_releaseVector (rx)
    
    ! Cubature done.
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
