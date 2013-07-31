!##############################################################################
!# Tutorial 023g: Calculate the L2 error on different meshes, use feevaluation1
!##############################################################################

module tutorial023g

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use basicgeometry
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
  
  use feevaluation
  
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  implicit none
  private
  
  public :: start_tutorial023g

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
    real(DP), dimension(:,:), pointer :: p_Dfunc1,p_Dfunc2
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    type(t_vectorBlock), pointer :: p_rx
    integer, dimension(:), pointer :: p_Ielements
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Initial integral value
    Dintvalues(1) = 0.0_DP

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Element numbers of the elements currently in process
    p_Ielements => rassemblyData%p_IelementList
    
    ! Get the reference function
    p_rx => rcollection%p_rvectorQuickAccess1

    ! Get the data array with the values of the FEM function
    ! in the cubature points. 1st vector is u_1, 2nd vector
    ! is temporary memory for u_2 which we have to fill with values here.
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dfunc2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)
    
    ! Evaluate u_2 in the cubature points, write to p_Dfunc1.
    ! This is rather expensive as the evaluation points in the fine mesh
    ! have to be searched for!
    do iel = 1,nelements
      call fevl_evaluate (DER_FUNC, p_Dfunc2(:,iel), &
          p_rx%RvectorBlock(1), p_Dpoints(:,:,p_Ielements(iel)))      
    end do

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the coordinates of the point
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Get the value of the error
        dval = p_Dfunc1(icubp,iel) - p_Dfunc2(icubp,iel)
        
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
    real(DP) :: dx,dy,dderivX,dderivY
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dderiv1X,p_Dderiv1Y
    real(DP), dimension(:,:,:), pointer :: p_Dderiv2
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    type(t_vectorBlock), pointer :: p_rx
    integer, dimension(:), pointer :: p_Ielements
    integer, dimension(2) :: Cderiv
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Initial integral value
    Dintvalues(1) = 0.0_DP

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Element numbers of the elements currently in process
    p_Ielements => rassemblyData%p_IelementList
    
    ! Get the reference function
    p_rx => rcollection%p_rvectorQuickAccess1

    ! Get the data array with the values of the FEM function
    ! in the cubature points. 1st vector is Du_1, 2nd vector
    ! is temporary memory for Du_2 which we have to fill with values here.
    p_Dderiv1X => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_X)
    p_Dderiv1Y => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_Y)
    
    ! Space for Du2 is given as a vector field
    p_Dderiv2 => revalVectors%p_RvectorData(2)%p_DdataVec(:,:,:,1)
    
    ! Evaluate u_2 in the cubature points, write to p_Dfunc1.
    ! This is rather expensive as the evaluation points in the fine mesh
    ! have to be searched for!
    ! p_Dderiv2(1,...) receives the X-derivative, p_Dderiv2(2,...) the Y-derivative.
    Cderiv = (/DER_DERIV2D_X,DER_DERIV2D_Y/)
    do iel = 1,nelements
      call fevl_evaluate (Cderiv, p_Dderiv2(:,:,iel), &
          p_rx%RvectorBlock(1), p_Dpoints(:,:,p_Ielements(iel)))      
    end do

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the coordinates of the point
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)
          
        ! Get the value of the error
        dderivX = p_Dderiv1X(icubp,iel) - p_Dderiv2(1,icubp,iel)
        dderivY = p_Dderiv1Y(icubp,iel) - p_Dderiv2(2,icubp,iel)
        
        ! Multiply the values by the cubature weight and sum up
        ! into the integral value
        Dintvalues(1) = Dintvalues(1) + p_DcubWeight(icubp,iel) * &
            ( dderivX**2 + dderivY**2 )
              
      end do ! icubp
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial023g

    ! Declare some variables.
    type(t_triangulation) :: rtria1,rtria2
    type(t_spatialDiscretisation) :: rspatialDiscr1,rspatialDiscr2
    type(t_blockDiscretisation) :: rblockDiscr1,rblockDiscr2
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_vectorBlock), target :: rx1,rx2
    type(t_fev2Vectors) :: rcoeffVectors
    type(t_collection) :: rcollection

    integer :: ivt
    real(DP) :: dx, dintvalue
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 023g")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create cparse and fine brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 9x9-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtria1, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 8, 8)
    call tria_initStandardMeshFromRaw (rtria1)

    ! And a "fine" mesh with 33x33 vertices
    call meshgen_rectangular2DQuadMesh (rtria2, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 32, 32)
    call tria_initStandardMeshFromRaw (rtria2)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rspatialDiscr
    ! which describes the discretisation.
    ! We generate a 1x1 system with Q1.
    ! =================================

    ! -----------
    ! COARSE MESH
    ! -----------

    ! Create a spatial discretisation with Q1
    call spdiscr_initDiscr_simple (rspatialDiscr1,EL_Q1_2D,rtria1)
    
    ! Create a block discretisation with 1 block Q1.
    call spdiscr_initBlockDiscr (rblockDiscr1,rtria1)
    call spdiscr_appendBlockComponent (rblockDiscr1,rspatialDiscr1)
    call spdiscr_commitBlockDiscr (rblockDiscr1)

    ! -----------
    ! FINE MESH
    ! -----------

    ! Create a spatial discretisation with Q1
    call spdiscr_initDiscr_simple (rspatialDiscr2,EL_Q1_2D,rtria2)
    
    ! Create a block discretisation with 1 block Q1.
    call spdiscr_initBlockDiscr (rblockDiscr2,rtria2)
    call spdiscr_appendBlockComponent (rblockDiscr2,rspatialDiscr2)
    call spdiscr_commitBlockDiscr (rblockDiscr2)

    ! =================================
    ! Create vectors.
    ! =================================

    call lsysbl_createVector (rblockDiscr1,rx1)
    call lsysbl_createVector (rblockDiscr2,rx2)
    
    ! =================================
    ! Fill the vectors with data.
    ! =================================
    
    ! -----------
    ! COARSE MESH
    ! -----------

    ! Get a pointer to the data.
    call lsyssc_getbase_double (rx1%RvectorBlock(1),p_Ddata)
    
    ! Get a pointer to the point coordinates.
    call storage_getbase_double2d (rtria1%h_DvertexCoords,p_DvertexCoords)
    
    ! Set the entries of the vector according to the function
    !    u(x,y) = 1/exp(x)
    do ivt=1,rx1%NEQ
      dx = p_DvertexCoords(1,ivt)
      p_Ddata(ivt) = 1.0_DP / exp(dx)
    end do

    ! -----------
    ! FINE MESH
    ! -----------

    ! Get a pointer to the data.
    call lsyssc_getbase_double (rx2%RvectorBlock(1),p_Ddata)
    
    ! Get a pointer to the point coordinates.
    call storage_getbase_double2d (rtria2%h_DvertexCoords,p_DvertexCoords)
    
    ! Set the entries of the vector according to the function
    !    u(x,y) = 1/exp(x)
    do ivt=1,rx2%NEQ
      dx = p_DvertexCoords(1,ivt)
      p_Ddata(ivt) = 1.0_DP / exp(dx)
    end do

    ! =================================
    ! Define cubature formula
    ! =================================

    ! Use a Gauss 3x3 formula for the discretisation.
    ! Coarse mesh.
    call spdiscr_createDefCubStructure (rspatialDiscr1,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Calculate the L2-error of u_1 (coarse grid)
    ! against u_2 (fine grid)
    ! =================================
    
    ! Pass u_1 via rcoeffVectors to bma_buildIntegral
    call fev2_addVectorToEvalList (rcoeffVectors,rx1%RvectorBlock(1),0)
    
    ! Add temporary space for the values of u_2 in the cubature points
    call fev2_addDummyVectorToEvalList (rcoeffVectors)
    
    ! Pass u_2 via rcollection. It "lives" on a different mesh, so we
    ! cannot pass it directly. The evaluation is much more complicated
    ! and expensive.
    rcollection%p_rvectorQuickAccess1 => rx2

    ! Calculate ||u_1-u_2||_L2^2
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,fcalc_L2error,&
        rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)

    call fev2_releaseVectorList(rcoeffVectors)
        
    ! Take the square root to get ||u_1-u_2||
    dintvalue = sqrt(dintvalue)
        
    call output_line ("L2-error = "//trim(sys_sdEL(dintvalue,10)))

    ! =================================
    ! Calculate the H1-error of u_1 (coarse grid)
    ! against u_2 (fine grid)
    ! =================================
    
    ! Pass u_1 and Du_1 via rcoeffVectors to bma_buildIntegral
    call fev2_addVectorToEvalList (rcoeffVectors,rx1%RvectorBlock(1),1)
    
    ! Add temporary space for the values of Du_2 in the cubature points.
    call fev2_addDummyVecFieldToEvalList (rcoeffVectors,NDIM2D)
    
    ! Pass u_2 via rcollection. It "lives" on a different mesh, so we
    ! cannot pass it directly. The evaluation is much more complicated
    ! and expensive.
    rcollection%p_rvectorQuickAccess1 => rx2

    ! Calculate ||D(u_1-u_2)||_L2^2
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,fcalc_H1error,&
        rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)

    call fev2_releaseVectorList(rcoeffVectors)
        
    ! Take the square root to get ||u_1-u_2||
    dintvalue = sqrt(dintvalue)
        
    call output_line ("H1-error = "//trim(sys_sdEL(dintvalue,10)))

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vectors
    call lsysbl_releaseVector (rx2)
    call lsysbl_releaseVector (rx1)
    
    ! Cubature done.
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the discretisations
    call spdiscr_releaseBlockDiscr (rblockDiscr2)
    call spdiscr_releaseBlockDiscr (rblockDiscr1)
    call spdiscr_releaseDiscr (rspatialDiscr2)
    call spdiscr_releaseDiscr (rspatialDiscr1)

    ! Release the triangulations
    call tria_done (rtria2)
    call tria_done (rtria1)
    
  end subroutine

end module
