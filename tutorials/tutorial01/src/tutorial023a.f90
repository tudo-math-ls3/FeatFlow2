!##############################################################################
!# Tutorial 023a: Calculate an L2 and H1 norm with the block integration
!##############################################################################

module tutorial023a

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
  
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  implicit none
  private
  
  public :: start_tutorial023a

contains

  ! ***************************************************************************

  subroutine start_tutorial023a

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
    call output_line ("This is FEAT-2. Tutorial 023a")
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
    ! Calculate the L2-norm of u
    ! =================================
    
    ! Pass u via rcoeffVectors to bma_buildIntegral
    call fev2_addVectorToEvalList (rcoeffVectors,rx%RvectorBlock(1),0)

    ! Calculate ||u||_L2^2
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
        bma_fcalc_L2norm,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)

    call fev2_releaseVectorList(rcoeffVectors)
        
    ! Take the square root to get ||u||
    dintvalue = sqrt(dintvalue)
        
    call output_line ("L2-norm = "//trim(sys_sdL(dintvalue,10)))

    ! =================================
    ! Calculate the H1-(semi)norm of u
    ! =================================
    
    ! Pass u and Du via rcoeffVectors to bma_buildIntegral
    call fev2_addVectorToEvalList (rcoeffVectors,rx%RvectorBlock(1),1)

    ! Calculate |u|_H1^2 = ||Du||_L2^2
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
        bma_fcalc_H1norm,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)

    call fev2_releaseVectorList(rcoeffVectors)
        
    ! Take the square root to get |u|_H1
    dintvalue = sqrt(dintvalue)
        
    call output_line ("H1-norm = "//trim(sys_sdL(dintvalue,10)))

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
