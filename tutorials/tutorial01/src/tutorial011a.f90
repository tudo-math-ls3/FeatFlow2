!##############################################################################
!# Tutorial 011a: Create a 2x2 block system, assemble a vector with block meth.
!##############################################################################

module tutorial011a

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
  
  public :: start_tutorial011a

contains

  ! ***************************************************************************

  subroutine start_tutorial011a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_vectorBlock) :: rrhs
    type(t_collection) :: rcollection

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 011a")
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
    ! Use a 3-point Gauss Formula for the assembly
    ! =================================
    
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Create a RHS vector.
    ! =================================
    
    ! Use the block discretisation to create a vector.
    call lsysbl_createVector (rblockDiscr,rrhs)
    
    ! =================================
    ! Assemble rrhs=(f,g) with
    !   f=1,
    !   g=32*y*(1-y)+32*x*(1-x)
    ! =================================
    
    ! Clear the vector.
    call lsysbl_clearVector (rrhs)

    ! Calculate f=1 to the 1st component.
    rcollection%DquickAccess(1) = 1.0_DP  ! Multiplier
    rcollection%IquickAccess(1) = 1       ! Component
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,bma_fcalc_rhsConst,rcollection,&
        rcubatureInfo=rcubatureInfo)
        
    ! Calculate g = 1 * 32*y*(1-y)+32*x*(1-x) to the 2nd component.
    rcollection%DquickAccess(1) = 1.0_DP  ! Multiplier
    rcollection%IquickAccess(1) = 2       ! Component
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,bma_fcalc_rhsBubble,rcollection,&
        rcubatureInfo=rcubatureInfo)

    ! =================================
    ! Output of the matrix structure
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rrhs, "vector", .true., 0, &
        "post/tutorial011a_vector.txt", "(E11.2)")

    ! Write the vector to a MATLAB file.
    call vecio_spyBlockVector(&
        "post/tutorial011a_vector","vector",rrhs,.true.)
    
    ! Write the vector to a MAPLE file
    call vecio_writeBlockVectorMaple (rrhs, "vector", .true., 0,&
        "post/tutorial011a_vector.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo)

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
