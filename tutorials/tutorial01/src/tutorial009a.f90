!##############################################################################
!# Tutorial 009a: Create a 2x2 block system / vector
!##############################################################################

module tutorial009a

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  
  use vectorio

  implicit none
  private
  
  public :: start_tutorial009a

contains

  ! ***************************************************************************

  subroutine start_tutorial009a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    
    type(t_vectorBlock) :: rx
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 009a")
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
    call lsysbl_createVector (rblockDiscr,rx)

    ! =================================
    ! Clear the vector
    ! =================================
    
    call lsysbl_clearVector (rx)

    ! =================================
    ! Fill the vector with data.
    ! 1st block: 1, 2, 3,...
    ! 2nd block: NEQ, NEQ-1, NEQ-2,...
    ! =================================
    
    ! Get a pointer to the data of the first block
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Ddata)
    
    ! Initialise the data: 1,2,3,...
    do i=1,rx%RvectorBlock(1)%NEQ
      p_Ddata(i) = real(i,DP)
    end do

    ! Get a pointer to the data of the second block
    call lsyssc_getbase_double (rx%RvectorBlock(2),p_Ddata)
    
    ! Initialise the data: NEQ, NEQ-1,....
    do i=rx%RvectorBlock(1)%NEQ,1,-1
      p_Ddata(i) = real(i,DP)
    end do

    ! =================================
    ! Output of the vector
    ! =================================
    call output_line ("Writing vectors to texts file...")
    
    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rx, "vector", .true., 0, &
        "post/tutorial009a_vector.txt", "(E11.2)")

    ! Write the vector to a MATLAB file.
    call vecio_spyBlockVector(&
        "post/tutorial009a_vector","vector",rx,.true.)
    
    ! Write the vector to a MAPLE file
    call vecio_writeBlockVectorMaple (rx, "vector", .true., 0,&
        "post/tutorial009a_vector.maple", "(E11.2)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsysbl_releaseVector (rx)
    
    ! Release the block discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
