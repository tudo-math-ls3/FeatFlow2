! This file describes the interface for a callback routine that is used
! as a hit-test routine for calculating the index arrays of a mesh region
! structure.
INTERFACE

!<subroutine>
  
  SUBROUTINE fmshregHitTest(inumCells,Dcoords,Ihit,rcollection)
  
  USE collection
  
!<description>
  ! This subroutine is called for hit-testing cells to decide whether a
  ! cell belongs to a specific mesh region or not.
!</description>
  
!<input>
  ! Number of cells (i.e. vertices,edges,faces,etc.) for which the
  ! hit-test is to be performed.
  INTEGER, INTENT(IN)                          :: inumCells
    
  ! Coordinates of the points for which the hit-test is to be performed.
  ! The dimension of the array is at least (1:idim,1:inumCells), where:
  ! -> idim is the dimension of the mesh, i.e. 1 for 1D, 2 for 2D, etc.
  ! -> inumCells is the parameter passed to this routine.
  REAL(DP), DIMENSION(:,:), INTENT(IN)         :: Dcoords
!</input>

!<output>
  ! An array that recieves the result of the hit-test.
  ! The dimension of the array is at least (1:inumCells).
  INTEGER, DIMENSION(:), INTENT(OUT)           :: Ihit
!</output>

!</inputoutput>
  ! OPTIONAL: A collection structure to provide additional information
  ! to the hit-test routine.
  TYPE(t_collection), INTENT(INOUT), OPTIONAL  :: rcollection
!</inputoutput>

!</subroutine>

  END SUBROUTINE

END INTERFACE
