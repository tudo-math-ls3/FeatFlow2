  ! This file describes the interface for a callback routine that is used
  ! as a hit-test routine for calculating the index arrays of a mesh region
  ! structure.

  interface

!<subroutine>

  subroutine fmshregHitTest(inumCells,Dcoords,Ihit,rcollection)

  use collection
  use fsystem

!<description>
  ! This subroutine is called for hit-testing cells to decide whether a
  ! cell belongs to a specific mesh region or not.
!</description>

!<input>
  ! Number of cells (i.e. vertices,edges,faces,etc.) for which the
  ! hit-test is to be performed.
  integer, intent(in)                          :: inumCells

  ! Coordinates of the points for which the hit-test is to be performed.
  ! The dimension of the array is at least (1:idim,1:inumCells), where:
  ! -> idim is the dimension of the mesh, i.e. 1 for 1D, 2 for 2D, etc.
  ! -> inumCells is the parameter passed to this routine.
  real(DP), dimension(:,:), intent(in)         :: Dcoords
!</input>

!<output>
  ! An array that recieves the result of the hit-test.
  ! The dimension of the array is at least (1:inumCells).
  integer, dimension(:), intent(out)           :: Ihit
!</output>

!</inputoutput>
  ! OPTIONAL: A collection structure to provide additional information
  ! to the hit-test routine.
  type(t_collection), intent(inout), optional  :: rcollection
!</inputoutput>

!</subroutine>

  end subroutine

end interface
