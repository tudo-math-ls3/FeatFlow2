  ! This file describes the interface for the subroutine that
  ! calculates the function values of the monitor function
  ! used in grid deformation

  interface

!<subroutine>
  subroutine def_monitorfct(DvertexCoords,Dentries)

  use triangulation
  use linearformevaluation
  use spdiscprojection
  use pprocerror
  use element
  use elementpreprocessing
  use pprocgradients
  use fsystem
  !<description>
    ! In this function we build the nodewise area distribution out
    ! of an elementwise distribution
  !</description>

  !<inputoutput>
   real(DP), dimension(:,:) :: DvertexCoords
   real(DP), dimension(:) :: Dentries
  !</inputoutput>

!</subroutine>


  end subroutine

end interface
