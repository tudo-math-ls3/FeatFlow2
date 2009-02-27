  ! -*- mode: f90; -*-
  
  ! This file describes the interface for the subroutine that
  ! calculates the function values of the monitor function
  ! used in grid deformation

  interface

!<subroutine>  
  SUBROUTINE def_monitorfct(DvertexCoords,Dentries)
  
  USE triangulation
  USE linearformevaluation
  USE spdiscprojection
  USE pprocerror
  USE element
  USE elementpreprocessing
  USE pprocgradients
  
  !<description>
    ! In this function we build the nodewise area distribution out 
    ! of an elementwise distribution
  !</description>

  !<inputoutput>
   REAL(DP), DIMENSION(:,:) :: DvertexCoords
   REAL(DP), DIMENSION(:) :: Dentries
  !</inputoutput>

!</subroutine>


  end subroutine

end interface