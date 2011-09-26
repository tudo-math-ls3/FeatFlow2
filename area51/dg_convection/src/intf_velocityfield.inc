  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar matrices.
  
  interface
 
  
  !<subroutine>

    subroutine velocityfield (dpoint,dvelocity)
    
    use fsystem
    
  !<description>
    ! This subroutine is called during the matrix assembly of the face terms
    ! The input is the coordinates of the points dpoints
    ! The output is the velocity at this point dvelocity
  !</description>
    
  !<input>
    ! The point, where to evaluate the velocity
    real(DP), dimension(:), intent(IN) :: dpoint
  !</input>
  
  !<output>
    ! The velocity
    real(DP), dimension(size(dpoint)), intent(OUT) :: dvelocity
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
