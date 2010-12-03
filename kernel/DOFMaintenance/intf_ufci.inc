  ! -*- mode: f90; -*-
  
  ! 
  
  interface
 
  !<subroutine>
  
    subroutine ffunction(rdata, rcollection)
    
    use ufci
    use collection
  
  !<inputoutput>
    ! The data structure of the callback function.
    type(t_ufciData), target, intent(inout) :: rdata
    
    ! Optional: A collection structure providing further information.
    type(t_collection), optional, intent(inout) :: rcollection
  !</inputoutput>
  
  !</subroutine>
  
      ! TODO: Add function code here
      
    end subroutine
  
  end interface
