!##############################################################################
!# ****************************************************************************
!# <name> bloodflow_app </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the main application which calls the subroutines.
!# </purpose>
!##############################################################################

program bloodflow_app

  use bloodflow
  use fsystem
  use paramlist
  
  implicit none 

  type(t_bloodflow) :: rbloodflow
  real(DP) :: dtime,dtimeStop,dtimeStep

  ! Initialization
  call bloodflow_init(rbloodflow, 'blood1.dat', 'data')

  ! Get time stepping values
  call parlst_getvalue_double(rbloodflow%rparlist, 'Timestepping', 'dtimeStart', dtime)
  call parlst_getvalue_double(rbloodflow%rparlist, 'Timestepping', 'dtimeStop',  dtimeStop)
  call parlst_getvalue_double(rbloodflow%rparlist, 'Timestepping', 'dtimeStep',  dtimeStep)

  ! Time stepping loop
  do while (dtime+dtimeStep .le. dtimeStop)
  
    ! Update simulation time
    dtime = dtime+dtimeStep

    ! Evaluate the object location at simulation time
    call bloodflow_evalObject(rbloodflow, dtime)
    
    ! Evaluate the bloodflow structure
    call bloodflow_evalIndicator(rbloodflow)
    
    ! Write the content of the bloodflow structure to GMV file
!    call bloodflow_outputStructure(rbloodflow)

  end do
  
  ! Finalization
  call bloodflow_done(rbloodflow)
  
end program bloodflow_app
