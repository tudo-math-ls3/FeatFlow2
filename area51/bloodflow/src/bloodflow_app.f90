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
  
  implicit none 

  type(t_bloodflow) :: rbloodflow
  real(DP) :: dtime

  ! Initialization
  call bloodflow_init(rbloodflow, 'blood1.dat', 'data')

  ! Set simulation time by hand
  dtime = 0.4

  ! Evaluate the object location at simulation time
  call bloodflow_evalObject(rbloodflow, dtime)
  
  ! Evaluate the bloodflow structure
  call bloodflow_evalIndicator(rbloodflow)

  ! Write the content of the bloodflow structure to GMV file
  call bloodflow_outputStructure(rbloodflow)

  ! Finalization
  call bloodflow_done(rbloodflow)
  
end program bloodflow_app
