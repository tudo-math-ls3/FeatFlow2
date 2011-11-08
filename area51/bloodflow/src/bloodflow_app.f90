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
  use genoutput
  use paramlist
  use triangulation
  
  implicit none

  type(t_bloodflow) :: rbloodflow
  character(LEN=SYS_STRLEN) :: sparamfile
  real(DP) :: dtime,dtimeStop,dtimeStep

  ! Initialization
  if (command_argument_count() .eq. 0) then
    call output_line('Missing parameter file!',&
        OU_CLASS_ERROR, OU_MODE_STD, 'bloodflow')
    call sys_halt()
  else
    call get_command_argument(command_argument_count(), sparamfile)
    sparamfile = trim(adjustl(sparamfile))
  end if
  call bloodflow_init(rbloodflow, sparamfile, 'data')

  ! Get time stepping values
  call parlst_getvalue_double(rbloodflow%rparlist, 'Timestepping', 'dtimeStart', dtime)
  call parlst_getvalue_double(rbloodflow%rparlist, 'Timestepping', 'dtimeStop',  dtimeStop)
  call parlst_getvalue_double(rbloodflow%rparlist, 'Timestepping', 'dtimeStep',  dtimeStep)

  ! Time stepping loop
  do while (dtime+dtimeStep .le. dtimeStop)
  
    ! Update simulation time
    dtime = dtime+dtimeStep

    call output_line('Simulation time: '//sys_sdE(dtime,3))

    ! Evaluate the object location at simulation time
    call bloodflow_evalObject(rbloodflow, dtime)
    
    ! Duplicate triangulation of base mesh
    call tria_duplicate(rbloodflow%rtriangulationBase, rbloodflow%rtriangulation, &
      TR_SHARE_NONE, .false.)

    ! Adapt the grid to the thin object
    call bloodflow_adaptObject(rbloodflow)
    
    ! Write the content of the bloodflow structure to GMV file
    call bloodflow_outputStructure(rbloodflow, dtime)

  end do
  
  ! Finalization
  call bloodflow_done(rbloodflow)
  
end program bloodflow_app
