!##############################################################################
!# ****************************************************************************
!# <name> ccinitgeneralparameters </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic initialisation routines for CC2D:
!# Initialisation of the main structures with parameters:
!#
!# 1.) cc_getLogFiles
!#     -> Get information for LOG files from DAT files.
!#
!# 2.) cc2d_getDAT
!#     -> Read all DAT files.
!#
!# 3.) cc_initParameters
!#     -> Init the problem structure with data from the INI/DAT files
!#
!# 4.) cc_doneParameters
!#     -> Clean up the problem structure
!#
!# </purpose>
!##############################################################################

module ccinitgeneralparameters

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  
  use collection
  use convection
  use geometry
  use ccbasic
  use ccnonstationary
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

subroutine cc_initParticleDescriptor(rPDescriptor)
  
!<description>
  ! here we initialise a particle descriptor
  ! with the values we want to use
!</description>
  
!<inputoutput>
  type(t_ParticleDescriptor), intent(inout) :: rPDescriptor
!</inputoutput>

!</subroutine>
  ! locals
  integer :: i
  real(dp) :: dx,dy,drho,drad
  
  rPDescriptor%iparticles=2
  
  allocate(rPDescriptor%pparameters(4,rPDescriptor%iparticles))
  
  drad = 0.125_dp
  drho = 1.05_dp
  dx = 1.00001_dp
  dy = 7.3_dp

  rPDescriptor%pparameters(1,1)= dx
  rPDescriptor%pparameters(2,1)= dy
  rPDescriptor%pparameters(3,1)= drad
  rPDescriptor%pparameters(4,1)= drho

  rPDescriptor%pparameters(1,2)= 0.9999_dp
  rPDescriptor%pparameters(2,2)= 6.8_dp
  rPDescriptor%pparameters(3,2)= drad
  rPDescriptor%pparameters(4,2)= drho


  
!  do i=1,rPDescriptor%iparticles
!    rPDescriptor%pparameters(1,i)= dx
!    rPDescriptor%pparameters(2,i)= dy
!    rPDescriptor%pparameters(3,i)= drad
!    rPDescriptor%pparameters(4,i)= drho
!    dx = dx + 2*drad + 0.2_dp
!    dy = dy ! + 0.5_dp
!  end do
  
end subroutine ! end cc_initParticleDescriptor

! ***************************************************************************


! ***************************************************************************
  
!<subroutine>
!subroutine cc_initParticleDescriptor(rPDescriptor)
!
!!<description>
!  ! here we initialise a particle descriptor
!  ! with the values we want to use
!!</description>
!
!!<inputoutput>
!  type(t_ParticleDescriptor), intent(inout) :: rPDescriptor
!!</inputoutput>
!
!!</subroutine>
!  ! locals
!  integer :: i,j,m,n
!  real(dp) :: dx,dy,drho,drad
!  real(dp) :: x,y,e,h,dia,rem,quot,alpha
!
!  rPDescriptor%iparticles=1
!
!  allocate(rPDescriptor%pparameters(4,rPDescriptor%iparticles))
!
!!   drad = 0.05_dp
!  drad = 0.025_dp
!  drho = 1.25_dp
!!   dx = 0.2_dp
!!   dy = 0.2_dp
!  dx = 1.0_dp
!  dy = 1.0_dp
!
!  n = rPDescriptor%iparticles
!  dia = 2*drad
!  e = 0.05_dp ! minimum distance between particles and wall
!  h = .05_dp ! height to start
!
!  ! increment of x and y
!  x = 0.0_dp + (drad + e) ! 0.0_dp is the left x-boundary value
!  y = 0.5_dp - (drad + e + h) ! 6.0_dp is the upper y-boundary value
!
!  rPDescriptor%pparameters(1,1) = x
!  rPDescriptor%pparameters(2,1) = y
!
!  m = 0
!  do i = 2,n
!    x = x + (dia+e)
!    rPDescriptor%pparameters(1,i) = x
!    rPDescriptor%pparameters(2,i) = y
!    m = m + 1;
!    if ((1.0_dp-(x+drad)) .lt. e) then ! 2.0_dp is the right x-boundary value
!        exit
!    end if
!  end do
!
!!  rem = mod(n,m) ! remainder
!!  quot = (n-rem)/m ! quotient
!!  alpha = 1.0_dp ! starting increment for alternate line
!!
!!  do i = 2,quot
!!    y = y -(dia+e)
!!      do j = 1,m
!!        rPDescriptor%pparameters(1,j + (i-1)*m) = &
!!        rPDescriptor%pparameters(1,j) - (1.0_dp + (-1.0_dp)**(i))/2*(-alpha*e)
!!!         rPDescriptor%pparameters(1,j + (i-1)*m) = &
!!!         rPDescriptor%pparameters(1,j)
!!        rPDescriptor%pparameters(2,j + (i-1)*m) = y
!!      end do
!!  end do
!!
!!  do i = 1,rem
!!    rPDescriptor%pparameters(1,m*quot + i) = &
!!    rPDescriptor%pparameters(1,i) - (1.0_dp + (-1.0_dp)**(quot+1))/2*(-alpha*e)
!!!     rPDescriptor%pparameters(1,m*quot + i) = &
!!!     rPDescriptor%pparameters(1,i);
!!    rPDescriptor%pparameters(2,m*quot + i) = y - (dia + e)
!!  end do
!!
!!  do i=1,rPDescriptor%iparticles
!!    rPDescriptor%pparameters(3,i)= drad
!!    rPDescriptor%pparameters(4,i)= drho
!!  end do
!
!
!
!end subroutine ! end cc_initParticleDescriptor

! ***************************************************************************

!<subroutine>

  subroutine cc_getLogFiles (slogfile,serrorfile,sbenchlogfile)
  
!<description>
  ! Temporarily reads the output DAT file to get the names of the output
  ! files.
!</description>

!<output>
  ! Name of the message log file.
  character(LEN=*), intent(out) :: slogfile
  
  ! Name of the error log file.
  character(LEN=*), intent(out) :: serrorfile

  ! Name of the benchmark log file.
  character(LEN=*), intent(out) :: sbenchlogfile
!</output>

!</subroutine>

    type(t_parlist) :: rparlist
    character(LEN=SYS_STRLEN) :: sstring,smaster
    logical :: bexists

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)
    
    ! Check if a command line parameter specifies the master.dat file.
    call sys_getcommandLineArg(1,smaster,sdefault='./data/master.dat')

    ! Read parameters that configure the output
    inquire(file=smaster, exist=bexists)
    
    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparlist, smaster)
    else
      call parlst_readfromfile (rparlist, './data/output.dat')
    end if
    
    ! Now the real initialisation of the output including log file stuff!
    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'smsgLog',sstring,'''''')
    read(sstring,*) slogfile

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'serrorLog',sstring,'''''')
    read(sstring,*) serrorfile

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'sbenchLog',sstring,'''''')
    read(sstring,*) sbenchlogfile
    
    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)
    
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2d_getDAT (rparamList)
  
!<description>
  ! Reads in all DAT files into the parameter list rparlist
!</description>

!<inputoutput>
  ! The parameter list where the values of the DAT files should be stored.
  ! The structure must have been initialised, the parameters are just added
  ! to the list.
  type(t_parlist), intent(inout) :: rparamList
!</inputoutput>

!</subroutine>

    logical :: bexists
    character(LEN=SYS_STRLEN) :: smaster
    
    ! Check if a command line parameter specifies the master.dat file.
    call sys_getcommandLineArg(1,smaster,sdefault='./data/master.dat')
    
    ! Read the file 'master.dat'.
    ! If that does not exist, try to manually read files with parameters from a
    ! couple of files.
    inquire(file=smaster, exist=bexists)
    
    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparamList, smaster)
    else
      ! Each 'readfromfile' command adds the parameter of the specified file
      ! to the parameter list.
      call parlst_readfromfile (rparamList, './data/discretisation.dat')
      call parlst_readfromfile (rparamList, './data/linsol_cc2d.dat')
      call parlst_readfromfile (rparamList, './data/nonlinsol_cc2d.dat')
      call parlst_readfromfile (rparamList, './data/output.dat')
      call parlst_readfromfile (rparamList, './data/paramtriang.dat')
      call parlst_readfromfile (rparamList, './data/bdconditions.dat')
      call parlst_readfromfile (rparamList, './data/timediscr.dat')
      call parlst_readfromfile (rparamList, './data/postprocessing.dat')
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initOutput (rproblem)
  
!<description>
  ! Initialises basic output settings based on the parameters in the DAT file.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initParameters (rproblem)
  
!<description>
  ! Initialises the structure rproblem with data from the initialisation
  ! files.
  !
  ! The parameters in rproblem\%rparameters are evaluated.
  ! Important parameters are written to the problem structure
  ! rproblem and/or the enclosed collection rproblem\%rcollection.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    real(DP) :: dnu
    integer :: ilvmin,ilvmax,i1

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

    ! Get the viscosity model
    ! Standard = 0 = constant viscosity
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'cviscoModel',rproblem%cviscoModel,0)
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'dviscoexponent',rproblem%dviscoexponent,2.0_DP)
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'dviscoEps',rproblem%dviscoEps,0.01_DP)

    ! Get the viscosity parameter, save it to the problem structure
    ! as well as into the collection.
    ! Note that the parameter in the DAT file is 1/nu !
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'RE',dnu,1000.0_DP)

    dnu = 1.0_DP/dnu
    rproblem%dnu = dnu
    
    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)

    ! Initialise the level in the problem structure
    if(ilvmin .le. 0) then
      rproblem%NLMIN = max(1,ilvmax+ilvmin)
    else
      rproblem%NLMIN = ilvmin
    end if
    
    rproblem%NLMAX = ilvmax
    
    ! Allocate memory for all the levels.
    allocate(rproblem%RlevelInfo(1:ilvmax))

    ! Which type of problem to discretise? (Stokes, Navier-Stokes,...)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iEquation',i1,0)
    rproblem%iequation = i1

    ! Type of subproblem (gradient tensor, deformation tensor,...)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isubEquation',i1,0)
    rproblem%isubEquation = i1

    ! Type of boundary conditions
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iBoundary',rproblem%iboundary,0)

    ! Time dependence
    call cc_initParTimeDependence (rproblem,'TIME-DISCRETISATION',&
        rproblem%rparamList)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneParameters (rproblem)
  
!<description>
  ! Cleans up parameters read from the DAT files. Removes all references to
  ! parameters from the collection rproblem\%rcollection that were
  ! set up in cc_initParameters.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! Deallocate memory
    deallocate(rproblem%RlevelInfo)

  end subroutine

end module
