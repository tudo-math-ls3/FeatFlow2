!##############################################################################
!# ****************************************************************************
!# <name> cc2dmedium_method2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module solves an optimal control problem for the stationary and
!# nonstationary Navier-Stokes optimal control problem 
!#
!#  $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alphga/2||u||^2 $$
!#
!#  $$- \nu Delta(y) + y*\Nabla(y) + \Nabla p = f $$
!#  $$ \Nabla \cdot y = 0$$
!#  $$- \nu Delta(\lambda) - y*\Nabla(\lambda) + \lambda\Nabla y + \Nabla \xi = y-z $$
!#  $$ \Nabla \cdot \lambda = 0$$
!#              
!#
!# on a 2D domain for a 2D function $y=(y_1,y_2)$, a pressure $p$,
!# a dual velocity $\lambda$ and a dual pressure $\xi$. $u$ is the control
!# and $z$ a desired flow field.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# For the nonlinearity, the nonlinear solver is invoked. The
!# defect that is setted up there is preconditioned by a linear Multigrid
!# solver with a simple-VANKA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type. As coarse grid solver,
!# UMFPACK is used.
!# </purpose>
!##############################################################################

module cc2dmedium_method2

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
    
  use cc2dmediumm2basic
  use cc2dmediumm2init
  use cc2dmediumm2boundary
  use cc2dmediumm2discretisation
  use cc2dmediumm2postprocessing
  use cc2dmediumm2stationary
  
  use externalstorage
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmedium2_getLogFiles (slogfile,serrorfile)
  
!<description>
  ! Temporarily reads the output DAT file to get the names of the output
  ! files.
!</description>

!<output>
  ! Name of the message log file.
  character(LEN=*), intent(OUT) :: slogfile
  
  ! Name of the error log file.
  character(LEN=*), intent(OUT) :: serrorfile
!</output>

!</subroutine>

    type(t_parlist) :: rparlist
    character(LEN=SYS_STRLEN) :: sstring

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)

    ! Read parameters that configure the output
    call parlst_readfromfile (rparlist, trim(DIR_DATA)//'/output.dat')
    
    ! Now the real initialisation of the output including log file stuff!
    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'smsgLog',sstring,'')
    read(sstring,*) slogfile

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'serrorLog',sstring,'')
    read(sstring,*) serrorfile
    
    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)
    
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmedium2_evalParameters ()
  
!<description>
  ! Evaluates command line parameters.
  ! In the current implementation, command line parameters are passed as
  ! text file. This routine searches in the main directory for a file
  ! "cmdline.dat". If this file is found, it's opened and evaluated.
  ! Every line may contain a command line parameter in the form of
  ! a DAT file (name=value pairs).
  !
  ! Supported command line parameters:
  !   "datdirectory = [Directory, where DAT files can be found]"
!</description>

!</subroutine>

    ! local variables
    type(t_parlist) :: rparamList
    logical :: bexists
    character(SYS_STRLEN) :: sdata

    ! Figure out if the file exists.
    inquire(file='./cmdline.dat', exist=bexists)
    
    if (bexists) then
      ! Read the file
      call parlst_init (rparamList)
      call parlst_readfromfile (rparamList, './cmdline.dat')
      
      ! Evaluate parameters
      call parlst_getvalue_string_direct ( &
          rparamList, '','datdirectory', sdata, DIR_DATA)
      read(sdata,*) DIR_DATA

      call parlst_done (rparamList)
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmedium2_getDAT (rparamList)
  
!<description>
  ! Reads in all DAT files into the parameter list rparlist
!</description>

!<inputoutput>
  ! The parameter list where the values of the DAT files should be stored.
  ! The structure must have been initialised, the parameters are just added
  ! to the list.
  type(t_parlist), intent(INOUT) :: rparamList
!</inputoutput>

!</subroutine>

    logical :: bexists
    
    ! Read the file 'master.dat'.
    ! If that does not exist, try to manually read files with parameters from a
    ! couple of files.
    inquire(file=trim(DIR_DATA)//'/master.dat', exist=bexists)
    
    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/master.dat',trim(DIR_DATA))
    else
      ! Each 'readfromfile' command adds the parameter of the specified file 
      ! to the parameter list.
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/discretisation.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'//linsol_cc2d.dat')
      ! CALL parlst_readfromfile (rparamList, TRIM(DIR_DATA)//'//nonlinsol_cc2d.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'//output.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'//paramtriang.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'//bdconditions.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'//timediscr.dat')
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmedium2optc
  
!<description>
  ! This is a 'separated' Navier-Stokes solver for solving a Navier-Stokes
  ! problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the 
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! (This is somehow a cleaner implementation than using a collection!).
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly), a collection is used.
  !
  ! The following tasks are performed by the subroutines:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! A problem structure for our problem
    type(t_problem), pointer :: p_rproblem
    
    ! A structure for the solution vector and the RHS vector of the problem.
    type(t_vectorBlock) :: rvector,rrhs
    
    integer :: i
    
    ! Ok, let's start. 
    !
    ! Initialise the external storage management.
    
    call exstor_init (999,100)
    !CALL exstor_attachDirectory('./ff2storage')
    
    ! Allocate memory for the problem; it's rather large.
    allocate (p_rproblem)
    
    ! Initialise the collection
    call collct_init (p_rproblem%rcollection)

    ! Initialise the parameter list object. This creates an empty parameter list.
    call parlst_init (p_rproblem%rparamList)
    
    ! Read parameters from the INI/DAT files into the parameter list. 
    call cc2dmedium2_getDAT (p_rproblem%rparamList)
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    call cc_initOutput (p_rproblem)
    OU_LINE_LENGTH = 132
    
    ! Print the parameters 
    call parlst_info (p_rproblem%rparamList)
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    call cc_initParameters (p_rproblem)
    
    do i=1,p_rproblem%NLMAX
      call collct_addlevel_all (p_rproblem%rcollection)
    end do
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call cc_initParamTriang (p_rproblem)
    call cc_initDiscretisation (p_rproblem)    
    call cc_allocMatVec (p_rproblem,rvector,rrhs)   
    
    ! Print information about the discretisation
    call output_line ('Discretisation statistics:')
    call output_line ('--------------------------')
    do i=p_rproblem%NLMIN,p_rproblem%NLMAX
      call output_lbrk ()
      call output_line ('Level '//sys_siL(i,10))
      call output_line ('---------')
      call dof_infoDiscrBlock (p_rproblem%RlevelInfo(i)%rdiscretisation,.false.)
    end do
     
    ! On all levels, generate the static matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    call cc_generateBasicMatrices (p_rproblem)

    ! Create the solution vector -- zero or read from file.
    call cc_initInitialSolution (p_rproblem,rvector)
    
    ! Now choose the algorithm. Stationary or time-dependent simulation?
    if (p_rproblem%itimedependence .eq. 0) then
    
      ! Stationary simulation

      ! Read the (stationary) target flow.
      call cc_initTargetFlow (p_rproblem)

      ! Generate the RHS vector.
      call cc_generateBasicRHS (p_rproblem,rrhs)
      
      ! Generate discrete boundary conditions
      call cc_initDiscreteBC (p_rproblem,rvector,rrhs)

      ! Implementation of boundary conditions
      call cc_implementBC (p_rproblem,rvector=rvector,rrhs=rrhs)
    
      ! Solve the problem
      call cc_solve (p_rproblem,rvector,rrhs)
    
      ! Postprocessing
      call cc_postprocessingStationary (p_rproblem,rvector)
      
      ! Release the target flow
      call cc_doneTargetFlow (p_rproblem)
      
    else
    
      ! Time dependent simulation with explicit time stepping.
      
      ! Initialise the boundary conditions for the 0th time step, but 
      ! don't implement any boundary conditions as the nonstationary solver
      ! doesn't like this.
      call cc_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Don't read the target flow, this is done in 
      ! cc_solveNonstationaryDirect!

      ! Call the nonstationary solver to solve the problem.
      call cc_solveNonstationaryDirect (p_rproblem,rvector)
      
    end if
    
    ! (Probably) write final solution vector
    call cc_writeSolution (p_rproblem,rvector)
    
    ! Cleanup
    call cc_doneMatVec (p_rproblem,rvector,rrhs)
    call cc_doneBC (p_rproblem)
    call cc_doneDiscretisation (p_rproblem)
    call cc_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    call cc_doneParameters (p_rproblem)

    ! Release the parameter list
    call parlst_done (p_rproblem%rparamList)
    
    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ('Remaining collection statistics:')
    call output_line ('--------------------------------')
    call output_lbrk ()
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    call collct_done (p_rproblem%rcollection)
    
    deallocate(p_rproblem)
    
    ! Information about external storage usage
    call output_lbrk ()
    call exstor_info ()
    
    ! Clean up the external storage management
    call exstor_done ()
    
  end subroutine

end module
