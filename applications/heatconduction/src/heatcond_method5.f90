!##############################################################################
!# ****************************************************************************
!# <name> heatcond_method5 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a simple
!# heat conduction problem with constant coefficients
!# on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# The different tasks of setting up the triangulation, discretisation, solver,
!# etc. are separated into different sub-modules as follows:
!#
!# - heatcond_basic
!#   -> Defines a basic problem structure that contains all parameters for
!#      the simulation
!#
!# - heatcond_callback
!#   -> Contains user-defined routines for setting up the RHS or for
!#      configuring values on the boundary.
!#
!# - heatcond_partridiscr
!#   -> Contains routines to read in the parametrisation, triangulation and
!#      to set up the discretisation (element, cubature rules,...=
!#
!# - heatcond_matvec
!#   -> Contains routines to initialise matrices and vectors.
!#
!# - heatcond_bounadrycondition
!#   -> Contains routines to set up analytical and discrete boundary conditions.
!#
!# - heatcond_solver
!#   -> Configures the linear solver for the problem.
!#
!# - heatcond_timeloop
!#   -> Contains the actual time stepping algorithm and nonstationary solver.
!# </purpose>
!##############################################################################

module heatcond_method5

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use sortstrategy
  use triangulation
  use element
  use spatialdiscretisation
  use coarsegridcorrection
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use linearsolver
  use discretebc
  use ucd
  
  use collection
  use paramlist
    
  use heatcond_callback
  
  use heatcond_basic
  use heatcond_matvec
  use heatcond_boundarycondition
  use heatcond_partridiscr
  use heatcond_solver
  use heatcond_timeloop
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_initparameters (rparams,rproblem)
  
!<description>
  ! Reads the DAT file from disc into the parameter list rparams and
  ! initialises basic variables (number of levels, time stepping technique)
  ! in rproblem according to these settings.
!</description>

!<inputoutput>
  ! A parameter list structure accepting the parameters from the DAT file.
  type(t_parlist), intent(inout) :: rparams

  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: cscheme,niterations
  real(DP) :: dtheta,dtstep,dtimemin,dtimemax
  character(len=SYS_STRLEN) :: smaster,sstring

    ! Read the parameters from disc and put a reference to it
    ! to the collection
    call sys_getcommandLineArg(1,smaster,sdefault="./data/heatcond.dat")
    call parlst_readfromfile (rparams, smaster)

    ! We want to solve our Laplace problem on level...
    call parlst_getvalue_int (rparams, "GENERAL", "NLMIN", rproblem%ilvmin, 7)
    call parlst_getvalue_int (rparams, "GENERAL", "NLMAX", rproblem%ilvmax, 7)
    
    ! Get the parameters for the time stepping scheme from the parameter list
    call parlst_getvalue_int (rparams, "TIMESTEPPING", "CSCHEME", cscheme, 0)
    call parlst_getvalue_int (rparams, "TIMESTEPPING", "NITERATIONS", niterations, 1000)
    call parlst_getvalue_double (rparams, "TIMESTEPPING", "DTHETA", dtheta, 1.0_DP)
    call parlst_getvalue_double (rparams, "TIMESTEPPING", "DTSTEP", dtstep, 0.1_DP)
    call parlst_getvalue_double (rparams, "TIMESTEPPING", "DTIMEMIN", dtimemin, 0.0_DP)
    call parlst_getvalue_double (rparams, "TIMESTEPPING", "DTIMEMAX", dtimemax, 1.0_DP)
    
    ! Get the path where to write gmv`s to.
    call parlst_getvalue_string (rparams, "GENERAL", &
                                 "sucddir", sstring,"""./pre/QUAD.prm""")
    read(sstring,*) rproblem%sucddir
    
    ! Get the path of the prm/tri files.
    call parlst_getvalue_string (rparams, "GENERAL", &
                                 "sprmfile", sstring,"""./pre/QUAD.prm""")
    read(sstring,*) rproblem%sprmfile

    call parlst_getvalue_string (rparams, "GENERAL", &
                                 "strifile", sstring,"""./pre/QUAD.tri""")
    read(sstring,*) rproblem%strifile

    ! Initialise the time stepping in the problem structure
    call timstp_init (rproblem%rtimedependence%rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)
                     
    rproblem%rtimedependence%niterations = niterations
    
    rproblem%rtimedependence%dtimemin = dtimemin
    rproblem%rtimedependence%dtime = dtimemin
    rproblem%rtimedependence%dtimemax = dtimemax

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine heatcond5
  
!<description>
  ! This is a "separated" heatcond solver for solving a nonstationary heat
  ! conduction problem. The different tasks of the problem are separated into
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
  ! 7.) Write solution to VTK file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! A paramlist structure with parameters from the dat file
    type(t_parlist) :: rparams

    ! A problem structure for our problem
    type(t_problem), target :: rproblem
    
    ! An initial RHS vector and a solution vector
    type(t_vectorBlock) :: rrhs,rvector
    
    integer :: i
    
    ! Initialise the parameter list
    call parlst_init(rparams)
    
    ! Initialise the collection.
    call collct_init (rproblem%rcollection)
    call collct_setvalue_parlst (rproblem%rcollection, "PARAMS", rparams, .true.)
    
    ! Read in the parameters from the DAT file and initialise the basic
    ! structures with these.
    call hc5_initparameters (rparams,rproblem)

    ! Add space for level information in the collection
    do i=1,rproblem%ilvmax
      call collct_addlevel (rproblem%rcollection)
    end do
    
    ! Allocate memory for the level information
    allocate (rproblem%RlevelInfo(rproblem%ilvmax))
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call hc5_initParamTriang (rproblem%ilvmin,rproblem%ilvmax,rproblem)
    call hc5_initDiscretisation (rproblem)
    call hc5_initMatVec (rproblem,rparams)

    ! Use the auxiliary RHS vector on the finest level to create an
    ! initial RHS and solution vector, which we pass later to the timeloop.
    call lsysbl_createVecBlockIndirect (rproblem%rrhs,rrhs,.false.)
    call lsysbl_createVecBlockIndirect (rproblem%rrhs,rvector,.true.)

    ! Calculate the initial RHS, do not incorporate any BC`s.
    call hc5_calcRHS (rproblem,rrhs)
    
    ! Discretise the boundary conditions
    call hc5_initDiscreteBC (rproblem)
    
    ! Implement them into the initial solution vector, as we have a zero
    ! vector as initial solution.
    rvector%p_rdiscreteBC => rproblem%rrhs%p_rdiscreteBC
    call vecfil_discreteBCsol (rvector)
    
    ! Initialise the solver
    call hc5_initSolver (rproblem)
    
    ! Call the timeloop to solve the problem
    call hc5_timeloop (rproblem,rvector,rrhs)
    
    ! Release the solver, we don't need it anymore
    call hc5_doneSolver (rproblem)
    
    ! Cleanup
    call hc5_doneMatVec (rproblem)
    call hc5_doneBC (rproblem)
    call hc5_doneDiscretisation (rproblem)
    call hc5_doneParamTriang (rproblem)
    
    ! Release memory for level information
    deallocate (rproblem%RlevelInfo)
    
    ! Release parameter list
    call collct_deletevalue (rproblem%rcollection,"PARAMS")
    call parlst_done (rparams)
    
    ! Release RHS and solution vector
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)

    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ("Remaining collection statistics:")
    call output_line ("--------------------------------")
    call output_lbrk ()
    call collct_printStatistics (rproblem%rcollection)
    
    ! Finally release the collection.
    call collct_done (rproblem%rcollection)
    
  end subroutine

end module
