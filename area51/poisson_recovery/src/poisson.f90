!##############################################################################
!# ****************************************************************************
!# <name> poisson </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#              - Laplace(u) = f
!#
!# for a scalar function u.
!#
!# </purpose>
!##############################################################################

program poisson

  use collection
  use convergencetable
  use fsystem
  use genoutput
  use linearsystemblock
  use paramlist
  use storage

  use poisson2d

  ! local variables
  type(t_parlist) :: rparlist
  type(t_problem) :: rproblem
  type(t_convergenceTable) :: rtable
  character(len=SYS_STRLEN) :: slogdir,slogfile
  integer :: NLMIN,NLMAX,i

  ! Initialise system-wide settings
  call sys_init()

  ! Initialise the output system.
  !
  ! Normally, we write all the output to the screen and to a file
  ! "./log/output.txt".
  ! In the case that environment variables "$logdir"/"$resultsfile" exists,
  ! we write all the output to that file. This can be used e.g. in
  ! regression tests to compare results to reference results.
  if (sys_getenv_string("LOGDIR",slogdir) .and. &
      sys_getenv_string("RESULTFILE",slogfile)) then
    call output_init(trim(slogdir)//"/"//trim(slogfile))
  else
    call output_init("./log/output.txt")
  end if
  
  ! Initialise the FEAT 2.0 storage management
  call storage_init(999, 100)

  ! Read parameter file
  call parlst_init(rparlist)
  call parlst_readfromfile(rparlist, 'data/master.dat', 'data', .true.)

  ! Get levels from parameter list
  call parlst_getvalue_int(rparlist, '', 'NLMIN', NLMIN)
  call parlst_getvalue_int(rparlist, '', 'NLMAX', NLMAX)

  ! Initialise the collection
  call collct_init (rproblem%rcollection)
  do i=1,NLMAX
    call collct_addlevel (rproblem%rcollection)
  end do

  ! Initialise the convergence table
  call ctab_init(rtable)
  
  do i=NLMIN,NLMAX
    ! Initialisation
    call poisson_initParamTriang(rparlist,NLMIN,i,rproblem)
    call poisson_initDiscretisation(rparlist,rproblem)
    call poisson_initMatVec(rparlist,rproblem)
    call poisson_initDiscreteBC(rproblem)
    
    ! Implementation of boundary conditions
    call poisson_implementBC(rproblem)
    
    ! Solve the problem
    call poisson_solve(rproblem)

    ! Post-processing
    call poisson_postprocessing(rproblem,rtable)

    ! Cleanup
    call poisson_doneMatVec(rproblem)
    call poisson_doneBC(rproblem)
    call poisson_doneDiscretisation(rproblem)
    call poisson_doneParamTriang(rproblem)
  end do
  
  ! Compute reduction rates
  call ctab_evalConvergenceRate(rtable,"L2-error",CTAB_REDUCTION_RATE)
  call ctab_evalConvergenceRate(rtable,"H1-error",CTAB_REDUCTION_RATE)
  call ctab_evalConvergenceRate(rtable,"L2-error DerivX",CTAB_REDUCTION_RATE)
  call ctab_evalConvergenceRate(rtable,"L2-error DerivY",CTAB_REDUCTION_RATE)
  call ctab_evalConvergenceRate(rtable,"H1-error DerivX",CTAB_REDUCTION_RATE)
  call ctab_evalConvergenceRate(rtable,"H1-error DerivY",CTAB_REDUCTION_RATE)

  ! Adjust format of convergence table
  call ctab_setPrecision(rtable,"L2-error",3)
  call ctab_setPrecision(rtable,"L2-error-convrate",3)
  call ctab_setPrecision(rtable,"L2-error DerivX",3)
  call ctab_setPrecision(rtable,"L2-error DerivY",3)
  call ctab_setPrecision(rtable,"L2-error DerivX-convrate",3)
  call ctab_setPrecision(rtable,"L2-error DerivY-convrate",3)

  call ctab_setPrecision(rtable,"H1-error",3)
  call ctab_setPrecision(rtable,"H1-error-convrate",3)
  call ctab_setPrecision(rtable,"H1-error DerivX",3)
  call ctab_setPrecision(rtable,"H1-error DerivY",3)
  call ctab_setPrecision(rtable,"H1-error DerivX-convrate",3)
  call ctab_setPrecision(rtable,"H1-error DerivY-convrate",3)

  call ctab_setScientific(rtable,"L2-error",.true.)
  call ctab_setScientific(rtable,"L2-error DerivX",.true.)
  call ctab_setScientific(rtable,"L2-error DerivY",.true.)

  call ctab_setScientific(rtable,"H1-error",.true.)
  call ctab_setScientific(rtable,"H1-error DerivX",.true.)
  call ctab_setScientific(rtable,"H1-error DerivY",.true.)

  call ctab_setTexCaption(rtable,"cells","\# cells")
  call ctab_setTexCaption(rtable,"dofs","\# dofs")

  call ctab_setTexCaption(rtable,"L2-error","$L^2$")
  call ctab_setTexCaption(rtable,"L2-error-convrate","$L^2$-rate")
  call ctab_setTexCaption(rtable,"L2-error DerivX","$L^2\partial_x$")
  call ctab_setTexCaption(rtable,"L2-error DerivY","$L^2\partial_y$")
  call ctab_setTexCaption(rtable,"L2-error DerivX-convrate","$L^2\partial_x$-rate")
  call ctab_setTexCaption(rtable,"L2-error DerivY-convrate","$L^2\partial_y$-rate")

  call ctab_setTexCaption(rtable,"H1-error","$H^1$")
  call ctab_setTexCaption(rtable,"H1-error-convrate","$H^1$-rate")
  call ctab_setTexCaption(rtable,"H1-error DerivX","$H^1\partial_x$")
  call ctab_setTexCaption(rtable,"H1-error DerivY","$H^1\partial_y$")
  call ctab_setTexCaption(rtable,"H1-error DerivX-convrate","$H^1\partial_x$-rate")
  call ctab_setTexCaption(rtable,"H1-error DerivY-convrate","$H^1\partial_y$-rate")

  call ctab_setTexFormat(rtable,"cells","r")
  call ctab_setTexFormat(rtable,"dofs","r")

  call ctab_setTexTableCaption(rtable,"$L^2$-Convergence table")
  call ctab_setTexTableLabel(rtable,"tab:l2_convergence_rate")

  call ctab_setHidden(rtable,"H1-error",.true.)
  call ctab_setHidden(rtable,"H1-error-convrate",.true.)
  call ctab_setHidden(rtable,"H1-error DerivX",.true.)
  call ctab_setHidden(rtable,"H1-error DerivY",.true.)
  call ctab_setHidden(rtable,"H1-error DerivX-convrate",.true.)
  call ctab_setHidden(rtable,"H1-error DerivY-convrate",.true.)
  
  ! Write convergence table to Tex file
  call ctab_outputTex(rtable,'./table_l2.tex')

  call ctab_setTexTableCaption(rtable,"$H^1$-Convergence table")
  call ctab_setTexTableLabel(rtable,"tab:h1_convergence_rate")

  call ctab_setHidden(rtable,"H1-error",.false.)
  call ctab_setHidden(rtable,"H1-error-convrate",.false.)
  call ctab_setHidden(rtable,"H1-error DerivX",.false.)
  call ctab_setHidden(rtable,"H1-error DerivY",.false.)
  call ctab_setHidden(rtable,"H1-error DerivX-convrate",.false.)
  call ctab_setHidden(rtable,"H1-error DerivY-convrate",.false.)

  call ctab_setHidden(rtable,"L2-error",.true.)
  call ctab_setHidden(rtable,"L2-error-convrate",.true.)
  call ctab_setHidden(rtable,"L2-error DerivX",.true.)
  call ctab_setHidden(rtable,"L2-error DerivY",.true.)
  call ctab_setHidden(rtable,"L2-error DerivX-convrate",.true.)
  call ctab_setHidden(rtable,"L2-error DerivY-convrate",.true.)
  
  ! Write convergence table to Tex file
  call ctab_outputTex(rtable,'./table_h1.tex')

  ! Clear convergence table
  call ctab_done(rtable)

  ! Clear parameter list
  call parlst_done(rparlist)
  
  ! Print out heap statistics
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program poisson
