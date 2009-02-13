program gendatfile

  use paramlist
  use storage
  use genoutput

  implicit none

  type(t_parlist) :: rparlist
  character(LEN=SYS_STRLEN) :: cbuffer, sparameterfile
  integer :: iunit,idata

  ! Initialize Feat2 subsystem
  call system_init()

  ! Get parameterfile name from command line arguments
  call get_command_argument(command_argument_count(), cbuffer)
  sparameterfile = adjustl(cbuffer)

  ! Initialize parameter list from file
  call parlst_init(rparlist)
  call parlst_readfromfile(rparlist, trim(sparameterfile))

  ! Open input file
  call io_openFileForReading(trim(sparameterfile), iunit)

  ! Open output file
  call output_init('output.dat')

  call output_line("# -*- mode: sh; -*-", OU_CLASS_MSG)
  call output_lbrk
  call output_line("simportdatafiles(1) =", OU_CLASS_MSG)
  call output_line("  'data/codire/default.dat'", OU_CLASS_MSG)
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("application = codire", OU_CLASS_MSG)
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[Codire]", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# file which contains the application specific data", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'InputOutput', 'indatfile', cbuffer)
  call output_line("indatfile = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of boundary condition for primal problem", OU_CLASS_MSG)
  call output_line("sprimalbdrcondname = boundary_conditions_primal", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of boundary condition for dual problem", OU_CLASS_MSG)
  call output_line("sdualbdrcondname = boundary_conditions_dual", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of spatial dimensions", OU_CLASS_MSG)
  call output_line("ndimension = 1", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# solution algorithm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'iflowtype', idata)
  select case(idata)
  case(0)
    call output_line("algorithm = transient_primal", OU_CLASS_MSG)
  case(1)
    call output_line("algorithm = stationary_primal", OU_CLASS_MSG)
  case(2)
    call output_line("algorithm = pseudotransient_primal", OU_CLASS_MSG)
  case default
    stop
  end select
  call output_lbrk
  call output_line("#-------------------------------------------------------------------------------", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of mass matrix", OU_CLASS_MSG)
  call output_line("# 0 = no mass matrix", OU_CLASS_MSG)
  call output_line("# 1 = consistent mass matrix", OU_CLASS_MSG)
  call output_line("# 2 = lumped mass matrix", OU_CLASS_MSG)
  if (idata == 1) then
    call output_line("imasstype = 0", OU_CLASS_MSG)
  else
    call output_line("imasstype = 2", OU_CLASS_MSG)
  end if
  call output_lbrk
  call output_line("# type of mass antidiffusion", OU_CLASS_MSG)
  call output_line("# 0 = no mass antidiffusion", OU_CLASS_MSG)
  call output_line("# 1 = consistent mass antidiffusion", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Convection', 'istabilisation', idata)
  select case(idata)
  case(0)
    call output_line("imassantidiffusion = 1", OU_CLASS_MSG)
  case(1,20)
    call output_line("imassantidiffusion = 0", OU_CLASS_MSG)
  case(10,11,12,21)
    call output_line("imassantidiffusion = 2", OU_CLASS_MSG)
  case default
    stop
  end select
  call output_lbrk
  call output_line("# type of flow velocity", OU_CLASS_MSG)
  call output_line("# 0 = zero velocity v=0", OU_CLASS_MSG)
  call output_line("# 1 = linear constant velocity v=v(x)", OU_CLASS_MSG)
  call output_line("# 2 = linear time-dependent velocity v=v(x,t)", OU_CLASS_MSG)
  call output_line("# 3 = Burgers equation in space-time", OU_CLASS_MSG)
  call output_line("# 4 = Buckley Leverett in space-time", OU_CLASS_MSG)
  call output_line("# 6 = Burgers' equation in 2D", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'ivelocitytype', idata)
  call output_line("ivelocitytype = "//trim(sys_siL(2,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of velocity field", OU_CLASS_MSG)
  call output_line("svelocityname = velocityfield", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of flow diffusion", OU_CLASS_MSG)
  call output_line("# 0 = zero diffusion", OU_CLASS_MSG)
  call output_line("# 1 = isotropic diffusion", OU_CLASS_MSG)
  call output_line("# 2 = anisotropic diffusion", OU_CLASS_MSG)
  call output_line("# 3 = variable diffusion", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'idiffusiontype', idata)
  call output_line("idiffusiontype = "//trim(sys_siL(2,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of diffusion tensor", OU_CLASS_MSG)
  call output_line("sdiffusionname = diffusiontensor", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of flow reaction", OU_CLASS_MSG)
  call output_line("# 0 = zero reactive term", OU_CLASS_MSG)
  call output_line("# 1 = analytical initial solution", OU_CLASS_MSG)
  call output_line("ireactiontype = 0", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of reactive term", OU_CLASS_MSG)
  call output_line("sreactionname = reactionterm", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of right-hand side vector", OU_CLASS_MSG)
  call output_line("# 0 = zero right-hand side", OU_CLASS_MSG)
  call output_line("# 1 = analytical right-hand side", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'irhstype', idata, 0)
  call output_line("irhstype = "//trim(sys_siL(2,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of right-hand side vector", OU_CLASS_MSG)
  call output_line("srhsname = rhs", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of initial solution profile", OU_CLASS_MSG)
  call output_line("# 0 = zero initial solution", OU_CLASS_MSG)
  call output_line("# 1 = analytical initial solution", OU_CLASS_MSG)
  call output_line("# 2 = PGM image", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'InputOutput', 'pgmfile', cbuffer, '')
  if (trim(cbuffer) .eq. '') then
    call output_line("isolutiontype = 1", OU_CLASS_MSG)
  else
    call output_line("isolutiontype = 2", OU_CLASS_MSG)
  end if
  call output_lbrk
  call output_line("# section name of initial solution", OU_CLASS_MSG)
  if (trim(cbuffer) .eq. '') then
    call output_line("ssolutionname = initial_solution", OU_CLASS_MSG)
  else
    call output_line("ssolutionname = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  end if
  call output_lbrk  
  call output_line("# type of exact solution profile", OU_CLASS_MSG)
  call output_line("# 0 = no exact solution available", OU_CLASS_MSG)
  call output_line("# 1 = analytical initial solution", OU_CLASS_MSG)
  call output_line("# 2 = PGM image", OU_CLASS_MSG)
  call output_line("iexactsolutiontype = 0", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of exact solution", OU_CLASS_MSG)
  call output_line("sexactsolutionname = exact_solution", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of target functional", OU_CLASS_MSG)
  call output_line("# 0 = no target variable", OU_CLASS_MSG)
  call output_line("# 1 = volume integral", OU_CLASS_MSG)
  call output_line("# 2 = surface integral", OU_CLASS_MSG)
  call output_line("itargetfunctype = 0", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of target functional", OU_CLASS_MSG)
  call output_line("stargetfuncname = targetfunctional", OU_CLASS_MSG)
  call output_lbrk
  call output_line("#-------------------------------------------------------------------------------", OU_CLASS_MSG)
  call output_lbrk
!!$  call output_line("# file which contains the boundary parametrisation", OU_CLASS_MSG)
!!$  call parlst_getvalue_string(rparlist, 'InputOutput', 'prmfile', cbuffer)
!!$  call output_line("prmfile = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
!!$  call output_lbrk
  call output_line("# file which contains the triangulation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'InputOutput', 'trifile', cbuffer)
  call output_line("trifile = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of finite element space", OU_CLASS_MSG)
  call output_line("#  1 = P1 finite elements", OU_CLASS_MSG)
  call output_line("#  2 = P2 finite elements", OU_CLASS_MSG)
  call output_line("# 11 = Q1 finite elements", OU_CLASS_MSG)
  call output_line("# 12 = Q2 finite elements", OU_CLASS_MSG)
  call output_line("# -1 = mixed P1/Q1 finite elements", OU_CLASS_MSG)
  call output_line("# -2 = mixed P2/Q2 finite elements", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'ieltype', idata)
  call output_line("ieltype = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# convert mesh to triangular mesh?", OU_CLASS_MSG)
!!$  call parlst_getvalue_int(rparlist, 'Benchmark', 'iconvtotria', idata)
!!$  call output_line("iconvToTria = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
!!$  call output_lbrk
  call output_line("# type of matrix format", OU_CLASS_MSG)
  call output_line("# 7 = matrix format 7", OU_CLASS_MSG)
  call output_line("# 9 = matrix format 9", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'imatrixformat', idata)
  call output_line("imatrixformat = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of Jacobian matrix", OU_CLASS_MSG)
  call output_line("# 0 = standard sparsity pattern", OU_CLASS_MSG)
  call output_line("# 1 = extended sparsity pattern", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Convection', 'iextendedJacobian', idata)
  call output_line("ijacobianformat = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("#-------------------------------------------------------------------------------", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of the convection stabilization", OU_CLASS_MSG)
  call output_line("convection = Convection", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of the diffusion stabilization", OU_CLASS_MSG)
  call output_line("diffusion = Diffusion", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of the time-stepping algorithm", OU_CLASS_MSG)
  call output_line("timestep = Timestepping", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of the top-level solver", OU_CLASS_MSG)
  call output_line("solver = FullMultigridSolver", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of the output configuration", OU_CLASS_MSG)
  call output_line("output = Output", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of the adaptation configuration", OU_CLASS_MSG)
  call output_line("adaptivity = Adaptivity", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# section name of the error estimation configuration", OU_CLASS_MSG)
  call output_line("errorestimator = ErrorEstimator", OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[Convection]", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of spatial stabilization:", OU_CLASS_MSG)
  call output_line("#  0 = no stabilization (Galerkin)", OU_CLASS_MSG)
  call output_line("#  1 = discrete upwinding", OU_CLASS_MSG)
  call output_line("# 10 = semi-impl. FEM-FCT", OU_CLASS_MSG)
  call output_line("# 11 = semi-expl. FEM-FCT", OU_CLASS_MSG)
  call output_line("# 12 = linearized FEM-FCT", OU_CLASS_MSG)
  call output_line("# 20 = FEM-TVD", OU_CLASS_MSG)
  call output_line("# 21 = FEM-GP", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Convection', 'istabilisation', idata)
  call output_line("istabilisation = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[Diffusion]", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk  
  call output_line("# type of spatial stabilization", OU_CLASS_MSG)
  call output_line("#  0 = no stabilization (Galerkin)", OU_CLASS_MSG)
  call output_line("#  2 = maximum principle preservation", OU_CLASS_MSG)
  call output_line("# 30 = symmetric limiting", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Diffusion', 'istabilisation', idata)
  call output_line("istabilisation = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[Output]", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of UCD output format", OU_CLASS_MSG)
  call output_line("# 0 = no output", OU_CLASS_MSG)
  call output_line("# 1 = GMV", OU_CLASS_MSG)
  call output_line("# 2 = AVS", OU_CLASS_MSG)
  call output_line("# 3 = VTK", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'InputOutput', 'iformatUCD', idata)
  call output_line("iformatUCD = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# time interval for UCD output", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'InputOutput', 'dstepUCD', cbuffer)
  call output_line("dstepUCD = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# file for UCD output of solution", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'InputOutput', 'sfilenameUCDsolution', cbuffer)
  call output_line("ucdsolution = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# file for UCD output of error", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'InputOutput', 'sfilenameUCDerror', cbuffer)
  call output_line("ucderror = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[Timestepping]", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# level of output information:", OU_CLASS_MSG)
  call output_line("# 0 = no output", OU_CLASS_MSG)
  call output_line("# 1 = errors", OU_CLASS_MSG)
  call output_line("# 2 = errors+warnings", OU_CLASS_MSG)
  call output_line("# 3 = information", OU_CLASS_MSG)
  call output_line("# 4 = verbose output", OU_CLASS_MSG)
  call output_line("# >9 = file output", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Timestepping', 'ioutputlevel', idata)
  call output_line("ioutputlevel = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of time-stepping algorithm", OU_CLASS_MSG)
  call output_line("# 5 = two-level theta-scheme", OU_CLASS_MSG)
  call output_line("# 6 = explicit Runge-Kutta scheme", OU_CLASS_MSG)
  call output_line("# 7 = full MG", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Timestepping', 'ctimestepType', idata)
  call output_line("ctimestepType = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# norm to use for solution variation checking", OU_CLASS_MSG) 
  call output_line("# 0 = Euclidian norm", OU_CLASS_MSG)
  call output_line("# 1 = L1-norm", OU_CLASS_MSG)
  call output_line("# 2 = L2-norm", OU_CLASS_MSG)
  call output_line("# 3 = MAX-norm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Timestepping', 'isolNorm', idata)
  call output_line("isolNorm = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# implicitness parameter: Valid values of THETA are in the range [0,1]", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'theta', cbuffer)
  call output_line("theta = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of steps for multisteps scheme", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Timestepping', 'multisteps', idata)
  call output_line("multisteps = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# initial time for the simulation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'dinitialTime', cbuffer)
  call output_line("dinitialTime = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# final time for the simulation.", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'dfinalTime', cbuffer)
  call output_line("dfinalTime = ", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# initial time step size", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'dinitialStep', cbuffer)
  call output_line("dinitialStep = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# lower bound for the admissible time step", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'dminStep', cbuffer)
  call output_line("dminStep = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# upper bound for the admissible time step", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'dmaxStep', cbuffer)
  call output_line("dmaxStep = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# adaptive time stepping algorithm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Timestepping', 'iadaptTimestep', idata, 0)
  call output_line("iadaptTimestep = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum factor by which time step may change", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'ddecreaseFactor', cbuffer)
  call output_line("ddecreaseFactor = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum factor by which time step may change", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'dincreaseFactor', cbuffer)
  call output_line("dincreaseFactor = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# reduction factor by which time step is reduced if simulation fails", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'dstepReductionFactor', cbuffer)
  call output_line("dstepReductionFactor = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# target tolerace for relative changes", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'depsRel', cbuffer, '0.0')
  call output_line("depsRel = ", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# absolute tolerance for relative changes", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'depsAbs', cbuffer, '0.0')
  call output_line("depsAbs = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# tolerance for steady state convergence", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Timestepping', 'depsSteady', cbuffer, '0.0')
  call output_line("depsSteady = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk

  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[FullMultigridSolver]", OU_CLASS_MSG)
  call output_line(" csolverType = 7", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# level of output information:", OU_CLASS_MSG)
  call output_line("# 0 = no output", OU_CLASS_MSG)
  call output_line("# 1 = errors", OU_CLASS_MSG)
  call output_line("# 2 = errors+warnings", OU_CLASS_MSG)
  call output_line("# 3 = information", OU_CLASS_MSG)
  call output_line("# 4 = verbose output", OU_CLASS_MSG)
  call output_line("# >9 = file output", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'FullMultigridSolver', 'ioutputlevel', idata)
  call output_line("ioutputlevel = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum multigrid level number", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'FullMultigridSolver', 'nlmin', idata)
  call output_line("nlmin = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum multigrid level number", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'FullMultigridSolver', 'nlmax', idata)
  call output_line("nlmax = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of multigrid cycle", OU_CLASS_MSG)
  call output_line("# 0 = F-cylce", OU_CLASS_MSG)
  call output_line("# 1 = V-cycle", OU_CLASS_MSG)
  call output_line("# 2 = W-cycle", OU_CLASS_MSG)
  call output_line("# 3,4,5 = corresponding saw tooth variants", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'FullMultigridSolver', 'icycle', idata)
  call output_line("icycle = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum number of multigrid steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'FullMultigridSolver', 'ilmin', idata)
  call output_line("ilmin = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum number of multigrid steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'FullMultigridSolver', 'ilmax', idata)
  call output_line("ilmax = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# nonlinear subsolver", OU_CLASS_MSG)
  call output_line("# 1 = nonlinear single-grid solver", OU_CLASS_MSG)
  call output_line("# 3 = nonlinear multigrid solver", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'FullMultigridSolver', 'isolver', idata)
  call output_line("isolver = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of the subsolver", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'FullMultigridSolver', 'ssolvername', cbuffer)
  call output_line("ssolvername = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[NonlinearMultigridSolver]", OU_CLASS_MSG)
  call output_line(" csolverType = 3", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# level of output information:", OU_CLASS_MSG)
  call output_line("# 0 = no output", OU_CLASS_MSG)
  call output_line("# 1 = errors", OU_CLASS_MSG)
  call output_line("# 2 = errors+warnings", OU_CLASS_MSG)
  call output_line("# 3 = information", OU_CLASS_MSG)
  call output_line("# 4 = verbose output", OU_CLASS_MSG)
  call output_line("# >9 = file output", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'ioutputlevel', idata)
  call output_line("ioutputlevel = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# norm to use for defect checking", OU_CLASS_MSG)
  call output_line("# 0 = Euclidian norm", OU_CLASS_MSG)
  call output_line("# 1 = L1-norm", OU_CLASS_MSG)
  call output_line("# 2 = L2-norm", OU_CLASS_MSG)
  call output_line("# 3 = MAX-norm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'iresNorm', idata)
  call output_line("iresNorm = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum multigrid level number", OU_CLASS_MSG)
  call output_line("# If NLMIN = NLMAX then the nonlinear problem is solved by a single-grid solver", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'nlmin', idata)
  call output_line("nlmin = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum multigrid level number (see NLMIN)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'nlmax', idata)
  call output_line("nlmax = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum number of linear mg steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'ilmin', idata)
  call output_line("ilmin = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum number of linear mg steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'ilmax', idata)
  call output_line("ilmax = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of multigrid cycle", OU_CLASS_MSG)
  call output_line("# 0 = F-cylce", OU_CLASS_MSG)
  call output_line("# 1 = V-cycle", OU_CLASS_MSG)
  call output_line("# 2 = W-cycle", OU_CLASS_MSG)
  call output_line("# 3,4,5 =  corresponding saw tooth variants", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'icycle', idata)
  call output_line("icycle = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# nonlinear coarse grid solver", OU_CLASS_MSG)
  call output_line("# 1 = nonlinear single-grid solver", OU_CLASS_MSG)
  call output_line("# 3 = nonlinear multigrid solver", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'isolver', idata)
  call output_line("isolver = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of nonlinear coarse grid solver", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearMultigridSolver', 'ssolvername', cbuffer)
  call output_line("ssolvername = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# smoother for the nonlinear multigrid solver", OU_CLASS_MSG)
  call output_line("# =1, nonlinear single-grid smoother, =3 nonlinear multigrid smoother", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'ismoother', idata)
  call output_line("ismoother = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of nonlinear smoother", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearMultigridSolver', 'ssmoothername', cbuffer)
  call output_line("ssmoothername = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of presmoothing steps (if saw tooth variant is used, then no", OU_CLASS_MSG)
  call output_line("# presmoothing steps are performed in the coarse-to-fine part of the cycle)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'npresmooth', idata)
  call output_line("npresmooth = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of postsmoothing steps (if saw tooth variant is used, then no", OU_CLASS_MSG)
  call output_line("# postsmoothing steps are performed in the coarse-to-fine part of the cycle)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'npostsmooth', idata)
  call output_line("npostsmooth = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# factor for pre/postsm. on coarser levels (On each level l, the number of", OU_CLASS_MSG)
  call output_line("# smoothing steps is computed from SM_l=SM_L*NSFAC**(L-l), where L stands", OU_CLASS_MSG)
  call output_line("# for the finest grid level and SM_L is the prescribed number of smoothing steps)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearMultigridSolver', 'nsfac', idata)
  call output_line("nsfac = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# nonlinear relaxation parameter", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearMultigridSolver', 'domega', cbuffer, '1.0')
  call output_line("domega = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# absolute tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearMultigridSolver', 'depsAbs', cbuffer, '0.0')
  call output_line("depsAbs = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearMultigridSolver', 'depsRel', cbuffer, '0.0')
  call output_line("depsRel = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for stagnation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearMultigridSolver', 'depsStag', cbuffer, '0.0')
  call output_line("depsStag = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[NonlinearSolver]", OU_CLASS_MSG)
  call output_line(" csolverType = 1", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# level of output information:", OU_CLASS_MSG)
  call output_line("# 0 = no output", OU_CLASS_MSG)
  call output_line("# 1 = errors", OU_CLASS_MSG)
  call output_line("# 2 = errors+warnings", OU_CLASS_MSG)
  call output_line("# 3 = information", OU_CLASS_MSG)
  call output_line("# 4 = verbose output", OU_CLASS_MSG)
  call output_line("# >9 = file output", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'ioutputlevel', idata)
  call output_line("ioutputlevel = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# norm to use for defect checking", OU_CLASS_MSG)
  call output_line("# 0 = Euclidian norm", OU_CLASS_MSG)
  call output_line("# 1 = L1-norm", OU_CLASS_MSG)
  call output_line("# 2 = L2-norm", OU_CLASS_MSG)
  call output_line("# 3 = MAX-norm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'iresNorm', idata)
  call output_line("iresNorm = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum number of nonlinear steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'nminIterations', idata)
  call output_line("nminIterations = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum number of nonlinear steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'nmaxIterations', idata)
  call output_line("nmaxIterations = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# preconditioner for the nonlinear solver", OU_CLASS_MSG)
  call output_line("# 1 = block-diagonal preconditioner", OU_CLASS_MSG)
  call output_line("# 2 = defect-correction algorithm", OU_CLASS_MSG)
  call output_line("# 3 = algebraic Newton algorithm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'iprecond', idata)
  call output_line("iprecond = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# nonlinear solver", OU_CLASS_MSG)
  call output_line("# 101 = fixed-point iteration", OU_CLASS_MSG)
  call output_line("# 102 = predictor-corrector fixed-point iteration", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'isolver', idata)
  call output_line("isolver = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# strategy for choosing the perturbation parameter in Newton's method", OU_CLASS_MSG)
  call output_line("# 1 = NITSOL", OU_CLASS_MSG)
  call output_line("# 2 = SQRT(EPS)", OU_CLASS_MSG)
  call output_line("# otherwise, user-defined value", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearSolver', 'dperturbationStrategy', cbuffer)
  call output_line("dperturbationStrategy = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# strategy for choosing the forcing term in Newton's method", OU_CLASS_MSG)
  call output_line("# 1 = choice 1 by Eisenstat/Walker", OU_CLASS_MSG)
  call output_line("# 2 = choice 2 by Eisenstat/Walker", OU_CLASS_MSG)
  call output_line("# 3 = choice by Brown/Saad", OU_CLASS_MSG)
  call output_line("# 4 = choice by Dembo/Steihaug", OU_CLASS_MSG)
  call output_line("# otherwise, user-defined fixed value", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearSolver', 'dforcingStrategy', cbuffer)
  call output_line("dforcingStrategy = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# check sufficient decrease condition in globalization", OU_CLASS_MSG)
  call output_line("# 0 = apply Newton increment without globalization", OU_CLASS_MSG)
  call output_line("# 1 = check sufficient decrease condition and perform backtracking if required", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'icheckSufficientDecrease', idata)
  call output_line("icheckSufficientDecrease = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum number of backtracking steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'nmaxBacktrackingSteps', idata)
  call output_line("nmaxBacktrackingSteps = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# update frequency of Jacobian", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'NonlinearSolver', 'iupdateFrequency', idata)
  call output_line("iupdateFrequency = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# nonlinear relaxation parameter", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearSolver', 'domega', cbuffer, '1.0')
  call output_line("domega = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# absolute tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearSolver', 'depsAbs', cbuffer, '0.0')
  call output_line("depsAbs = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearSolver', 'depsRel', cbuffer, '0.0')
  call output_line("depsRel = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for stagnation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearSolver', 'depsStag', cbuffer, '0.0')
  call output_line("depsStag = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of linear sub-solver", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'NonlinearSolver', 'ssolvername', cbuffer)
  call output_line("ssolvername = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[LinearMultigridSolver]", OU_CLASS_MSG)
  call output_line(" csolverType = 4", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# level of output information:", OU_CLASS_MSG)
  call output_line("# 0 = no output", OU_CLASS_MSG)
  call output_line("# 1 = errors", OU_CLASS_MSG)
  call output_line("# 2 = errors+warnings", OU_CLASS_MSG)
  call output_line("# 3 = information", OU_CLASS_MSG)
  call output_line("# 4 = verbose output", OU_CLASS_MSG)
  call output_line("# >9 = file output", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'ioutputlevel', idata)
  call output_line("ioutputlevel = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# norm to use for defect checking", OU_CLASS_MSG)
  call output_line("# 0 = Euclidian norm", OU_CLASS_MSG)
  call output_line("# 1 = L1-norm", OU_CLASS_MSG)
  call output_line("# 2 = L2-norm", OU_CLASS_MSG)
  call output_line("# 3 = MAX-norm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'iresNorm', idata)
  call output_line("iresNorm = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum multigrid level number", OU_CLASS_MSG)
  call output_line("# If NLMIN = NLMAX then the nonlinear problem is solved by a single-grid solver", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'nlmin', idata)
  call output_line("nlmin = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum multigrid level number (see NLMIN)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'nlmax', idata)
  call output_line("nlmax = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum number of linear mg steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'ilmin', idata)
  call output_line("ilmin = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum number of linear mg steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'ilmax', idata)
  call output_line("ilmax = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of multigrid cycle", OU_CLASS_MSG)
  call output_line("# 0 = F-cylce", OU_CLASS_MSG)
  call output_line("# 1 = V-cycle", OU_CLASS_MSG)
  call output_line("# 2 = W-cycle", OU_CLASS_MSG)
  call output_line("# 3,4,5 =  corresponding saw tooth variants", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'icycle', idata)
  call output_line("icycle = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# linear coarse grid solver", OU_CLASS_MSG)
  call output_line("# 2 = linear single-grid solver", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'isolver', idata)
  call output_line("isolver = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of linear coarse grid solver", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearMultigridSolver', 'ssolvername', cbuffer)
  call output_line("ssolvername = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# smoother for the linear multigrid solver", OU_CLASS_MSG)
  call output_line("# 2 = linear single-grid solver", OU_CLASS_MSG)
  call output_line("# 4 = linear multigrid solver", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'ismoother', idata)
  call output_line("ismoother = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of smoother", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearMultigridSolver', 'ssmoothername', cbuffer)
  call output_line("ssmoothername = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of presmoothing steps (if saw tooth variant is used, then no", OU_CLASS_MSG)
  call output_line("# presmoothing steps are performed in the coarse-to-fine part of the cycle)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'npresmooth', idata)
  call output_line("npresmooth = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of postsmoothing steps (if saw tooth variant is used, then no", OU_CLASS_MSG)
  call output_line("# postsmoothing steps are performed in the coarse-to-fine part of the cycle)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'npostsmooth', idata)
  call output_line("npostsmooth = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# factor for pre/postsm. on coarser levels (On each level l, the number of", OU_CLASS_MSG)
  call output_line("# smoothing steps is computed from SM_l=SM_L*NSFAC**(L-l), where L stands", OU_CLASS_MSG)
  call output_line("# for the finest grid level and SM_L is the prescribed number of smoothing steps)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearMultigridSolver', 'nsmoothfactor', idata)
  call output_line("nsmoothfactor = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# linear relaxation parameter", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearMultigridSolver', 'domega', cbuffer, '1.0')
  call output_line("domega = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# absolute tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearMultigridSolver', 'depsAbs', cbuffer, '0.0')
  call output_line("depsAbs = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearMultigridSolver', 'depsRel', cbuffer, '0.0')
  call output_line("depsRel = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for stagnation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearMultigridSolver', 'depsStag', cbuffer, '0.0')
  call output_line("depsStag = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[LinearSmoother]", OU_CLASS_MSG)
  call output_line(" csolverType = 2", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# smoother for the linear multigrid solver", OU_CLASS_MSG)
  call output_line("#  2 = Jacobi", OU_CLASS_MSG)
  call output_line("#  4 = SOR", OU_CLASS_MSG)
  call output_line("#  5 = SSOR", OU_CLASS_MSG)
  call output_line("#  7 = BiCGSTAB", OU_CLASS_MSG)
  call output_line("#  8 = FGMRES", OU_CLASS_MSG)
  call output_line("# 50 = ILU", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSmoother', 'isolver', idata)
  call output_line("isolver = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# dimension of the Krylov subspace for FGMRES method", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSmoother', 'nkrylov', idata)
  call output_line("nkrylov = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# tolerance of (M)ILU-preconditioner", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSmoother', 'depsILU', cbuffer, '0.0')
  call output_line("depsILU = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# size of fill-in for (M)ILU(s)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSmoother', 'ifill', idata)
  call output_line("ifill = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# (M)ILU(s) relaxation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSmoother', 'domega', cbuffer, '1.0')
  call output_line("domega = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of preconditioner", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSmoother', 'sprecondName', cbuffer)
  call output_line("sprecondName = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[LinearSolver]", OU_CLASS_MSG)
  call output_line(" csolverType = 2", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# level of output information:", OU_CLASS_MSG)
  call output_line("# 0 = no output", OU_CLASS_MSG)
  call output_line("# 1 = errors", OU_CLASS_MSG)
  call output_line("# 2 = errors+warnings", OU_CLASS_MSG)
  call output_line("# 3 = information", OU_CLASS_MSG)
  call output_line("# 4 = verbose output", OU_CLASS_MSG)
  call output_line("# >9 = file output", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSolver', 'ioutputlevel', idata)
  call output_line("ioutputlevel = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# norm to use for defect checking", OU_CLASS_MSG)
  call output_line("# 0 = Euclidian norm", OU_CLASS_MSG)
  call output_line("# 1 = L1-norm", OU_CLASS_MSG)
  call output_line("# 2 = L2-norm", OU_CLASS_MSG)
  call output_line("# 3 = MAX-norm", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSolver', 'iresNorm', idata)
  call output_line("iresNorm = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# minimum number of linear steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSolver', 'nminIterations', idata)
  call output_line("nminIterations = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum number of linear steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSolver', 'nmaxIterations', idata)
  call output_line("nmaxIterations = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of linear solver", OU_CLASS_MSG)
  call output_line("#  2 = Jacobi", OU_CLASS_MSG)
  call output_line("#  4 = SOR", OU_CLASS_MSG)
  call output_line("#  5 = SSOR", OU_CLASS_MSG)
  call output_line("#  7 = BiCGSTAB", OU_CLASS_MSG)
  call output_line("#  8 = FGMRES", OU_CLASS_MSG)
  call output_line("# 11 = UMFPACK4", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSolver', 'isolver', idata)
  call output_line("isolver = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# dimension of the Krylov subspace for FGMRES method", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearSolver', 'nkrylov', idata)
  call output_line("nkrylov = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# linear relaxation parameter", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSolver', 'domega', cbuffer, '1.0')
  call output_line("domega = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# absolute tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSolver', 'depsAbs', cbuffer, '0.0')
  call output_line("depsAbs = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for residual", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSolver', 'depsRel', cbuffer, '0.0')
  call output_line("depsRel = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# relative tolerance for stagnation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSolver', 'depsStag', cbuffer, '0.0')
  call output_line("depsStag = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# name of preconditioner", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearSolver', 'sprecondName', cbuffer)
  call output_line("sprecondName = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[LinearPrecond]", OU_CLASS_MSG)
  call output_line(" csolverType = 2", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# preconditioner for the linear single-grid solver", OU_CLASS_MSG)
  call output_line("#  2 = Jacobi", OU_CLASS_MSG)
  call output_line("#  4 = SOR", OU_CLASS_MSG)
  call output_line("#  5 = SSOR", OU_CLASS_MSG)
  call output_line("# 50 = ILU", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearPrecond', 'isolver', idata)
  call output_line("isolver = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# tolerance of (M)ILU-preconditioner", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearPrecond', 'depsILU', cbuffer, '0.0')
  call output_line("depsILU = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# size of fill-in for (M)ILU(s)", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'LinearPrecond', 'ifill', idata)
  call output_line("ifill = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# (M)ILU(s) relaxation", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'LinearPrecond', 'domega', cbuffer)
  call output_line("domega = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[Adaptivity]", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of pre-adaptation steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'npreadapt', idata)
  call output_line("npreadapt = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of adaptation steps", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Benchmark', 'nadapt', idata)
  call output_line("nadapt = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# time for grid adaptivity", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Adaptivity', 'dtimeadapt', cbuffer)
  call output_line("dtimeadapt = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# time step for grid adaptivity", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Adaptivity', 'dstepadapt', cbuffer)
  call output_line("dstepadapt = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# maximum number of refinement levels", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Adaptivity', 'nsubdividemax', idata)
  call output_line("nsubdividemax = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# adaptation strategy: 1 = red-green", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'Adaptivity', 'iadaptationStrategy', idata)
  call output_line("iadaptationStrategy = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# refinement tolerance", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Adaptivity', 'drefinementTolerance', cbuffer)
  call output_line("drefinementTolerance = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# coarsening tolerance", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'Adaptivity', 'dcoarseningTolerance', cbuffer)
  call output_line("dcoarseningTolerance = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_lbrk
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_line("[ErrorEstimator]", OU_CLASS_MSG)
  call output_line("################################################################################", OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type of error estimator", OU_CLASS_MSG)
  call output_line("# 1 = L2-projection", OU_CLASS_MSG)
  call output_line("# 2 = node-based SPR", OU_CLASS_MSG)
  call output_line("# 3 = element-based SPR", OU_CLASS_MSG)
  call output_line("# 4 = face-based SPR", OU_CLASS_MSG)
  call output_line("# 5 = limited gradient averaging", OU_CLASS_MSG)
  call output_line("# 6 = second-difference indicator", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'ErrorEstimator', 'ierrorestimator', idata)
  call output_line("ierrorestimator = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# type or grid indicator", OU_CLASS_MSG)
  call output_line("# 0 = as is", OU_CLASS_MSG)
  call output_line("# 1 = equidistribution", OU_CLASS_MSG)
  call output_line("# 2 = logarithmic equidistribution", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'ErrorEstimator', 'igridindicator', idata)
  call output_line("igridindicator = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# noise filter ", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'ErrorEstimator', 'dnoisefilter', cbuffer)
  call output_line("dnoisefilter = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# absolute filter", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'ErrorEstimator', 'dabsfilter', cbuffer)
  call output_line("dabsfilter = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# refinement tolerance", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'ErrorEstimator', 'drefinementTolerance', cbuffer)
  call output_line("drefinementTolerance = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# coarsening tolerance", OU_CLASS_MSG)
  call parlst_getvalue_string(rparlist, 'ErrorEstimator', 'dcoarseningTolerance', cbuffer)
  call output_line("dcoarseningTolerance = "//trim(adjustl(cbuffer)), OU_CLASS_MSG)
  call output_lbrk
  call output_line("# number of protection layers", OU_CLASS_MSG)
  call parlst_getvalue_int(rparlist, 'ErrorEstimator', 'nprotectLayers', idata)
  call output_line("nprotectLayers = "//trim(sys_siL(idata,5)), OU_CLASS_MSG)
  
  ! Close output file
  call output_done()


  print *, "!!!"
  print *, "!!! Successfully converted parameter file"
  print *, "!!!"

  ! Release parameter list
  call parlst_done(rparlist)

  ! Release storage
  call storage_info(.true.)
  call storage_done()

end program gendatfile
