! SVSTD_IRESNORM - this parameter defines the norm in which residual/rhs should be measured
integer, parameter :: SVSTD_IRESNORM = 2

! SVSTD_IRESNORM - this parameter defines the norm in which solution changes should be measured
integer, parameter :: SVSTD_ISOLNORM = 2

! SVSTD_ICYCLE - this parameter defines the type of multigrid cycle
integer, parameter :: SVSTD_ICYCLE = 0

! SVSTD_NMINITERATIONS - this parameter defines the minimum number of iterations
integer, parameter :: SVSTD_NMINITERATIONS = 1

! SVSTD_NMAXITERATIONS - this parameter defines the maximum number of iterations
integer, parameter :: SVSTD_NMAXITERATIONS = 100

! SVSTD_NSMOOTHFACTOR - this parameter defines the factor for smoothing
!                       steps on different levels
integer, parameter :: SVSTD_NSMOOTHFACTOR = 1

! SVSTD_DOMEGA - this parameter defines the relaxation factor DOMEGA
real(DP), parameter :: SVSTD_DOMEGA = 1._DP

! SVSTD_DEPSREL - this parameter defines the relative tolerance for iterative solver
real(DP), parameter :: SVSTD_DEPSREL = 1E-5_DP

! SVSTD_DMAXREL - this parameter defines the maximum tolerance for relative changes
real(DP), parameter :: SVSTD_DMAXREL = 1E-1_DP

! SVSTD_DEPSABS - this parameter defines the absolute tolerance for iterative solver
real(DP), parameter :: SVSTD_DEPSABS = 1E-5_DP

! SVSTD_DDIVREL - this parameter defines the relative tolerance for divergence checking
real(DP), parameter :: SVSTD_DDIVREL = 1E6_DP

! SVSTD_DDIVABS - this parameter defines the absolute tolerance for divergence checking
real(DP), parameter :: SVSTD_DDIVABS = SYS_INFINITY

! SVSTD_DRHSZERO - this parameter defines the tolerance for zero r.h.s checking
real(DP), parameter :: SVSTD_DRHSZERO = 1E-90_DP

! SVSTD_DRHSZERO - this parameter defines the tolerance for zero defect checking
real(DP), parameter :: SVSTD_DDEFZERO = 1E-90_DP

! SVSTD_DEPSSTAG - this parameter defines the relative tolerance for stagnation checking
real(DP), parameter :: SVSTD_DEPSSTAG = 1E-5_DP

! SVSTD_DEPSILU - this parameter defines the tolerance for incomplete LU-decomposition
real(DP), parameter :: SVSTD_DEPSILU = 1E-10_DP

! SVSTD_IOLEVEL - this parameter defines the standard level for information output
integer, parameter :: SVSTD_IOLEVEL = 2

! SVSTD_IADAPTTIMESTEP - this parameter defines if adaptive time-stepping should be performed
integer, parameter :: SVSTD_IADAPTTIMESTEP = 0

! SVSTD_NPREADAPTSTEP - this parameter defines how many time steps are performed before adaptive time stepping
integer, parameter :: SVSTD_NPREADAPTSTEPS = 5

! SVSTD_DPROPORTIONALEXPONENT - this parameter is used as exponent of the proportional term in PID controller
real(DP), parameter :: SVSTD_DPROPORTIONALEXPONENT = 0.075_DP

! SVSTD_DINTEGRALEXPONENT - this parameter is used as exponent of the integral term in PID controller
real(DP), parameter :: SVSTD_DINTEGRALEXPONENT = 0.175_DP

! SVSTD_DDERIVATIVEEXPONENT - this parameter is used as exponent of the derivative term in PID controller
real(DP), parameter :: SVSTD_DDERIVATIVEEXPONENT = 0.01_DP

! SVSTD_DINCREASEFACTOR - this parameter defines the maximum factor by which the
!                         relative changes may vary in the PID controller
real(DP), parameter :: SVSTD_DINCREASEFACTOR = 2._DP

! SVSTD_DDECREASEFACTOR - this parameter defines the minimum factor by which the
!                         relative changes may vary in the time step controller
real(DP), parameter :: SVSTD_DDECREASEFACTOR = 0.5_DP

! SVSTD_STEPREDUCTIONFACTOR - this parameter defines the step reduction factor by which the
!                             current time step is scaled of the solution algorithm diverged
real(DP), parameter :: SVSTD_DSTEPREDUCTIONFACTOR = 0.5_DP

! SVSTD_DEPSSTEADY - this parameter defines the tolerance for steady state convergence
real(DP), parameter :: SVSTD_DEPSSTEADY = 0.0_DP

